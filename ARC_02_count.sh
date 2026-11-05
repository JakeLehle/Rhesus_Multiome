#!/usr/bin/env bash
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
# Job name and log paths are set per run on the sbatch command line, which
# overrides anything written here.
# =============================================================================
# ARC_02_count.sh
#
# One self-contained script that both plans and runs the 20 cellranger-arc
# count jobs. No array jobs: each library gets its own independent SLURM job,
# resubmittable on its own.
#
#   ./ARC_02_count.sh                  print the run table and verify paths
#   ./ARC_02_count.sh submit first     submit one brain library, for validation
#   ./ARC_02_count.sh submit brain      submit the 8 brain libraries
#   ./ARC_02_count.sh submit colon      submit the 12 colon libraries
#   ./ARC_02_count.sh submit all        submit all 20
#   ./ARC_02_count.sh status            which runs are done, running, failed
#   ./ARC_02_count.sh libraries         write every libraries.csv and stop
#
# `submit` re-submits THIS script with `run <run_id>`, so there is one file to
# maintain and the worker logic cannot drift from the plan.
#
# -----------------------------------------------------------------------------
# WHY 20 COUNT RUNS, AND WHY THE RUN IDS LOOK LIKE THIS
# -----------------------------------------------------------------------------
#
# count operates on ONE library: one GEX and one ATAC FASTQ pair from the same
# nuclei. So 8 brain plus 12 colon is 20 runs however the analysis is grouped.
#
# The run IDs are the vendor's own sample labels from the loom Samples
# attribute. ARC_00 recovered their library ordering and it is simply the
# alphabetical sort of those labels:
#
#   Brain   _1.._8   40872 .. 40879
#   Colon   _1.._12  40872N 40873N 40874N 40874P 40875N 40875P
#                    40876N 40877N 40877P 40878N 40878P 40879N
#
# That matters because the vendor's CellID format is BARCODE-1_<index>, which
# is Seurat merge() appending a list position, NOT aggr's BARCODE-<n>. So they
# ran count per library, read each filtered_feature_bc_matrix, and merged per
# tissue. Reading these outputs in the order below reproduces their cell names
# exactly, and the annotation transfer becomes a plain join rather than a
# barcode remap.
#
# -----------------------------------------------------------------------------
# RESUME BEHAVIOUR
# -----------------------------------------------------------------------------
#
# count resumes a partial run if invoked again with the same --id in the same
# working directory. This script distinguishes three states rather than
# refusing on any existing directory:
#
#   outs/ present and complete   skip, already done
#   _lock present                refuse, another job may be running
#   partial, no lock             resume, which is what count is designed to do
#
# Jake Lehle, 2026-10
# =============================================================================

set -uo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

FASTQ_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome/Raw Data_Colon_Brain"
REFERENCE="/master/jlehle/WORKING/REFS/Mmul8/Mmul8"

ANALYSIS_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
LIB_DIR="${ANALYSIS_ROOT}/02_arc/01_libraries"
COUNT_ROOT="${ANALYSIS_ROOT}/02_arc/02_count"

CONDA_BASE="/master/jlehle/anaconda3"
CONDA_ENV="sc_pre_R"
ARC_DIR="/master/jlehle/cellranger-arc-2.0.2"
LOG_DIR="/master/jlehle/WORKING/LOGS"

SLURM_PARTITION="normal"
SLURM_CPUS=16
SLURM_MEM_GB=128
SLURM_TIME="48:00:00"
# Kept under the SLURM allocation so the pipeline has headroom.
ARC_LOCALMEM_GB=120

# FASTQ directory layout, from the delivered tree.
BRAIN_ATAC_DIR="${FASTQ_ROOT}/Brain data/Mohan_M_ATAC_03132023"
BRAIN_GEX_DIR="${FASTQ_ROOT}/Brain data/Mohan_M_GEX_03132023"
COLON_ATAC_DIR="${FASTQ_ROOT}/Colon data/Mohan_M_Colon_ATAC_04032023"
COLON_GEX_DIR="${FASTQ_ROOT}/Colon data/Mohan_M_Colon_GEX_04032023"

# -----------------------------------------------------------------------------
# The run table
#
# run_id | tissue | library_index | atac_sample | gex_sample
#
# library_index is the vendor's Seurat merge() position, recovered by ARC_00.
# Listed in that order so sorting is never needed again.
# -----------------------------------------------------------------------------

RUNS=(
  # --- Brain, 8 libraries, single timepoint ---
  "Brain_40872|Brain|1|BG_40872_arcATAC|BG_40872_arcGEX"
  "Brain_40873|Brain|2|BG_40873_arcATAC|BG_40873_arcGEX"
  "Brain_40874|Brain|3|BG_40874_arcATAC|BG_40874_arcGEX"
  "Brain_40875|Brain|4|BG_40875_arcATAC|BG_40875_arcGEX"
  "Brain_40876|Brain|5|BG_40876_arcATAC|BG_40876_arcGEX"
  "Brain_40877|Brain|6|BG_40877_arcATAC|BG_40877_arcGEX"
  "Brain_40878|Brain|7|BG_40878_arcATAC|BG_40878_arcGEX"
  "Brain_40879|Brain|8|BG_40879_arcATAC|BG_40879_arcGEX"
  # --- Colon, 12 libraries. N = necropsy, P = pre-infection ---
  "Colon_40872N|Colon|1|C40872_NA|C40872_NG"
  "Colon_40873N|Colon|2|C40873_NA|C40873_NG"
  "Colon_40874N|Colon|3|C40874_NA|C40874_NG"
  "Colon_40874P|Colon|4|C40874_PA|C40874_PG"
  "Colon_40875N|Colon|5|C40875_NA|C40875_NG"
  "Colon_40875P|Colon|6|C40875_PA|C40875_PG"
  "Colon_40876N|Colon|7|C40876_NA|C40876_NG"
  "Colon_40877N|Colon|8|C40877_NA|C40877_NG"
  "Colon_40877P|Colon|9|C40877_PA|C40877_PG"
  "Colon_40878N|Colon|10|C40878_NA|C40878_NG"
  "Colon_40878P|Colon|11|C40878_PA|C40878_PG"
  "Colon_40879N|Colon|12|C40879_NA|C40879_NG"
)

FIRST_RUN="Brain_40872"

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

say()  { echo "$*"; }
hdr()  { say ""; say "=============================================================================="; say "$*"; say "=============================================================================="; }
ok()   { say "  [ok]      $*"; }
bad()  { say "  [MISSING] $*"; }
warnx(){ say "  [warn]    $*"; }
die()  { say ""; say "FATAL: $*"; exit 1; }

run_field() {  # run_field <run_id> <1=tissue 2=index 3=atac 4=gex>
  local id="$1" n="$2" r
  for r in "${RUNS[@]}"; do
    if [[ "${r%%|*}" == "${id}" ]]; then
      IFS='|' read -r _ tis idx atac gex <<< "${r}"
      case "${n}" in
        1) echo "${tis}" ;;  2) echo "${idx}" ;;
        3) echo "${atac}" ;; 4) echo "${gex}" ;;
      esac
      return 0
    fi
  done
  return 1
}

atac_dir_for() { [[ "$1" == "Brain" ]] && echo "${BRAIN_ATAC_DIR}" || echo "${COLON_ATAC_DIR}"; }
gex_dir_for()  { [[ "$1" == "Brain" ]] && echo "${BRAIN_GEX_DIR}"  || echo "${COLON_GEX_DIR}";  }

write_libraries() {  # write_libraries <run_id> -> prints the csv path
  local id="$1"
  local tis atac gex ad gd out
  tis="$(run_field "${id}" 1)" || return 1
  atac="$(run_field "${id}" 3)"
  gex="$(run_field "${id}" 4)"
  ad="$(atac_dir_for "${tis}")/${atac}"
  gd="$(gex_dir_for  "${tis}")/${gex}"
  out="${LIB_DIR}/${id}_libraries.csv"
  mkdir -p "${LIB_DIR}"
  # Unix line endings deliberately: a \r\n here breaks downstream bash parsing.
  {
    printf 'fastqs,sample,library_type\n'
    printf '%s,%s,Chromatin Accessibility\n' "${ad}" "${atac}"
    printf '%s,%s,Gene Expression\n'          "${gd}" "${gex}"
  } > "${out}"
  echo "${out}"
}

run_state() {  # done | locked | partial | absent
  local d="${COUNT_ROOT}/$1"
  [[ -d "${d}" ]] || { echo absent; return; }
  [[ -e "${d}/_lock" ]] && { echo locked; return; }
  [[ -s "${d}/outs/summary.csv" && -s "${d}/outs/atac_fragments.tsv.gz" ]] \
    && { echo done; return; }
  echo partial
}

select_runs() {  # select_runs <first|brain|colon|all>
  local sel="$1" r
  for r in "${RUNS[@]}"; do
    local id="${r%%|*}" tis
    IFS='|' read -r _ tis _ _ _ <<< "${r}"
    case "${sel}" in
      first) [[ "${id}" == "${FIRST_RUN}" ]] && echo "${id}" ;;
      brain) [[ "${tis}" == "Brain" ]] && echo "${id}" ;;
      colon) [[ "${tis}" == "Colon" ]] && echo "${id}" ;;
      all)   echo "${id}" ;;
    esac
  done
}

# -----------------------------------------------------------------------------
# MODE: run   (the SLURM worker)
# -----------------------------------------------------------------------------

do_run() {
  local id="$1"
  local tis idx atac gex libs

  tis="$(run_field "${id}" 1)" || die "unknown run id '${id}'. Try: $0"
  idx="$(run_field "${id}" 2)"
  atac="$(run_field "${id}" 3)"
  gex="$(run_field "${id}" 4)"

  hdr "cellranger-arc count : ${id}"
  say "host:          $(hostname)"
  say "job:           ${SLURM_JOB_ID:-interactive}"
  say "tissue:        ${tis}"
  say "library index: ${idx}   (vendor merge() position, so cells become BARCODE-1_${idx})"
  say "ATAC sample:   ${atac}"
  say "GEX sample:    ${gex}"
  say "started:       $(date '+%Y-%m-%d %H:%M:%S')"

  # cellranger-arc is self-contained and sits outside conda. The explicit
  # CONDA_PREFIX/bin prepend is the established fix for the broken system
  # samtools at /usr/local/bin/samtools.
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/bin/activate" "${CONDA_ENV}" \
    || die "could not activate ${CONDA_ENV}"
  export PATH="${CONDA_PREFIX}/bin:${PATH}"
  export PATH="${ARC_DIR}:${PATH}"

  command -v cellranger-arc >/dev/null 2>&1 \
    || die "cellranger-arc not on PATH after adding ${ARC_DIR}"
  say ""
  cellranger-arc --version 2>&1 | head -2 | sed 's/^/  /'

  [[ -s "${REFERENCE}/reference.json" ]] \
    || die "reference not found at ${REFERENCE}"
  ok "reference: ${REFERENCE}"

  libs="$(write_libraries "${id}")" || die "could not write libraries.csv"
  say ""
  say "  libraries.csv:"
  sed 's/^/    /' "${libs}"

  # Verify every FASTQ directory the sheet names, and that it holds files for
  # the sample it claims.
  local n_bad=0
  while IFS=, read -r d s _; do
    [[ "${d}" == "fastqs" ]] && continue
    if [[ -d "${d}" ]]; then
      local n
      n="$(find "${d}" -maxdepth 1 -name "${s}_S*_L*_*.fastq.gz" | wc -l)"
      if [[ "${n}" -gt 0 ]]; then
        ok "${s}: ${n} FASTQ files"
      else
        bad "${s}: directory exists but holds no matching FASTQs"
        n_bad=$((n_bad+1))
      fi
    else
      bad "${s}: directory not found: ${d}"
      n_bad=$((n_bad+1))
    fi
  done < "${libs}"
  [[ ${n_bad} -eq 0 ]] || die "${n_bad} FASTQ path problem(s); nothing submitted"

  mkdir -p "${COUNT_ROOT}"
  cd "${COUNT_ROOT}" || die "cannot cd to ${COUNT_ROOT}"

  case "$(run_state "${id}")" in
    done)
      say ""
      ok "${id} already complete; nothing to do"
      say "  outs: ${COUNT_ROOT}/${id}/outs"
      return 0
      ;;
    locked)
      die "${COUNT_ROOT}/${id}/_lock exists. Another job may be running this library. If you are certain it is dead, remove the lock and resubmit."
      ;;
    partial)
      say ""
      warnx "${id} exists but is incomplete. count resumes from where it stopped, which is the designed behaviour, so letting it continue."
      ;;
  esac

  local cores="${SLURM_CPUS_PER_TASK:-${SLURM_CPUS}}"
  say ""
  say "  cellranger-arc count --id=${id} --reference=${REFERENCE} \\"
  say "      --libraries=${libs} --jobmode=local \\"
  say "      --localcores=${cores} --localmem=${ARC_LOCALMEM_GB}"
  say ""

  local t0 t1 rc
  t0=$(date +%s)
  cellranger-arc count \
    --id="${id}" \
    --reference="${REFERENCE}" \
    --libraries="${libs}" \
    --jobmode=local \
    --localcores="${cores}" \
    --localmem="${ARC_LOCALMEM_GB}"
  rc=$?
  t1=$(date +%s)

  say ""
  say "  finished: $(date '+%Y-%m-%d %H:%M:%S')"
  say "  elapsed:  $(( (t1-t0)/3600 ))h $(( ((t1-t0)%3600)/60 ))m"
  say "  exit:     ${rc}"

  if [[ ${rc} -ne 0 ]]; then
    say ""
    say "  FAILED. The stage that broke is named in:"
    say "    ${COUNT_ROOT}/${id}/_log"
    say "  Resubmitting this same job resumes rather than restarting."
    exit ${rc}
  fi

  local outs="${COUNT_ROOT}/${id}/outs"
  say ""
  say "--- outs ---"
  ls -la "${outs}" 2>/dev/null | sed 's/^/  /'

  say ""
  say "--- files the rest of the project needs ---"
  for f in atac_fragments.tsv.gz atac_fragments.tsv.gz.tbi atac_peaks.bed \
           atac_possorted_bam.bam gex_possorted_bam.bam \
           filtered_feature_bc_matrix.h5 filtered_feature_bc_matrix \
           per_barcode_metrics.csv summary.csv; do
    if [[ -e "${outs}/${f}" ]]; then
      printf '  [ok]      %-36s %s\n' "${f}" "$(du -sh "${outs}/${f}" | cut -f1)"
    else
      printf '  [MISSING] %s\n' "${f}"
    fi
  done
  say ""
  say "  atac_fragments.tsv.gz is the file the collaborators never received and"
  say "  that every coverage figure depends on."
  say ""
  say "  filtered_feature_bc_matrix/ holds the matrix, features and barcodes"
  say "  that the Seurat rebuild reads. Its barcodes end in -1, and merge()"
  say "  appends _${idx}, which reproduces the vendor's CellID exactly."

  say ""
  say "--- outs size ---"
  du -sh "${outs}" 2>/dev/null | sed 's/^/  /'

  say ""
  say "--- summary.csv ---"
  cat "${outs}/summary.csv" 2>/dev/null | sed 's/^/  /'

  return 0
}

# -----------------------------------------------------------------------------
# MODE: plan
# -----------------------------------------------------------------------------

do_plan() {
  hdr "ARC_02_count : plan"
  say "reference:   ${REFERENCE}"
  say "count root:  ${COUNT_ROOT}"
  say "libraries:   ${LIB_DIR}"
  say "resources:   ${SLURM_CPUS} cores, ${SLURM_MEM_GB} GB, ${SLURM_TIME}, --localmem=${ARC_LOCALMEM_GB}"

  say ""
  if [[ -s "${REFERENCE}/reference.json" ]]; then
    ok "reference present"
    [[ -e "${REFERENCE}/regions/motifs.pfm" ]] \
      && ok "regions/motifs.pfm present, so JASPAR was accepted at mkref"
  else
    bad "reference not found at ${REFERENCE}; run SETUP_04_mkref.sh"
  fi

  say ""
  say "  SETUP_03 reproduced 100.00% of the vendor's 43,128 TSS distances"
  say "  against this reference, so the assembly and annotation are confirmed."

  hdr "Run table"
  say "  Order is the vendor's Seurat merge() position, recovered by ARC_00 and"
  say "  equal to the alphabetical sort of the loom Samples labels."
  say ""
  printf '  %-14s %-7s %-4s %-18s %-18s %-9s %s\n' \
    "run_id" "tissue" "idx" "atac_sample" "gex_sample" "state" "input"
  say "  $(printf '%.0s-' {1..92})"

  local total=0 n_absent=0 n_done=0 n_partial=0 n_bad=0
  local r id tis idx atac gex ad gd state sz
  for r in "${RUNS[@]}"; do
    IFS='|' read -r id tis idx atac gex <<< "${r}"
    ad="$(atac_dir_for "${tis}")/${atac}"
    gd="$(gex_dir_for  "${tis}")/${gex}"
    state="$(run_state "${id}")"

    local bytes=0 missing=""
    if [[ -d "${ad}" ]]; then
      bytes=$(( bytes + $(du -sb "${ad}" 2>/dev/null | cut -f1) ))
    else
      missing="ATAC"
    fi
    if [[ -d "${gd}" ]]; then
      bytes=$(( bytes + $(du -sb "${gd}" 2>/dev/null | cut -f1) ))
    else
      missing="${missing:+${missing}+}GEX"
    fi
    total=$(( total + bytes ))
    sz="$(awk -v b="${bytes}" 'BEGIN{printf "%.0fG", b/1e9}')"

    printf '  %-14s %-7s %-4s %-18s %-18s %-9s %s%s\n' \
      "${id}" "${tis}" "${idx}" "${atac}" "${gex}" "${state}" "${sz}" \
      "${missing:+   <-- MISSING ${missing}}"

    [[ -n "${missing}" ]] && n_bad=$((n_bad+1))
    case "${state}" in
      done)    n_done=$((n_done+1)) ;;
      partial) n_partial=$((n_partial+1)) ;;
      absent)  n_absent=$((n_absent+1)) ;;
    esac
  done

  say ""
  say "  runs:    ${#RUNS[@]}  (done ${n_done}, partial ${n_partial}, not started ${n_absent})"
  say "  input:   $(awk -v b="${total}" 'BEGIN{printf "%.2f TB", b/1e12}')"
  if [[ ${n_bad} -gt 0 ]]; then
    say ""
    bad "${n_bad} run(s) have a missing FASTQ directory. Fix the paths at the top of this script before submitting."
  fi

  say ""
  say "  Output estimate: outs/ is dominated by the two position-sorted BAMs,"
  say "  which run comparable to the input FASTQ size, plus pipeline working"
  say "  space. So roughly 2 to 3 TB across all 20. Against 1.1 PB free on"
  say "  /master that is not a constraint, but run 1 gives a measured number."

  if [[ -d "${COUNT_ROOT}" ]]; then
    say ""
    say "  disk at the count root:"
    df -h "${COUNT_ROOT}" 2>/dev/null | sed 's/^/    /'
  fi

  hdr "Submit"
  say "  Start with one library. A single brain run validates the reference, the"
  say "  libraries.csv format, the resource request and the outs/ layout before"
  say "  19 more are committed, and replaces every estimate above with a"
  say "  measured number."
  say ""
  say "      $0 submit first      # ${FIRST_RUN} alone"
  say "      $0 submit brain      # the 8 brain libraries"
  say "      $0 submit colon      # the 12 colon libraries"
  say "      $0 submit all        # all ${#RUNS[@]}"
  say ""
  say "      $0 status            # progress"
  say "      $0 libraries         # write every libraries.csv without submitting"
  say ""
  say "  Each library becomes its own SLURM job, so any one can be resubmitted"
  say "  alone and a failure resumes rather than restarting."
}

# -----------------------------------------------------------------------------
# MODE: submit
# -----------------------------------------------------------------------------

do_submit() {
  local sel="${1:-first}"
  case "${sel}" in first|brain|colon|all) ;; *) die "submit takes first, brain, colon or all" ;; esac

  mapfile -t ids < <(select_runs "${sel}")
  [[ ${#ids[@]} -gt 0 ]] || die "no runs matched '${sel}'"

  local self
  self="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
  mkdir -p "${LOG_DIR}" "${COUNT_ROOT}" "${LIB_DIR}"

  hdr "Submitting ${#ids[@]} job(s) : ${sel}"
  local id state skip=()
  for id in "${ids[@]}"; do
    state="$(run_state "${id}")"
    printf '  %-14s %s\n' "${id}" "${state}"
    [[ "${state}" == "done" || "${state}" == "locked" ]] && skip+=("${id}")
  done

  if [[ ${#skip[@]} -gt 0 ]]; then
    say ""
    warnx "skipping ${#skip[@]} run(s) already complete or locked: ${skip[*]}"
  fi

  mapfile -t todo < <(for id in "${ids[@]}"; do
    s="$(run_state "${id}")"
    [[ "${s}" == "done" || "${s}" == "locked" ]] || echo "${id}"
  done)
  [[ ${#todo[@]} -gt 0 ]] || { say ""; ok "nothing to submit"; return 0; }

  say ""
  read -r -p "  submit ${#todo[@]} job(s)? [y/N] " ans
  [[ "${ans}" == "y" || "${ans}" == "Y" ]] || { say "  aborted"; return 0; }

  say ""
  for id in "${todo[@]}"; do
    local jid
    jid="$(sbatch --parsable \
            --job-name="arc_${id}" \
            --partition="${SLURM_PARTITION}" \
            --cpus-per-task="${SLURM_CPUS}" \
            --mem="${SLURM_MEM_GB}G" \
            --time="${SLURM_TIME}" \
            --output="${LOG_DIR}/arc_count_${id}_%j.out" \
            --error="${LOG_DIR}/arc_count_${id}_%j.err" \
            "${self}" run "${id}")"
    printf '  %-10s %s\n' "${jid}" "${id}"
  done

  say ""
  say "  squeue -u ${USER}"
  say "  logs: ${LOG_DIR}/arc_count_*"
  say ""
  say "  then: $0 status"
}

# -----------------------------------------------------------------------------
# MODE: status
# -----------------------------------------------------------------------------

do_status() {
  hdr "Status"
  printf '  %-14s %-9s %-10s %s\n' "run_id" "state" "outs" "cells (filtered)"
  say "  $(printf '%.0s-' {1..60})"

  local r id state outs sz cells n_done=0
  for r in "${RUNS[@]}"; do
    id="${r%%|*}"
    state="$(run_state "${id}")"
    outs="${COUNT_ROOT}/${id}/outs"
    sz="-"; cells="-"
    if [[ -d "${outs}" ]]; then
      sz="$(du -sh "${outs}" 2>/dev/null | cut -f1)"
      if [[ -s "${outs}/summary.csv" ]]; then
        cells="$(awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i ~ /[Ee]stimated number of cells/) c=i}
                          NR==2 && c{gsub(/"/,"",$c); print $c}' \
                 "${outs}/summary.csv" 2>/dev/null)"
        [[ -z "${cells}" ]] && cells="see summary.csv"
      fi
    fi
    printf '  %-14s %-9s %-10s %s\n' "${id}" "${state}" "${sz}" "${cells}"
    [[ "${state}" == "done" ]] && n_done=$((n_done+1))
  done

  say ""
  say "  complete: ${n_done} of ${#RUNS[@]}"

  if [[ ${n_done} -eq ${#RUNS[@]} ]]; then
    say ""
    ok "all runs complete"
    say ""
    say "  Next, and the decision to make first:"
    say ""
    say "  The vendor's CellID is BARCODE-1_<index>, which is Seurat merge(),"
    say "  not aggr's BARCODE-<n>. So the RNA objects can be built by reading"
    say "  each filtered_feature_bc_matrix and merging per tissue in the index"
    say "  order above, and the cell names come out identical to theirs with no"
    say "  remapping. The annotation transfer is then a plain join on CellID."
    say ""
    say "  The ATAC side still needs a consensus peak set and fragments per"
    say "  tissue. Two routes:"
    say "    a) cellranger-arc aggr per tissue. Gives one peak set and one"
    say "       merged fragments file, but its barcodes are BARCODE-<n> and"
    say "       need mapping onto BARCODE-1_<index>."
    say "    b) Union the per-library atac_peaks.bed in Signac, then"
    say "       FeatureMatrix against per-library Fragment objects. More R, no"
    say "       barcode gymnastics, and keeps per-library fragments."
    say ""
    say "  Either reproduces a per-tissue peak universe, which is what the"
    say "  delivered 204,136 brain and 345,553 colon figures imply."
  else
    say ""
    say "  squeue -u ${USER}"
    say "  A partial run resumes on resubmission; it does not restart."
  fi
}

# -----------------------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------------------

case "${1:-plan}" in
  plan|"")   do_plan ;;
  run)       [[ -n "${2:-}" ]] || die "run needs a run_id"; do_run "$2" ;;
  submit)    do_submit "${2:-first}" ;;
  status)    do_status ;;
  libraries)
    hdr "Writing libraries.csv for all ${#RUNS[@]} runs"
    for r in "${RUNS[@]}"; do
      p="$(write_libraries "${r%%|*}")" && ok "$(basename "${p}")"
    done
    say ""
    say "  in ${LIB_DIR}"
    ;;
  -h|--help|help)
    sed -n '/^# ARC_02_count.sh/,/^# ----/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    ;;
  *) die "unknown mode '${1}'. Try: $0 --help" ;;
esac

exit 0
