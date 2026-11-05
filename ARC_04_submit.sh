#!/usr/bin/env bash
# =============================================================================
# ARC_04_submit.sh
#
# Submits ARC_04_build_object.R, one SLURM job per tissue. Same self-dispatching
# pattern as ARC_02_count.sh: `submit` re-submits this file with `run <tissue>`,
# so the plan and the worker cannot drift apart.
#
#   ./ARC_04_submit.sh                 verify every input and print the plan
#   ./ARC_04_submit.sh check           R package check only
#   ./ARC_04_submit.sh submit brain    submit the brain object
#   ./ARC_04_submit.sh submit colon    submit the colon object
#   ./ARC_04_submit.sh submit all      submit both
#   ./ARC_04_submit.sh status          what exists, and how far FeatureMatrix got
#   ./ARC_04_submit.sh run brain       run in the foreground (what SLURM calls)
#
# -----------------------------------------------------------------------------
# WHY THE TWO TISSUES GET DIFFERENT RESOURCES
# -----------------------------------------------------------------------------
#
# The expensive step is FeatureMatrix: counting every consensus peak against
# every fragment file. ARC_03 measured both sides of that product.
#
#   Brain   8 fragment files    ~203,000 consensus peaks    ~57,000 cells
#   Colon  12 fragment files    ~340,000 consensus peaks   ~34,000 cells
#
# So colon is 1.5x the files against 1.7x the peaks, and it gets the larger
# allocation and the longer wall clock. Brain carries more cells, which costs
# memory in the GEX cbind but little in the peak counting.
#
# Neither job threads. Signac's FeatureMatrix is single-process here, so the
# cpus-per-task request is about the memory that comes attached to cores on
# this partition, not about parallelism. If the scheduler is busy, dropping
# cpus-per-task costs nothing as long as --mem holds.
#
# -----------------------------------------------------------------------------
# WHY A TIMEOUT IS NOT A DISASTER
# -----------------------------------------------------------------------------
#
# ARC_04 checkpoints each library's FeatureMatrix to
# 04_objects/<tissue>/featurematrix_checkpoints/ and validates the checkpoint
# against the consensus peak set before reusing it. So a job that dies partway
# resumes from the next library. Resubmitting the same tissue is always safe.
#
# Two knobs live in ARC_04_build_object.R rather than here, because they change
# what the object is rather than how it is scheduled:
#
#   VALIDATE_FRAGMENTS   TRUE reads through every fragments file once to confirm
#                        the requested barcodes are present. Worth it on the
#                        first build. Set FALSE on a rerun to save 20-40 min.
#   FM_CHUNK             peaks per FeatureMatrix chunk. Lower it if a job is
#                        killed for memory rather than time.
#
# Jake Lehle, 2026-10
# =============================================================================

set -uo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

ANALYSIS_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
COUNT_ROOT="${ANALYSIS_ROOT}/02_arc/02_count"
LOOM_META_DIR="${ANALYSIS_ROOT}/02_arc/00_loom_metadata"
OUT_ROOT="${ANALYSIS_ROOT}/02_arc/04_objects"

GTF="/master/jlehle/WORKING/REFS/Mmul8/genes.gtf"

CONDA_BASE="/master/jlehle/anaconda3"
# multiome_R is the env that actually carries Seurat, Signac, GenomicRanges and
# rtracklayer. sc_pre_R has r-base and little else, so pointing at it made the
# package check report everything missing. Override with ARC_CONDA_ENV=... to
# test another env without editing this file.
CONDA_ENV="${ARC_CONDA_ENV:-multiome_R}"
LOG_DIR="/master/jlehle/WORKING/LOGS"

SLURM_PARTITION="normal"

# Per-tissue allocation: cpus | mem GB | wall clock
declare -A TISSUE_CPUS=( [Brain]=16 [Colon]=16 )
declare -A TISSUE_MEM=(  [Brain]=160 [Colon]=256 )
declare -A TISSUE_TIME=( [Brain]="24:00:00" [Colon]="48:00:00" )

# Library order per tissue. Must match LIB_ORDER in ARC_04_build_object.R; the
# plan mode checks that the two agree rather than trusting this copy.
BRAIN_LIBS=(40872 40873 40874 40875 40876 40877 40878 40879)
COLON_LIBS=(40872N 40873N 40874N 40874P 40875N 40875P \
            40876N 40877N 40877P 40878N 40878P 40879N)

# R packages the build needs. Checked before anything is submitted, because a
# missing package should fail in two seconds at the terminal, not six hours
# into a colon job.
R_PACKAGES=(Matrix Seurat Signac GenomicRanges rtracklayer data.table)

# Not required. ARC_04 reads filtered_feature_bc_matrix/ through Read10X, which
# is pure R, and only falls back to the .h5 if a MEX directory is incomplete.
# hdf5r is what that fallback would need, so its absence is reported but never
# fatal.
R_PACKAGES_OPTIONAL=(hdf5r)

# -----------------------------------------------------------------------------
# Plumbing
# -----------------------------------------------------------------------------

self="$(readlink -f "${BASH_SOURCE[0]}")"
self_dir="$(dirname "${self}")"

# WHY THIS IS NOT JUST "${self_dir}/ARC_04_build_object.R"
#
# sbatch copies the batch script into the compute node's spool directory and
# runs it from there, so inside a job BASH_SOURCE is
# /var/spool/slurmd/job<id>/<name> and the R script is not beside it. Resolving
# the path from self_dir therefore works interactively and fails under SLURM.
#
# So: submit resolves the path while it is still in the real directory and
# passes it to the job as the third argument to `run`. The resolution order
# below tries that argument, then the env override, then self_dir when we are
# demonstrably not in the spool, then the fallback constant.
IN_SPOOL=0
case "${self_dir}" in /var/spool/*|/tmp/slurm*) IN_SPOOL=1 ;; esac

# Last resort, for a job submitted by an older copy of this script that does not
# pass the path through. Edit if you move the scripts.
R_SCRIPT_FALLBACK="/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome/ARC_04_build_object.R"

resolve_r_script() {
  local from_arg="${1:-}" c
  local -a cands=()
  [[ -n "${ARC_04_R:-}" ]] && cands+=("${ARC_04_R}")
  [[ -n "${from_arg}" ]]   && cands+=("${from_arg}")
  [[ ${IN_SPOOL} -eq 0 ]]  && cands+=("${self_dir}/ARC_04_build_object.R")
  cands+=("${R_SCRIPT_FALLBACK}")
  for c in "${cands[@]}"; do
    [[ -s "${c}" ]] && { readlink -f "${c}"; return 0; }
  done
  return 1
}

R_SCRIPT="$(resolve_r_script || true)"

# $0 inside a SLURM job is the spool copy's path, so printing it in a "run this
# next" hint hands you a command that does not exist. Use a plain name there.
if [[ ${IN_SPOOL} -eq 1 ]]; then INVOKE="./ARC_04_submit.sh"; else INVOKE="$0"; fi

say()   { printf '%s\n' "$*"; }
hdr()   { say ""; say "============================================================================"; say "$*"; say "============================================================================"; }
subh()  { say ""; say "--- $* ---"; }
ok()    { say "  [ok]    $*"; }
warnx() { say "  [warn]  $*"; }
bad()   { say "  [FAIL]  $*"; }
die()   { say ""; say "  [FAIL]  $*"; say ""; exit 1; }

tissue_libs() {
  case "$1" in
    Brain) printf '%s\n' "${BRAIN_LIBS[@]}" ;;
    Colon) printf '%s\n' "${COLON_LIBS[@]}" ;;
    *) return 1 ;;
  esac
}

# brain -> Brain, COLON -> Colon
canon_tissue() {
  local t="${1,,}"
  case "${t}" in
    brain) echo Brain ;;
    colon) echo Colon ;;
    *) return 1 ;;
  esac
}

# Must match the filename ARC_04_build_object.R writes in phase 8.
rds_path() { echo "${OUT_ROOT}/$1/seurat_$1_rebuilt.rds"; }
ckpt_dir() { echo "${OUT_ROOT}/$1/featurematrix_checkpoints"; }

# -----------------------------------------------------------------------------
# Input verification
#
# Every check here is a precondition of ARC_04, in the order the R script hits
# them. The point is that `plan` either prints an all-clear or names the one
# thing to fix.
# -----------------------------------------------------------------------------

verify_tissue() {
  local tis="$1" fails=0 libs frag tbi mtx n_h5=0
  mapfile -t libs < <(tissue_libs "${tis}")

  subh "${tis}  (${#libs[@]} libraries)"

  # The loom metadata is the cell list and the annotation source, so it is the
  # first thing the R script reads and the first thing worth checking.
  # Counted with a CSV parser, not wc -l. The vendor's cluster labels contain
  # literal newlines ("Entero\ncyte-3"), so those rows are quoted and span two
  # physical lines: wc -l reported 29,957 for colon where the real row count is
  # 27,296. A count that disagrees with what the R script reads is worse than no
  # count, because it looks like cells went missing.
  local meta="${LOOM_META_DIR}/cell_metadata_${tis}.csv"
  if [[ -s "${meta}" ]]; then
    local n_cells
    if command -v python3 >/dev/null 2>&1; then
      n_cells="$(python3 -I -c '
import csv, sys
with open(sys.argv[1], newline="") as fh:
    print(max(sum(1 for _ in csv.reader(fh)) - 1, 0))
' "${meta}" 2>/dev/null)"
    fi
    if [[ -n "${n_cells:-}" ]]; then
      ok "cell_metadata_${tis}.csv  ${n_cells} cells"
    else
      n_cells="$(( $(wc -l < "${meta}") - 1 ))"
      ok "cell_metadata_${tis}.csv  ~${n_cells} lines (no python3 for an exact row count)"
    fi
  else
    bad "missing ${meta}  (run ARC_00_loom_metadata.py)"
    fails=$((fails + 1))
  fi

  for id in "${libs[@]}"; do
    # The count directory is <Tissue>_<sample>, the run_id ARC_02 created and
    # the one ARC_04_build_object.R rebuilds as paste0(TISSUE, "_", libs).
    # NOT the bare sample id.
    local rid="${tis}_${id}"
    local outs="${COUNT_ROOT}/${rid}/outs"
    frag="${outs}/atac_fragments.tsv.gz"
    tbi="${frag}.tbi"
    # Either input satisfies ARC_04: it prefers the MEX directory through
    # Read10X and falls back to the .h5 only if that directory is incomplete.
    local mex="${outs}/filtered_feature_bc_matrix"
    mtx="${outs}/filtered_feature_bc_matrix.h5"
    local peaks="${outs}/atac_peaks.bed"

    # Distinguish "the run never happened" from "the run is incomplete", since
    # the fix is different: one is a resubmit of ARC_02, the other a look at
    # why the pipestance stopped.
    if [[ ! -d "${COUNT_ROOT}/${rid}" ]]; then
      bad "${rid}: no count directory at ${COUNT_ROOT}/${rid}"
      fails=$((fails + 1))
      continue
    fi
    if [[ ! -d "${outs}" ]]; then
      bad "${rid}: count directory exists but has no outs/ (run did not finish)"
      fails=$((fails + 1))
      continue
    fi

    local mex_ok=1
    for f in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
      [[ -s "${mex}/${f}" ]] || mex_ok=0
    done

    local miss=()
    [[ ${mex_ok} -eq 1 || -s "${mtx}" ]] || miss+=("filtered matrix (MEX dir and h5 both absent)")
    [[ ${mex_ok} -eq 1 ]] || n_h5=$((n_h5 + 1))
    [[ -s "${frag}" ]]  || miss+=("fragments")
    [[ -s "${tbi}" ]]   || miss+=("fragments index")
    [[ -s "${peaks}" ]] || miss+=("atac_peaks.bed")

    if [[ ${#miss[@]} -ne 0 ]]; then
      bad "${rid}: missing ${miss[*]}"
      fails=$((fails + 1))
    fi
  done
  [[ ${fails} -eq 0 ]] && ok "all ${#libs[@]} count outputs complete (matrix, fragments, .tbi, peaks)"
  if [[ ${n_h5} -gt 0 ]]; then
    warnx "${n_h5} library/libraries have no complete MEX directory, so those will"
    warnx "be read from the .h5 instead, which needs hdf5r installed"
  fi

  # State of a previous build, so `plan` doubles as "where did I leave off".
  local rds ck
  rds="$(rds_path "${tis}")"
  ck="$(ckpt_dir "${tis}")"
  if [[ -s "${rds}" ]]; then
    warnx "object already exists: ${rds}  ($(du -h "${rds}" | cut -f1)) - resubmitting overwrites it"
  fi
  if [[ -d "${ck}" ]]; then
    local n_ck
    n_ck="$(find "${ck}" -name '*.rds' -type f 2>/dev/null | wc -l)"
    [[ "${n_ck}" -gt 0 ]] && ok "${n_ck}/${#libs[@]} FeatureMatrix checkpoints present, will be reused"
  fi

  return "${fails}"
}

# The one consistency check that actually matters: the library order in this
# file against LIB_ORDER in the R script. They encode the vendor's merge()
# position, so a mismatch silently renames every cell.
verify_lib_order() {
  subh "library order agreement"
  local r_brain r_colon mine_brain mine_colon
  # The Brain entry runs from its own line down to the line before Colon's.
  # A /Brain/,/)/ range will not do: sed never closes a range on the line that
  # opened it, so it would swallow the Colon entry too.
  r_brain="$(sed -n '/^LIB_ORDER/,/^)/p' "${R_SCRIPT}" \
              | sed -n '/Brain *=/,/Colon *=/p' | sed '$d' \
              | grep -o '"[0-9A-Z]*"' | tr -d '"' | tr '\n' ' ')"
  r_colon="$(sed -n '/^LIB_ORDER/,/^)/p' "${R_SCRIPT}" \
              | sed -n '/Colon *=/,$p' | grep -o '"[0-9A-Z]*"' | tr -d '"' | tr '\n' ' ')"
  mine_brain="$(printf '%s ' "${BRAIN_LIBS[@]}")"
  mine_colon="$(printf '%s ' "${COLON_LIBS[@]}")"

  if [[ "${r_brain}" == "${mine_brain}" ]]; then
    ok "Brain order matches ARC_04_build_object.R"
  else
    bad "Brain order differs"
    say "         here: ${mine_brain}"
    say "         R:    ${r_brain}"
    return 1
  fi
  if [[ "${r_colon}" == "${mine_colon}" ]]; then
    ok "Colon order matches ARC_04_build_object.R"
  else
    bad "Colon order differs"
    say "         here: ${mine_colon}"
    say "         R:    ${r_colon}"
    return 1
  fi
  return 0
}

activate_env() {
  # The SLURM job activates CONDA_ENV explicitly rather than inheriting whatever
  # was active at submit time, so the build is reproducible from any shell. If
  # you are sitting in a different env, that is worth saying out loud, because
  # the packages the check reports are the job's, not your shell's.
  if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "${CONDA_ENV}" ]]; then
    warnx "your shell has '${CONDA_DEFAULT_ENV}' active; this script uses '${CONDA_ENV}'"
    warnx "override with: ARC_CONDA_ENV=${CONDA_DEFAULT_ENV} ${INVOKE} ..."
  fi
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/bin/activate" "${CONDA_ENV}" \
    || die "could not activate ${CONDA_ENV}"
  export PATH="${CONDA_PREFIX}/bin:${PATH}"
  command -v Rscript >/dev/null 2>&1 \
    || die "Rscript not on PATH after activating ${CONDA_ENV}"
}

# -----------------------------------------------------------------------------
# MODE: check
# -----------------------------------------------------------------------------

mode_check() {
  hdr "R environment"
  activate_env
  say "  env:      ${CONDA_PREFIX}"
  say "  Rscript:  $(command -v Rscript)"
  say "  R:        $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
  say ""

  local pkgs opt
  pkgs="$(printf '"%s",' "${R_PACKAGES[@]}")";          pkgs="${pkgs%,}"
  opt="$(printf '"%s",' "${R_PACKAGES_OPTIONAL[@]}")";  opt="${opt%,}"

  # R reports each package and emits a parseable MISSING: line. The conda
  # package names are composed in bash below, because the r- / bioconductor-
  # prefix is not derivable from the R name. Optional packages are reported but
  # never make it into MISSING:, so they cannot fail the check.
  local out rc
  out="$(Rscript -e "
pkgs <- c(${pkgs})
optional <- c(${opt})
missing <- character(0)
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) { missing <- c(missing, p); cat(sprintf('  [FAIL]  %-16s not installed\n', p)) }
  else cat(sprintf('  [ok]    %-16s %s\n', p, v))
}
for (p in optional) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) cat(sprintf('  [--]    %-16s not installed (optional, .h5 fallback only)\n', p))
  else cat(sprintf('  [ok]    %-16s %s   (optional)\n', p, v))
}
if (length(missing)) cat(paste0('MISSING:', paste(missing, collapse = ','), '\n'))
" 2>&1)"
  rc=$?

  printf '%s\n' "${out}" | grep -v '^MISSING:'

  local missing_line
  missing_line="$(printf '%s\n' "${out}" | sed -n 's/^MISSING://p')"

  if [[ ${rc} -ne 0 ]]; then
    say ""
    bad "Rscript failed while checking packages"
    return 1
  fi

  if [[ -z "${missing_line}" ]]; then
    say ""
    ok "all ${#R_PACKAGES[@]} packages present"
    return 0
  fi

  # R package name -> conda package name. CRAN takes r-, Bioconductor takes
  # bioconductor-, and nothing about the R name tells you which.
  declare -A CONDA_PKG=(
    [Matrix]=r-matrix
    [Seurat]=r-seurat
    [Signac]=r-signac
    [data.table]=r-data.table
    [hdf5r]=r-hdf5r
    [GenomicRanges]=bioconductor-genomicranges
    [rtracklayer]=bioconductor-rtracklayer
  )

  local -a want=()
  local IFS=','
  for p in ${missing_line}; do
    want+=("${CONDA_PKG[$p]:-r-${p,,}}")
  done
  unset IFS

  say ""
  say "  install into ${CONDA_ENV}:"
  say ""
  say "    conda install -n ${CONDA_ENV} --solver=libmamba \\"
  say "      -c conda-forge -c bioconda --strict-channel-priority \\"
  say "      ${want[*]}"
  say ""
  say "  Add --dry-run first if you want to see the solve before it commits."
  return 1
}

# -----------------------------------------------------------------------------
# MODE: plan  (default)
# -----------------------------------------------------------------------------

mode_plan() {
  hdr "ARC_04 build plan"
  say "R script:   ${R_SCRIPT:-<not found>}"
  say "count root: ${COUNT_ROOT}"
  say "loom meta:  ${LOOM_META_DIR}"
  say "output:     ${OUT_ROOT}"
  say "GTF:        ${GTF}"

  local fails=0

  subh "script and reference"
  if [[ -n "${R_SCRIPT}" && -s "${R_SCRIPT}" ]]; then ok "ARC_04_build_object.R found"
  else
    bad "ARC_04_build_object.R not found beside this script or at ${R_SCRIPT_FALLBACK}"
    say "          set ARC_04_R=/abs/path/ARC_04_build_object.R"
    fails=$((fails+1))
  fi
  if [[ -s "${GTF}" ]]; then ok "GTF: $(du -h "${GTF}" | cut -f1)"
  else bad "GTF not found at ${GTF}"; fails=$((fails+1)); fi
  if [[ -d "${LOG_DIR}" ]]; then ok "log dir: ${LOG_DIR}"
  else warnx "log dir ${LOG_DIR} does not exist, creating"; mkdir -p "${LOG_DIR}"; fi

  [[ -s "${R_SCRIPT}" ]] && { verify_lib_order || fails=$((fails+1)); }

  verify_tissue Brain || fails=$((fails + $?))
  verify_tissue Colon || fails=$((fails + $?))

  subh "SLURM allocation"
  printf '  %-8s %6s %8s %12s\n' tissue cpus mem time
  for tis in Brain Colon; do
    printf '  %-8s %6s %8s %12s\n' \
      "${tis}" "${TISSUE_CPUS[$tis]}" "${TISSUE_MEM[$tis]}G" "${TISSUE_TIME[$tis]}"
  done
  say ""
  say "  Colon gets more because FeatureMatrix runs 12 fragment files against"
  say "  ~340K peaks against brain's 8 files and ~203K peaks. Both checkpoint"
  say "  per library, so a timeout resumes."

  # Subshell, because activate_env calls die on failure and die exits. In a
  # subshell that only ends the check, so `plan` can report it and carry on.
  # The output is captured rather than discarded: a package failure should name
  # the package here, not send you off to run another command.
  subh "R packages"
  local chk_out chk_rc
  chk_out="$( ( mode_check ) 2>&1 )"; chk_rc=$?
  if [[ ${chk_rc} -ne 0 ]]; then
    printf '%s\n' "${chk_out}" | sed -n '/^  \[/,$p'
    bad "package check failed"
    fails=$((fails + 1))
  else
    ok "all ${#R_PACKAGES[@]} packages present in ${CONDA_ENV}"
  fi

  hdr "Verdict"
  if [[ ${fails} -eq 0 ]]; then
    say "  READY"
    say ""
    say "    ${INVOKE} submit brain      # 8 libraries, the faster one, run this first"
    say "    ${INVOKE} submit colon"
    say "    ${INVOKE} submit all"
    say ""
    say "  Brain first is worth it: if the barcode match gate trips, it trips"
    say "  there in under an hour rather than overnight in colon."
  else
    say "  NOT READY: ${fails} problem(s) above"
  fi
  say ""
  [[ ${fails} -eq 0 ]]
}

# -----------------------------------------------------------------------------
# MODE: submit
# -----------------------------------------------------------------------------

mode_submit() {
  local what="${1:-}"
  local -a tissues=()
  case "${what,,}" in
    brain) tissues=(Brain) ;;
    colon) tissues=(Colon) ;;
    all|both|"") tissues=(Brain Colon) ;;
    *) die "submit takes brain, colon, or all" ;;
  esac

  hdr "Submitting ${#tissues[@]} job(s)"

  local fails=0
  for tis in "${tissues[@]}"; do
    verify_tissue "${tis}" || fails=$((fails + $?))
  done
  [[ ${fails} -eq 0 ]] || die "${fails} input problem(s), see above. Run: ${INVOKE}"

  # Resolved here, in the real directory, because the job will run from the
  # spool copy where the R script is not reachable by relative means.
  [[ -n "${R_SCRIPT}" && -s "${R_SCRIPT}" ]] \
    || die "ARC_04_build_object.R not found. Put it beside this script, or pass ARC_04_R=/abs/path"
  ok "R script: ${R_SCRIPT}"

  mkdir -p "${LOG_DIR}"

  say ""
  for tis in "${tissues[@]}"; do
    printf '  %-8s %s cpus  %sG  %s\n' \
      "${tis}" "${TISSUE_CPUS[$tis]}" "${TISSUE_MEM[$tis]}" "${TISSUE_TIME[$tis]}"
  done
  say ""
  read -r -p "  submit? [y/N] " ans
  [[ "${ans}" == "y" || "${ans}" == "Y" ]] || { say "  aborted"; return 0; }

  say ""
  for tis in "${tissues[@]}"; do
    local jid
    jid="$(sbatch --parsable \
            --job-name="arc04_${tis,,}" \
            --partition="${SLURM_PARTITION}" \
            --cpus-per-task="${TISSUE_CPUS[$tis]}" \
            --mem="${TISSUE_MEM[$tis]}G" \
            --time="${TISSUE_TIME[$tis]}" \
            --output="${LOG_DIR}/arc04_${tis,,}_%j.out" \
            --error="${LOG_DIR}/arc04_${tis,,}_%j.err" \
            "${self}" run "${tis}" "${R_SCRIPT}")" \
      || { bad "sbatch failed for ${tis}"; continue; }
    printf '  %-10s %s\n' "${jid}" "${tis}"
  done

  say ""
  say "  squeue -u ${USER}"
  say "  logs: ${LOG_DIR}/arc04_*"
  say ""
  say "  then: ${INVOKE} status"
  say ""
}

# -----------------------------------------------------------------------------
# MODE: run   (what SLURM executes)
# -----------------------------------------------------------------------------

mode_run() {
  local tis
  tis="$(canon_tissue "${1:-}")" || die "run takes Brain or Colon"

  # The second argument is the absolute R script path that submit resolved. We
  # are running from the spool copy, so this is the only reliable source.
  R_SCRIPT="$(resolve_r_script "${2:-}")" || R_SCRIPT=""
  if [[ -z "${R_SCRIPT}" ]]; then
    hdr "ARC_04 run : ${tis}"
    say "  Tried, in order:"
    [[ -n "${ARC_04_R:-}" ]] && say "    ARC_04_R=${ARC_04_R}"
    [[ -n "${2:-}" ]]        && say "    argument: $2"
    [[ ${IN_SPOOL} -eq 0 ]]  && say "    ${self_dir}/ARC_04_build_object.R"
    say "    ${R_SCRIPT_FALLBACK}"
    say ""
    say "  This job is running from ${self_dir}, which is SLURM's spool copy of"
    say "  the batch script, so the R script cannot be found relative to it."
    die "ARC_04_build_object.R not found; resubmit with ARC_04_R=/abs/path exported"
  fi

  hdr "ARC_04 run : ${tis}"
  say "host:     $(hostname)"
  say "job:      ${SLURM_JOB_ID:-interactive}"
  say "cpus:     ${SLURM_CPUS_PER_TASK:-?}"
  say "mem:      ${SLURM_MEM_PER_NODE:-?} MB"
  say "R script: ${R_SCRIPT}"
  say "started:  $(date '+%Y-%m-%d %H:%M:%S')"

  activate_env
  say ""
  ok "env: ${CONDA_PREFIX}"
  ok "R:   $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"

  # Signac and data.table both try to grab every core they can see. On a shared
  # partition that oversubscribes the allocation, so pin them to what SLURM gave.
  local nthreads="${SLURM_CPUS_PER_TASK:-4}"
  export OMP_NUM_THREADS="${nthreads}"
  export OPENBLAS_NUM_THREADS="${nthreads}"
  export MKL_NUM_THREADS="${nthreads}"
  export R_DATATABLE_NUM_THREADS="${nthreads}"
  ok "threads pinned to ${nthreads}"

  # Already resolved and confirmed at the top of this function.

  say ""
  say "============================================================================"
  say ""
  Rscript --vanilla "${R_SCRIPT}" "${tis}"
  local rc=$?
  say ""
  say "============================================================================"

  local rds
  rds="$(rds_path "${tis}")"
  if [[ ${rc} -eq 0 && -s "${rds}" ]]; then
    hdr "Done : ${tis}"
    ok "object:   ${rds}  ($(du -h "${rds}" | cut -f1))"
    [[ -s "${OUT_ROOT}/${tis}/MANIFEST_${tis}.txt" ]] && {
      say ""
      sed 's/^/  /' "${OUT_ROOT}/${tis}/MANIFEST_${tis}.txt"
    }
    say ""
    say "  finished: $(date '+%Y-%m-%d %H:%M:%S')"
  else
    hdr "FAILED : ${tis}  (exit ${rc})"
    say "  Checkpoints are kept, so resubmitting resumes rather than restarting:"
    say "    ${INVOKE} submit ${tis,,}"
    say ""
    say "  If the job was killed for memory, lower FM_CHUNK in"
    say "  ARC_04_build_object.R before resubmitting."
  fi
  say ""
  return "${rc}"
}

# -----------------------------------------------------------------------------
# MODE: status
# -----------------------------------------------------------------------------

mode_status() {
  hdr "ARC_04 status"
  printf '  %-8s %-10s %-12s %-14s %s\n' tissue state size checkpoints cells
  for tis in Brain Colon; do
    local rds ck n_ck n_libs state size cells
    rds="$(rds_path "${tis}")"
    ck="$(ckpt_dir "${tis}")"
    mapfile -t _libs < <(tissue_libs "${tis}")
    n_libs="${#_libs[@]}"
    n_ck=0
    [[ -d "${ck}" ]] && n_ck="$(find "${ck}" -name '*.rds' -type f 2>/dev/null | wc -l)"

    if [[ -s "${rds}" ]]; then
      state="done"
      size="$(du -h "${rds}" | cut -f1)"
    elif [[ "${n_ck}" -gt 0 ]]; then
      state="partial"; size="-"
    else
      state="not built"; size="-"
    fi

    cells="-"
    local cm="${OUT_ROOT}/${tis}/cell_metadata_rebuilt_${tis}.csv"
    [[ -s "${cm}" ]] && cells="$(( $(wc -l < "${cm}") - 1 ))"

    printf '  %-8s %-10s %-12s %-14s %s\n' \
      "${tis}" "${state}" "${size}" "${n_ck}/${n_libs}" "${cells}"
  done

  say ""
  if command -v squeue >/dev/null 2>&1; then
    subh "queue"
    squeue -u "${USER}" -o '%.12i %.14j %.8T %.10M %.6C %.8m %R' 2>/dev/null \
      | grep -E 'JOBID|arc04' || say "  no arc04 jobs in the queue"
  fi

  subh "recent logs"
  ls -1t "${LOG_DIR}"/arc04_*.out 2>/dev/null | head -4 | sed 's/^/  /' \
    || say "  none yet"
  say ""
}

# -----------------------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------------------

case "${1:-plan}" in
  plan|"")  mode_plan ;;
  check)    mode_check ;;
  submit)   mode_submit "${2:-all}" ;;
  # The third argument is the absolute R script path that submit resolved. It
  # must be forwarded, because the job runs from SLURM's spool copy.
  run)      mode_run "${2:-}" "${3:-}" ;;
  status)   mode_status ;;
  -h|--help|help)
    sed -n '3,20p' "${self}" | sed 's/^# \{0,1\}//'
    ;;
  *) die "unknown mode: $1   (plan | check | submit | run | status)" ;;
esac
