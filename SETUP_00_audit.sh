#!/usr/bin/env bash
# =============================================================================
# SETUP_00_audit.sh
#
# Read-only audit of what the THC/CBD multiome project needs and what is already
# present on titan. Changes NOTHING: creates no environments, installs no
# packages, downloads nothing, writes only into OUT_DIR.
#
# Run this, read the report, then run SETUP_01_build_env.sh.
#
# WHAT IT CHECKS
#   1. conda / mamba, and which environments exist
#   2. whether sc_pre looks healthy enough to clone
#   3. every command-line binary the two tracks need, with versions
#   4. every R package needed, with versions
#   5. every Python package needed, with versions
#   6. cellranger-arc presence and version
#   7. reference files, and whether they are complete and indexed
#   8. the tidy tables from MULTIOME_01, which the reference validation needs
#   9. disk headroom against a real estimate of what the 20 ARC runs will write
#  10. the conda solve for the new environment, with --dry-run
#
# Exit code is 0 even when things are missing. The point is the report, not a
# pass/fail gate.
#
# Jake Lehle, 2026-10
# =============================================================================

set -uo pipefail

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

PROJECT_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome"
ANALYSIS_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
OUT_DIR="${ANALYSIS_ROOT}/setup_audit"

CONDA_BASE="${HOME}/anaconda3"
SOURCE_ENV="sc_pre"
TARGET_ENV="sc_pre_R"

# Where cellranger-arc is or will be unpacked. Checked, not created.
ARC_EXPECTED_DIR="/master/jlehle/WORKING/SOFTWARE/cellranger-arc-2.0.2"
ARC_EXPECTED_VERSION="2.0.2"

# Where the Mmul8 reference will live.
REF_DIR="/master/jlehle/WORKING/REFS/Mmul8"

# Filesystem to check for headroom, and the estimate to check it against.
# 1.7 TB of input across 20 library pairs. ARC output is dominated by the two
# position-sorted BAMs, which run comparable to input FASTQ size, plus pipeline
# working space. 5 TB is a deliberately conservative floor.
DISK_CHECK_PATH="/master/jlehle/WORKING"
DISK_NEEDED_TB=5

# Try the solve for the new environment. Slow (minutes) but catches a failure
# before it costs an afternoon. Set to 0 to skip.
TRY_SOLVE=1

# -----------------------------------------------------------------------------
# Requirements
# -----------------------------------------------------------------------------

# Command-line binaries: name|conda package|why
BINARIES=(
  "R|r-base|R interpreter for the Seurat/Signac/muscat stack"
  "Rscript|r-base|non-interactive R"
  "samtools|samtools|faidx, BAM indexing. NOTE: system binary is broken, must come from conda"
  "bedtools|bedtools|getfasta, nuc, shuffle for motif backgrounds"
  "sea|meme|MEME known-motif enrichment (Track B primary)"
  "xstreme|meme|MEME de novo discovery (optional)"
  "tomtom|meme|map de novo motifs back to JASPAR"
  "fimo|meme|locate individual motif sites"
  "meme-chip|meme|alternative combined pipeline"
  "bgzip|htslib|compress fragments files"
  "tabix|htslib|index fragments files"
  "wget|wget|reference downloads"
  "curl|curl|reference downloads"
  "md5sum|coreutils|checksum verification"
)

# R packages: name|conda package|why
R_PACKAGES=(
  "Seurat|r-seurat|object reconstruction"
  "Signac|r-signac|ChromatinAssay, CoveragePlot, reference validation"
  "SeuratObject|r-seuratobject|Seurat dependency"
  "muscat|bioconductor-muscat|per-animal pseudobulk DA, matches the vendor RNA method"
  "edgeR|bioconductor-edger|muscat backend"
  "SingleCellExperiment|bioconductor-singlecellexperiment|muscat input class"
  "GenomicRanges|bioconductor-genomicranges|interval arithmetic"
  "rtracklayer|bioconductor-rtracklayer|GTF and bigWig IO"
  "Matrix|r-matrix|sparse matrices"
  "data.table|r-data.table|fast IO"
)

# Python packages: import name|conda or pip package|why
PY_PACKAGES=(
  "pandas|pandas|tidy tables"
  "numpy|numpy|everything"
  "h5py|h5py|read loom column attributes without SeuratDisk"
  "scanpy|scanpy|already in sc_pre"
  "anndata|anndata|already in sc_pre"
  "pysam|pysam|fragment and BAM access"
  "snapatac2|snapatac2 (pip)|fragment ingestion, per-cell-type bigWig export"
  "pygenometracks|pygenometracks|the composite figure"
  "matplotlib|matplotlib|figures"
  "Bio|biopython|FASTA handling for motif input"
  "logomaker|logomaker (pip)|motif logos, optional"
)

# Reference files expected under REF_DIR once SETUP_02 has run.
REF_FILES=(
  "genome.fa|Mmul8 FASTA with chr-prefixed contigs"
  "genome.fa.fai|samtools faidx index"
  "genes.gtf|Mmul8 annotation, contigs renamed to match"
  "jaspar_core_vertebrates.pfm|JASPAR motifs for mkref input_motifs and MEME"
)

# Tidy tables the reference validation test depends on.
TIDY_FILES=(
  "01_tidy/peak_gene_map.csv|268,436 vendor (peak, gene, TSS distance) triples"
  "01_tidy/peaks_da.csv|DA peaks with accessibility statistics"
  "01_tidy/genes_rna.csv|tested gene universe with expression statistics"
)

# -----------------------------------------------------------------------------
# Plumbing
# -----------------------------------------------------------------------------

mkdir -p "${OUT_DIR}" 2>/dev/null || {
  echo "FATAL: cannot create ${OUT_DIR}"; exit 1; }

REPORT="${OUT_DIR}/audit_report.txt"
MANIFEST="${OUT_DIR}/missing_manifest.tsv"
: > "${REPORT}"
: > "${MANIFEST}"
printf "category\titem\tconda_package\tstatus\tnote\n" >> "${MANIFEST}"

N_OK=0; N_MISSING=0; N_WARN=0

say()  { echo "$*" | tee -a "${REPORT}"; }
hdr()  { say ""; say "=============================================================================="; say "$*"; say "=============================================================================="; }
sub()  { say ""; say "--- $* ---"; }

ok()   { say "  [ok]      $1"; N_OK=$((N_OK+1)); }
miss() { say "  [MISSING] $1"; N_MISSING=$((N_MISSING+1));
         printf "%s\t%s\t%s\tmissing\t%s\n" "$2" "$3" "$4" "${5:-}" >> "${MANIFEST}"; }
warnx(){ say "  [warn]    $1"; N_WARN=$((N_WARN+1));
         printf "%s\t%s\t%s\twarn\t%s\n" "${2:-general}" "${3:-}" "${4:-}" "$1" >> "${MANIFEST}"; }

hdr "SETUP_00_audit  —  $(date '+%Y-%m-%d %H:%M:%S')"
say "host:          $(hostname 2>/dev/null || echo unknown)"
say "user:          ${USER:-unknown}"
say "project root:  ${PROJECT_ROOT}"
say "analysis root: ${ANALYSIS_ROOT}"
say "report:        ${REPORT}"
say ""
say "This script changes nothing. It only looks and reports."

# -----------------------------------------------------------------------------
# 1. conda / mamba and the environment list
# -----------------------------------------------------------------------------

hdr "1. Conda and environments"

CONDA_OK=0
if [[ -x "${CONDA_BASE}/bin/conda" ]]; then
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/bin/activate" base 2>/dev/null || true
  CONDA_OK=1
  ok "conda found: $(conda --version 2>&1 | head -1)"
else
  miss "conda not found at ${CONDA_BASE}/bin/conda" "conda" "conda" "n/a" \
       "adjust CONDA_BASE at the top of this script"
fi

if command -v mamba >/dev/null 2>&1; then
  ok "mamba found: $(mamba --version 2>&1 | head -1)"
else
  warnx "mamba not found. Strongly recommended: a conda solve across this many Bioconductor packages can run for hours or fail. Install with: conda install -n base -c conda-forge mamba" \
        "conda" "mamba" "mamba"
fi

sub "Existing environments"
if [[ ${CONDA_OK} -eq 1 ]]; then
  conda env list 2>/dev/null | tee -a "${REPORT}"

  if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "${SOURCE_ENV}"; then
    ok "source env '${SOURCE_ENV}' exists and can be cloned"
  else
    miss "source env '${SOURCE_ENV}' not found" "conda" "${SOURCE_ENV}" "n/a" \
         "nothing to clone from; SETUP_01 will need to build from scratch"
  fi

  if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "${TARGET_ENV}"; then
    warnx "target env '${TARGET_ENV}' ALREADY EXISTS. SETUP_01 must not clobber it. Rename, remove, or pick another name." \
          "conda" "${TARGET_ENV}" "n/a"
    TARGET_EXISTS=1
  else
    say "  [note]    target env '${TARGET_ENV}' does not exist yet, as expected"
    TARGET_EXISTS=0
  fi
else
  TARGET_EXISTS=0
fi

# -----------------------------------------------------------------------------
# 2. Activate the best available environment for probing
#
# If sc_pre_R already exists, probe that. Otherwise probe sc_pre, which tells us
# what the clone will inherit.
# -----------------------------------------------------------------------------

hdr "2. Environment used for the probes below"

PROBE_ENV=""
if [[ ${CONDA_OK} -eq 1 ]]; then
  if [[ ${TARGET_EXISTS} -eq 1 ]]; then
    PROBE_ENV="${TARGET_ENV}"
  elif conda env list 2>/dev/null | awk '{print $1}' | grep -qx "${SOURCE_ENV}"; then
    PROBE_ENV="${SOURCE_ENV}"
  fi
fi

if [[ -n "${PROBE_ENV}" ]]; then
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/bin/activate" "${PROBE_ENV}" 2>/dev/null || true
  export PATH="${CONDA_PREFIX:-}/bin:${PATH}"
  say "probing inside: ${PROBE_ENV}"
  say "CONDA_PREFIX:   ${CONDA_PREFIX:-<unset>}"
  say ""
  say "Anything reported missing below is what SETUP_01 needs to install into"
  say "${TARGET_ENV}."
else
  warnx "no conda environment available to probe; binary and package checks will reflect the bare system only"
  say "probing inside: <system PATH only>"
fi

# -----------------------------------------------------------------------------
# 3. Command-line binaries
# -----------------------------------------------------------------------------

hdr "3. Command-line binaries"

for spec in "${BINARIES[@]}"; do
  IFS='|' read -r bin pkg why <<< "${spec}"
  if command -v "${bin}" >/dev/null 2>&1; then
    path="$(command -v "${bin}")"
    case "${bin}" in
      R|Rscript) ver="$("${bin}" --version 2>&1 | head -1)" ;;
      samtools|bedtools) ver="$("${bin}" --version 2>&1 | head -1)" ;;
      sea|xstreme|tomtom|fimo|meme-chip) ver="$("${bin}" -version 2>&1 | head -1)" ;;
      *) ver="$("${bin}" --version 2>&1 | head -1)" ;;
    esac
    ok "$(printf '%-12s %s' "${bin}" "${ver}")"
    say "              path: ${path}"

    # The known hazard: a samtools from outside conda.
    if [[ "${bin}" == "samtools" && -n "${CONDA_PREFIX:-}" ]]; then
      if [[ "${path}" != "${CONDA_PREFIX}"* ]]; then
        warnx "samtools resolves to ${path}, OUTSIDE the conda env. The system binary at /usr/local/bin/samtools is broken. Every script must 'export PATH=\"\${CONDA_PREFIX}/bin:\${PATH}\"' after activation." \
              "binary" "samtools" "samtools"
      fi
    fi
  else
    miss "$(printf '%-12s %s' "${bin}" "(${why})")" "binary" "${bin}" "${pkg}" "${why}"
  fi
done

# -----------------------------------------------------------------------------
# 4. MEME motif databases
#
# A MEME install from bioconda or the tarball usually ships db/motif_databases,
# which includes JASPAR CORE vertebrates. Finding it saves a download.
# -----------------------------------------------------------------------------

hdr "4. MEME motif databases"

MEME_DB_FOUND=""
if command -v sea >/dev/null 2>&1; then
  sea_path="$(command -v sea)"
  sea_prefix="$(dirname "$(dirname "${sea_path}")")"
  for cand in \
      "${sea_prefix}/share/meme"*"/db/motif_databases" \
      "${sea_prefix}/db/motif_databases" \
      "${sea_prefix}/share/meme/motif_databases" \
      "${CONDA_PREFIX:-/nonexistent}/share/meme"*"/db/motif_databases" ; do
    for d in ${cand}; do
      [[ -d "${d}" ]] && MEME_DB_FOUND="${d}" && break 2
    done
  done
fi

if [[ -n "${MEME_DB_FOUND}" ]]; then
  ok "motif_databases found: ${MEME_DB_FOUND}"
  jaspar_hits="$(find "${MEME_DB_FOUND}" -iname '*JASPAR*' -maxdepth 3 2>/dev/null | head -20)"
  if [[ -n "${jaspar_hits}" ]]; then
    say "  JASPAR files present:"
    echo "${jaspar_hits}" | sed 's/^/    /' | tee -a "${REPORT}"
    vert="$(echo "${jaspar_hits}" | grep -i 'vertebrate' | head -1)"
    if [[ -n "${vert}" ]]; then
      ok "JASPAR CORE vertebrates available: ${vert}"
    else
      warnx "JASPAR present but no vertebrates file matched; check the listing above" \
            "motifs" "jaspar_vertebrates" "download"
    fi
  else
    warnx "motif_databases exists but contains no JASPAR files; download from jaspar.elixir.no" \
          "motifs" "jaspar" "download"
  fi
else
  miss "MEME motif_databases not found (JASPAR will need downloading separately)" \
       "motifs" "motif_databases" "meme / download" \
       "bioconda meme sometimes omits the databases; grab them from meme-suite.org or JASPAR directly"
fi

# -----------------------------------------------------------------------------
# 5. R packages
# -----------------------------------------------------------------------------

hdr "5. R packages"

if command -v Rscript >/dev/null 2>&1; then
  say "R library paths:"
  Rscript -e 'cat(paste("   ", .libPaths(), collapse="\n"), "\n")' 2>/dev/null | tee -a "${REPORT}"
  say ""
  for spec in "${R_PACKAGES[@]}"; do
    IFS='|' read -r rpkg pkg why <<< "${spec}"
    ver="$(Rscript -e "suppressWarnings(suppressMessages({ if (requireNamespace('${rpkg}', quietly=TRUE)) cat(as.character(utils::packageVersion('${rpkg}'))) else cat('NA') }))" 2>/dev/null)"
    if [[ -n "${ver}" && "${ver}" != "NA" ]]; then
      ok "$(printf '%-24s %s' "${rpkg}" "${ver}")"
    else
      miss "$(printf '%-24s %s' "${rpkg}" "(${why})")" "r_package" "${rpkg}" "${pkg}" "${why}"
    fi
  done
else
  say "  Rscript not available, so no R packages can be checked."
  say "  All of the following will be needed in ${TARGET_ENV}:"
  for spec in "${R_PACKAGES[@]}"; do
    IFS='|' read -r rpkg pkg why <<< "${spec}"
    miss "$(printf '%-24s %s' "${rpkg}" "(${why})")" "r_package" "${rpkg}" "${pkg}" "${why}"
  done
fi

# -----------------------------------------------------------------------------
# 6. Python packages
# -----------------------------------------------------------------------------

hdr "6. Python packages"

if command -v python >/dev/null 2>&1; then
  ok "python: $(python --version 2>&1)"
  say "              path: $(command -v python)"
  say ""
  for spec in "${PY_PACKAGES[@]}"; do
    IFS='|' read -r mod pkg why <<< "${spec}"
    ver="$(python - <<PYEOF 2>/dev/null
import importlib, sys
try:
    m = importlib.import_module("${mod}")
    print(getattr(m, "__version__", "present"))
except Exception:
    sys.exit(1)
PYEOF
)"
    if [[ -n "${ver}" ]]; then
      ok "$(printf '%-16s %s' "${mod}" "${ver}")"
    else
      miss "$(printf '%-16s %s' "${mod}" "(${why})")" "py_package" "${mod}" "${pkg}" "${why}"
    fi
  done
else
  miss "python not on PATH" "py_package" "python" "python" "no environment active"
fi

# -----------------------------------------------------------------------------
# 7. cellranger-arc
# -----------------------------------------------------------------------------

hdr "7. cellranger-arc ${ARC_EXPECTED_VERSION}"

ARC_BIN=""
if command -v cellranger-arc >/dev/null 2>&1; then
  ARC_BIN="$(command -v cellranger-arc)"
elif [[ -x "${ARC_EXPECTED_DIR}/cellranger-arc" ]]; then
  ARC_BIN="${ARC_EXPECTED_DIR}/cellranger-arc"
fi

if [[ -n "${ARC_BIN}" ]]; then
  ok "cellranger-arc found: ${ARC_BIN}"
  arc_ver="$("${ARC_BIN}" --version 2>&1 | head -2 | tr '\n' ' ')"
  say "              version: ${arc_ver}"
  if echo "${arc_ver}" | grep -q "${ARC_EXPECTED_VERSION}"; then
    ok "version matches the ${ARC_EXPECTED_VERSION} Genohub used"
  else
    warnx "version does NOT look like ${ARC_EXPECTED_VERSION}. ARC peak calling changed between 1.0 and 2.0, so a mismatch will not reproduce the vendor's 204,136 / 345,553 peak universes." \
          "software" "cellranger-arc" "10x download"
  fi
  say ""
  say "  Next step once the environment is built, and do not skip it:"
  say "    cellranger-arc testrun --id=tiny"
  say "  This is 10x's own smoke test over the full pipeline. A broken install"
  say "  found 18 runs into a queue is expensive."
  if command -v "${ARC_BIN}" >/dev/null 2>&1; then
    say ""
    say "  Also worth capturing for the methods record:"
    say "    ${ARC_BIN} sitecheck > ${OUT_DIR}/arc_sitecheck.txt"
  fi
else
  miss "cellranger-arc not found (looked on PATH and at ${ARC_EXPECTED_DIR})" \
       "software" "cellranger-arc" "10x download" \
       "self-contained tarball, Linux x86_64, no conda. Download URLs are signed and expire."
  say ""
  say "  Note: CellRanger v10.0.0 is the gene expression pipeline and cannot do"
  say "  multiome. cellranger-arc is separate software and is not bundled with it."
fi

# -----------------------------------------------------------------------------
# 8. Reference files
# -----------------------------------------------------------------------------

hdr "8. Mmul8 reference"

say "expected location: ${REF_DIR}"
say ""
say "Recipe we believe Genohub followed, from the evidence in the data"
say "(Note 03, Section 3.3):"
say "  Peaks are chr1..chr20, chrX, chrY. 22 contigs, no scaffolds, no chrM."
say "  Gene symbols include C1H1orf35 and C1H1orf52, which is Ensembl/VGNC rhesus"
say "  ortholog naming. Ensembl symbols with UCSC-style contig names is what you"
say "  get by following 10x's own reference recipe, which renames Ensembl contigs"
say "  to chr*. No chrM peaks is consistent with chrM declared in"
say "  non_nuclear_contigs, which excludes it from peak calling."
say ""

if [[ -d "${REF_DIR}" ]]; then
  ok "reference directory exists"
  for spec in "${REF_FILES[@]}"; do
    IFS='|' read -r f why <<< "${spec}"
    p="${REF_DIR}/${f}"
    if [[ -s "${p}" ]]; then
      sz="$(du -h "${p}" 2>/dev/null | cut -f1)"
      ok "$(printf '%-34s %8s' "${f}" "${sz}")"
      if [[ "${f}" == "genome.fa" ]]; then
        say "  contig names (first 5):"
        grep '^>' "${p}" 2>/dev/null | head -5 | sed 's/^/    /' | tee -a "${REPORT}"
        n_contigs="$(grep -c '^>' "${p}" 2>/dev/null || echo 0)"
        say "  total contigs: ${n_contigs}"
        if grep -q '^>chr1\b' "${p}" 2>/dev/null; then
          ok "contigs are chr-prefixed, matching the peak coordinates"
        else
          warnx "contigs are NOT chr-prefixed. They must be renamed to chr1..chr20, chrX, chrY or every join against the existing peak coordinates fails silently." \
                "reference" "genome.fa" "rename"
        fi
      fi
      if [[ "${f}" == "genes.gtf" ]]; then
        say "  first non-comment contig name:"
        grep -v '^#' "${p}" 2>/dev/null | head -1 | cut -f1 | sed 's/^/    /' | tee -a "${REPORT}"
      fi
    else
      miss "$(printf '%-34s %s' "${f}" "(${why})")" "reference" "${f}" "download" "${why}"
    fi
  done
else
  miss "reference directory ${REF_DIR} does not exist" "reference" "${REF_DIR}" "n/a" \
       "SETUP_02 will create and populate it"
  for spec in "${REF_FILES[@]}"; do
    IFS='|' read -r f why <<< "${spec}"
    miss "$(printf '%-34s %s' "${f}" "(${why})")" "reference" "${f}" "download" "${why}"
  done
fi

sub "Existing ARC reference"
if [[ -d "${REF_DIR}/Mmul8_arc" ]]; then
  ok "an ARC reference directory exists at ${REF_DIR}/Mmul8_arc"
  [[ -s "${REF_DIR}/Mmul8_arc/reference.json" ]] && \
    { say "  reference.json:"; sed 's/^/    /' "${REF_DIR}/Mmul8_arc/reference.json" | head -30 | tee -a "${REPORT}"; }
else
  say "  [note]    no ARC reference built yet, as expected"
fi

# -----------------------------------------------------------------------------
# 9. Tidy tables, which the reference validation test depends on
# -----------------------------------------------------------------------------

hdr "9. Tidy tables from MULTIOME_01"

say "The reference acceptance test recomputes TSS distances for the vendor's"
say "268,436 (peak, gene, distance) triples from our own GTF and compares."
say "A wrong reference does not announce itself: every join still runs and every"
say "figure still renders, pointing at the wrong places. This test is the only"
say "cheap way to catch it, and it must run BEFORE mkref and before any"
say "alignment job is queued."
say ""

for spec in "${TIDY_FILES[@]}"; do
  IFS='|' read -r f why <<< "${spec}"
  p="${ANALYSIS_ROOT}/${f}"
  if [[ -s "${p}" ]]; then
    n="$(wc -l < "${p}" 2>/dev/null)"
    sz="$(du -h "${p}" 2>/dev/null | cut -f1)"
    ok "$(printf '%-28s %10s lines  %8s' "${f}" "${n}" "${sz}")"
  else
    miss "$(printf '%-28s %s' "${f}" "(${why})")" "tidy" "${f}" "MULTIOME_01" "${why}"
  fi
done

# -----------------------------------------------------------------------------
# 10. Disk headroom
# -----------------------------------------------------------------------------

hdr "10. Disk headroom"

say "Input is 1.7 TB across 20 library pairs, roughly 85 GB per pair."
say "ARC output is dominated by atac_possorted_bam.bam and gex_possorted_bam.bam,"
say "which run comparable to the input FASTQ size, plus pipeline working space on"
say "top. Expect on the order of 85 to 120 GB of outs/ per run."
say ""
say "Checking ${DISK_CHECK_PATH} against a conservative ${DISK_NEEDED_TB} TB floor."
say ""

if [[ -d "${DISK_CHECK_PATH}" ]]; then
  df -h "${DISK_CHECK_PATH}" 2>/dev/null | tee -a "${REPORT}"
  avail_kb="$(df -Pk "${DISK_CHECK_PATH}" 2>/dev/null | awk 'NR==2{print $4}')"
  if [[ -n "${avail_kb}" ]]; then
    avail_tb="$(awk -v k="${avail_kb}" 'BEGIN{printf "%.2f", k/1024/1024/1024}')"
    say ""
    say "  available: ${avail_tb} TB"
    if awk -v a="${avail_tb}" -v n="${DISK_NEEDED_TB}" 'BEGIN{exit !(a >= n)}'; then
      ok "headroom looks sufficient for all 20 runs"
    else
      warnx "only ${avail_tb} TB available against a ${DISK_NEEDED_TB} TB estimate. Do not queue 20 runs against this. Run one library, measure what it actually writes, then decide whether to stage in batches and move BAMs off as they complete." \
            "disk" "${DISK_CHECK_PATH}" "n/a"
    fi
  fi
else
  warnx "${DISK_CHECK_PATH} not found; adjust DISK_CHECK_PATH at the top of this script" \
        "disk" "${DISK_CHECK_PATH}" "n/a"
fi

sub "Input size, measured"
if [[ -d "${PROJECT_ROOT}" ]]; then
  say "  (du across 1.7 TB takes a moment)"
  du -sh "${PROJECT_ROOT}" 2>/dev/null | sed 's/^/  /' | tee -a "${REPORT}"
  for d in "${PROJECT_ROOT}/Raw Data_Colon_Brain"/*; do
    [[ -d "${d}" ]] && du -sh "${d}" 2>/dev/null | sed 's/^/  /' | tee -a "${REPORT}"
  done
else
  warnx "project root ${PROJECT_ROOT} not found" "disk" "${PROJECT_ROOT}" "n/a"
fi

# -----------------------------------------------------------------------------
# 11. Try the solve
#
# Dry run only. Nothing is installed. Catches an unsolvable combination before
# it costs an afternoon.
# -----------------------------------------------------------------------------

hdr "11. Dependency solve, dry run"

if [[ ${TRY_SOLVE} -eq 1 && ${CONDA_OK} -eq 1 ]]; then
  SOLVER="conda"
  command -v mamba >/dev/null 2>&1 && SOLVER="mamba"

  CONDA_PKGS=(
    r-base r-seurat r-signac r-seuratobject r-matrix r-data.table
    bioconductor-muscat bioconductor-edger bioconductor-singlecellexperiment
    bioconductor-genomicranges bioconductor-rtracklayer
    meme bedtools samtools htslib pysam
    pygenometracks matplotlib biopython
  )

  say "solver:   ${SOLVER}"
  say "channels: conda-forge, bioconda"
  say "packages: ${CONDA_PKGS[*]}"
  say ""
  say "Running a dry-run solve. This takes minutes and installs nothing."
  say ""

  solve_log="${OUT_DIR}/solve_dryrun.log"
  if ${SOLVER} create --dry-run -n "__audit_probe_delete_me" \
       -c conda-forge -c bioconda "${CONDA_PKGS[@]}" \
       > "${solve_log}" 2>&1; then
    ok "the solve succeeds; a single ${TARGET_ENV} can hold the whole stack"
    say "  full log: ${solve_log}"
    say ""
    say "  snapatac2 and logomaker are pip-only and are not in this solve."
    say "  Install them with pip after the env is created."
  else
    warnx "the combined solve FAILED. See ${solve_log}. Recommended fallback, following the precedent of isolating viral_id_genomad rather than fighting the solver: split into sc_pre_R (R stack, samtools, bedtools, MEME) and sc_viz (snapatac2, pyGenomeTracks, matplotlib)." \
          "conda" "solve" "n/a"
    say ""
    say "  last 40 lines of the solver output:"
    tail -40 "${solve_log}" 2>/dev/null | sed 's/^/    /' | tee -a "${REPORT}"
  fi
else
  say "  skipped (TRY_SOLVE=0 or conda unavailable)"
fi

# -----------------------------------------------------------------------------
# 12. Summary
# -----------------------------------------------------------------------------

hdr "12. Summary"

say "  ok:       ${N_OK}"
say "  missing:  ${N_MISSING}"
say "  warnings: ${N_WARN}"
say ""
say "  report:   ${REPORT}"
say "  manifest: ${MANIFEST}"
say ""

if [[ ${N_MISSING} -gt 0 ]]; then
  say "Missing items by category:"
  awk -F'\t' 'NR>1 && $4=="missing" {print "  " $1 "\t" $2 "\t" $3}' "${MANIFEST}" \
    | sort | tee -a "${REPORT}"
  say ""
fi

if [[ ${N_WARN} -gt 0 ]]; then
  say "Warnings:"
  awk -F'\t' 'NR>1 && $4=="warn" {print "  " $5}' "${MANIFEST}" | tee -a "${REPORT}"
  say ""
fi

say "Next steps, in order:"
say "  1. Read ${REPORT} and the manifest."
say "  2. Confirm the disk headroom decision in Section 10."
say "  3. Pull cellranger-arc 2.0.2 from 10x if Section 7 reported it missing."
say "  4. Run SETUP_01_build_env.sh, which will be written against this manifest."
say "  5. cellranger-arc testrun --id=tiny"
say "  6. Download and prepare the reference, then run the TSS validation test"
say "     from Section 9 BEFORE mkref and before queueing any alignment."
say ""
say "Done. Nothing was changed."

exit 0
