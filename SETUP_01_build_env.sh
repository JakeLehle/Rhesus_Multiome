#!/usr/bin/env bash
# =============================================================================
# SETUP_01_plan_env.sh   (revision 2)
#
# Reports what is missing and prints the commands to install it. Installs
# nothing. Runs in seconds.
#
#   bash SETUP_01_plan_env.sh          what is missing, and the commands
#   bash SETUP_01_plan_env.sh verify   did the installs work
#
# -----------------------------------------------------------------------------
# WHAT CHANGED, AND WHY
# -----------------------------------------------------------------------------
#
# FLAG FIX. Revision 1 emitted `--channel-priority strict`. That is the
# CONFIG FILE syntax. On the command line `--channel-priority` is a BOOLEAN
# flag, so conda parsed `strict` as a package name and libmamba reported
# "nothing provides requested strict". The correct flag is
# `--strict-channel-priority`.
#
# DROPPED --override-channels. Combined with that stray positional it broke
# argparse before the solve was reached, producing a bare `usage: conda` line
# instead of a useful error. It also buys little: strict priority already stops
# the plotly and hcc channels from winning anything.
#
# THREE ENVIRONMENTS INSTEAD OF ONE. The reason for cloning sc_pre was to put R
# in it. That premise is gone, and the same argument applies to MEME: nothing
# calls MEME from Python in-process, since the motif scripts invoke sea, fimo
# and bedtools as subprocesses. So there is no coupling to gain by co-locating
# them with scanpy, and a working scanpy env to risk perturbing.
#
#     meme_suite   meme, bedtools            motif track, no Python or R
#     multiome_R   R and Bioconductor        object rebuild, pseudobulk DA
#     sc_pre_R     the sc_pre clone + pysam  parsing, tidy tables, figures
#
# Each solves in seconds. None constrains the others.
#
# PIN REPORTING. sc_pre_R inherits a `python 3.11.*` pin from sc_pre, which
# surfaced in the failed solve. It was almost certainly collateral from the
# bogus `strict` spec, since bioconda's meme has no Python dependency, but the
# plan now reports pins per environment so a real pin conflict is visible before
# a solve rather than during one.
#
# DRY-RUN PROBES. Every command is paired with a --dry-run form. A tractable
# solve returns in under a minute; one that would hang shows it for free.
#
# Jake Lehle, 2026-10
# =============================================================================

set -uo pipefail

# -----------------------------------------------------------------------------
# Cell 1: parameters
# -----------------------------------------------------------------------------

CONDA_BASE="/master/jlehle/anaconda3"

MOTIF_ENV="meme_suite"     # meme, bedtools
R_ENV="multiome_R"         # R and Bioconductor
PY_ENV="sc_pre_R"          # the sc_pre clone, Python side

ANALYSIS_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
OUT_DIR="${ANALYSIS_ROOT}/setup_audit"
CMD_FILE="${OUT_DIR}/install_commands.sh"
ARC_DIR="/master/jlehle/cellranger-arc-2.0.2"

MODE="${1:-plan}"

NEED_MOTIF=(meme bedtools)
NEED_R=(
  r-base r-seurat r-signac r-seuratobject r-matrix r-data.table
  bioconductor-muscat bioconductor-edger bioconductor-singlecellexperiment
  bioconductor-genomicranges bioconductor-rtracklayer
)
NEED_PY=(pysam)
NEED_PIP=(snapatac2)
# Deferred: pygenometracks pins matplotlib tightly and sc_pre_R has a working
# matplotlib 3.10.9 with scanpy on it. Own env, when the coverage figures come
# up, which is blocked on the rebuild anyway.
NEED_DEFERRED=(pygenometracks)

# -----------------------------------------------------------------------------
# Cell 2: plumbing
# -----------------------------------------------------------------------------

mkdir -p "${OUT_DIR}" 2>/dev/null || { echo "FATAL: cannot create ${OUT_DIR}"; exit 1; }

say()  { echo "$*"; }
hdr()  { say ""; say "=============================================================================="; say "$*"; say "=============================================================================="; }
sub()  { say ""; say "--- $* ---"; }
ok()   { say "  [ok]      $*"; }
miss() { say "  [missing] $*"; }
warnx(){ say "  [warn]    $*"; }
die()  { say ""; say "FATAL: $*"; exit 1; }

hdr "SETUP_01_plan_env rev2  —  mode=${MODE}  —  $(date '+%Y-%m-%d %H:%M:%S')"
say "host: $(hostname)"

[[ -x "${CONDA_BASE}/bin/conda" ]] || die "conda not found at ${CONDA_BASE}/bin/conda"
# shellcheck disable=SC1091
source "${CONDA_BASE}/bin/activate" base >/dev/null 2>&1 || die "could not activate base"

env_exists() { conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"; }
env_prefix() { conda env list 2>/dev/null | awk -v e="$1" '$1==e{print $NF}'; }

# -----------------------------------------------------------------------------
# Cell 3: solver and flags
# -----------------------------------------------------------------------------

hdr "Solver and flags"

CONDA_VER="$(conda --version 2>&1 | awk '{print $2}')"
say "  conda ${CONDA_VER}"

LIBMAMBA_VER="$(conda list -n base 2>/dev/null | awk '$1=="conda-libmamba-solver"{print $2}' | head -1)"
if [[ -n "${LIBMAMBA_VER}" ]]; then
  ok "conda-libmamba-solver ${LIBMAMBA_VER}; --solver=libmamba is accepted"
else
  miss "conda-libmamba-solver absent. Install it first or every solve stays slow:"
  say  "        conda install -n base -c conda-forge conda-libmamba-solver"
fi

SOLVER_FLAG="--solver=libmamba"
# --strict-channel-priority, NOT --channel-priority strict. The latter is config
# file syntax; on the CLI --channel-priority is boolean and `strict` is then
# parsed as a package name.
CHAN_FLAGS="--strict-channel-priority -c conda-forge -c bioconda"

say ""
say "  flags used below: ${SOLVER_FLAG} ${CHAN_FLAGS}"

sub "channels configured"
conda config --show channels 2>/dev/null | sed 's/^/    /'
CP="$(conda config --show channel_priority 2>/dev/null | awk '{print $2}')"
say ""
say "  channel_priority: ${CP:-<unset, meaning flexible>}"
[[ "${CP}" != "strict" ]] && say "  Each command passes --strict-channel-priority rather than changing your global config."

# -----------------------------------------------------------------------------
# Cell 4: current state, including pins
# -----------------------------------------------------------------------------

hdr "Current state"

have() { grep -qx "$2" "${OUT_DIR}/.pkg_$1" 2>/dev/null; }

for e in "${MOTIF_ENV}" "${R_ENV}" "${PY_ENV}"; do
  if env_exists "${e}"; then
    conda list -n "${e}" 2>/dev/null | awk '!/^#/{print $1}' > "${OUT_DIR}/.pkg_${e}"
    ok "${e} exists, $(wc -l < "${OUT_DIR}/.pkg_${e}") packages"
    pf="$(env_prefix "${e}")/conda-meta/pinned"
    if [[ -s "${pf}" ]]; then
      warnx "  pinned specs in ${e}:"
      sed 's/^/        /' "${pf}"
      warnx "  A pin that conflicts with a requested package makes a solve fail or crawl. This is what surfaced in the failed meme solve."
    fi
  else
    : > "${OUT_DIR}/.pkg_${e}"
    miss "${e} does not exist"
  fi
done

MISS_MOTIF=(); MISS_R=(); MISS_PY=()

sub "motif track, ${MOTIF_ENV}"
for p in "${NEED_MOTIF[@]}"; do
  if have "${MOTIF_ENV}" "${p}"; then ok "${p}"; else miss "${p}"; MISS_MOTIF+=("${p}"); fi
done

sub "R stack, ${R_ENV}"
for p in "${NEED_R[@]}"; do
  if have "${R_ENV}" "${p}"; then ok "${p}"; else miss "${p}"; MISS_R+=("${p}"); fi
done

sub "python extras, ${PY_ENV}"
for p in "${NEED_PY[@]}"; do
  if have "${PY_ENV}" "${p}"; then ok "${p}"; else miss "${p}"; MISS_PY+=("${p}"); fi
done

sub "already present, not reinstalled"
for p in samtools htslib python pandas numpy h5py scanpy anndata matplotlib biopython; do
  have "${PY_ENV}" "${p}" && ok "${p} in ${PY_ENV}"
done

# -----------------------------------------------------------------------------
# Cell 5: plan
# -----------------------------------------------------------------------------

if [[ "${MODE}" == "plan" ]]; then

hdr "Commands"

: > "${CMD_FILE}"
emit() { say "$*"; echo "$*" >> "${CMD_FILE}"; }

emit "#!/usr/bin/env bash"
emit "# Generated by SETUP_01_plan_env.sh on $(date '+%Y-%m-%d %H:%M:%S')"
emit "# Run one at a time, not as a script, so a slow solve is visible."
emit ""

if [[ -z "${LIBMAMBA_VER}" ]]; then
  say ""
  say "### 0. Solver first, or everything below crawls."
  emit "conda install -y -n base -c conda-forge conda-libmamba-solver"
  emit ""
fi

if [[ ${#MISS_MOTIF[@]} -gt 0 ]]; then
  say ""
  say "### 1. Motif track, own env. Highest priority: this is all that blocks"
  say "###    Track B once genome.fa exists. No Python or R entanglement, so"
  say "###    the solve is fast and nothing else is at risk."
  say "###    Probe first:"
  say "###      conda create --dry-run -n _probe ${SOLVER_FLAG} ${CHAN_FLAGS} ${MISS_MOTIF[*]}"
  emit "conda create -y -n ${MOTIF_ENV} ${SOLVER_FLAG} ${CHAN_FLAGS} ${MISS_MOTIF[*]}"
  emit ""
fi

if [[ ${#MISS_R[@]} -gt 0 ]]; then
  say ""
  say "### 2. R stack, fresh create. No pins inherited from sc_pre, which is"
  say "###    what a create buys over an install."
  say "###    Probe first:"
  say "###      conda create --dry-run -n _probe ${SOLVER_FLAG} ${CHAN_FLAGS} ${MISS_R[*]}"
  emit "conda create -y -n ${R_ENV} ${SOLVER_FLAG} ${CHAN_FLAGS} ${MISS_R[*]}"
  emit ""
fi

if [[ ${#MISS_PY[@]} -gt 0 ]]; then
  say ""
  say "### 3. Python extras. Not urgent: these serve coverage figures, which"
  say "###    are blocked on the rebuild."
  emit "source ${CONDA_BASE}/bin/activate ${PY_ENV}"
  emit "conda install -y ${SOLVER_FLAG} ${CHAN_FLAGS} ${MISS_PY[*]}"
  emit "python -m pip install ${NEED_PIP[*]}"
  emit ""
fi

say ""
say "### Deferred: ${NEED_DEFERRED[*]}"
say "###   conda create -n sc_viz ${SOLVER_FLAG} ${CHAN_FLAGS} ${NEED_DEFERRED[*]}"

chmod +x "${CMD_FILE}" 2>/dev/null
say ""
ok "also written to ${CMD_FILE}"

hdr "How the scripts find these tools"

say "  Three environments means the analysis scripts need the right one on PATH."
say "  MEME and bedtools are called as subprocesses, so a PATH prepend is all it"
say "  takes. The motif script will do this:"
say ""
say "      MOTIF_BIN=\"${CONDA_BASE}/envs/${MOTIF_ENV}/bin\""
say "      export PATH=\"\${MOTIF_BIN}:\${PATH}\""
say ""
say "  and the R steps run through:"
say ""
say "      conda run -n ${R_ENV} Rscript <script>"
say ""
say "  so no environment has to hold everything."

hdr "If a solve is still slow"

say "  Probe one package at a time to find the culprit. meme is the likely one"
say "  of the pair, because of its Perl dependency tree:"
say ""
say "      conda create --dry-run -n _probe ${SOLVER_FLAG} ${CHAN_FLAGS} bedtools"
say "      conda create --dry-run -n _probe ${SOLVER_FLAG} ${CHAN_FLAGS} meme"
say ""
say "  And add -v to watch a real install rather than waiting blind."

fi

# -----------------------------------------------------------------------------
# Cell 6: verify
# -----------------------------------------------------------------------------

if [[ "${MODE}" == "verify" ]]; then

hdr "Verification"

N_FAIL=0; N_MOTIF_FAIL=0

sub "motif binaries in ${MOTIF_ENV}"
if env_exists "${MOTIF_ENV}"; then
  BIN="$(env_prefix "${MOTIF_ENV}")/bin"
  for b in sea xstreme tomtom fimo meme-chip bedtools; do
    if [[ -x "${BIN}/${b}" ]]; then ok "$(printf '%-12s %s' "${b}" "${BIN}/${b}")"
    else miss "$(printf '%-12s not in %s' "${b}" "${MOTIF_ENV}")"; N_FAIL=$((N_FAIL+1)); N_MOTIF_FAIL=$((N_MOTIF_FAIL+1)); fi
  done
  say ""
  [[ -x "${BIN}/sea" ]] && { say "  sea:"; "${BIN}/sea" --version 2>&1 | head -1 | sed 's/^/    /'; }
  [[ -x "${BIN}/bedtools" ]] && { say "  bedtools:"; "${BIN}/bedtools" --version 2>&1 | head -1 | sed 's/^/    /'; }

  sub "MEME motif databases"
  DB=""
  for c in "$(env_prefix "${MOTIF_ENV}")"/share/meme*/db/motif_databases \
           "$(env_prefix "${MOTIF_ENV}")"/share/meme/motif_databases; do
    for d in ${c}; do [[ -d "${d}" ]] && DB="${d}" && break 2; done
  done
  if [[ -n "${DB}" ]]; then ok "${DB}"
  else warnx "bioconda meme ships without motif_databases, as the audit predicted. SETUP_02_reference.sh downloads JASPAR directly, so this is expected."; fi
else
  miss "${MOTIF_ENV} does not exist"
  N_FAIL=$((N_FAIL+1)); N_MOTIF_FAIL=$((N_MOTIF_FAIL+1))
fi

sub "R stack in ${R_ENV}"
if env_exists "${R_ENV}"; then
  for p in Seurat Signac SeuratObject muscat edgeR SingleCellExperiment GenomicRanges rtracklayer Matrix data.table; do
    v="$(conda run -n "${R_ENV}" Rscript -e "suppressWarnings(suppressMessages({ if (requireNamespace('${p}', quietly=TRUE)) cat(as.character(utils::packageVersion('${p}'))) else cat('NA') }))" 2>/dev/null)"
    if [[ -n "${v}" && "${v}" != "NA" ]]; then ok "$(printf '%-24s %s' "${p}" "${v}")"
    else miss "$(printf '%-24s not importable' "${p}")"; N_FAIL=$((N_FAIL+1)); fi
  done

  sub ".libPaths(), which must not include the system personal library"
  conda run -n "${R_ENV}" Rscript -e 'cat(paste("   ", .libPaths(), collapse="\n"), "\n")' 2>&1
  if conda run -n "${R_ENV}" Rscript -e 'q(status = as.integer(any(grepl("x86_64-redhat-linux-gnu-library", .libPaths()))))' 2>/dev/null; then
    ok "clean"
  else
    miss "the system R personal library is on .libPaths(). Conda R loading packages built against system R 4.4.3 segfaults in ways that look like package bugs. Fix with an activate.d hook, printed below."
    N_FAIL=$((N_FAIL+1))
  fi

  sub "Signac smoke test, real peak coordinates from this project"
  if conda run -n "${R_ENV}" Rscript -e '
suppressMessages({library(Signac); library(GenomicRanges); library(Matrix)})
gr <- GRanges(seqnames = c("chr1","chr1"),
              ranges   = IRanges(start = c(56999, 64098), end = c(57879, 65061)))
m  <- Matrix(matrix(c(1,0,2,3), nrow = 2), sparse = TRUE)
rownames(m) <- c("chr1-56999-57879","chr1-64098-65061")
colnames(m) <- c("cellA","cellB")
cat("ChromatinAssay built with",
    length(granges(CreateChromatinAssay(counts = m, ranges = gr))), "ranges\n")
' 2>&1 | tail -3 | sed 's/^/    /'; then
    ok "ran"
  else
    miss "Signac smoke test failed, which is the symptom of the .libPaths() problem above"
    N_FAIL=$((N_FAIL+1))
  fi

  RP="$(env_prefix "${R_ENV}")"
  if [[ -f "${RP}/etc/conda/activate.d/zz_rhesus_multiome.sh" ]]; then
    ok "activate.d hook present"
  else
    sub "install the R_LIBS hook"
    say "    conda activate ${R_ENV}"
    say "    mkdir -p \${CONDA_PREFIX}/etc/conda/activate.d"
    say "    cat > \${CONDA_PREFIX}/etc/conda/activate.d/zz_rhesus_multiome.sh <<'EOF'"
    say "    export R_LIBS=\"\${CONDA_PREFIX}/lib/R/library\""
    say "    export R_LIBS_USER=\"\""
    say "    export PATH=\"\${CONDA_PREFIX}/bin:\${PATH}\""
    say "    EOF"
  fi
else
  miss "${R_ENV} does not exist yet"
fi

sub "python in ${PY_ENV}"
if env_exists "${PY_ENV}"; then
  for m in pandas numpy h5py scanpy anndata matplotlib Bio logomaker pysam snapatac2; do
    v="$(conda run -n "${PY_ENV}" python -c "import importlib,sys
try:
    mm=importlib.import_module('${m}'); print(getattr(mm,'__version__','present'))
except Exception: sys.exit(1)" 2>/dev/null)"
    if [[ -n "${v}" ]]; then ok "$(printf '%-16s %s' "${m}" "${v}")"
    else warnx "$(printf '%-16s not importable' "${m}")"; fi
  done
fi

sub "cellranger-arc"
if [[ -x "${ARC_DIR}/cellranger-arc" ]]; then
  ok "${ARC_DIR}"
  "${ARC_DIR}/cellranger-arc" --version 2>&1 | head -2 | sed 's/^/    /'
else
  warnx "not found at ${ARC_DIR}"
fi

hdr "Summary"
say "  failures: ${N_FAIL}"
say ""
if [[ ${N_MOTIF_FAIL} -eq 0 ]]; then
  say "  MOTIF TRACK UNBLOCKED. Once SETUP_02_reference.sh lands genome.fa,"
  say "  Track B runs regardless of the R stack."
else
  say "  Motif track blocked: ${N_MOTIF_FAIL} tool(s) missing. Nearest deadline"
  say "  sits behind this one."
fi
say ""
say "  Next:"
say "    bash SETUP_02_reference.sh"
say "    python SETUP_03_validate_gtf.py"
say "    sbatch SETUP_04_mkref.sh"

fi

rm -f "${OUT_DIR}"/.pkg_* 2>/dev/null
say ""
exit 0
