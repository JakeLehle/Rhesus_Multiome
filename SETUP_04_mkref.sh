#!/usr/bin/env bash
#SBATCH --job-name=arc_mkref_Mmul8
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/arc_mkref_Mmul8_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/arc_mkref_Mmul8_%j.err
# =============================================================================
# SETUP_04_mkref.sh
#
# Build the cellranger-arc reference for Mmul8. Submit with sbatch, or run
# directly with bash for an interactive build.
#
# This is a STAR genome index build under the hood, so it is memory-hungry and
# slow: expect roughly 1 to 3 hours on a ~3 Gb genome at 16 threads, and give it
# real headroom on RAM. The SLURM header above asks for 128 GB, comfortably more
# than STAR needs for a genome this size.
#
# GATE
#   Refuses to run unless SETUP_03_validate_gtf.py has written a PASS verdict.
#   The reason is in Note 03, Section 3.4: a wrong reference does not announce
#   itself. Every join still runs, every figure still renders, and the peaks
#   point at the wrong places. Burning three hours and then 20 alignment runs on
#   an unvalidated reference is the expensive mistake this gate prevents.
#
#   To proceed on an INVESTIGATE verdict, set ALLOW_INVESTIGATE=1 and record why
#   in the manifest. Never override a FAIL.
#
# Jake Lehle, 2026-10
# =============================================================================

set -uo pipefail

# -----------------------------------------------------------------------------
# Cell 1: parameters
# -----------------------------------------------------------------------------

CONDA_BASE="/master/jlehle/anaconda3"
TARGET_ENV="sc_pre_R"
ARC_DIR="/master/jlehle/cellranger-arc-2.0.2"

REF_DIR="/master/jlehle/WORKING/REFS/Mmul8"
CONFIG="${REF_DIR}/Mmul8_arc.config"

ANALYSIS_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
OUT_DIR="${ANALYSIS_ROOT}/setup_audit"
VERDICT="${OUT_DIR}/gtf_validation_verdict.json"
LOG="${OUT_DIR}/setup_04_mkref.log"

# mkref writes the reference into the current working directory, named after the
# `genome` field in the config. Build in REF_DIR so it lands beside its inputs.
BUILD_DIR="${REF_DIR}"
EXPECTED_REF="${BUILD_DIR}/Mmul8"

NTHREADS="${SLURM_CPUS_PER_TASK:-16}"
MEMGB=120

# Only for a documented INVESTIGATE verdict. Never for FAIL.
ALLOW_INVESTIGATE=0

# -----------------------------------------------------------------------------
# Cell 2: plumbing
# -----------------------------------------------------------------------------

mkdir -p "${OUT_DIR}" "$(dirname "${LOG}")" 2>/dev/null
: > "${LOG}"

say()  { echo "$*" | tee -a "${LOG}"; }
hdr()  { say ""; say "=============================================================================="; say "$*"; say "=============================================================================="; }
ok()   { say "  [ok]    $*"; }
warnx(){ say "  [warn]  $*"; }
die()  { say ""; say "FATAL: $*"; say "log: ${LOG}"; exit 1; }

hdr "SETUP_04_mkref  —  $(date '+%Y-%m-%d %H:%M:%S')"
say "host:     $(hostname)"
say "job:      ${SLURM_JOB_ID:-interactive}"
say "threads:  ${NTHREADS}"
say "memgb:    ${MEMGB}"
say "log:      ${LOG}"

# -----------------------------------------------------------------------------
# Cell 3: the gate
# -----------------------------------------------------------------------------

hdr "Gate: GTF validation verdict"

[[ -f "${VERDICT}" ]] || die "no validation verdict at ${VERDICT}. Run SETUP_03_validate_gtf.py first. This gate exists because a wrong reference fails silently."

V="$(python - <<'PY' 2>/dev/null
import json, sys
try:
    d = json.load(open("/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS/setup_audit/gtf_validation_verdict.json"))
    print(d.get("verdict", "UNKNOWN"))
except Exception:
    print("UNREADABLE")
PY
)"

say "  verdict: ${V}"
say ""
python - <<'PY' 2>/dev/null | tee -a "${LOG}"
import json
d = json.load(open("/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS/setup_audit/gtf_validation_verdict.json"))
for k in ("best_tss_source", "best_peak_anchor", "frac_exact",
          "systematic_offset_bp", "n_triples_tested", "gene_label_coverage"):
    print(f"    {k:<24} {d.get(k)}")
print()
print("    " + str(d.get("message", "")))
PY

case "${V}" in
  PASS)
    ok "validation passed, proceeding"
    ;;
  INVESTIGATE)
    if [[ ${ALLOW_INVESTIGATE} -eq 1 ]]; then
      warnx "verdict is INVESTIGATE and ALLOW_INVESTIGATE=1. Proceeding deliberately."
      warnx "Write the reason into ${REF_DIR}/REFERENCE_MANIFEST.txt so the methods section can state it."
    else
      die "verdict is INVESTIGATE. Resolve it, or set ALLOW_INVESTIGATE=1 having decided the discrepancy is acceptable and documented why."
    fi
    ;;
  FAIL)
    die "verdict is FAIL. The reference does not reproduce the vendor's TSS distances. Do not build it and do not queue alignments. Check the assembly version and contig naming first."
    ;;
  *)
    die "verdict is '${V}', which is not a recognized value. Re-run SETUP_03_validate_gtf.py."
    ;;
esac

# -----------------------------------------------------------------------------
# Cell 4: environment and inputs
# -----------------------------------------------------------------------------

hdr "Environment"

# shellcheck disable=SC1091
source "${CONDA_BASE}/bin/activate" "${TARGET_ENV}" 2>/dev/null \
  || die "could not activate ${TARGET_ENV}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"
export PATH="${ARC_DIR}:${PATH}"
ok "env:  ${CONDA_PREFIX}"
ok "arc:  ${ARC_DIR}"

command -v cellranger-arc >/dev/null 2>&1 \
  || die "cellranger-arc not on PATH after adding ${ARC_DIR}"
say ""
cellranger-arc --version 2>&1 | head -2 | sed 's/^/    /' | tee -a "${LOG}"

# The known hazard: a samtools from outside the env.
sam="$(command -v samtools 2>/dev/null || echo none)"
if [[ "${sam}" == "${CONDA_PREFIX}"* ]]; then
  ok "samtools from the env: ${sam}"
else
  warnx "samtools resolves to ${sam}, outside the env. The system binary is broken. mkref does not call it, but later scripts do."
fi

hdr "Inputs"

[[ -f "${CONFIG}" ]] || die "config not found at ${CONFIG}. Run SETUP_02_reference.sh first."
ok "config: ${CONFIG}"
say ""
sed 's/^/    /' "${CONFIG}" | tee -a "${LOG}"

# Verify every path the config names actually exists.
say ""
say "--- checking the files the config points at ---"
while IFS= read -r p; do
  [[ -z "${p}" ]] && continue
  if [[ -s "${p}" ]]; then
    ok "$(printf '%-52s %8s' "$(basename "${p}")" "$(du -h "${p}" | cut -f1)")"
  else
    die "config references a missing or empty file: ${p}"
  fi
done < <(grep -oP '"/[^"]+"' "${CONFIG}" | tr -d '"')

[[ -s "${REF_DIR}/genome.fa.fai" ]] && ok "genome.fa.fai present" \
  || warnx "genome.fa.fai missing; mkref does not need it but later scripts do"

if [[ -f "${REF_DIR}/REFERENCE_MANIFEST.txt" ]]; then
  say ""
  say "--- reference manifest ---"
  sed 's/^/    /' "${REF_DIR}/REFERENCE_MANIFEST.txt" | tee -a "${LOG}"
fi

# -----------------------------------------------------------------------------
# Cell 5: pre-flight
# -----------------------------------------------------------------------------

hdr "Pre-flight"

if [[ -d "${EXPECTED_REF}" ]]; then
  die "a reference already exists at ${EXPECTED_REF}. Move or remove it first; mkref will not overwrite and a half-built reference is worse than none."
fi

avail_kb="$(df -Pk "${BUILD_DIR}" 2>/dev/null | awk 'NR==2{print $4}')"
if [[ -n "${avail_kb}" ]]; then
  avail_gb="$(awk -v k="${avail_kb}" 'BEGIN{printf "%.0f", k/1024/1024}')"
  say "  available on $(df -P "${BUILD_DIR}" | awk 'NR==2{print $6}'): ${avail_gb} GB"
  # A STAR index on a ~3 Gb genome runs ~30-40 GB; 100 GB is comfortable room.
  if [[ "${avail_gb}" -lt 100 ]]; then
    warnx "under 100 GB free. The STAR index on a 3 Gb genome runs roughly 30 to 40 GB and mkref needs working space on top."
  else
    ok "disk headroom fine"
  fi
fi

say ""
say "  mkref builds a STAR genome index, so this is the slow, memory-hungry"
say "  step: expect roughly 1 to 3 hours at ${NTHREADS} threads."

# -----------------------------------------------------------------------------
# Cell 6: build
# -----------------------------------------------------------------------------

hdr "Building"

cd "${BUILD_DIR}" || die "cannot cd to ${BUILD_DIR}"
say "  cwd: $(pwd)"
say "  mkref writes the reference here, named after the config's genome field."
say ""
say "  command:"
say "    cellranger-arc mkref --config=${CONFIG} --nthreads=${NTHREADS} --memgb=${MEMGB}"
say ""
say "  started: $(date '+%Y-%m-%d %H:%M:%S')"
say ""

START_TS=$(date +%s)
if cellranger-arc mkref \
      --config="${CONFIG}" \
      --nthreads="${NTHREADS}" \
      --memgb="${MEMGB}" 2>&1 | tee -a "${LOG}"; then
  BUILD_OK=1
else
  BUILD_OK=0
fi
END_TS=$(date +%s)
ELAPSED=$(( END_TS - START_TS ))

say ""
say "  finished: $(date '+%Y-%m-%d %H:%M:%S')"
say "  elapsed:  $(( ELAPSED / 3600 ))h $(( (ELAPSED % 3600) / 60 ))m"

if [[ ${BUILD_OK} -ne 1 ]]; then
  say ""
  say "  mkref FAILED. Common causes, in rough order of likelihood:"
  say ""
  say "   1. JASPAR format rejected. 10x's own references were built against"
  say "      JASPAR2020 and the format changed slightly afterwards. Drop the"
  say "      input_motifs line from the config and rebuild; ARC's own motif"
  say "      scoring is lost but the MEME work in Track B is unaffected."
  say ""
  say "   2. GTF problem, usually a transcript whose exons disagree with its"
  say "      bounds, or a gene_id that appears on two contigs. The error names"
  say "      the offending record."
  say ""
  say "   3. non_nuclear_contigs names a contig absent from the FASTA. Check"
  say "      that chrM survived the Phase 3 filter in SETUP_02."
  say ""
  say "   4. Out of memory. Raise --mem in the SLURM header and --memgb here."
  die "see the log above"
fi

# -----------------------------------------------------------------------------
# Cell 7: verify
# -----------------------------------------------------------------------------

hdr "Verification"

[[ -d "${EXPECTED_REF}" ]] || die "mkref reported success but ${EXPECTED_REF} does not exist"
ok "reference directory created"

say ""
say "--- layout ---"
ls -la "${EXPECTED_REF}" | sed 's/^/    /' | tee -a "${LOG}"

N_MISSING=0
for f in reference.json fasta/genome.fa fasta/genome.fa.fai genes/genes.gtf.gz star/Genome star/SA; do
  if [[ -e "${EXPECTED_REF}/${f}" ]]; then
    ok "$(printf '%-28s %8s' "${f}" "$(du -h "${EXPECTED_REF}/${f}" 2>/dev/null | cut -f1)")"
  else
    warnx "$(printf '%-28s absent' "${f}")"
    N_MISSING=$((N_MISSING+1))
  fi
done

if [[ -f "${EXPECTED_REF}/reference.json" ]]; then
  say ""
  say "--- reference.json ---"
  sed 's/^/    /' "${EXPECTED_REF}/reference.json" | tee -a "${LOG}"
fi

if [[ -f "${EXPECTED_REF}/regions/motifs.pfm" ]]; then
  n="$(grep -c '^>' "${EXPECTED_REF}/regions/motifs.pfm" 2>/dev/null || echo 0)"
  ok "motifs baked into the reference: ${n}"
else
  warnx "no regions/motifs.pfm; the reference was built without input_motifs. ARC's own motif scoring is unavailable, which does not affect the MEME work in Track B."
fi

{
  echo ""
  echo "ARC reference built"
  echo "  when:      $(date '+%Y-%m-%d %H:%M:%S')"
  echo "  host:      $(hostname)"
  echo "  job:       ${SLURM_JOB_ID:-interactive}"
  echo "  elapsed:   $(( ELAPSED / 3600 ))h $(( (ELAPSED % 3600) / 60 ))m"
  echo "  threads:   ${NTHREADS}"
  echo "  path:      ${EXPECTED_REF}"
  echo "  arc:       $(cellranger-arc --version 2>&1 | head -1)"
  echo "  verdict at build time: ${V}"
} >> "${REF_DIR}/REFERENCE_MANIFEST.txt"
ok "manifest appended"

# -----------------------------------------------------------------------------
# Cell 8: next
# -----------------------------------------------------------------------------

hdr "Next"

say "Reference: ${EXPECTED_REF}"
say ""
if [[ ${N_MISSING} -gt 0 ]]; then
  warnx "${N_MISSING} expected file(s) were absent. Read the listing above before using this reference."
  say ""
fi
say "  1. cellranger-arc testrun --id=tiny"
say "     Run it from a scratch directory if it has not been run yet. 10x's own"
say "     smoke test over the full pipeline. A broken install found 18 runs into"
say "     a queue is expensive."
say ""
say "  2. Regenerate the libraries.csv sheets."
say "     The drafts in 00_inventory/libraries_draft/ are all named Brain_*,"
say "     including the 12 colon ones, from the tissue-inference bug that is"
say "     fixed in multiome_common.py. 20 runs: 8 brain, 12 colon."
say ""
say "  3. Run ONE brain library end to end before queueing the rest."
say "     Brain is the smaller tissue with the cleaner contrast, so a single run"
say "     validates the reference, the sheet format, the resource request and the"
say "     output layout, and gives a measured number for disk and wall time"
say "     instead of an estimate."
say ""
say "  Track B needs nothing from here beyond genome.fa, so the motif analysis"
say "  can run in parallel with the 20 alignment jobs."
say ""
say "log: ${LOG}"

exit 0
