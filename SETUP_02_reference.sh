#!/usr/bin/env bash
# =============================================================================
# SETUP_02_reference.sh
#
# Download and prepare the Mmul8 reference for cellranger-arc, plus JASPAR
# motifs. Stops deliberately BEFORE mkref so the TSS validation test can run
# first.
#
# WHY THIS RECIPE
#
#   The delivered peak coordinates are chr1..chr20, chrX, chrY. Twenty-two
#   contigs, no scaffolds, no chrM. Gene symbols include C1H1orf35 and
#   C1H1orf52, which is Ensembl/VGNC rhesus ortholog naming.
#
#   Ensembl gene symbols with UCSC-style contig names is not a combination
#   either source ships. It is what you get by following 10x's own reference
#   recipe, which takes an Ensembl FASTA and GTF and renames contigs to chr* to
#   match their human and mouse references. The absence of chrM peaks is
#   consistent with chrM declared in non_nuclear_contigs, which excludes it from
#   peak calling.
#
#   So: Ensembl Mmul_8.0.1 FASTA and GTF, contigs renamed, restricted to the 22
#   primary contigs plus chrM, GTF biotype-filtered.
#
#   This is a hypothesis, not a fact. SETUP_03_validate_gtf.py tests it against
#   the 268,483 (peak, gene, TSS distance) triples Genohub computed, which is the
#   only cheap way to find out before 20 alignment runs are queued.
#
# WHICH ENSEMBL RELEASE
#
#   Mmul_8.0.1 was superseded by Mmul_10 partway through the Ensembl series, and
#   I am not certain of the exact cutover. Rather than hardcode a guess, Phase 1
#   probes candidate releases newest-first and takes the highest one that
#   actually serves Mmul_8.0.1 files.
#
# A BIOTYPE DETAIL THAT MATTERS
#
#   Ensembl merged lincRNA and antisense into a single lncRNA biotype partway
#   through the same era. An older release will have the old names. The
#   allowlist below carries both, and Phase 5 reports which were actually found
#   so the filter can never silently drop everything.
#
# Jake Lehle, 2026-10
# =============================================================================

set -uo pipefail

# -----------------------------------------------------------------------------
# Cell 1: parameters
# -----------------------------------------------------------------------------

CONDA_BASE="/master/jlehle/anaconda3"
TARGET_ENV="sc_pre_R"

REF_DIR="/master/jlehle/WORKING/REFS/Mmul8"
DOWNLOAD_DIR="${REF_DIR}/downloads"

ANALYSIS_ROOT="/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
OUT_DIR="${ANALYSIS_ROOT}/setup_audit"
LOG="${OUT_DIR}/setup_02_reference.log"

ENSEMBL_SPECIES="macaca_mulatta"
ENSEMBL_ASSEMBLY="Mmul_8.0.1"
# Probed newest-first. The first release that serves Mmul_8.0.1 files wins.
ENSEMBL_CANDIDATE_RELEASES=(97 96 95 94 93 92 91 90)
ENSEMBL_BASE="https://ftp.ensembl.org/pub"

# Contigs to keep, in output order. chrM is kept so it can be declared
# non-nuclear in the mkref config, which is what excludes it from peak calling.
KEEP_CONTIGS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
              chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
              chrX chrY chrM)

# GTF gene_biotype allowlist. Mirrors what 10x retains, with the older Ensembl
# lncRNA spellings included.
KEEP_BIOTYPES=(
  protein_coding
  lncRNA lincRNA antisense
  IG_C_gene IG_D_gene IG_J_gene IG_LV_gene IG_V_gene
  IG_V_pseudogene IG_J_pseudogene IG_C_pseudogene
  TR_C_gene TR_D_gene TR_J_gene TR_V_gene
  TR_V_pseudogene TR_J_pseudogene
)

# JASPAR. The 2024 CORE vertebrates non-redundant set. jaspar format for the
# mkref input_motifs field, meme format for Track B.
JASPAR_MIRRORS=(
  "https://jaspar.elixir.no/download/data/2024/CORE"
  "https://jaspar.genereg.net/download/data/2024/CORE"
)
JASPAR_JASPAR_FMT="JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
JASPAR_MEME_FMT="JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"

# Re-download even if present.
FORCE_DOWNLOAD=0

# -----------------------------------------------------------------------------
# Cell 2: plumbing
# -----------------------------------------------------------------------------

mkdir -p "${OUT_DIR}" "${REF_DIR}" "${DOWNLOAD_DIR}" || { echo "FATAL: cannot create directories"; exit 1; }
: > "${LOG}"

say()  { echo "$*" | tee -a "${LOG}"; }
hdr()  { say ""; say "=============================================================================="; say "$*"; say "=============================================================================="; }
ok()   { say "  [ok]    $*"; }
warnx(){ say "  [warn]  $*"; }
fail() { say "  [FAIL]  $*"; }
die()  { say ""; say "FATAL: $*"; say "log: ${LOG}"; exit 1; }

hdr "SETUP_02_reference  —  $(date '+%Y-%m-%d %H:%M:%S')"
say "ref dir:  ${REF_DIR}"
say "log:      ${LOG}"

# shellcheck disable=SC1091
source "${CONDA_BASE}/bin/activate" "${TARGET_ENV}" 2>/dev/null \
  || die "could not activate ${TARGET_ENV}; run SETUP_01_build_env.sh first"
export PATH="${CONDA_PREFIX}/bin:${PATH}"
ok "env: ${CONDA_PREFIX}"

for b in curl wget samtools awk gzip sum; do
  command -v "${b}" >/dev/null 2>&1 || die "${b} not available"
done
ok "required tools present"

url_exists() { curl -sfI --max-time 30 "$1" >/dev/null 2>&1; }

fetch() {
  # fetch <url> <destination>
  local url="$1" dest="$2"
  if [[ -s "${dest}" && ${FORCE_DOWNLOAD} -eq 0 ]]; then
    ok "already present: $(basename "${dest}") ($(du -h "${dest}" | cut -f1))"
    return 0
  fi
  say "  downloading $(basename "${dest}")"
  if wget -q --show-progress -O "${dest}.partial" "${url}" 2>>"${LOG}"; then
    mv "${dest}.partial" "${dest}"
    ok "$(basename "${dest}") ($(du -h "${dest}" | cut -f1))"
    return 0
  fi
  rm -f "${dest}.partial"
  fail "download failed: ${url}"
  return 1
}

# -----------------------------------------------------------------------------
# Cell 3: find the right Ensembl release
# -----------------------------------------------------------------------------

hdr "Phase 1: locate ${ENSEMBL_ASSEMBLY} on Ensembl"

REL=""
FASTA_URL=""
FASTA_KIND=""

for r in "${ENSEMBL_CANDIDATE_RELEASES[@]}"; do
  say "  probing release ${r}"
  base_fa="${ENSEMBL_BASE}/release-${r}/fasta/${ENSEMBL_SPECIES}/dna"
  # Species without alt loci ship toplevel only, so try primary_assembly first
  # and fall back.
  for kind in primary_assembly toplevel; do
    cand="${base_fa}/Macaca_mulatta.${ENSEMBL_ASSEMBLY}.dna.${kind}.fa.gz"
    if url_exists "${cand}"; then
      REL="${r}"; FASTA_URL="${cand}"; FASTA_KIND="${kind}"
      ok "release ${r} serves ${ENSEMBL_ASSEMBLY} dna.${kind}"
      break 2
    fi
  done
done

[[ -n "${REL}" ]] || die "no candidate Ensembl release serves ${ENSEMBL_ASSEMBLY}. Check network access from titan, or extend ENSEMBL_CANDIDATE_RELEASES."

GTF_URL="${ENSEMBL_BASE}/release-${REL}/gtf/${ENSEMBL_SPECIES}/Macaca_mulatta.${ENSEMBL_ASSEMBLY}.${REL}.gtf.gz"
url_exists "${GTF_URL}" || die "FASTA found for release ${REL} but GTF missing at ${GTF_URL}"
ok "matching GTF found"

say ""
say "  release:      ${REL}"
say "  assembly:     ${ENSEMBL_ASSEMBLY}"
say "  FASTA kind:   dna.${FASTA_KIND}"
say ""
say "  Recording these in the manifest. The Genohub PDFs should be checked"
say "  against them: if Genohub used a different release, the GTF differs and"
say "  the TSS validation in SETUP_03 will show a systematic offset."

# -----------------------------------------------------------------------------
# Cell 4: download
# -----------------------------------------------------------------------------

hdr "Phase 2: download"

FASTA_GZ="${DOWNLOAD_DIR}/$(basename "${FASTA_URL}")"
GTF_GZ="${DOWNLOAD_DIR}/$(basename "${GTF_URL}")"

fetch "${FASTA_URL}" "${FASTA_GZ}" || die "FASTA download failed"
fetch "${GTF_URL}"   "${GTF_GZ}"   || die "GTF download failed"

# Ensembl CHECKSUMS use the old Unix `sum` format, not md5.
for pair in "${ENSEMBL_BASE}/release-${REL}/fasta/${ENSEMBL_SPECIES}/dna/CHECKSUMS|${DOWNLOAD_DIR}/CHECKSUMS.fasta|${FASTA_GZ}" \
            "${ENSEMBL_BASE}/release-${REL}/gtf/${ENSEMBL_SPECIES}/CHECKSUMS|${DOWNLOAD_DIR}/CHECKSUMS.gtf|${GTF_GZ}"; do
  IFS='|' read -r curl_url cks_file target <<< "${pair}"
  if fetch "${curl_url}" "${cks_file}"; then
    want="$(awk -v f="$(basename "${target}")" '$3==f {print $1" "$2}' "${cks_file}")"
    if [[ -n "${want}" ]]; then
      got="$(sum "${target}" | awk '{print $1" "$2}')"
      if [[ "${want}" == "${got}" ]]; then
        ok "checksum verified: $(basename "${target}")"
      else
        fail "CHECKSUM MISMATCH on $(basename "${target}"): want '${want}', got '${got}'"
        die "the download is corrupt or truncated; delete it and rerun"
      fi
    else
      warnx "no CHECKSUMS entry for $(basename "${target}"); falling back to gzip integrity"
    fi
  fi
done

say ""
say "--- gzip integrity, which catches truncation the checksum might not ---"
for f in "${FASTA_GZ}" "${GTF_GZ}"; do
  if gzip -t "${f}" 2>>"${LOG}"; then ok "$(basename "${f}") intact"
  else fail "$(basename "${f}") is TRUNCATED or corrupt"; die "delete it and rerun"; fi
done

# -----------------------------------------------------------------------------
# Cell 5: FASTA, rename and restrict
#
# Streamed in awk so a multi-GB genome never lands in memory. Ensembl names
# contigs 1..20, X, Y, MT; we emit chr1..chr20, chrX, chrY, chrM and drop
# everything else, which removes the scaffolds.
# -----------------------------------------------------------------------------

hdr "Phase 3: FASTA"

GENOME_FA="${REF_DIR}/genome.fa"
KEEP_RE="$(IFS='|'; echo "${KEEP_CONTIGS[*]}")"

if [[ -s "${GENOME_FA}" && ${FORCE_DOWNLOAD} -eq 0 ]]; then
  ok "genome.fa already present, skipping rebuild (set FORCE_DOWNLOAD=1 to redo)"
else
  say "  renaming contigs and dropping scaffolds. Streams, so this takes a few minutes."
  gzip -dc "${FASTA_GZ}" | awk -v keep="${KEEP_RE}" '
    BEGIN { n = split(keep, K, "|"); for (i = 1; i <= n; i++) want[K[i]] = 1 }
    /^>/ {
      # Ensembl header: ">1 dna:chromosome chromosome:Mmul_8.0.1:1:1:..."
      name = substr($1, 2)
      newname = (name == "MT") ? "chrM" : "chr" name
      emit = (newname in want)
      if (emit) print ">" newname
      next
    }
    emit { print }
  ' > "${GENOME_FA}.partial" 2>>"${LOG}"

  [[ -s "${GENOME_FA}.partial" ]] || die "FASTA rewrite produced nothing; check the Ensembl header format in ${FASTA_GZ}"
  mv "${GENOME_FA}.partial" "${GENOME_FA}"
  ok "genome.fa written ($(du -h "${GENOME_FA}" | cut -f1))"
fi

say ""
say "  indexing with samtools faidx"
samtools faidx "${GENOME_FA}" 2>>"${LOG}" || die "faidx failed"
ok "genome.fa.fai written"

say ""
say "--- contigs kept, with lengths ---"
awk '{printf "    %-8s %12d\n", $1, $2}' "${GENOME_FA}.fai" | tee -a "${LOG}"

N_CONTIG="$(wc -l < "${GENOME_FA}.fai")"
say ""
say "  contigs: ${N_CONTIG} (expected ${#KEEP_CONTIGS[@]})"
if [[ "${N_CONTIG}" -ne "${#KEEP_CONTIGS[@]}" ]]; then
  warnx "contig count differs from expectation. Any contig in KEEP_CONTIGS absent from the assembly is simply dropped, which is usually chrM or chrY. Check the listing above."
fi

for c in chr1 chrX chrY; do
  grep -qP "^${c}\t" "${GENOME_FA}.fai" && ok "${c} present" || warnx "${c} MISSING"
done

# -----------------------------------------------------------------------------
# Cell 6: GTF, rename, restrict, biotype filter
# -----------------------------------------------------------------------------

hdr "Phase 4: GTF"

GENES_GTF="${REF_DIR}/genes.gtf"
BIOTYPE_RE="$(IFS='|'; echo "${KEEP_BIOTYPES[*]}")"

say "  biotype allowlist: ${KEEP_BIOTYPES[*]}"
say ""
say "  Note on lncRNA: Ensembl merged lincRNA and antisense into lncRNA partway"
say "  through this era. All three spellings are in the allowlist, and the"
say "  survey below reports which this release actually uses."
say ""

say "--- gene_biotype values present in the source GTF ---"
gzip -dc "${GTF_GZ}" \
  | awk -F'\t' '$3=="gene"' \
  | grep -oP 'gene_biotype "[^"]+"' \
  | sed 's/gene_biotype "//; s/"//' \
  | sort | uniq -c | sort -rn \
  | head -30 | awk '{printf "    %8d  %s\n", $1, $2}' | tee -a "${LOG}"

if [[ -s "${GENES_GTF}" && ${FORCE_DOWNLOAD} -eq 0 ]]; then
  ok "genes.gtf already present, skipping rebuild"
else
  say ""
  say "  filtering and renaming"
  gzip -dc "${GTF_GZ}" | awk -F'\t' -v OFS='\t' \
      -v keep="${KEEP_RE}" -v bio="${BIOTYPE_RE}" '
    BEGIN {
      n = split(keep, K, "|"); for (i = 1; i <= n; i++) want[K[i]] = 1
      m = split(bio,  B, "|"); for (i = 1; i <= m; i++) wantbio[B[i]] = 1
    }
    /^#/ { print; next }
    {
      newname = ($1 == "MT") ? "chrM" : "chr" $1
      if (!(newname in want)) next
      if (match($9, /gene_biotype "[^"]+"/)) {
        bt = substr($9, RSTART + 14, RLENGTH - 15)
        if (!(bt in wantbio)) next
      }
      $1 = newname
      print
    }
  ' > "${GENES_GTF}.partial" 2>>"${LOG}"

  [[ -s "${GENES_GTF}.partial" ]] || die "GTF filter produced nothing; check the biotype survey above against KEEP_BIOTYPES"
  mv "${GENES_GTF}.partial" "${GENES_GTF}"
  ok "genes.gtf written ($(du -h "${GENES_GTF}" | cut -f1))"
fi

say ""
say "--- what survived ---"
N_GENE="$(awk -F'\t' '$3=="gene"' "${GENES_GTF}" | wc -l)"
N_TX="$(awk -F'\t' '$3=="transcript"' "${GENES_GTF}" | wc -l)"
N_EXON="$(awk -F'\t' '$3=="exon"' "${GENES_GTF}" | wc -l)"
say "    genes:       ${N_GENE}"
say "    transcripts: ${N_TX}"
say "    exons:       ${N_EXON}"

say ""
say "  The delivered wide gene-peak matrices carry 20,596 genes. A filtered GTF"
say "  gene count in that neighbourhood is a good sign; far from it suggests a"
say "  different biotype filter or a different release."
say ""
say "--- contigs in the filtered GTF ---"
awk -F'\t' '!/^#/{print $1}' "${GENES_GTF}" | sort -u | tr '\n' ' ' | fold -w 76 \
  | sed 's/^/    /' | tee -a "${LOG}"
say ""

if [[ "${N_GENE}" -lt 10000 ]]; then
  fail "only ${N_GENE} genes survived. That is too few. The biotype filter is probably wrong for this release; compare KEEP_BIOTYPES against the survey above."
fi

# -----------------------------------------------------------------------------
# Cell 7: JASPAR
# -----------------------------------------------------------------------------

hdr "Phase 5: JASPAR motifs"

JASPAR_FOR_ARC="${REF_DIR}/jaspar_core_vertebrates.pfm"
JASPAR_FOR_MEME="${REF_DIR}/jaspar_core_vertebrates.meme"

got_jaspar=0
got_meme=0
for mirror in "${JASPAR_MIRRORS[@]}"; do
  say "  trying ${mirror}"
  if [[ ${got_jaspar} -eq 0 ]] && url_exists "${mirror}/${JASPAR_JASPAR_FMT}"; then
    fetch "${mirror}/${JASPAR_JASPAR_FMT}" "${JASPAR_FOR_ARC}" && got_jaspar=1
  fi
  if [[ ${got_meme} -eq 0 ]] && url_exists "${mirror}/${JASPAR_MEME_FMT}"; then
    fetch "${mirror}/${JASPAR_MEME_FMT}" "${JASPAR_FOR_MEME}" && got_meme=1
  fi
  [[ ${got_jaspar} -eq 1 && ${got_meme} -eq 1 ]] && break
done

if [[ ${got_jaspar} -eq 1 ]]; then
  n="$(grep -c '^>' "${JASPAR_FOR_ARC}" 2>/dev/null || echo 0)"
  ok "JASPAR format: ${n} motifs, for the mkref input_motifs field"
else
  warnx "could not fetch the JASPAR-format file. mkref's input_motifs is optional, so SETUP_04 will proceed without it. Only ARC's own motif scoring is lost; the MEME work in Track B is unaffected."
fi

if [[ ${got_meme} -eq 1 ]]; then
  n="$(grep -c '^MOTIF' "${JASPAR_FOR_MEME}" 2>/dev/null || echo 0)"
  ok "MEME format: ${n} motifs, for Track B"
else
  warnx "could not fetch the MEME-format file. Convert from the JASPAR format instead:"
  say  "      jaspar2meme -pfm ${REF_DIR} > ${JASPAR_FOR_MEME}"
fi

say ""
say "  If mkref later rejects the JASPAR file, try the JASPAR2020 release"
say "  instead. 10x's own references were built against JASPAR2020 and the file"
say "  format changed slightly afterwards."

# -----------------------------------------------------------------------------
# Cell 8: mkref config, written but not run
# -----------------------------------------------------------------------------

hdr "Phase 6: mkref config"

CONFIG="${REF_DIR}/Mmul8_arc.config"
{
  echo "{"
  echo "    organism: \"Macaca mulatta\""
  echo "    genome: [\"Mmul8\"]"
  echo "    input_fasta: [\"${GENOME_FA}\"]"
  echo "    input_gtf: [\"${GENES_GTF}\"]"
  echo "    non_nuclear_contigs: [\"chrM\"]"
  [[ ${got_jaspar} -eq 1 ]] && echo "    input_motifs: \"${JASPAR_FOR_ARC}\""
  echo "}"
} > "${CONFIG}"

ok "written: ${CONFIG}"
say ""
sed 's/^/    /' "${CONFIG}" | tee -a "${LOG}"
say ""
say "  non_nuclear_contigs declares chrM, which excludes it from peak calling."
say "  That is what we think produced the absence of chrM peaks in the delivered"
say "  data, so it is part of reproducing the vendor's peak universe."

# -----------------------------------------------------------------------------
# Cell 9: manifest
# -----------------------------------------------------------------------------

hdr "Phase 7: manifest"

MANIFEST="${REF_DIR}/REFERENCE_MANIFEST.txt"
{
  echo "Mmul8 reference manifest"
  echo "built:            $(date '+%Y-%m-%d %H:%M:%S')"
  echo "host:             $(hostname)"
  echo "by:               SETUP_02_reference.sh"
  echo ""
  echo "Ensembl release:  ${REL}"
  echo "assembly:         ${ENSEMBL_ASSEMBLY}"
  echo "FASTA source:     ${FASTA_URL}"
  echo "GTF source:       ${GTF_URL}"
  echo "FASTA kind:       dna.${FASTA_KIND}"
  echo ""
  echo "contigs kept:     ${KEEP_CONTIGS[*]}"
  echo "contig count:     ${N_CONTIG}"
  echo "biotypes kept:    ${KEEP_BIOTYPES[*]}"
  echo ""
  echo "genes:            ${N_GENE}"
  echo "transcripts:      ${N_TX}"
  echo "exons:            ${N_EXON}"
  echo ""
  echo "genome.fa md5:    $(md5sum "${GENOME_FA}" | awk '{print $1}')"
  echo "genes.gtf md5:    $(md5sum "${GENES_GTF}" | awk '{print $1}')"
  if [[ ${got_jaspar} -eq 1 ]]; then
    echo "jaspar pfm md5:   $(md5sum "${JASPAR_FOR_ARC}" | awk '{print $1}')"
  fi
  echo ""
  echo "NOT YET VALIDATED. The recipe is inferred from the delivered peak"
  echo "coordinates and gene symbols, not read from Genohub documentation."
  echo "SETUP_03_validate_gtf.py must pass before mkref is run."
} > "${MANIFEST}"

ok "written: ${MANIFEST}"
say ""
sed 's/^/    /' "${MANIFEST}" | tee -a "${LOG}"

# -----------------------------------------------------------------------------
# Cell 10: stop here, on purpose
# -----------------------------------------------------------------------------

hdr "Stopping before mkref, deliberately"

say "The reference files are prepared but NOT yet built into an ARC reference."
say ""
say "mkref takes hours and a lot of RAM, and a wrong reference does not announce"
say "itself: every join still runs, every figure still renders, and the peaks"
say "point at the wrong places. The validation test is free and uses data we"
say "already have."
say ""
say "Next, in this order:"
say ""
say "  1. python SETUP_03_validate_gtf.py"
say ""
say "     Recomputes TSS distances for the 268,483 (peak, gene, distance)"
say "     triples Genohub produced, using this GTF, across several candidate"
say "     conventions. Reports which reproduces their numbers."
say ""
say "       near-exact match      reference reproduced, proceed"
say "       constant offset       same assembly, different GTF version"
say "       per-chromosome split  contig or coordinate problem"
say "       scattered mismatch    wrong assembly, stop"
say ""
say "  2. bash SETUP_04_mkref.sh"
say ""
say "     Refuses to run until SETUP_03 has written a passing report."
say ""
say "log: ${LOG}"

exit 0
