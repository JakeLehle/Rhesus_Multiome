# %%
# =============================================================================
# SETUP_03_validate_gtf.py
#
# The reference acceptance test. Read-only.
#
# WHAT IT DOES
#   Genohub computed 268,483 (peak, gene, signed TSS distance) triples against
#   their reference. Those triples are in peak_gene_map.csv. This script
#   recomputes the distance for every one of them from OUR GTF, under several
#   candidate conventions, and reports which convention reproduces their numbers.
#
# WHY IT MATTERS
#   A wrong reference does not announce itself. Every join still runs, every
#   figure still renders, and the peaks point at the wrong places. This is the
#   only cheap way to find out, and it uses data already in hand.
#
# WHAT THE OUTCOMES MEAN
#   near-exact match on one convention   reference and annotation reproduced
#   constant offset                      same assembly, different GTF version
#                                        or TSS definition, correctable
#   agreement on some chromosomes only   contig naming or coordinate system
#   scattered disagreement               wrong assembly, stop
#
# THE CONVENTION GRID
#   Two unknowns, tested as a grid rather than guessed:
#
#   TSS source      gene-level (one TSS per gene, by strand) or
#                   transcript-level (nearest transcript TSS)
#
#   Peak anchor     midpoint, nearest edge, or peak start
#
#   Six combinations. The winner tells us the reference is right AND fixes the
#   exact convention for all downstream region classification, which removes a
#   guess from every figure that follows.
#
# SIGN CONVENTION
#   The vendor note says: "If the peak is upstream, the distance is negative,
#   otherwise positive." Upstream is strand-aware, so for a minus-strand gene
#   upstream is to the right in genomic coordinates.
#
# WRITES
#   Only into OUT_DIR. A PASS/FAIL verdict file that SETUP_04_mkref.sh checks.
#
# Jake Lehle, 2026-10
# =============================================================================

# %%
# -----------------------------------------------------------------------------
# Cell 1: parameters
# -----------------------------------------------------------------------------

REF_DIR = "/master/jlehle/WORKING/REFS/Mmul8"
GTF_PATH = f"{REF_DIR}/genes.gtf"

ANALYSIS_ROOT = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
PEAK_GENE_MAP = f"{ANALYSIS_ROOT}/01_tidy/peak_gene_map.csv"
OUT_DIR = f"{ANALYSIS_ROOT}/setup_audit"

# A distance is "exact" if it matches the vendor's to within this many bp. Zero
# is the honest test; a small tolerance absorbs rounding in a midpoint
# definition (odd-width peaks divide unevenly).
EXACT_TOLERANCE_BP = 1

# Verdict thresholds on the best-performing convention.
PASS_EXACT_FRACTION = 0.95      # at or above this, PASS
INVESTIGATE_FRACTION = 0.50     # between the two, INVESTIGATE; below, FAIL

# A systematic offset is worth reporting separately from a scatter: if the
# modal non-zero difference accounts for this share of mismatches, it is an
# offset rather than noise.
OFFSET_MODE_FRACTION = 0.50

VERBOSE = True


# %%
# -----------------------------------------------------------------------------
# Cell 2: imports
# -----------------------------------------------------------------------------

import os
import re
import sys
import json
from collections import defaultdict

import numpy as np
import pandas as pd

os.makedirs(OUT_DIR, exist_ok=True)
NOTES = []


def say(msg=""):
    if VERBOSE:
        print(msg)


def note(msg):
    NOTES.append(msg)
    say("  [note] " + msg)


say("=" * 78)
say("SETUP_03_validate_gtf")
say("=" * 78)
say(f"GTF:           {GTF_PATH}")
say(f"vendor triples: {PEAK_GENE_MAP}")
say("")

for p in (GTF_PATH, PEAK_GENE_MAP):
    if not os.path.isfile(p):
        sys.exit(f"FATAL: not found: {p}")


# %%
# -----------------------------------------------------------------------------
# Cell 3: parse the GTF
#
# Gene-level and transcript-level TSS, both keyed on the same gene label the
# vendor used. Their GeneID column holds symbols where a gene has one and
# ENSMMUG IDs where it does not, which is the gene_name-else-gene_id convention.
# -----------------------------------------------------------------------------

say("-" * 78)
say("Parsing the GTF")
say("-" * 78)

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def parse_attrs(field):
    return dict(ATTR_RE.findall(field))


gene_rows, tx_rows = [], []
n_lines = 0

with open(GTF_PATH, "r", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        n_lines += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        feat = f[2]
        if feat not in ("gene", "transcript"):
            continue
        chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
        a = parse_attrs(f[8])
        gid = a.get("gene_id", "")
        label = a.get("gene_name") or gid          # the vendor's convention
        if not label:
            continue
        # GTF is 1-based inclusive. TSS is the 5' end.
        tss = start if strand == "+" else end
        rec = {"gene": label, "gene_id": gid, "chrom": chrom,
               "strand": strand, "tss": tss}
        (gene_rows if feat == "gene" else tx_rows).append(rec)

say(f"  GTF data lines read: {n_lines:,}")
say(f"  gene features:       {len(gene_rows):,}")
say(f"  transcript features: {len(tx_rows):,}")

genes = pd.DataFrame(gene_rows)
txs = pd.DataFrame(tx_rows)

if not len(genes):
    sys.exit("FATAL: no gene features parsed. Check the GTF format.")

dup = genes["gene"].duplicated().sum()
if dup:
    note(f"{dup:,} duplicate gene labels in the GTF (same symbol on several "
         f"loci). Keeping all; the match step takes the best-fitting locus per "
         f"peak, which is also what a closest-feature annotation would do.")

say("")
say("  contigs in the GTF:")
say("    " + " ".join(sorted(genes["chrom"].unique())))


# %%
# -----------------------------------------------------------------------------
# Cell 4: load the vendor triples
# -----------------------------------------------------------------------------

say("")
say("-" * 78)
say("Loading the vendor's triples")
say("-" * 78)

vend = pd.read_csv(PEAK_GENE_MAP)
need = {"peak_id", "gene", "distance_to_tss"}
missing = need - set(vend.columns)
if missing:
    sys.exit(f"FATAL: {PEAK_GENE_MAP} lacks columns {missing}")

vend = vend[vend["gene"].notna() & vend["distance_to_tss"].notna()].copy()
vend["distance_to_tss"] = vend["distance_to_tss"].astype(int)

# Deduplicate: the same (peak, gene, distance) appears once per cell type and
# contrast. The geometry is identical, so test each unique triple once.
trip = vend[["peak_id", "gene", "distance_to_tss"]].drop_duplicates()
say(f"  rows with a distance:  {len(vend):,}")
say(f"  unique triples:        {len(trip):,}")

coords = trip["peak_id"].str.rsplit("-", n=2, expand=True)
trip = trip.assign(
    chrom=coords[0],
    pstart=pd.to_numeric(coords[1], errors="coerce"),
    pend=pd.to_numeric(coords[2], errors="coerce"),
)
bad = trip["pstart"].isna() | trip["pend"].isna()
if bad.any():
    note(f"{bad.sum()} peak ids did not parse; dropped")
    trip = trip[~bad]
trip["pstart"] = trip["pstart"].astype(int)
trip["pend"] = trip["pend"].astype(int)

say("")
say("  vendor distance range: "
    f"{trip['distance_to_tss'].min():,} to {trip['distance_to_tss'].max():,}")
say(f"  peak width: min {int((trip['pend'] - trip['pstart']).min())}, "
    f"median {int((trip['pend'] - trip['pstart']).median())}, "
    f"max {int((trip['pend'] - trip['pstart']).max())}")

# Gene label coverage is the first real signal. A low rate here means the GTF
# annotation source is wrong, independent of any coordinate question.
gtf_labels = set(genes["gene"])
in_gtf = trip["gene"].isin(gtf_labels)
cover = in_gtf.mean()
say("")
say(f"  vendor genes found in our GTF: {in_gtf.sum():,} of {len(trip):,} "
    f"({100*cover:.1f}%)")

if cover < 0.80:
    note(f"only {100*cover:.1f}% of the vendor's gene labels appear in this "
         f"GTF. That points at a different annotation source or release, not a "
         f"coordinate problem. Examples we cannot find:")
    for g in trip.loc[~in_gtf, "gene"].drop_duplicates().head(15):
        say(f"      {g}")

trip = trip[in_gtf].copy()
if not len(trip):
    sys.exit("FATAL: no vendor genes found in this GTF. Wrong annotation source.")


# %%
# -----------------------------------------------------------------------------
# Cell 5: the convention grid
#
# For each (TSS source, peak anchor) pair, compute a signed distance and score
# it against the vendor's.
#
# Sign: negative upstream, strand-aware. For a plus-strand gene, a peak left of
# the TSS is upstream and negative. For a minus-strand gene, a peak right of
# the TSS is upstream and negative.
# -----------------------------------------------------------------------------

say("")
say("-" * 78)
say("Testing the convention grid")
say("-" * 78)

TSS_SOURCES = {
    "gene_level": genes,
    "transcript_level": txs if len(txs) else None,
}
ANCHORS = ["midpoint", "nearest_edge", "peak_start"]


def anchor_position(pstart, pend, tss, how):
    """Return the peak position to measure from, as an int array."""
    if how == "midpoint":
        return (pstart + pend) // 2
    if how == "peak_start":
        return pstart
    if how == "nearest_edge":
        # nearest edge of the peak to the TSS, 0 gap when overlapping
        return np.where(pend < tss, pend, np.where(pstart > tss, pstart, tss))
    raise ValueError(how)


def signed_distance(anchor, tss, strand):
    """Negative upstream, strand-aware."""
    raw = anchor - tss
    return np.where(strand == "+", raw, -raw)


results = []
detail = {}

for tss_name, tss_df in TSS_SOURCES.items():
    if tss_df is None or not len(tss_df):
        note(f"{tss_name}: no features of that kind in the GTF, skipped")
        continue

    # Candidate (gene, chrom, strand, tss) rows joined onto the triples. A gene
    # with several transcripts or loci yields several candidates; we keep the
    # one that best matches the vendor for each convention, which is what a
    # closest-feature annotation effectively does.
    cand = trip.merge(
        tss_df[["gene", "chrom", "strand", "tss"]],
        on="gene", how="left", suffixes=("", "_gtf"))

    # Guard against a gene label that lives on a different contig in our GTF.
    same_chrom = cand["chrom"] == cand["chrom_gtf"]
    n_cross = (~same_chrom).sum()
    if n_cross:
        note(f"{tss_name}: {n_cross:,} candidate pairings sit on a different "
             f"contig than the peak and were dropped")
    cand = cand[same_chrom].copy()
    if not len(cand):
        note(f"{tss_name}: nothing left after the contig check")
        continue

    for how in ANCHORS:
        a = anchor_position(cand["pstart"].to_numpy(),
                            cand["pend"].to_numpy(),
                            cand["tss"].to_numpy(), how)
        d = signed_distance(a, cand["tss"].to_numpy(),
                            cand["strand"].to_numpy())
        diff = d - cand["distance_to_tss"].to_numpy()

        tmp = pd.DataFrame({
            "peak_id": cand["peak_id"].to_numpy(),
            "gene": cand["gene"].to_numpy(),
            "chrom": cand["chrom"].to_numpy(),
            "vendor": cand["distance_to_tss"].to_numpy(),
            "ours": d,
            "diff": diff,
            "absdiff": np.abs(diff),
        })
        # best candidate per triple
        best = (tmp.sort_values("absdiff")
                   .drop_duplicates(subset=["peak_id", "gene"], keep="first"))

        n = len(best)
        exact = int((best["absdiff"] <= EXACT_TOLERANCE_BP).sum())
        key = f"{tss_name}|{how}"
        detail[key] = best

        results.append({
            "tss_source": tss_name,
            "peak_anchor": how,
            "n_tested": n,
            "n_exact": exact,
            "frac_exact": exact / n if n else np.nan,
            "median_abs_diff": float(best["absdiff"].median()) if n else np.nan,
            "mean_abs_diff": float(best["absdiff"].mean()) if n else np.nan,
            "p95_abs_diff": float(best["absdiff"].quantile(0.95)) if n else np.nan,
        })
        say(f"  {key:<34} exact {exact:>7,} / {n:>7,} "
            f"({100*exact/n if n else 0:5.1f}%)   median |diff| "
            f"{best['absdiff'].median() if n else float('nan'):.0f}")

res = pd.DataFrame(results).sort_values("frac_exact", ascending=False)
if not len(res):
    sys.exit("FATAL: no convention could be tested.")


# %%
# -----------------------------------------------------------------------------
# Cell 6: verdict
# -----------------------------------------------------------------------------

say("")
say("-" * 78)
say("Verdict")
say("-" * 78)

say("")
say(res.to_string(index=False))

best_row = res.iloc[0]
best_key = f"{best_row['tss_source']}|{best_row['peak_anchor']}"
best = detail[best_key]
frac = float(best_row["frac_exact"])

say("")
say(f"  best convention: {best_key}")
say(f"  exact agreement: {100*frac:.2f}%")

# Is the residual a systematic offset or scatter?
offset_bp = None
mm = best[best["absdiff"] > EXACT_TOLERANCE_BP]
if len(mm):
    modal = mm["diff"].value_counts()
    top_val, top_n = int(modal.index[0]), int(modal.iloc[0])
    share = top_n / len(mm)
    say("")
    say(f"  mismatches: {len(mm):,}")
    say(f"  modal difference: {top_val:+,} bp, in {100*share:.1f}% of them")
    if share >= OFFSET_MODE_FRACTION:
        offset_bp = top_val
        say(f"  that is a SYSTEMATIC OFFSET, not scatter")

    say("")
    say("  mismatch rate by contig, which separates a global problem from a")
    say("  contig-specific one:")
    by_chrom = (best.assign(bad=best["absdiff"] > EXACT_TOLERANCE_BP)
                    .groupby("chrom")
                    .agg(n=("bad", "size"), n_bad=("bad", "sum")))
    by_chrom["frac_bad"] = by_chrom["n_bad"] / by_chrom["n"]
    say(by_chrom.sort_values("frac_bad", ascending=False).to_string())

    worst = by_chrom["frac_bad"].max()
    bestc = by_chrom["frac_bad"].min()
    if worst - bestc > 0.25:
        note(f"mismatch rate varies a lot by contig ({100*bestc:.0f}% to "
             f"{100*worst:.0f}%). That points at a contig naming or coordinate "
             f"problem rather than a GTF version difference.")

if frac >= PASS_EXACT_FRACTION:
    verdict = "PASS"
    msg = (f"{100*frac:.2f}% of the vendor's TSS distances reproduce exactly "
           f"under {best_key}. The reference and annotation are reproduced. "
           f"Proceed to mkref.")
elif offset_bp is not None:
    verdict = "INVESTIGATE"
    msg = (f"only {100*frac:.2f}% match exactly, but the residual is a "
           f"systematic {offset_bp:+,} bp offset rather than scatter. That is "
           f"the signature of the same assembly with a different GTF version "
           f"or TSS definition. It is correctable, but decide how before "
           f"building: either accept the offset and document it, or try the "
           f"neighbouring Ensembl release.")
elif frac >= INVESTIGATE_FRACTION:
    verdict = "INVESTIGATE"
    msg = (f"{100*frac:.2f}% match exactly with no single dominant offset. "
           f"Partial agreement with scatter usually means a different "
           f"annotation release. Check the Genohub PDFs for the release they "
           f"used before building.")
else:
    verdict = "FAIL"
    msg = (f"only {100*frac:.2f}% match. This does not look like the reference "
           f"Genohub used. Do NOT run mkref or queue alignments. Check the "
           f"assembly version and contig naming first.")

say("")
say("=" * 78)
say(f"  {verdict}")
say("=" * 78)
say("")
for line in msg.split(". "):
    if line.strip():
        say(f"  {line.strip().rstrip('.')}.")

if verdict != "PASS":
    say("")
    say("  Worst mismatches, for inspection:")
    cols = ["peak_id", "gene", "chrom", "vendor", "ours", "diff"]
    say(best.sort_values("absdiff", ascending=False)[cols].head(20).to_string(index=False))


# %%
# -----------------------------------------------------------------------------
# Cell 7: write
# -----------------------------------------------------------------------------

say("")
say("-" * 78)
say(f"Writing to {OUT_DIR}")
say("-" * 78)


def dump(df, name):
    if df is None or not len(df):
        return
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=False, lineterminator="\n")
    say(f"  wrote {name}  ({len(df):,} rows)")


dump(res, "gtf_validation_conventions.csv")
dump(best.sort_values("absdiff", ascending=False).head(5000),
     "gtf_validation_worst_mismatches.csv")

verdict_path = os.path.join(OUT_DIR, "gtf_validation_verdict.json")
payload = {
    "verdict": verdict,
    "generated": pd.Timestamp.now().isoformat(),
    "gtf": GTF_PATH,
    "gtf_md5_note": "see REFERENCE_MANIFEST.txt in the reference directory",
    "n_triples_tested": int(best_row["n_tested"]),
    "best_tss_source": best_row["tss_source"],
    "best_peak_anchor": best_row["peak_anchor"],
    "frac_exact": frac,
    "systematic_offset_bp": offset_bp,
    "exact_tolerance_bp": EXACT_TOLERANCE_BP,
    "gene_label_coverage": float(cover),
    "message": msg,
    "notes": NOTES,
}
with open(verdict_path, "w") as fh:
    json.dump(payload, fh, indent=2, default=str)
say(f"  wrote gtf_validation_verdict.json  [{verdict}]")

say("")
if verdict == "PASS":
    say("  The winning convention is now a known quantity, not a guess:")
    say(f"    TSS source:  {best_row['tss_source']}")
    say(f"    peak anchor: {best_row['peak_anchor']}")
    say("  Use it for every region classification downstream so our numbers and")
    say("  the vendor's stay comparable.")
    say("")
    say("  Next: bash SETUP_04_mkref.sh")
else:
    say("  SETUP_04_mkref.sh will refuse to run on this verdict. Resolve the")
    say("  question above first.")

say("")
say("Done. Nothing outside the output directory was changed.")
