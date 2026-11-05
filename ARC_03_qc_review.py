# %%
# =============================================================================
# ARC_03_qc_review.py
#
# Read-only review of the 20 cellranger-arc count runs, before any object is
# built.
#
# THE FOUR QUESTIONS IT ANSWERS
#
#   1. Did every run actually finish, with every file the rest of the project
#      needs? A pipestance can report success and still be missing something.
#
#   2. How many cells did count keep versus how many the vendor annotated?
#      Colon_40878P gave 9,192 estimated cells against 5,023 in the loom, so
#      the vendor filtered roughly 45% out AFTER count. That ratio decides the
#      join: the rebuilt object must be subset to the cells they annotated, or
#      thousands of unannotated cells end up in the coverage plots. A library
#      whose retention is far from the others is worth knowing about before it
#      becomes a figure.
#
#   3. Are the ATAC metrics good enough to plot? TSS enrichment, fraction of
#      fragments in cells, fraction of transposition events in peaks, and the
#      non-nuclear fraction, against 10x's own thresholds.
#
#   4. Does a reduced union of the per-library peaks reproduce the vendor's
#      peak universe? They have 204,136 for brain and 345,553 for colon. If a
#      union lands near those, we have reproduced their peak set without
#      needing `cellranger-arc aggr` at all, and the whole aggr question goes
#      away. If it lands far off, aggr re-called peaks on the combined
#      fragments and we need to do the same.
#
# Nothing here is written outside OUT_DIR.
#
# Jake Lehle, 2026-10
# =============================================================================

# %%
# -----------------------------------------------------------------------------
# Cell 1: parameters
# -----------------------------------------------------------------------------

ANALYSIS_ROOT = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
COUNT_ROOT = f"{ANALYSIS_ROOT}/02_arc/02_count"
LOOM_META_DIR = f"{ANALYSIS_ROOT}/02_arc/00_loom_metadata"
OUT_DIR = f"{ANALYSIS_ROOT}/02_arc/03_qc"

# The vendor's per-tissue peak universes, recovered from the Bonferroni ratio
# in their FindMarkers output. The union computed below is compared to these.
VENDOR_PEAK_UNIVERSE = {"Brain": 204136, "Colon": 345553}

# 10x's published guidance for ARC. Used to flag, not to fail.
THRESHOLDS = {
    "ATAC TSS enrichment score": (5.0, "greater"),
    "ATAC Fraction of high-quality fragments in cells": (0.40, "greater"),
    "ATAC Fraction of transposition events in peaks in cells": (0.25, "greater"),
    "ATAC Non-nuclear read pairs": (0.20, "less"),
    "ATAC Valid barcodes": (0.85, "greater"),
    "GEX Valid barcodes": (0.85, "greater"),
    "GEX Reads mapped confidently to transcriptome": (0.30, "greater"),
    "GEX Median genes per cell": (300, "greater"),
}

# Flag a library whose cell retention deviates from its tissue median by more
# than this, in absolute proportion.
RETENTION_DEVIATION = 0.15

# Computing the reduced peak union reads 20 BED files. A couple of minutes.
DO_PEAK_UNION = True

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
WARNINGS = []


def say(m=""):
    if VERBOSE:
        print(m)


def warn(m):
    WARNINGS.append(m)
    say("  [warn] " + m)


say("=" * 78)
say("ARC_03_qc_review")
say("=" * 78)
say(f"count root: {COUNT_ROOT}")

if not os.path.isdir(COUNT_ROOT):
    sys.exit(f"FATAL: {COUNT_ROOT} not found")


# %%
# -----------------------------------------------------------------------------
# Cell 3: did every run finish, with everything we need
# -----------------------------------------------------------------------------

REQUIRED = [
    "summary.csv",
    "atac_fragments.tsv.gz",
    "atac_fragments.tsv.gz.tbi",
    "atac_peaks.bed",
    "filtered_feature_bc_matrix.h5",
    "filtered_feature_bc_matrix/barcodes.tsv.gz",
    "filtered_feature_bc_matrix/features.tsv.gz",
    "filtered_feature_bc_matrix/matrix.mtx.gz",
    "per_barcode_metrics.csv",
    "atac_possorted_bam.bam",
    "gex_possorted_bam.bam",
]

runs = sorted(d for d in os.listdir(COUNT_ROOT)
              if os.path.isdir(os.path.join(COUNT_ROOT, d))
              and not d.startswith("."))

say("")
say("-" * 78)
say(f"Completeness, {len(runs)} run directories")
say("-" * 78)

rows = []
for rid in runs:
    base = os.path.join(COUNT_ROOT, rid)
    outs = os.path.join(base, "outs")
    rec = {"run_id": rid,
           "tissue": "Brain" if rid.startswith("Brain") else "Colon",
           "sample": rid.split("_", 1)[1] if "_" in rid else rid}

    # _finalstate is written by martian when a pipestance completes.
    rec["finalstate"] = os.path.exists(os.path.join(base, "_finalstate"))
    rec["locked"] = os.path.exists(os.path.join(base, "_lock"))

    missing = [f for f in REQUIRED if not os.path.exists(os.path.join(outs, f))]
    rec["n_missing"] = len(missing)
    rec["missing"] = ";".join(missing)
    rec["outs_exists"] = os.path.isdir(outs)

    try:
        rec["frag_gb"] = os.path.getsize(
            os.path.join(outs, "atac_fragments.tsv.gz")) / 1e9
    except Exception:
        rec["frag_gb"] = np.nan

    rows.append(rec)

comp = pd.DataFrame(rows)
say("")
say(f"  {'run_id':<14} {'final':<6} {'lock':<5} {'missing':<8} {'fragments':>10}")
say("  " + "-" * 50)
for r in comp.itertuples():
    say(f"  {r.run_id:<14} {str(r.finalstate):<6} {str(r.locked):<5} "
        f"{r.n_missing:<8} {r.frag_gb:>9.1f}G")

for r in comp.itertuples():
    if r.locked:
        warn(f"{r.run_id}: _lock still present, so a job may be live or died "
             f"mid-run")
    if not r.finalstate:
        warn(f"{r.run_id}: no _finalstate, so the pipestance did not complete")
    if r.n_missing:
        warn(f"{r.run_id}: missing {r.n_missing} required file(s): {r.missing}")

n_ok = int(((comp["n_missing"] == 0) & comp["finalstate"] & ~comp["locked"]).sum())
say("")
say(f"  complete and intact: {n_ok} of {len(comp)}")
if n_ok == 20:
    say("  [ok] all 20 runs finished with every required file present")


# %%
# -----------------------------------------------------------------------------
# Cell 4: the metrics
# -----------------------------------------------------------------------------

say("")
say("-" * 78)
say("summary.csv metrics")
say("-" * 78)

mets = []
for rid in runs:
    p = os.path.join(COUNT_ROOT, rid, "outs", "summary.csv")
    if not os.path.isfile(p):
        continue
    try:
        d = pd.read_csv(p)
    except Exception as exc:
        warn(f"{rid}: could not read summary.csv ({exc})")
        continue
    if not len(d):
        warn(f"{rid}: summary.csv is empty")
        continue
    rec = d.iloc[0].to_dict()
    rec["run_id"] = rid
    rec["tissue"] = "Brain" if rid.startswith("Brain") else "Colon"
    mets.append(rec)

m = pd.DataFrame(mets)
if not len(m):
    sys.exit("FATAL: no summary.csv could be read")

for c in m.columns:
    if c not in ("run_id", "tissue", "Sample ID", "Genome", "Pipeline version"):
        m[c] = pd.to_numeric(m[c], errors="coerce")

say("")
say("  reference and version consistency:")
for c in ("Genome", "Pipeline version"):
    if c in m.columns:
        v = m[c].dropna().unique()
        if len(v) == 1:
            say(f"    [ok] {c}: {v[0]} across all runs")
        else:
            warn(f"{c} is not consistent across runs: {list(v)}. Every library "
                 f"must come from the same reference and pipeline version or "
                 f"the peak coordinates are not comparable.")

# Column lookup that tolerates a metric being absent or renamed in a future
# cellranger-arc version. An earlier revision read metrics positionally off
# itertuples (getattr(r, '_3')), which assumed summary.csv's own column order;
# that index landed on a string column and crashed on a numeric format spec.
# Nothing here depends on column position.
def col(df, name):
    """The column as a numeric Series, or an all-NaN Series if it is absent."""
    if name in df.columns:
        return pd.to_numeric(df[name], errors="coerce")
    return pd.Series(np.nan, index=df.index, dtype="float64")


METRIC_COLS = {
    "cells":     "Estimated number of cells",
    "peaks":     "ATAC Number of peaks",
    "TSS":       "ATAC TSS enrichment score",
    "frag/cell": "ATAC Median high-quality fragments per cell",
    "in_cells":  "ATAC Fraction of high-quality fragments in cells",
    "FRiP":      "ATAC Fraction of transposition events in peaks in cells",
    "nonnuc":    "ATAC Non-nuclear read pairs",
    "UMI":       "GEX Median UMI counts per cell",
    "genes":     "GEX Median genes per cell",
    "txome":     "GEX Reads mapped confidently to transcriptome",
}

for short, full in METRIC_COLS.items():
    if full not in m.columns:
        warn(f"summary.csv has no column '{full}'; it shows as '-' below. "
             f"10x renames metrics between versions, so check the header if "
             f"several are missing.")

disp = pd.DataFrame({"run_id": m["run_id"], "_tissue": m["tissue"]})
for short, full in METRIC_COLS.items():
    disp[short] = col(m, full)
disp = disp.sort_values(["_tissue", "run_id"]).drop(columns="_tissue")

# Nullable integers so a missing value prints as <NA> rather than 0 or nan.
for c in ("cells", "peaks", "frag/cell", "UMI", "genes"):
    disp[c] = disp[c].round(0).astype("Int64")
disp["TSS"] = disp["TSS"].round(2)
for c in ("in_cells", "FRiP", "nonnuc", "txome"):
    disp[c] = disp[c].round(3)

say("")
say("  " + disp.to_string(index=False, na_rep="-").replace("\n", "\n  "))
say("")
say("  in_cells = fraction of high-quality fragments in cells")
say("  FRiP     = fraction of transposition events in peaks, in cells")
say("  nonnuc   = non-nuclear read pairs, which is chrM")
say("  txome    = GEX reads mapped confidently to the transcriptome")

say("")
say("  against 10x's thresholds:")
for metric, (thr, direction) in THRESHOLDS.items():
    if metric not in m.columns:
        continue
    bad = m[m[metric] < thr] if direction == "greater" else m[m[metric] > thr]
    if len(bad):
        for r in bad.itertuples():
            warn(f"{r.run_id}: {metric} = {m.loc[r.Index, metric]:.4g}, "
                 f"{'below' if direction == 'greater' else 'above'} {thr}")
    else:
        say(f"    [ok] {metric}")

if "ATAC Non-nuclear read pairs" in m.columns:
    say("")
    say(f"  non-nuclear read pairs: median "
        f"{m['ATAC Non-nuclear read pairs'].median():.3f}")
    say("  That is chrM, which the mkref config declared in non_nuclear_contigs")
    say("  so it is excluded from peak calling. A figure in this range is")
    say("  expected and is part of reproducing the absence of chrM peaks in the")
    say("  delivered data.")


# %%
# -----------------------------------------------------------------------------
# Cell 5: cells kept by count versus cells the vendor annotated
#
# The number that decides how the annotation join has to work.
# -----------------------------------------------------------------------------

say("")
say("-" * 78)
say("Cell retention against the loom")
say("-" * 78)

loom_counts = {}
for tissue in ("Brain", "Colon"):
    p = os.path.join(LOOM_META_DIR, f"cell_metadata_{tissue}.csv")
    if not os.path.isfile(p):
        warn(f"{p} not found; run ARC_00_loom_metadata.py to enable this check")
        continue
    d = pd.read_csv(p, usecols=["Samples"])
    for s, n in d["Samples"].astype(str).value_counts().items():
        loom_counts[(tissue, s)] = int(n)

if loom_counts:
    cells_col = col(m, "Estimated number of cells")
    ret = []
    for idx in m.index:
        tissue = m.at[idx, "tissue"]
        run_id = m.at[idx, "run_id"]
        sample = run_id.split("_", 1)[1] if "_" in run_id else run_id
        n_count = cells_col.at[idx]
        n_loom = loom_counts.get((tissue, sample), np.nan)
        ok = (not pd.isna(n_count)) and n_count > 0 and (not pd.isna(n_loom))
        ret.append({
            "run_id": run_id, "tissue": tissue, "sample": sample,
            "cells_count": n_count, "cells_loom": n_loom,
            "retention": (n_loom / n_count) if ok else np.nan,
        })
    rt = pd.DataFrame(ret)

    rdisp = rt.sort_values(["tissue", "run_id"])[
        ["run_id", "cells_count", "cells_loom", "retention"]].copy()
    rdisp["cells_count"] = rdisp["cells_count"].round(0).astype("Int64")
    rdisp["cells_loom"] = rdisp["cells_loom"].round(0).astype("Int64")
    rdisp["retained_%"] = (100 * rdisp["retention"]).round(1)
    rdisp = rdisp.drop(columns="retention")

    say("")
    say("  " + rdisp.to_string(index=False, na_rep="-").replace("\n", "\n  "))

    say("")
    for tissue, g in rt.groupby("tissue"):
        med = g["retention"].median()
        say(f"  {tissue}: median retention {100*med:.1f}%  "
            f"({g['cells_loom'].sum():,.0f} of {g['cells_count'].sum():,.0f} cells)")
        off = g[(g["retention"] - med).abs() > RETENTION_DEVIATION]
        for r in off.itertuples():
            warn(f"{r.run_id}: retention {100*r.retention:.1f}% deviates from "
                 f"the {tissue} median of {100*med:.1f}% by more than "
                 f"{100*RETENTION_DEVIATION:.0f} points. Worth checking its "
                 f"web_summary before it contributes to a figure.")

    say("")
    say("  The vendor applied their own QC after count. The loom floors are")
    say("  roughly nUMI >= 1,000 and nGene >= 300, which accounts for most of")
    say("  the gap.")
    say("")
    say("  CONSEQUENCE FOR THE REBUILD: the object must be SUBSET to the cells")
    say("  present in the loom, with an inner join on CellID. An outer join")
    say("  would put thousands of unannotated cells into the coverage plots,")
    say("  and they would not be the cells behind the collaborators' figures.")

    rt.to_csv(os.path.join(OUT_DIR, "cell_retention.csv"),
              index=False, lineterminator="\n")
else:
    rt = pd.DataFrame()


# %%
# -----------------------------------------------------------------------------
# Cell 6: does a reduced peak union reproduce the vendor's universe
#
# If it does, aggr is unnecessary: a union of the per-library peaks IS their
# peak set, and the whole aggr-versus-Signac question is settled.
# -----------------------------------------------------------------------------

say("")
say("-" * 78)
say("Peak union against the vendor's universe")
say("-" * 78)

union_result = {}

if DO_PEAK_UNION:
    def read_bed(path):
        iv = defaultdict(list)
        try:
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    f = line.split("\t")
                    if len(f) < 3:
                        continue
                    try:
                        iv[f[0]].append((int(f[1]), int(f[2])))
                    except ValueError:
                        continue
        except Exception as exc:
            warn(f"could not read {path} ({exc})")
        return iv

    def merge_intervals(by_chrom):
        """Standard sweep merge. Returns the merged interval count."""
        total = 0
        merged = {}
        for c, lst in by_chrom.items():
            lst.sort()
            out = []
            cs, ce = lst[0]
            for s, e in lst[1:]:
                if s <= ce:
                    ce = max(ce, e)
                else:
                    out.append((cs, ce))
                    cs, ce = s, e
            out.append((cs, ce))
            merged[c] = out
            total += len(out)
        return total, merged

    for tissue in ("Brain", "Colon"):
        tissue_runs = [r for r in runs if r.startswith(tissue)]
        if not tissue_runs:
            continue
        say("")
        say(f"  {tissue}: reading {len(tissue_runs)} peak BEDs")

        pooled = defaultdict(list)
        per_lib = {}
        for rid in tissue_runs:
            bed = os.path.join(COUNT_ROOT, rid, "outs", "atac_peaks.bed")
            iv = read_bed(bed)
            n = sum(len(v) for v in iv.values())
            per_lib[rid] = n
            for c, lst in iv.items():
                pooled[c].extend(lst)

        say("")
        say("    per-library peak counts:")
        for rid, n in per_lib.items():
            say(f"      {rid:<14} {n:>8,}")

        if not pooled:
            warn(f"{tissue}: no peaks read")
            continue

        n_union, merged = merge_intervals(pooled)
        vendor = VENDOR_PEAK_UNIVERSE.get(tissue)
        say("")
        say(f"    reduced union:      {n_union:>8,}")
        say(f"    vendor's universe:  {vendor:>8,}")
        ratio = n_union / vendor if vendor else np.nan
        say(f"    ratio:              {ratio:>8.3f}")

        widths = [e - s for lst in merged.values() for s, e in lst]
        say(f"    merged peak width:  min {min(widths)}, "
            f"median {int(np.median(widths))}, max {max(widths)}")
        say(f"    contigs:            {len(merged)}")
        if "chrM" in merged:
            warn(f"{tissue}: the union contains chrM peaks, which the vendor's "
                 f"data does not. non_nuclear_contigs should have excluded it; "
                 f"drop chrM before building the consensus set.")

        union_result[tissue] = {
            "n_union": int(n_union), "vendor": vendor, "ratio": float(ratio),
            "per_library": per_lib, "n_contigs": len(merged),
        }

        say("")
        if 0.85 <= ratio <= 1.15:
            say("    [ok] within 15% of the vendor's universe. A reduced union")
            say("    of the per-library peaks reproduces their peak set closely")
            say("    enough, so `cellranger-arc aggr` is NOT needed. Build the")
            say("    consensus set this way and run FeatureMatrix against the")
            say("    per-library fragments.")
        elif ratio > 1.15:
            say("    The union is larger than the vendor's universe. aggr")
            say("    re-calls peaks on the pooled fragments, which merges")
            say("    library-specific calls and yields fewer, wider peaks. Two")
            say("    options: run aggr per tissue to match them, or proceed with")
            say("    the union and state in the methods that our peak set is")
            say("    broader. The union is the more permissive choice and loses")
            say("    nothing, since every vendor peak should still be inside it.")
        else:
            say("    The union is SMALLER than the vendor's universe, which is")
            say("    unexpected: a union cannot lose peaks. Check that all")
            say("    libraries were read and that the universe estimate from")
            say("    the Bonferroni ratio is right.")

    # The decisive test: are the vendor's actual peak coordinates inside our union?
    peaks_da = os.path.join(ANALYSIS_ROOT, "01_tidy", "peaks_da.csv")
    if os.path.isfile(peaks_da) and union_result:
        say("")
        say("  Containment: are the vendor's own DA peaks inside our union?")
        try:
            pk = pd.read_csv(peaks_da, usecols=["tissue", "peak_id", "chrom",
                                                "start", "end"])
            pk = pk.drop_duplicates(subset=["tissue", "peak_id"])
            for tissue, g in pk.groupby("tissue"):
                if tissue not in union_result:
                    continue
                tissue_runs = [r for r in runs if r.startswith(tissue)]
                pooled = defaultdict(list)
                for rid in tissue_runs:
                    for c, lst in read_bed(os.path.join(
                            COUNT_ROOT, rid, "outs", "atac_peaks.bed")).items():
                        pooled[c].extend(lst)
                _, merged = merge_intervals(pooled)

                n_in = 0
                for r in g.itertuples():
                    lst = merged.get(r.chrom)
                    if not lst:
                        continue
                    mid = (r.start + r.end) // 2
                    lo, hi = 0, len(lst) - 1
                    while lo <= hi:
                        k = (lo + hi) // 2
                        s, e = lst[k]
                        if mid < s:
                            hi = k - 1
                        elif mid > e:
                            lo = k + 1
                        else:
                            n_in += 1
                            break
                frac = n_in / len(g)
                say(f"    {tissue}: {n_in:,} of {len(g):,} vendor DA peak "
                    f"midpoints fall inside our union ({100*frac:.2f}%)")
                union_result[tissue]["vendor_peak_containment"] = float(frac)
                if frac >= 0.98:
                    say("    [ok] our peak set contains theirs, so every DA peak")
                    say("    they reported can be quantified and plotted.")
                else:
                    warn(f"{tissue}: only {100*frac:.1f}% of the vendor's DA "
                         f"peaks fall inside our union. Peaks outside it cannot "
                         f"be quantified against our consensus set, so check "
                         f"whether the missing ones cluster anywhere.")
        except Exception as exc:
            warn(f"containment test failed ({exc})")
else:
    say("  skipped (DO_PEAK_UNION=0)")


# %%
# -----------------------------------------------------------------------------
# Cell 7: write
# -----------------------------------------------------------------------------

say("")
say("=" * 78)
say(f"Writing to {OUT_DIR}")
say("=" * 78)

comp.to_csv(os.path.join(OUT_DIR, "run_completeness.csv"),
            index=False, lineterminator="\n")
say(f"  wrote run_completeness.csv  ({len(comp)} runs)")

m.to_csv(os.path.join(OUT_DIR, "arc_metrics.csv"), index=False, lineterminator="\n")
say(f"  wrote arc_metrics.csv  ({len(m)} runs, {m.shape[1]} metrics)")

if len(rt):
    say("  wrote cell_retention.csv")

if WARNINGS:
    pd.DataFrame({"warning": WARNINGS}).to_csv(
        os.path.join(OUT_DIR, "warnings.csv"), index=False, lineterminator="\n")
    say(f"  wrote warnings.csv  ({len(WARNINGS)})")

payload = {
    "generated": pd.Timestamp.now().isoformat(),
    "count_root": COUNT_ROOT,
    "n_runs": len(runs),
    "n_complete": n_ok,
    "genome": sorted(m["Genome"].dropna().unique().tolist()) if "Genome" in m else [],
    "pipeline": sorted(m["Pipeline version"].dropna().unique().tolist())
                if "Pipeline version" in m else [],
    "cells_count_total": int(m["Estimated number of cells"].sum()),
    "cells_loom_total": int(rt["cells_loom"].sum()) if len(rt) else None,
    "retention_median_by_tissue": (
        rt.groupby("tissue")["retention"].median().round(4).to_dict()
        if len(rt) else {}),
    "peak_union": union_result,
    "n_warnings": len(WARNINGS),
}
with open(os.path.join(OUT_DIR, "qc_summary.json"), "w") as fh:
    json.dump(payload, fh, indent=2, default=str)
say("  wrote qc_summary.json")

say("")
say(f"  {len(WARNINGS)} warning(s)")
for w in WARNINGS:
    say(f"    - {w}")

say("")
say("Next:")
say("")
say("  The two things this run settles for ARC_04:")
say("")
say("  1. The annotation join is an INNER join on CellID, subsetting to the")
say("     cells the vendor annotated. See the retention table above.")
say("")
say("  2. Whether the consensus peak set comes from a reduced union of the")
say("     per-library BEDs or needs aggr. See the ratio and containment")
say("     figures above.")
say("")
say("  Then: sbatch ARC_04_submit.sh")
say("")
say("Done. Nothing outside the output directory was changed.")
