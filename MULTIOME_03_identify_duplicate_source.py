#!/usr/bin/env python3
# =============================================================================
# MULTIOME_03_identify_duplicate_source.py
#
# MULTIOME_02 proved that two pairs of brain cell types carry identical RNA
# statistics, and that both files are real and correctly labelled in their own
# headers:
#
#   ATAC_results_Brain/ATAC-Results/GlutamatergicNeurons/ == .../Interneuron/
#   ATAC_results_Brain/ATAC-Results/Microglia-1/         == .../Oligodendrocytes/
#
# So one label in each pair is attached to the wrong numbers. This decides WHICH
# one, from the data rather than from a guess.
#
# WHY IT CAN BE DECIDED AT ALL
#
# The vendor's files carry logCPM, an expression level. Cell types have marker
# genes whose expression is not subtle: P2RY12 and CSF1R are high in microglia
# and near-absent in oligodendrocytes, while MBP and PLP1 are the reverse by an
# enormous margin. So the numbers themselves say which population they came
# from, whatever the directory is called.
#
# The sorted-adjacency argument predicts that the SECOND member of each pair is
# the corrupted one (Interneuron and Oligodendrocytes), because both pairs sit
# next to each other in a case-insensitive sort and in each case the second
# carries the first's numbers, which is what a loop that reuses the previous
# iteration's result after a failure looks like. This script tests that
# prediction instead of assuming it.
#
# Read-only.
#
#   python3 -I /master/jlehle/WORKING/SC/fastq/Rhesus_Multiome/MULTIOME_03_identify_duplicate_source.py
#
# Jake Lehle, 2026-11
# =============================================================================

import argparse
import os
import sys

import pandas as pd

AUDIT_CSV = ("/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
             "/01_tidy/duplicate_audit.csv")
ROOT = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome"
OUT_CSV = ("/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
           "/01_tidy/duplicate_source_identity.csv")

# Marker panels. Chosen to be unambiguous and well separated, not exhaustive.
# These are the symbols as they appear in the Mmul8 annotation the vendor used,
# which carries human-style names (NRP1, VSIG4 and SIGLEC9 all resolved there).
MARKERS = {
    "Microglia_1": ["P2RY12", "CX3CR1", "CSF1R", "C1QA", "C1QB", "TMEM119",
                    "AIF1", "ITGAM", "TYROBP", "FCER1G"],
    "Oligodendrocytes": ["MBP", "PLP1", "MOG", "MAG", "CNP", "MOBP", "CLDN11",
                         "ST18", "MYRF", "ERMN"],
    # Our earlier finding: the directory called GlutamatergicNeurons is in fact
    # the CCK+ GABAergic cluster, confirmed independently by ATAC_all_clusters
    # and by the loom. So it is separated from the TAC3+ interneuron cluster on
    # CCK versus TAC3, not on glutamatergic versus GABAergic markers.
    "Glutamatergic_FILENAME": ["CCK", "SLC32A1", "GAD1", "GAD2"],
    "Interneuron": ["TAC3", "NPY", "SST", "PVALB"],
}

# Cross-checks that should be high in every brain cell type, to confirm the
# marker lookup is working at all rather than silently matching nothing.
HOUSEKEEPING = ["ACTB", "GAPDH", "EEF1A1", "RPL18"]


def say(*a):
    print(*a, flush=True)


def hdr(x):
    say("")
    say("=" * 78)
    say(x)
    say("=" * 78)


def subh(x):
    say("")
    say("--- %s ---" % x)


def load_stats(path):
    df = pd.read_csv(path, sep="\t")
    gene_col = df.columns[0]
    if "logCPM" not in df.columns:
        return None, "no logCPM column"
    out = df[[gene_col, "logCPM"]].copy()
    out.columns = ["gene", "logCPM"]
    out["gene"] = out["gene"].astype(str).str.upper()
    return out, None


def panel_score(stats, genes):
    """Mean logCPM over the markers present, plus how many were found."""
    hit = stats[stats["gene"].isin([g.upper() for g in genes])]
    if not len(hit):
        return float("nan"), 0, []
    return (float(hit["logCPM"].mean()), len(hit),
            sorted(hit["gene"].tolist()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", default=AUDIT_CSV)
    ap.add_argument("--root", default=ROOT)
    ap.add_argument("--out", default=OUT_CSV)
    args = ap.parse_args()

    if not os.path.isfile(args.audit):
        sys.exit("audit CSV not found: %s\nRun MULTIOME_02_duplicate_audit.py first."
                 % args.audit)

    hdr("MULTIOME_03_identify_duplicate_source")
    say("audit: %s" % args.audit)
    say("root:  %s" % args.root)

    aud = pd.read_csv(args.audit)
    aud = aud[aud["stats_fingerprint"].notna()]
    groups = (aud.groupby(["tissue", "stats_fingerprint"])
                 .filter(lambda g: len(g) > 1)
                 .groupby(["tissue", "stats_fingerprint"]))

    if not len(groups):
        say("")
        say("  no duplicated fingerprints in the audit. Nothing to identify.")
        return

    rows = []
    for (tis, fp), g in groups:
        labels = sorted(g["header_celltype"].astype(str).tolist())
        subh("%s   %s" % (tis, " == ".join(labels)))

        path = os.path.join(args.root, g["rel_path"].iloc[0])
        stats, err = load_stats(path)
        if stats is None:
            say("  cannot score: %s" % err)
            continue
        say("  scoring: %s" % path)
        say("  %d genes in the file" % len(stats))

        hk, n_hk, _ = panel_score(stats, HOUSEKEEPING)
        say("  housekeeping control: mean logCPM %.2f over %d/%d markers"
            % (hk, n_hk, len(HOUSEKEEPING)))
        if n_hk == 0:
            say("  [FAIL] no housekeeping marker matched. Gene symbols in this")
            say("         file do not match the panels; identification aborted.")
            continue

        say("")
        say("  %-26s %10s %8s  %s" % ("candidate panel", "mean logCPM",
                                      "markers", "found"))
        scored = []
        for lab in labels:
            genes = MARKERS.get(lab)
            if not genes:
                say("  %-26s %10s %8s  (no marker panel defined)" % (lab, "-", "-"))
                continue
            m, n, found = panel_score(stats, genes)
            scored.append((lab, m, n))
            say("  %-26s %10.2f %5d/%-3d %s"
                % (lab, m, n, len(genes), ", ".join(found) if found else "none"))

        # Also score the OTHER pair's panels, as a negative control: a
        # microglial file should score low on oligodendrocyte markers AND low
        # on neuronal markers.
        others = [l for l in MARKERS if l not in labels]
        if others:
            say("")
            say("  negative controls (panels for cell types not in this pair)")
            for lab in others:
                m, n, _ = panel_score(stats, MARKERS[lab])
                say("  %-26s %10.2f %5d/%-3d" % (lab, m, n, len(MARKERS[lab])))

        valid = [s for s in scored if s[2] >= 3 and s[1] == s[1]]
        verdict = "undetermined"
        if len(valid) >= 2:
            valid.sort(key=lambda s: -s[1])
            top, second = valid[0], valid[1]
            gap = top[1] - second[1]
            say("")
            say("  %s scores %.2f against %s at %.2f, a gap of %.2f logCPM"
                % (top[0], top[1], second[0], second[1], gap))
            if gap >= 1.0:
                verdict = top[0]
                say("  VERDICT: these numbers are %s." % top[0])
                say("           The other label in the pair is attached to data")
                say("           that was never produced for it.")
            else:
                say("  VERDICT: too close to call on these markers. Widen the")
                say("           panels before relying on either label.")
        else:
            say("")
            say("  fewer than two panels had 3+ markers present; cannot compare.")

        rows.append({
            "tissue": tis,
            "pair": " == ".join(labels),
            "scored_file": path,
            "genes_in_file": len(stats),
            "housekeeping_logCPM": hk,
            "verdict_data_belongs_to": verdict,
            "corrupted_label": (" / ".join(l for l in labels if l != verdict)
                                if verdict != "undetermined" else ""),
        })

    if rows:
        out = pd.DataFrame(rows)
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        out.to_csv(args.out, index=False)
        hdr("Summary")
        say(out[["tissue", "pair", "verdict_data_belongs_to",
                 "corrupted_label"]].to_string(index=False))
        say("")
        say("  wrote %s" % args.out)
        say("")
        say("  The sorted-adjacency prediction was that the SECOND label in")
        say("  each pair is the corrupted one, so Interneuron and")
        say("  Oligodendrocytes. If the verdicts above say Glutamatergic and")
        say("  Microglia_1, the prediction held and every ARC_06 hit stands,")
        say("  because all of them were in Microglia-1.")
        say("")
        say("  If a verdict comes back the other way, the Microglia-1 hits are")
        say("  oligodendrocyte data and the whole brain selection has to be")
        say("  redone against the correct label.")


if __name__ == "__main__":
    main()
