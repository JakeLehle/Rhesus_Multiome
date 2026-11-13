#!/usr/bin/env python3
# =============================================================================
# MULTIOME_02_duplicate_audit.py
#
# ARC_06's Phase 0 gate found two pairs of brain cell types carrying identical
# RNA statistics:
#
#     Interneuron      == Glutamatergic_FILENAME   (5,767 genes compared)
#     Oligodendrocytes == Microglia_1              (6,723 genes compared)
#
# Two independent muscat fits cannot produce that. This script finds out which
# of two possible causes it is, because they need completely different fixes.
#
#   CAUSE A, our side: MULTIOME_01 registered one source file under two
#   cell-type tokens. Fix is in MULTIOME_01_build_tidy.py.
#
#   CAUSE B, the vendor's side: two distinct files, each self-describing as a
#   different cell type, carry the same numbers. Fix is a question to Genohub,
#   and nothing we write can repair it.
#
# Reading peek_gpm_header() already argues for B: shape is "long" only when the
# header has exactly one DA_ column, so celltypes[0] cannot be ambiguous, and a
# dict assignment keyed on it cannot fan one path out to two tokens. But that is
# an argument, not evidence, so this measures it.
#
# HOW IT TELLS THEM APART
#
# For every gene-peak-matrix-peaks.txt it computes a content fingerprint from
# the (gene, logFC, p_val) triples, sorted, with the DA_ column excluded. The
# fingerprint ignores column naming and row order, so:
#
#   same path appearing twice          -> cause A
#   two paths, same fingerprint        -> cause B
#   two paths, different fingerprints  -> the duplication is downstream of here
#
# Read-only. Writes one CSV next to the tidy tables and touches nothing else.
#
#   python3 -I MULTIOME_02_duplicate_audit.py
#   python3 -I MULTIOME_02_duplicate_audit.py --root /some/other/tree
#
# Jake Lehle, 2026-11
# =============================================================================

import argparse
import hashlib
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiome_common import (            # noqa: E402
    peek_gpm_header, infer_tissue, norm_celltype,
)

ROOT_DEFAULT = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome"
OUT_DEFAULT = ("/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
               "/01_tidy/duplicate_audit.csv")

TARGET = "gene-peak-matrix-peaks.txt"


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


def fingerprint(path, da_col):
    """
    A content hash of the statistics only.

    Deliberately excludes the DA_ column, because that column is the ONLY thing
    that differs between two files shipped as different cell types with the same
    numbers. Including it would give every file a unique hash and hide exactly
    the thing we are looking for.

    Rows are sorted so row order cannot matter, and the numbers are formatted to
    17 significant digits so a float that round-trips through text still matches.
    """
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        return None, 0, "unreadable: %s" % exc

    gene_col = df.columns[0]
    stat_cols = [c for c in df.columns
                 if c in ("logFC", "logCPM", "F", "p_val", "p_adj.loc",
                          "p_adj.glb")]
    if not stat_cols:
        return None, len(df), "no recognised statistics columns"

    use = [gene_col] + stat_cols
    sub = df[use].copy()
    sub[gene_col] = sub[gene_col].astype(str)
    sub = sub.sort_values(by=use, kind="mergesort")

    h = hashlib.md5()
    for row in sub.itertuples(index=False):
        h.update(("\t".join("%.17g" % v if isinstance(v, float) else str(v)
                            for v in row) + "\n").encode())
    return h.hexdigest(), len(df), ",".join(stat_cols)


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=ROOT_DEFAULT)
    ap.add_argument("--out", default=OUT_DEFAULT)
    args = ap.parse_args()

    if not os.path.isdir(args.root):
        sys.exit("root not found: %s" % args.root)

    hdr("MULTIOME_02_duplicate_audit")
    say("root: %s" % args.root)

    rows = []
    seen_paths = {}
    dup_paths = []

    for dirpath, dirnames, filenames in os.walk(args.root):
        dirnames.sort()
        rel_dir = os.path.relpath(dirpath, args.root)
        parts = [] if rel_dir == "." else rel_dir.split(os.sep)
        for fn in sorted(filenames):
            if fn.lower() != TARGET:
                continue
            path = os.path.join(dirpath, fn)
            rel = os.path.relpath(path, args.root)

            real = os.path.realpath(path)
            if real in seen_paths:
                dup_paths.append((rel, seen_paths[real]))
            else:
                seen_paths[real] = rel

            shape, contrast, celltypes = peek_gpm_header(path)
            tissue = infer_tissue(parts + [fn])
            if tissue is None:
                tissue = ("Brain" if any(c in ("Astrocytes", "D1M_M", "OPC")
                                         for c in celltypes) else "Colon")
            da_col = celltypes[0] if celltypes else None
            fp, nrow, note = fingerprint(path, da_col)

            rows.append({
                "rel_path": rel,
                "dir_token": norm_celltype(parts[-1]) if parts else "",
                "tissue": tissue,
                "contrast": contrast,
                "shape": shape,
                "header_celltype": da_col,
                "n_celltypes_in_header": len(celltypes),
                "rows": nrow,
                "stat_cols": note,
                "file_md5": md5_file(path),
                "stats_fingerprint": fp,
            })

    if not rows:
        sys.exit("no %s found under %s" % (TARGET, args.root))

    df = pd.DataFrame(rows)

    subh("every gene-peak-matrix-peaks.txt found")
    say(df[["tissue", "shape", "header_celltype", "n_celltypes_in_header",
            "rows", "rel_path"]].to_string(index=False))

    # ---- cause A: the same file reached the walker twice -------------------
    subh("cause A: one path registered twice")
    if dup_paths:
        for a, b in dup_paths:
            say("  %s  is the same file as  %s" % (a, b))
        say("")
        say("  CAUSE A CONFIRMED. A symlink or a second copy resolves to one")
        say("  inode. Fix the discovery in MULTIOME_01_build_tidy.py.")
    else:
        say("  no path appears twice; every file is distinct on disk")

    # ---- the ambiguity that peek_gpm_header rules out ---------------------
    subh("long files whose header names more than one cell type")
    bad = df[(df["shape"] == "long") & (df["n_celltypes_in_header"] > 1)]
    if len(bad):
        say(bad[["rel_path", "n_celltypes_in_header"]].to_string(index=False))
        say("")
        say("  celltypes[0] is a coin flip for these. That IS a MULTIOME_01 bug.")
    else:
        say("  none. Every long file names exactly one cell type, so the key")
        say("  MULTIOME_01 builds from celltypes[0] is unambiguous.")

    # ---- cause B: distinct files, identical numbers ------------------------
    subh("cause B: distinct files carrying identical statistics")
    ok = df[df["stats_fingerprint"].notna()]
    groups = (ok.groupby(["tissue", "stats_fingerprint"])
                .agg(n_files=("rel_path", "size"),
                     celltypes=("header_celltype", lambda s: " == ".join(sorted(map(str, s)))),
                     paths=("rel_path", lambda s: " | ".join(sorted(s))),
                     rows=("rows", "first"))
                .reset_index())
    dups = groups[groups["n_files"] > 1].sort_values(["tissue", "celltypes"])

    if len(dups):
        for r in dups.itertuples(index=False):
            say("  %s  %s" % (r.tissue, r.celltypes))
            say("      %d rows, identical in every statistics column" % r.rows)
            for p in r.paths.split(" | "):
                say("      %s" % p)
            say("")
        say("  CAUSE B CONFIRMED. These are different files on disk, each")
        say("  self-describing as a different cell type, carrying the same")
        say("  numbers. Nothing in MULTIOME_01 can repair that: the data for")
        say("  one of each pair was never delivered.")
        say("")
        say("  Also compare file_md5 in the CSV. If the full-file hashes DIFFER")
        say("  while the fingerprints match, only the DA_ column header was")
        say("  changed, which is the signature of a copied file whose label was")
        say("  updated and whose contents were not.")
    else:
        say("  no two distinct files share a statistics fingerprint.")
        say("")
        say("  That means the duplication ARC_06 found is not present in these")
        say("  source files and enters somewhere between here and")
        say("  genes_rna.csv. Re-check Cell 7 of MULTIOME_01_build_tidy.py.")

    # ---- row counts, as a sanity cross-check ------------------------------
    subh("rows per cell type, by tissue")
    for tis in sorted(df["tissue"].dropna().unique()):
        say("  %s" % tis)
        t = df[df["tissue"] == tis].sort_values("rows", ascending=False)
        for r in t.itertuples(index=False):
            say("    %-28s %6d rows" % (str(r.header_celltype), r.rows))
    say("")
    say("  A duplicated pair has identical row counts. Two genuinely different")
    say("  cell types almost never do, because filterByExpr keeps a different")
    say("  gene universe in each.")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, index=False)
    hdr("Done")
    say("  %d file(s) audited" % len(df))
    say("  wrote %s" % args.out)
    say("")
    say("  If cause B is confirmed, the question for Genohub is specific:")
    say("  the per-cell-type RNA tables for the pairs listed above contain")
    say("  the same numbers under different cell-type labels. Ask which label")
    say("  is correct and for the missing table to be regenerated.")


if __name__ == "__main__":
    main()
