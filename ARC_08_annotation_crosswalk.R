#!/usr/bin/env Rscript
# =============================================================================
# ARC_08_annotation_crosswalk.R
#
# Recover gene symbols the Mmul8 release 97 GTF never assigned, by mapping
# Ensembl gene ID to symbol from a later Ensembl release.
#
# WHY THIS EXISTS
#
# ARC_07 measured the real problem. Of 23,296 distinct gene_name values on the
# ATAC annotation, 7,773 (33.4%) are bare ENSMMUG identifiers. They are not
# corrupt: ARC_04 does this on purpose, because Signac needs a non-empty
# gene_name and the GTF supplies none --
#
#     gtf$gene_name <- ifelse(is.na(gtf$gene_name) | gtf$gene_name == "",
#                             gtf$gene_id, gtf$gene_name)
#
# -- so a third of the annotation cannot be reached by symbol or by any alias.
# MRC2 is one casualty. It is absent from the ATAC annotation, from the RNA
# reference and from the vendor's genes_rna.csv, which are all the same
# annotation, so that is one fact rather than three. It is NOT absent from the
# macaque genome: NCBI annotates it as Gene ID 717576, and the human gene is
# ENSG00000011028.
#
# Ensembl gene IDs are stable across releases and assemblies. Release 97 had not
# yet assigned symbols to those 7,773 genes; later releases named many of them.
# So the coordinates we already have are fine and only the labels are missing,
# which makes this a join rather than a realignment.
#
# WHAT IT DOES AND DOES NOT TOUCH
#
# It rewrites gene_name on the annotation attached to the ATAC assay. It does
# not touch peaks, counts, fragments, cell calls or coordinates. Nothing is
# recomputed and no object is rebuilt. The patched annotation is written as its
# own .rds so the change is reversible by deleting one file, and the original
# object on disk is left exactly as it was unless --patch-object is passed.
#
# HOW TO GET THE CROSSWALK
#
# Ensembl BioMart, or the REST API, or the current GTF. Two attributes only:
#
#   Ensembl BioMart -> Ensembl Genes -> Macaca mulatta genes
#   Attributes: "Gene stable ID" and "Gene name"
#   Export as TSV
#
# Or from a shell with network access:
#
# Or from a shell with network access. Note that curl does NOT expand a '*' in
# a remote path, so the filename has to be discovered first; an earlier version
# of these instructions used a glob, fetched nothing, and left a zero-byte file:
#
#   cd /master/jlehle/WORKING/REFS/Mmul8
#   BASE=https://ftp.ensembl.org/pub/current_gtf/macaca_mulatta/
#   FN=$(curl -fsS "$BASE" | grep -oE 'Macaca_mulatta\.[A-Za-z0-9_.]+\.chr\.gtf\.gz' | sort -u | head -1)
#   curl -fSL -o "$FN" "${BASE}${FN}"
#   gzip -t "$FN" && echo "gzip OK"
#   zcat "$FN" | awk -F'\t' '$3=="gene"' \
#     | sed -nE 's/.*gene_id "([^"]+)".*gene_name "([^"]+)".*/\1\t\2/p' \
#     | sort -u > ensmmug_to_symbol.tsv
#   wc -l ensmmug_to_symbol.tsv
#
# Any two-column file works: Ensembl ID in one column, symbol in another, with
# or without a header. The script finds the columns itself.
#
# READ-ONLY BY DEFAULT.
#
#   source /master/jlehle/anaconda3/bin/activate multiome_R
#
#   # 1. measure the gap, no crosswalk needed
#   Rscript ARC_08_annotation_crosswalk.R Colon
#
#   # 2. apply a crosswalk and see what it recovers
#   Rscript ARC_08_annotation_crosswalk.R Colon \
#       --crosswalk=/master/jlehle/WORKING/REFS/Mmul8/ensmmug_to_symbol.tsv
#
#   # 3. same, and also rewrite the object in place (makes a .bak first)
#   Rscript ARC_08_annotation_crosswalk.R Colon \
#       --crosswalk=... --patch-object
#
# Jake Lehle, 2026-11
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(data.table)
})

ANALYSIS_ROOT <- "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
OBJ_ROOT      <- file.path(ANALYSIS_ROOT, "02_arc", "04_objects")
GTF           <- "/master/jlehle/WORKING/REFS/Mmul8/genes.gtf"

# The genes we actually care about resolving, so the report answers the
# question that prompted this rather than a generic count.
TARGET_GENES <- c(
  # colon panel
  "HMOX1","ISG15","MERTK","CD36","COLEC12","MSR1","MRC2","CD14",
  "CEMIP","PDGFRB","TNC","MMP1","COL5A3","PLAU","F13A1",
  "STAT1","CXCL9","DDX58","MX1","OAS2",
  "NLRP3","CASP9","ULK1","LACC1","SERPINB1",
  "IL10RA","HSD11B1","HIF1A","ITGA6","PKNOX2","FCAR","IGHM",
  # brain panel
  "NRP1","VSIG4","FOLR2","SIGLEC9","SIGLEC5","COX4I1","NDUFB4","ATP5MC3",
  "ATP6V0A1","LAMTOR2","RNASET2","RPL18","RPL29","RPS6","RPS9","EEF1A1",
  "UBA52","NEDD8","PSMB6","HSP90AB1","AIM2")

# -----------------------------------------------------------------------------
# plumbing
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: ARC_08_annotation_crosswalk.R <Brain|Colon> ",
       "[--crosswalk=FILE] [--patch-object]")
}
TISSUE <- args[1]
if (!TISSUE %in% c("Brain", "Colon")) stop("tissue must be Brain or Colon")

.rest  <- args[-1]
.xw    <- grep("^--crosswalk=", .rest, value = TRUE)
PATCH  <- "--patch-object" %in% .rest
.unk   <- setdiff(.rest[grepl("^--", .rest)],
                  c(.xw, "--patch-object"))
if (length(.unk)) stop("unrecognised option(s): ", paste(.unk, collapse = ", "))
if (length(.rest[!grepl("^--", .rest)])) {
  stop("unexpected bare argument(s). Options are --crosswalk=FILE and --patch-object")
}
XW_PATH <- if (length(.xw)) base::sub("^--crosswalk=", "", .xw[1]) else NA_character_

say  <- function(...) cat(paste0(..., "\n"))
hdr  <- function(x) { say(""); say(strrep("=", 76)); say(x); say(strrep("=", 76)) }
subh <- function(x) { say(""); say(paste0("--- ", x, " ---")) }
ok   <- function(x) say(paste0("  [ok]    ", x))
wrn  <- function(x) say(paste0("  [warn]  ", x))

for (fn in c("say","hdr","subh","ok","wrn")) {
  if (exists(fn, envir = baseenv(), inherits = FALSE)) {
    stop("helper '", fn, "' collides with base::", fn, ". Rename it.")
  }
}

OBJ_DIR <- file.path(OBJ_ROOT, TISSUE)
OBJ_RDS <- file.path(OBJ_DIR, paste0("seurat_", TISSUE, "_rebuilt.rds"))
ANN_OUT <- file.path(OBJ_DIR, paste0("annotation_", TISSUE, "_crosswalked.rds"))
REP_OUT <- file.path(OBJ_DIR, paste0("annotation_crosswalk_report_", TISSUE, ".csv"))

hdr(paste0("ARC_08_annotation_crosswalk : ", TISSUE))
say(paste0("object:    ", OBJ_RDS))
say(paste0("crosswalk: ", if (is.na(XW_PATH)) "(none supplied, measure-only run)" else XW_PATH))
say(paste0("patch:     ", if (PATCH) "YES, the object on disk will be rewritten" else "no, read-only"))
if (!file.exists(OBJ_RDS)) stop("object not found: ", OBJ_RDS)

# -----------------------------------------------------------------------------
# Phase 1: measure the gap
# -----------------------------------------------------------------------------

hdr("Phase 1: how much of the annotation has no symbol")

say("  reading the object")
obj <- readRDS(OBJ_RDS)
ann <- Annotation(obj[["ATAC"]])
if (!"gene_id" %in% names(mcols(ann))) {
  stop("the annotation has no gene_id column, so an ID crosswalk cannot be ",
       "applied. Rebuild the annotation in ARC_04 keeping gene_id.")
}

gname <- as.character(ann$gene_name)
gid   <- as.character(ann$gene_id)
is_ens <- grepl("^ENSMMUG", gname)

n_rec <- length(ann)
u_all <- unique(gname[!is.na(gname) & nzchar(gname)])
u_ens <- unique(gname[is_ens])

say(sprintf("  annotation records:            %s", format(n_rec, big.mark = ",")))
say(sprintf("  distinct gene_name values:     %s", format(length(u_all), big.mark = ",")))
say(sprintf("  unnamed (bare ENSMMUG):        %s  (%.1f%% of genes)",
            format(length(u_ens), big.mark = ","),
            100 * length(u_ens) / length(u_all)))
say(sprintf("  records carrying no symbol:    %s  (%.1f%% of records)",
            format(sum(is_ens), big.mark = ","), 100 * mean(is_ens)))
say("")
say("  Those genes are unreachable by symbol and by every alias, which is why")
say("  adding aliases never recovered MRC2. The gene_name is literally the")
say("  Ensembl ID, substituted by ARC_04 because the GTF supplied none.")

# Sanity check that gene_id is usable as a join key at all.
if (!any(grepl("^ENSMMUG", gid))) {
  wrn("gene_id does not look like Ensembl IDs; a crosswalk keyed on ENSMMUG will not join")
}

# -----------------------------------------------------------------------------
# Phase 2: where the target genes stand right now
# -----------------------------------------------------------------------------

hdr("Phase 2: the panel genes, before any crosswalk")

present_now <- TARGET_GENES[TARGET_GENES %in% gname]
missing_now <- setdiff(TARGET_GENES, gname)
say(sprintf("  panel genes resolved:   %d of %d", length(present_now), length(TARGET_GENES)))
if (length(missing_now)) {
  wrn(sprintf("unresolved: %s", paste(missing_now, collapse = ", ")))
} else {
  ok("every panel gene already resolves; a crosswalk would add nothing here")
}

if (is.na(XW_PATH)) {
  hdr("No crosswalk supplied")
  say("  This was a measure-only run and nothing was changed.")
  say("")
  say("  To recover the unnamed genes, build a two-column Ensembl ID to symbol")
  say("  table from a CURRENT Ensembl release and pass it with --crosswalk=.")
  say("  Gene IDs are stable across releases, so a later release's symbols")
  say("  attach to the Mmul8 coordinates you already have. Coordinates do not")
  say("  change; only labels do.")
  say("")
  say("    BioMart -> Ensembl Genes -> Macaca mulatta genes")
  say("    Attributes: Gene stable ID, Gene name   ->  export TSV")
  say("")
  say("  or, from a host with network access. curl does NOT expand a '*' in a")
  say("  remote path, so discover the filename first and verify each step:")
  say("")
  say("    cd /master/jlehle/WORKING/REFS/Mmul8")
  say("    BASE=https://ftp.ensembl.org/pub/current_gtf/macaca_mulatta/")
  say("    FN=$(curl -fsS \"$BASE\" \\")
  say("         | grep -oE 'Macaca_mulatta\\.[A-Za-z0-9_.]+\\.chr\\.gtf\\.gz' \\")
  say("         | sort -u | head -1)")
  say("    echo \"filename: ${FN:-NONE FOUND}\"")
  say("    curl -fSL -o \"$FN\" \"${BASE}${FN}\"")
  say("    gzip -t \"$FN\" && echo 'gzip OK'")
  say("    zcat \"$FN\" | awk -F'\\t' '$3==\"gene\"' \\")
  say("      | sed -nE 's/.*gene_id \"([^\"]+)\".*gene_name \"([^\"]+)\".*/\\1\\t\\2/p' \\")
  say("      | sort -u > ensmmug_to_symbol.tsv")
  say("    wc -l ensmmug_to_symbol.tsv")
  say("    awk -F'\\t' '$2==\"MRC2\"' ensmmug_to_symbol.tsv \\")
  say("      || echo 'MRC2 is not in the current Ensembl macaque annotation'")
  say("")
  say("  The -f flag matters: without it curl saves an error page or nothing at")
  say("  all and the redirect still creates the file, so the failure only shows")
  say("  up later as a parse error.")
  say("")
  say("  A caveat worth stating up front: if MRC2 is still absent after the")
  say("  crosswalk, then Ensembl does not annotate it in macaque at all and the")
  say("  gene is only in NCBI RefSeq (Gene ID 717576). That is a real finding")
  say("  for the collaborators, not something more aliasing can fix.")
  q(status = 0)
}

# -----------------------------------------------------------------------------
# Phase 3: load and validate the crosswalk
# -----------------------------------------------------------------------------

hdr("Phase 3: the crosswalk")

if (!file.exists(XW_PATH)) stop("crosswalk file not found: ", XW_PATH)

# A FAILED DOWNLOAD LOOKS LIKE A MALFORMED FILE, SO SAY WHICH IT IS.
#
# Added 2026-11-13. The first attempt at this reported "the crosswalk needs at
# least two columns; found 0", which is true and useless: the actual problem was
# a zero-byte file left behind by a curl that fetched nothing. A shell pipeline
# that ends in a redirect creates the output file before the pipeline runs, so a
# failed download leaves an empty file rather than no file, and every downstream
# error then describes the symptom instead of the cause.
.sz <- file.info(XW_PATH)$size
if (is.na(.sz) || .sz == 0) {
  stop("the crosswalk file is EMPTY (0 bytes): ", XW_PATH, "\n",
       "  Nothing is wrong with the file format; the download that was meant to\n",
       "  create it produced no output. A pipeline ending in '> file' makes the\n",
       "  file before it runs, so a failed fetch leaves it empty.\n",
       "  Re-fetch it, confirm with:\n",
       "      ls -l ", XW_PATH, "\n",
       "      wc -l ", XW_PATH, "\n",
       "  and only then rerun this script.")
}
.pk <- tryCatch(readLines(XW_PATH, n = 5, warn = FALSE),
                error = function(e) character(0))
if (!length(.pk) || all(!nzchar(.pk))) {
  stop("the crosswalk file has bytes but no readable lines: ", XW_PATH)
}
if (any(grepl("<!DOCTYPE|<html|^404|Not Found|Access denied", .pk, ignore.case = TRUE))) {
  stop("the crosswalk file contains an HTML or error page, not data: ", XW_PATH, "\n",
       "  first line: ", substr(.pk[1], 1, 120), "\n",
       "  The download returned a web page. Use curl -fSL so it fails loudly on\n",
       "  a bad URL instead of saving the error page.")
}

xw_raw <- fread(XW_PATH, header = "auto", sep = "auto", quote = "")
if (!is.data.frame(xw_raw) || ncol(xw_raw) < 2) {
  stop("the crosswalk needs at least two columns; found ",
       if (is.data.frame(xw_raw)) ncol(xw_raw) else 0, ".\n",
       "  First line of the file: ", substr(.pk[1], 1, 160), "\n",
       "  Expected two tab-separated columns: Ensembl gene ID, then symbol.")
}

# Find the ID column by content rather than by name, because BioMart exports,
# REST dumps and awk one-liners all head their columns differently.
id_col <- NA_integer_
for (j in seq_len(ncol(xw_raw))) {
  v <- as.character(xw_raw[[j]])
  if (mean(grepl("^ENSMMUG[0-9]+", v), na.rm = TRUE) > 0.8) { id_col <- j; break }
}
if (is.na(id_col)) {
  stop("no column in ", XW_PATH, " looks like Ensembl macaque gene IDs ",
       "(expected most values to match ^ENSMMUG[0-9]+)")
}
sym_candidates <- setdiff(seq_len(ncol(xw_raw)), id_col)
sym_col <- sym_candidates[which.max(vapply(sym_candidates, function(j) {
  v <- as.character(xw_raw[[j]])
  mean(nzchar(v) & !grepl("^ENSMMUG", v), na.rm = TRUE)
}, numeric(1)))]

xw <- data.table(gene_id = as.character(xw_raw[[id_col]]),
                 symbol  = as.character(xw_raw[[sym_col]]))
say(sprintf("  columns used: ID = column %d (%s), symbol = column %d (%s)",
            id_col, names(xw_raw)[id_col], sym_col, names(xw_raw)[sym_col]))

# Strip version suffixes on both sides. Ensembl IDs sometimes carry ".3".
xw[, gene_id := base::sub("\\..*$", "", gene_id)]
xw <- xw[!is.na(symbol) & nzchar(symbol) & !grepl("^ENSMMUG", symbol)]
xw <- unique(xw, by = c("gene_id", "symbol"))

# A gene ID mapping to two different symbols is ambiguous; keep neither
# silently. Report and drop, because guessing here mislabels a locus.
dup <- xw[, .N, by = gene_id][N > 1]
if (nrow(dup)) {
  wrn(sprintf("%d gene ID(s) map to more than one symbol; dropping those",
              nrow(dup)))
  say(sprintf("          e.g. %s", paste(head(dup$gene_id, 5), collapse = ", ")))
  xw <- xw[!gene_id %in% dup$gene_id]
}
ok(sprintf("%s usable ID-to-symbol pairs", format(nrow(xw), big.mark = ",")))

# -----------------------------------------------------------------------------
# Phase 4: apply it, only to records that have no symbol
# -----------------------------------------------------------------------------

hdr("Phase 4: apply")

say("  Only records whose gene_name is a bare ENSMMUG ID are rewritten. A gene")
say("  the GTF already named keeps its name, so the crosswalk cannot silently")
say("  rename anything the collaborators have already seen in a figure.")

gid_bare <- base::sub("\\..*$", "", gid)
m <- xw$symbol[match(gid_bare, xw$gene_id)]
fill <- is_ens & !is.na(m) & nzchar(m)

new_gname <- gname
new_gname[fill] <- m[fill]

n_genes_named <- length(unique(gid_bare[fill]))
say("")
say(sprintf("  records relabelled:        %s of %s unnamed",
            format(sum(fill), big.mark = ","), format(sum(is_ens), big.mark = ",")))
say(sprintf("  distinct genes recovered:  %s of %s unnamed genes (%.1f%%)",
            format(n_genes_named, big.mark = ","),
            format(length(u_ens), big.mark = ","),
            100 * n_genes_named / max(1L, length(u_ens))))

u_new <- unique(new_gname[!is.na(new_gname) & nzchar(new_gname)])
still_ens <- unique(new_gname[grepl("^ENSMMUG", new_gname)])
say(sprintf("  symbol coverage:           %.1f%% -> %.1f%%",
            100 * (1 - length(u_ens) / length(u_all)),
            100 * (1 - length(still_ens) / length(u_new))))

# Collisions: a recovered symbol that now equals a symbol already in use is a
# real risk, because gene_region() in ARC_05 spans every record sharing a
# gene_name and would merge two loci into one window.
recovered <- setdiff(unique(m[fill]), NA_character_)
collide <- intersect(recovered, u_all)
if (length(collide)) {
  wrn(sprintf("%d recovered symbol(s) already exist in the annotation:", length(collide)))
  say(sprintf("          %s", paste(head(collide, 20), collapse = ", ")))
  say("          ARC_05 spans every record sharing a gene_name, so a duplicate")
  say("          symbol would merge two loci into one window. ARC_05 rev 14's")
  say("          multi-contig guard catches the cross-contig case and warns;")
  say("          check any of these you intend to plot.")
}

subh("the panel genes, after the crosswalk")
present_after <- TARGET_GENES[TARGET_GENES %in% new_gname]
gained <- setdiff(present_after, present_now)
still_missing <- setdiff(TARGET_GENES, new_gname)
say(sprintf("  panel genes resolved: %d -> %d of %d",
            length(present_now), length(present_after), length(TARGET_GENES)))
if (length(gained)) {
  ok(sprintf("RECOVERED: %s", paste(gained, collapse = ", ")))
  for (g in gained) {
    gg <- ann[which(new_gname == g)]
    sq <- unique(as.character(seqnames(gg)))
    say(sprintf("    %-10s %s:%s-%s (%s bp, %d records, strand %s)", g,
                paste(sq, collapse = "/"),
                format(min(start(gg)), big.mark = ","),
                format(max(end(gg)), big.mark = ","),
                format(max(end(gg)) - min(start(gg)) + 1, big.mark = ","),
                length(gg),
                paste(unique(as.character(strand(gg))), collapse = "/")))
  }
}
if (length(still_missing)) {
  wrn(sprintf("STILL MISSING: %s", paste(still_missing, collapse = ", ")))
  say("          For these, Ensembl does not name the gene in macaque in any")
  say("          release. That is an annotation finding to report, not a bug.")
  say("          NCBI RefSeq may still annotate them (MRC2 is Gene ID 717576),")
  say("          in which case the only routes are a RefSeq-derived annotation")
  say("          or dropping the gene from the panel.")
}

# -----------------------------------------------------------------------------
# Phase 5: write
# -----------------------------------------------------------------------------

hdr("Phase 5: outputs")

rep <- data.table(
  gene = TARGET_GENES,
  resolved_before = TARGET_GENES %in% present_now,
  resolved_after  = TARGET_GENES %in% present_after,
  recovered_by_crosswalk = TARGET_GENES %in% gained)
fwrite(rep, REP_OUT)
ok(paste0("wrote ", REP_OUT))

ann_new <- ann
ann_new$gene_name <- new_gname
saveRDS(ann_new, ANN_OUT)
ok(paste0("wrote ", ANN_OUT))
say("")
say("  That .rds is the patched annotation on its own. Delete it and nothing")
say("  has changed. To use it without rewriting the object, add this to ARC_05")
say("  right after the object is loaded:")
say("")
say(sprintf('    ANN_PATCH <- "%s"', ANN_OUT))
say('    if (file.exists(ANN_PATCH)) {')
say('      Annotation(obj[["ATAC"]]) <- readRDS(ANN_PATCH)')
say('      ok("using the crosswalked annotation")')
say('    }')

if (PATCH) {
  bak <- paste0(OBJ_RDS, ".bak")
  if (!file.exists(bak)) {
    say("")
    say("  backing up the object before rewriting it (this is slow)")
    file.copy(OBJ_RDS, bak)
    ok(paste0("backup at ", bak))
  } else {
    wrn(paste0("a backup already exists and was left alone: ", bak))
  }
  Annotation(obj[["ATAC"]]) <- ann_new
  say("  rewriting the object")
  saveRDS(obj, OBJ_RDS)
  ok(paste0("patched ", OBJ_RDS))
  say("          Peaks, counts, fragments and cell calls are untouched; only")
  say("          the gene_name field on the annotation changed.")
} else {
  say("")
  say("  The object on disk was NOT modified. Pass --patch-object to rewrite it")
  say("  in place, which makes a .bak first.")
}

hdr("Done")
say("  Symbol coverage is the ceiling on every gene-name-driven request on this")
say("  project, including the motif work in Aim 3. It is worth fixing once here")
say("  rather than per gene, which is what this does.")
