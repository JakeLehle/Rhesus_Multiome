#!/usr/bin/env Rscript
# =============================================================================
# ARC_07_symbol_audit.R
#
# Why a panel gene produces no coverage figure, answered per gene and per layer.
#
# ARC_05 Phase 3 resolves a requested symbol against the gene_name field of the
# annotation attached to the ATAC assay. When that lookup misses, the gene is
# skipped and nothing is drawn. The miss has at least five distinct causes and
# they need different fixes, so this audit separates them instead of leaving
# "not in the annotation" as the whole answer:
#
#   1. UNNAMED IN THE GTF. ARC_04 does this, deliberately:
#
#        gtf$gene_name <- ifelse(is.na(gtf$gene_name) | gtf$gene_name == "",
#                                gtf$gene_id, gtf$gene_name)
#
#      Mmul8 release 97 leaves a large fraction of genes with no gene_name, so
#      those genes exist in the annotation ONLY as ENSMMUG00000xxxxx. No symbol
#      and no alias can ever match them. The fix is a coordinate override, not
#      another alias.
#
#   2. NAMED UNDER A DIFFERENT SYMBOL. An older or primate-specific symbol that
#      is not in GENE_ALIASES yet. The fix is one alias entry.
#
#   3. IN THE RNA REFERENCE BUT NOT THE ATAC ANNOTATION. The RNA assay comes
#      from the cellranger-arc reference; the ATAC annotation is parsed from
#      /master/jlehle/WORKING/REFS/Mmul8/genes.gtf. These are different files
#      and can disagree. The fix is a coordinate override.
#
#   4. GENUINELY ABSENT FROM THE BUILD. Some genes are pseudogenes or absent in
#      Old World monkeys (FCAR/CD89), and some loci are annotated too poorly to
#      use (the immunoglobulin locus, IGHM). That is a finding to report back to
#      the collaborators, not a bug to fix.
#
#   5. RESOLVED BUT UNDRAWN. The symbol was found and something later failed:
#      no region, no cells, a plotting error. ARC_05's own warnings cover this,
#      and this audit says so rather than guessing.
#
# For causes 1 and 3 the locus can usually still be recovered WITHOUT the GTF,
# because the vendor's own tables carry coordinates for the gene: genes_rna.csv
# has peak_string, and peak_gene_map.csv maps peaks to genes. This script emits
# a ready-to-paste GENE_REGION_OVERRIDE block from those coordinates for every
# gene it can recover, which is the actual deliverable here.
#
# READ-ONLY. Opens the object, the GTF and the tidy tables; writes one CSV and
# one R snippet into the 05_tracks output directory and nothing else.
#
#   source /master/jlehle/anaconda3/bin/activate multiome_R
#   Rscript /master/jlehle/WORKING/SC/fastq/Rhesus_Multiome/ARC_07_symbol_audit.R Colon
#   Rscript /master/jlehle/WORKING/SC/fastq/Rhesus_Multiome/ARC_07_symbol_audit.R Colon --only=MRC2
#
# Jake Lehle, 2026-11
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(data.table)
})

# -----------------------------------------------------------------------------
# Cell 1: paths and the panels. Kept identical to ARC_05 on purpose.
# -----------------------------------------------------------------------------

ANALYSIS_ROOT <- "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
OBJ_ROOT      <- file.path(ANALYSIS_ROOT, "02_arc", "04_objects")
OUT_ROOT      <- file.path(ANALYSIS_ROOT, "02_arc", "05_tracks")
RNA_CSV       <- file.path(ANALYSIS_ROOT, "01_tidy", "genes_rna.csv")
PEAK_MAP_CSV  <- file.path(ANALYSIS_ROOT, "01_tidy", "peak_gene_map.csv")
DA_PEAKS_CSV  <- file.path(ANALYSIS_ROOT, "01_tidy", "peaks_da.csv")
GTF           <- "/master/jlehle/WORKING/REFS/Mmul8/genes.gtf"

# Same gene lists and aliases as ARC_05 rev 14. If you change one, change both.
PANEL <- list(
  Brain = list(
    celltype_regex = "^Microglia-1$",
    celltype_label = "Microglia-1",
    vendor_token   = "Microglia_1",
    genes = c("NRP1","VSIG4","FOLR2","SIGLEC9","SIGLEC5",
              "COX4I1","NDUFB4","ATP5MC3",
              "ATP6V0A1","LAMTOR2","RNASET2",
              "RPL18","RPL29","RPS6","RPS9","EEF1A1","UBA52","NEDD8",
              "PSMB6","HSP90AB1","AIM2")),
  Colon = list(
    celltype_regex = "^Macrophage",
    celltype_label = "Macrophages",
    vendor_token   = "Macrophage",
    genes = c("HMOX1","ISG15",
              "MERTK","CD36","COLEC12","MSR1","MRC2","CD14",
              "CEMIP","PDGFRB","TNC","MMP1","COL5A3","PLAU","F13A1",
              "STAT1","CXCL9","DDX58","MX1","OAS2",
              "NLRP3","CASP9","ULK1","LACC1","SERPINB1",
              "IL10RA","HSD11B1","HIF1A",
              "ITGA6","PKNOX2",
              "FCAR","IGHM"))
)

GENE_ALIASES <- list(
  HMOX1    = c("HMOX1","HO-1","HSP32"),
  ISG15    = c("ISG15","G1P2","UCRP","IFI15"),
  SIGLEC9  = c("SIGLEC9","CD329"),
  SIGLEC5  = c("SIGLEC5","CD170","OBBP2"),
  ATP5MC3  = c("ATP5MC3","ATP5G3"),
  ATP6V0A1 = c("ATP6V0A1","ATP6V0A1","ATP6N1","VPP1"),
  NRP1     = c("NRP1","NRP","VEGF165R","CD304"),
  VSIG4    = c("VSIG4","CRIg","Z39IG"),
  FOLR2    = c("FOLR2","FBP"),
  DDX58    = c("DDX58","RIGI","RIG-I","RIGI1","RLR1"),
  CEMIP    = c("KIAA1199","CEMIP","TMEM2L"),
  MMP1     = c("MMP1","MMP1A","MMP1L1","CLG"),
  MERTK    = c("MERTK","MER","EYK","NYK"),
  CD36     = c("CD36","GP4","GP3B","SCARB3","FAT"),
  COLEC12  = c("COLEC12","CL-P1","SCARA4","NSR2"),
  MSR1     = c("MSR1","CD204","SCARA1","MSR"),
  MRC2     = c("MRC2","ENDO180","CD280","KIAA0709","UPARAP"),
  CD14     = c("CD14"),
  PDGFRB   = c("PDGFRB","PDGFR","PDGFR1","CD140B"),
  TNC      = c("TNC","HXB","TN"),
  COL5A3   = c("COL5A3"),
  PLAU     = c("PLAU","UPA","URK"),
  F13A1    = c("F13A1","F13A"),
  STAT1    = c("STAT1"),
  CXCL9    = c("CXCL9","MIG","SCYB9","CMK"),
  MX1      = c("MX1","MXA","IFI78"),
  OAS2     = c("OAS2"),
  NLRP3    = c("NLRP3","NALP3","CIAS1","PYPAF1"),
  CASP9    = c("CASP9","MCH6","APAF3"),
  ULK1     = c("ULK1","ATG1","ATG1A","UNC51"),
  LACC1    = c("LACC1","C13orf31","FAMIN"),
  SERPINB1 = c("SERPINB1","ELANH2","MNEI","PI2"),
  IL10RA   = c("IL10RA","IL10R","CD210","CDW210A"),
  HSD11B1  = c("HSD11B1","HSD11","HSD11L","HSD11B"),
  HIF1A    = c("HIF1A","HIF1","MOP1","PASD8","BHLHE78"),
  ITGA6    = c("ITGA6","CD49F","VLA6"),
  PKNOX2   = c("PKNOX2"),
  FCAR     = c("FCAR","CD89","FCAR1"),
  IGHM     = c("IGHM","IGHM1","MU","VH")
)

GENE_FAMILY_HINT <- list(
  MRC2 = "MRC", CEMIP = "KIAA11", DDX58 = "DDX5", MMP1 = "MMP",
  FCAR = "FC", IGHM = "IGH", COLEC12 = "COLEC", MSR1 = "MSR",
  MERTK = "MERT", PDGFRB = "PDGFR", COL5A3 = "COL5", HSD11B1 = "HSD11",
  SERPINB1 = "SERPINB", IL10RA = "IL10R", PKNOX2 = "PKNOX", ITGA6 = "ITGA",
  LACC1 = "LACC", OAS2 = "OAS", NLRP3 = "NLRP", ULK1 = "ULK", CASP9 = "CASP"
)

# -----------------------------------------------------------------------------
# Cell 2: arguments and plumbing
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: ARC_07_symbol_audit.R <Brain|Colon> [--only=GENE1,GENE2,...]")
}
TISSUE <- args[1]
if (!TISSUE %in% names(PANEL)) stop("tissue must be Brain or Colon")

PANEL_GENES <- PANEL[[TISSUE]]$genes
.only <- grep("^--only=", args[-1], value = TRUE)
GENES <- if (length(.only)) {
  g <- unlist(strsplit(base::sub("^--only=", "", .only), "[,[:space:]]+"))
  trimws(g[nzchar(g)])
} else PANEL_GENES

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

OBJ_RDS <- file.path(OBJ_ROOT, TISSUE, paste0("seurat_", TISSUE, "_rebuilt.rds"))
OUT_DIR <- file.path(OUT_ROOT, TISSUE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

hdr(paste0("ARC_07_symbol_audit : ", TISSUE))
say(paste0("object: ", OBJ_RDS))
say(paste0("gtf:    ", GTF))
say(paste0("genes:  ", length(GENES), " (", paste(GENES, collapse = ", "), ")"))
say(paste0("out:    ", OUT_DIR))

for (p in c(OBJ_RDS, GTF)) {
  if (!file.exists(p)) stop("required input not found: ", p)
}

# -----------------------------------------------------------------------------
# Cell 3: the three layers a symbol can live in
# -----------------------------------------------------------------------------

hdr("Phase 1: load the layers")

say("  reading the object (this is the slow part)")
obj <- readRDS(OBJ_RDS)
ann <- Annotation(obj[["ATAC"]])
ok(sprintf("%s annotation records, %s RNA features",
           format(length(ann), big.mark = ","),
           format(nrow(obj[["RNA"]]), big.mark = ",")))

ann_names <- unique(as.character(ann$gene_name))
ann_names <- ann_names[!is.na(ann_names) & nzchar(ann_names)]
rna_names <- rownames(obj[["RNA"]])

# HOW BAD IS THE UNNAMED PROBLEM. ARC_04 substitutes gene_id for a missing
# gene_name, so every unnamed gene shows up as ENSMMUG00000xxxxx. Counting them
# says how many of the panel's misses are structurally unfixable by aliasing.
n_ens <- sum(grepl("^ENSMMUG", ann_names))
n_loc <- sum(grepl("^LOC[0-9]+$", ann_names))
subh("how many annotation symbols are not real symbols")
say(sprintf("  distinct gene_name values:        %s", format(length(ann_names), big.mark = ",")))
say(sprintf("  bare Ensembl IDs (ENSMMUG...):    %s  (%.1f%%)",
            format(n_ens, big.mark = ","), 100 * n_ens / length(ann_names)))
say(sprintf("  LOC identifiers:                  %s  (%.1f%%)",
            format(n_loc, big.mark = ","), 100 * n_loc / length(ann_names)))
if (n_ens > 0) {
  say("")
  say("  Every ENSMMUG entry is a gene the GTF did not name, carried over by")
  say("  ARC_04's gene_name <- gene_id fallback. Those genes cannot be found")
  say("  by symbol or by any alias. A coordinate override is the only route.")
}

# The vendor tables, which are a second and independent source of coordinates.
rna_tab <- if (file.exists(RNA_CSV)) fread(RNA_CSV) else NULL
map_tab <- if (file.exists(PEAK_MAP_CSV)) fread(PEAK_MAP_CSV) else NULL
da_tab  <- if (file.exists(DA_PEAKS_CSV)) fread(DA_PEAKS_CSV) else NULL

subh("vendor tables")
for (nm in c("genes_rna.csv", "peak_gene_map.csv", "peaks_da.csv")) {
  p <- switch(nm, "genes_rna.csv" = RNA_CSV,
                  "peak_gene_map.csv" = PEAK_MAP_CSV,
                  "peaks_da.csv" = DA_PEAKS_CSV)
  t <- switch(nm, "genes_rna.csv" = rna_tab,
                  "peak_gene_map.csv" = map_tab,
                  "peaks_da.csv" = da_tab)
  if (is.null(t)) { wrn(sprintf("%-18s NOT FOUND at %s", nm, p)) } else {
    ok(sprintf("%-18s %s rows", nm, format(nrow(t), big.mark = ",")))
    say(sprintf("            columns: %s", paste(names(t), collapse = ", ")))
  }
}

# -----------------------------------------------------------------------------
# Cell 4: does the raw GTF mention the symbol at all
#
# The annotation in the object has been through rtracklayer and ARC_04's
# gene_name fallback. Grepping the GTF text asks the upstream question: is the
# string in the file at all, perhaps on a feature type or attribute that did not
# survive into the object.
# -----------------------------------------------------------------------------

gtf_hits <- function(sym) {
  out <- tryCatch(
    system2("grep", c("-m", "5", "-w", shQuote(sym), shQuote(GTF)),
            stdout = TRUE, stderr = FALSE),
    warning = function(w) character(0), error = function(e) character(0))
  out[nzchar(out)]
}

# -----------------------------------------------------------------------------
# Cell 5: the per-gene audit
# -----------------------------------------------------------------------------

hdr("Phase 2: per-gene audit")

resolve_try <- function(g) unique(c(g, GENE_ALIASES[[g]]))

rows <- list()
override_lines <- character(0)

for (g in GENES) {
  subh(g)
  tries <- resolve_try(g)
  say(sprintf("  aliases tried: %s", paste(tries, collapse = ", ")))

  # --- layer 1: the ATAC annotation, exactly as ARC_05 does it --------------
  hit_exact <- tries[tries %in% ann_names]
  hit_ci <- character(0)
  if (!length(hit_exact)) {
    for (t in tries) {
      i <- which(toupper(ann_names) == toupper(t))
      if (length(i)) { hit_ci <- ann_names[i[1]]; break }
    }
  }
  ann_symbol <- if (length(hit_exact)) hit_exact[1]
                else if (length(hit_ci)) hit_ci[1] else NA_character_

  # --- layer 2: the RNA assay ----------------------------------------------
  rna_hit <- tries[tries %in% rna_names]
  in_rna <- length(rna_hit) > 0

  # --- layer 3: the vendor RNA table ---------------------------------------
  vendor_hit <- NA_character_
  if (!is.null(rna_tab) && "gene" %in% names(rna_tab)) {
    vh <- tries[tries %in% unique(rna_tab$gene)]
    if (length(vh)) vendor_hit <- vh[1]
  }

  # --- layer 4: the raw GTF text -------------------------------------------
  gtf_found <- character(0)
  for (t in tries) {
    h <- gtf_hits(t)
    if (length(h)) { gtf_found <- h; break }
  }

  # --- report ---------------------------------------------------------------
  if (!is.na(ann_symbol)) {
    ok(sprintf("ATAC annotation: FOUND as '%s'", ann_symbol))
    gg <- ann[!is.na(ann$gene_name) & ann$gene_name == ann_symbol]
    sq <- unique(as.character(seqnames(gg)))
    say(sprintf("    %d records on %s, span %s-%s (%s bp), strand %s",
                length(gg), paste(sq, collapse = "/"),
                format(min(start(gg)), big.mark = ","),
                format(max(end(gg)), big.mark = ","),
                format(max(end(gg)) - min(start(gg)) + 1, big.mark = ","),
                paste(unique(as.character(strand(gg))), collapse = "/")))
    if (length(sq) > 1) {
      wrn("on more than one contig; ARC_05 rev 14 keeps the best-supported one")
    }
  } else {
    wrn("ATAC annotation: NOT FOUND under any alias")
    hint <- GENE_FAMILY_HINT[[g]]
    if (!is.null(hint)) {
      cand <- sort(grep(paste0("^", hint), ann_names, value = TRUE, ignore.case = TRUE))
      say(sprintf("    symbols starting '%s': %s", hint,
                  if (length(cand)) paste(head(cand, 30), collapse = ", ") else "NONE"))
    }
  }

  say(sprintf("  RNA assay:       %s", if (in_rna) paste0("present as '", rna_hit[1], "'") else "absent"))
  say(sprintf("  vendor RNA table:%s",
              if (!is.na(vendor_hit)) paste0(" present as '", vendor_hit, "'") else " absent"))
  say(sprintf("  raw GTF text:    %s",
              if (length(gtf_found)) sprintf("%d matching line(s)", length(gtf_found)) else "no match"))
  if (length(gtf_found)) {
    say(sprintf("    %s", substr(gtf_found[1], 1, 150)))
  }

  # --- coordinate recovery for the unfixable-by-alias cases ----------------
  rec_chr <- NA_character_; rec_start <- NA_integer_; rec_end <- NA_integer_
  rec_src <- NA_character_
  if (is.na(ann_symbol)) {
    sym_for_map <- if (!is.na(vendor_hit)) vendor_hit else g
    cand_peaks <- NULL
    if (!is.null(map_tab)) {
      gcol <- intersect(c("gene","gene_name","symbol"), names(map_tab))
      if (length(gcol)) {
        cand_peaks <- map_tab[get(gcol[1]) == sym_for_map]
        if (nrow(cand_peaks)) rec_src <- "peak_gene_map.csv"
      }
    }
    if ((is.null(cand_peaks) || !nrow(cand_peaks)) &&
        !is.null(rna_tab) && "peak_string" %in% names(rna_tab)) {
      pr <- rna_tab[gene == sym_for_map & !is.na(peak_string) & nzchar(peak_string)]
      if (nrow(pr)) {
        ps <- unlist(strsplit(paste(unique(pr$peak_string), collapse = ";"), "[;,[:space:]]+"))
        ps <- ps[grepl("^chr[0-9XYM]+[-:_][0-9]+[-:_][0-9]+$", ps)]
        if (length(ps)) {
          parts <- tstrsplit(ps, "[-:_]")
          cand_peaks <- data.table(chrom = parts[[1]],
                                   start = as.numeric(parts[[2]]),
                                   end   = as.numeric(parts[[3]]))
          rec_src <- "genes_rna.csv peak_string"
        }
      }
    }
    if (!is.null(cand_peaks) && nrow(cand_peaks)) {
      cc <- intersect(c("chrom","chr","seqnames"), names(cand_peaks))
      sc <- intersect(c("start","peak_start"), names(cand_peaks))
      ec <- intersect(c("end","peak_end"), names(cand_peaks))
      if (length(cc) && length(sc) && length(ec)) {
        tabc <- sort(table(as.character(cand_peaks[[cc[1]]])), decreasing = TRUE)
        rec_chr <- names(tabc)[1]
        sel <- cand_peaks[as.character(get(cc[1])) == rec_chr]
        rec_start <- as.integer(min(as.numeric(sel[[sc[1]]])))
        rec_end   <- as.integer(max(as.numeric(sel[[ec[1]]])))
        ok(sprintf("RECOVERED from %s: %s:%s-%s (%s bp, %d peaks)",
                   rec_src, rec_chr,
                   format(rec_start, big.mark = ","),
                   format(rec_end, big.mark = ","),
                   format(rec_end - rec_start + 1, big.mark = ","),
                   nrow(sel)))
        override_lines <- c(override_lines, sprintf(
          '  %-9s = list(chrom = "%s", start = %d, end = %d),   # from %s',
          g, rec_chr, rec_start, rec_end, rec_src))
      }
    } else {
      wrn("no coordinates available from the vendor tables either")
    }
  }

  cause <- if (!is.na(ann_symbol)) {
    "RESOLVED (look at ARC_05's own warnings for this gene)"
  } else if (!is.na(rec_chr)) {
    "UNNAMED OR ABSENT IN GTF, COORDINATES RECOVERABLE -> use an override"
  } else if (in_rna || !is.na(vendor_hit)) {
    "IN RNA BUT NOT IN THE ATAC GTF, no peaks to localise it -> ask the vendor"
  } else {
    "ABSENT EVERYWHERE -> likely not in this build; report to collaborators"
  }
  say(sprintf("  VERDICT: %s", cause))

  rows[[g]] <- data.table(
    gene = g, aliases_tried = paste(tries, collapse = "|"),
    ann_symbol = ann_symbol, in_atac_annotation = !is.na(ann_symbol),
    in_rna_assay = in_rna, in_vendor_rna_table = !is.na(vendor_hit),
    in_raw_gtf = length(gtf_found) > 0,
    rec_chrom = rec_chr, rec_start = rec_start, rec_end = rec_end,
    rec_source = rec_src, verdict = cause)
}

# -----------------------------------------------------------------------------
# Cell 6: summary and the paste-ready override block
# -----------------------------------------------------------------------------

aud <- rbindlist(rows)
hdr("Summary")
print(aud[, .(gene, in_atac_annotation, in_rna_assay, in_vendor_rna_table,
              in_raw_gtf, rec_chrom)])

out_csv <- file.path(OUT_DIR, paste0("symbol_audit_", TISSUE, ".csv"))
fwrite(aud, out_csv)
say("")
ok(paste0("wrote ", out_csv))

say("")
say(sprintf("  resolved in the ATAC annotation: %d of %d",
            sum(aud$in_atac_annotation), nrow(aud)))
say(sprintf("  recoverable by coordinate override: %d",
            sum(!aud$in_atac_annotation & !is.na(aud$rec_chrom))))
say(sprintf("  no coordinates anywhere: %d",
            sum(!aud$in_atac_annotation & is.na(aud$rec_chrom))))

if (length(override_lines)) {
  snippet <- file.path(OUT_DIR, paste0("GENE_REGION_OVERRIDE_", TISSUE, ".R"))
  writeLines(c(
    "# Paste into ARC_05_coverage_tracks.R next to GENE_WINDOW.",
    "# Generated by ARC_07_symbol_audit.R from the vendor's own peak coordinates,",
    "# for genes the Mmul8 GTF does not name. These are the SPAN OF THE VENDOR'S",
    "# PEAKS for the gene, not the gene body, so the window is approximately the",
    "# regulatory neighbourhood rather than an exact transcript model. Treat the",
    "# gene-model track for these as absent, and say so in any figure legend.",
    "GENE_REGION_OVERRIDE <- list(",
    override_lines,
    ")"), snippet)
  say("")
  ok(paste0("wrote ", snippet))
  say("")
  say("  That block is the fix for every gene the GTF does not name. Send it to")
  say("  me, or paste it in, and ARC_05 can draw those loci from coordinates.")
}

hdr("Done")
say("  Read the per-gene VERDICT lines above. They separate four different")
say("  problems that all look identical in ARC_05's log:")
say("")
say("    RESOLVED            the symbol is fine; the failure is later")
say("    COORDINATES RECOVERABLE  the GTF does not name it; use the override")
say("    IN RNA NOT IN GTF   the two references disagree; ask the vendor")
say("    ABSENT EVERYWHERE   not in this build; a finding, not a bug")
