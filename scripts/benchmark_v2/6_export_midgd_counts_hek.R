# ==============================================================================
# HEK293 benchmark v2, step 6: export raw EXONIC UMI counts for miDGD.
#
# miDGD needs raw integer counts for ALL genes, exonic where there is a choice
# (miDGD setup notes). The scran matrix is normalised, protein-coding only and
# built from intron+exon counts, so it cannot be used. This reads the zUMIs
# output (umicount$exon$all) for plates 101/103/105, names cells <plate>_<barcode>
# and keeps exactly the 753 cells of the scran matrix.
#
# Input : datasets/HEK_SS3/raw/P32156_<plate>.dgecounts.rds   (on Dardel, not in git)
# Output: analysis/HEK_SS3/benchmark_v2/midgd/hek_exon_umi.tsv     genes x cells, gene_id column
# Then run step 4b (see the run order) with the midgd environment.
# ==============================================================================

suppressPackageStartupMessages({ library(Matrix); library(data.table) })

BASE_DIR <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
RAW_DIR  <- file.path(BASE_DIR, "datasets/HEK_SS3/raw")
OUT_DIR  <- file.path(BASE_DIR, "analysis/HEK_SS3/benchmark_v2/midgd")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
PLATES <- c("101", "103", "105")

cells <- colnames(readRDS(file.path(BASE_DIR, "datasets/HEK_SS3/processed/scran_normalized_linear.rds")))

mats <- lapply(PLATES, function(p) {
  f <- file.path(RAW_DIR, sprintf("P32156_%s.dgecounts.rds", p))
  if (!file.exists(f)) stop("Not found: ", f)
  m <- readRDS(f)$umicount$exon$all
  colnames(m) <- paste0(p, "_", colnames(m))
  m[, colnames(m) %in% cells, drop = FALSE]
})
genes <- Reduce(union, lapply(mats, rownames))
m <- do.call(cbind, lapply(mats, function(x) {
  y <- Matrix(0, length(genes), ncol(x), dimnames = list(genes, colnames(x)), sparse = TRUE)
  y[rownames(x), ] <- x
  y
}))

missing <- setdiff(cells, colnames(m))
cat(sprintf("Exonic UMI matrix: %d genes x %d cells; %d of %d scran cells found\n",
            nrow(m), ncol(m), length(cells) - length(missing), length(cells)))
if (length(missing) > 0.05 * length(cells))
  stop(length(missing), " scran cells not found in the dgecounts files -- check the <plate>_<barcode> naming.")
if (length(missing) > 0) cat("Cells without exonic counts (will be absent from miDGD):", head(missing), "\n")
stopifnot(all(m@x == round(m@x)))

m <- m[, intersect(cells, colnames(m))]
out <- file.path(OUT_DIR, "hek_exon_umi.tsv")
fwrite(as.data.table(as.matrix(m), keep.rownames = "gene_id"), out, sep = "\t")
cat("Wrote", out, "\n")
