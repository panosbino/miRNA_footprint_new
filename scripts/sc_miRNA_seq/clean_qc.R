#!/usr/bin/env Rscript
# =============================================================================
# 01_clean_qc.R
# Clean a single-cell miRNA count table (MirGeneDB names, one column per well)
# and flag failed wells.
#
# Usage: Rscript 01_clean_qc.R <count_table.tsv> [out_dir]
#
# What it does, and why:
#   1. Collapses EXACT duplicate rows. Multi-mapping reads are counted once per
#      locus (e.g. Mir-124 P1/P2/P3), which inflates library sizes and counts
#      the same molecules several times.
#   2. Collapses NEAR duplicates: isoforms/paralogs of the SAME family and arm
#      whose counts differ by <= NEAR_DUP_MAXDIFF (e.g. Mir-124-P1-v1 vs -v2).
#      Collapsed by per-cell MAXIMUM, not sum: the reads are largely shared, so
#      summing would double count. This is an approximation; re-quantifying with
#      unique/fractional assignment would be the proper fix.
#      Highly correlated features from DIFFERENT families are NOT merged, only
#      reported for manual inspection (can be cross-mapping or coincidence).
#   3. Per-well QC. A well FAILS only if it has almost no reads or a single
#      feature dominates it. Low depth is flagged (low_depth) but kept.
#      Composition similarity to a leave-one-out consensus is tested against a
#      depth-matched downsampling null and reported as composition_outlier, but
#      NOT used to exclude: without UMIs reads come in PCR clusters, so read-level
#      downsampling underestimates real well-to-well variation and the test
#      over-flags ordinary jackpot noise. Inspect these wells; don't auto-drop.
#   4. Jackpot screen: features with >= JACKPOT_SHARE of their reads in one well
#      (one molecule amplified many times; typical of no-UMI data).
#      Reported only; such features are excluded from multi-miRNA analyses by the
#      detection filter in script 02.
# =============================================================================

suppressPackageStartupMessages(library(stats))

## ---- Parameters -------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
in_file  <- if (length(args) >= 1) args[1] else "mirna_count_table_v2_2nd_run.tsv"
out_dir  <- if (length(args) >= 2) args[2] else "results"

COL_OFFSET        <- 6     # sequenced column + offset = physical plate column
                           # (+6 for run 2 = physical columns 7-12; UNCONFIRMED,
                           #  taken from the preprocessing script; set 0 for run 1)
TARGET_REGEX      <- "^Hsa-Mir-124-.*_3p$"  # experimental miRNA; excluded from the
                                            # consensus profile so high-miR-124 cells
                                            # are not penalised as "abnormal"
NEAR_DUP_MAXDIFF  <- 0.05  # max relative count difference to call a near duplicate
NEAR_DUP_MINCOUNT <- 10    # both features need >= this many total reads
MIN_READS_FAIL    <- 100   # below this a well is failed regardless
TOP_SHARE_FAIL    <- 0.90  # one feature >= 90% of reads -> failed
LOW_DEPTH         <- 2000  # flagged (not removed) below this many reads
PROFILE_MIN_PROP  <- 0.002 # features used for the similarity score
N_NULL            <- 200   # downsampling draws per well for the depth-matched null
SIM_P_FLAG        <- 0.01  # empirical p below this -> composition_outlier (flag only)
JACKPOT_SHARE     <- 0.90  # jackpot screen: share of a feature's reads in one well
JACKPOT_MIN       <- 50    # ... for features with at least this many reads
COR_MIN_WELLS     <- 3     # cross-family correlation check needs >= this many non-zero wells
SEED              <- 1

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)
log_file <- file.path(out_dir, "01_qc_log.txt")
sink(log_file, split = TRUE)
cat("01_clean_qc.R |", format(Sys.time()), "\ninput:", in_file, "\n\n")

## ---- Read ------------------------------------------------------------------
raw  <- read.delim(in_file, check.names = FALSE, stringsAsFactors = FALSE)
feat <- sub("^>", "", raw[[1]])
if (anyDuplicated(feat)) stop("Duplicated feature names in input.")
cnt  <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(cnt) <- "numeric"
rownames(cnt) <- feat
colnames(cnt) <- sub("_.*$", "", colnames(cnt))          # H01_S8_L001_R1_001 -> H01
if (!all(grepl("^[A-P][0-9]{2}$", colnames(cnt))))
  stop("Could not parse well IDs from column names: ",
       paste(head(colnames(cnt)[!grepl("^[A-P][0-9]{2}$", colnames(cnt))]), collapse = ", "))
if (anyDuplicated(colnames(cnt))) stop("Duplicated well IDs after parsing.")
if (any(cnt < 0) || any(cnt != round(cnt))) stop("Counts must be non-negative integers.")
cat(sprintf("raw: %d features x %d wells, %d features with zero counts\n",
            nrow(cnt), ncol(cnt), sum(rowSums(cnt) == 0)))

cnt <- cnt[rowSums(cnt) > 0, , drop = FALSE]

family_of <- function(x) sub("^([^-]+-[^-_]+-[^-_]+).*$", "\\1", x)  # Hsa-Mir-10-P1c-v2_5p -> Hsa-Mir-10
arm_of    <- function(x) ifelse(grepl("_[35]p", x), sub("^.*_([35]p).*$", "\\1", x), "NA")

## ---- Union-find helper -----------------------------------------------------
group_pairs <- function(n, pairs) {
  parent <- seq_len(n)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (k in seq_len(nrow(pairs))) {
    a <- find(pairs[k, 1]); b <- find(pairs[k, 2])
    if (a != b) parent[b] <- a
  }
  vapply(seq_len(n), find, integer(1))
}

## ---- 1. Exact duplicates -----------------------------------------------------
key   <- apply(cnt, 1, paste, collapse = ",")
grp1  <- match(key, unique(key))
members1 <- split(rownames(cnt), grp1)
cross_fam <- Filter(function(g) length(unique(family_of(g))) > 1, members1)
if (length(cross_fam))
  cat("WARNING: exact duplicates spanning different families:\n",
      paste(" ", vapply(cross_fam, paste, "", collapse = " | "), collapse = "\n"), "\n")

rep_of <- function(g, totals) g[order(-totals[g], g)][1]   # representative name
tot0   <- rowSums(cnt)
reps1  <- vapply(members1, rep_of, "", totals = tot0)
cnt1   <- cnt[reps1, , drop = FALSE]
map1   <- data.frame(feature = unlist(members1),
                     step1_rep = rep(reps1, lengths(members1)), stringsAsFactors = FALSE)
cat(sprintf("exact duplicates: %d non-zero rows -> %d unique (%.1f%% of raw counts were in duplicate rows)\n",
            nrow(cnt), nrow(cnt1),
            100 * sum(cnt[duplicated(key) | duplicated(key, fromLast = TRUE), ]) / sum(cnt)))

## ---- 2. Near duplicates (same family + arm only) ----------------------------
fam <- family_of(rownames(cnt1)); arm <- arm_of(rownames(cnt1)); tot1 <- rowSums(cnt1)
pairs <- matrix(integer(0), ncol = 2)
for (fa in unique(paste(fam, arm))) {
  idx <- which(paste(fam, arm) == fa & tot1 >= NEAR_DUP_MINCOUNT)
  if (length(idx) < 2) next
  for (i in idx) for (j in idx) if (i < j) {
    d <- sum(abs(cnt1[i, ] - cnt1[j, ])) / mean(c(tot1[i], tot1[j]))
    if (d <= NEAR_DUP_MAXDIFF) pairs <- rbind(pairs, c(i, j))
  }
}
grp2     <- group_pairs(nrow(cnt1), pairs)
members2 <- split(rownames(cnt1), grp2)
reps2    <- vapply(members2, rep_of, "", totals = tot1)
clean    <- t(vapply(members2, function(g) apply(cnt1[g, , drop = FALSE], 2, max),
                     numeric(ncol(cnt1))))
rownames(clean) <- reps2
colnames(clean) <- colnames(cnt1)
lookup <- setNames(rep(reps2, lengths(members2)), unlist(members2))
map1$final_feature <- lookup[map1$step1_rep]
cat(sprintf("near duplicates (same family+arm, <= %.0f%% different): %d -> %d features\n",
            100 * NEAR_DUP_MAXDIFF, nrow(cnt1), nrow(clean)))
for (g in members2[lengths(members2) > 1])
  cat("  merged:", paste(g, collapse = " + "), "\n")

# Report (do not merge) suspiciously correlated features from different families
big <- clean[rowSums(clean) >= 50 & rowSums(clean > 0) >= COR_MIN_WELLS, , drop = FALSE]
if (nrow(big) > 1) {
  cc <- cor(t(log1p(big))); cc[lower.tri(cc, diag = TRUE)] <- NA
  hit <- which(cc > 0.99, arr.ind = TRUE)
  hit <- hit[family_of(rownames(big)[hit[, 1]]) != family_of(rownames(big)[hit[, 2]]), , drop = FALSE]
  if (nrow(hit)) {
    cat("NOTE: r > 0.99 between different families (NOT merged; inspect for cross-mapping):\n")
    for (k in seq_len(nrow(hit)))
      cat(sprintf("  %s ~ %s (r = %.3f)\n", rownames(big)[hit[k, 1]], colnames(cc)[hit[k, 2]],
                  cc[hit[k, 1], hit[k, 2]]))
  }
}

## ---- 3. Well QC --------------------------------------------------------------
is_target <- grepl(TARGET_REGEX, rownames(clean))
cat(sprintf("\ntarget feature(s) matching '%s': %s\n", TARGET_REGEX,
            if (any(is_target)) paste(rownames(clean)[is_target], collapse = ", ") else "NONE"))
nt   <- clean[!is_target, , drop = FALSE]
lib  <- colSums(clean)
prop <- sweep(nt, 2, pmax(colSums(nt), 1), "/")

basic_fail <- lib < MIN_READS_FAIL | apply(prop, 2, max) >= TOP_SHARE_FAIL
# Consensus profile from deep, non-failed wells (top half by depth)
ref_cells <- names(lib)[!basic_fail & lib >= median(lib[!basic_fail])]
# Leave-one-out consensus: a well is never compared with a profile containing itself
loo_profile <- function(excl) rowMeans(prop[, setdiff(ref_cells, excl), drop = FALSE])
sim_to <- function(v, prof) {
  use <- prof >= PROFILE_MIN_PROP
  suppressWarnings(cor(v[use], prof[use], method = "spearman"))
}
sim <- vapply(colnames(prop), function(w) sim_to(prop[, w], loo_profile(w)), numeric(1))
loo_cache <- lapply(setNames(ref_cells, ref_cells), loo_profile)

downsample <- function(v, n) {           # sampling reads without replacement
  if (sum(v) <= n) return(v)
  tabulate(sample(rep.int(seq_along(v), v), n), nbins = length(v))
}
# Depth-matched null: a random reference well, downsampled to this well's depth,
# compared with the consensus that excludes that reference well.
# Caveat: read-level downsampling ignores PCR clustering, so this null is too
# narrow for no-UMI data (see header) -> used as a flag only.
sim_p <- sim_null5 <- setNames(rep(NA_real_, ncol(clean)), colnames(clean))
for (w in colnames(clean)) {
  n_w <- sum(nt[, w]); if (n_w < 1) next
  draws <- vapply(seq_len(N_NULL), function(i) {
    r <- sample(ref_cells, 1)
    s <- downsample(nt[, r], n_w)
    sim_to(s / sum(s), loo_cache[[r]])
  }, numeric(1))
  sim_p[w]     <- (sum(draws <= sim[w], na.rm = TRUE) + 1) / (sum(!is.na(draws)) + 1)
  sim_null5[w] <- quantile(draws, 0.05, na.rm = TRUE)
}

qc <- data.frame(
  well        = colnames(clean),
  row         = substr(colnames(clean), 1, 1),
  col_seq     = as.integer(substr(colnames(clean), 2, 3)),
  col_phys    = as.integer(substr(colnames(clean), 2, 3)) + COL_OFFSET,
  lib_raw     = colSums(cnt)[colnames(clean)],
  lib_clean   = lib,
  lib_nontarget = colSums(nt),
  target_reads  = if (any(is_target)) colSums(clean[is_target, , drop = FALSE]) else NA,
  n_detected  = colSums(clean > 0),
  n_ge5       = colSums(clean >= 5),
  top_feature = rownames(prop)[apply(prop, 2, which.max)],
  top_share   = round(apply(prop, 2, max), 3),
  similarity  = round(sim, 3),
  sim_null_q05 = round(sim_null5, 3),
  sim_p       = signif(sim_p, 3),
  stringsAsFactors = FALSE)
qc$fail_reason <- with(qc, ifelse(lib_clean < MIN_READS_FAIL, "too_few_reads",
                        ifelse(top_share >= TOP_SHARE_FAIL, "single_feature_dominates", "")))
qc$failed    <- qc$fail_reason != ""
qc$composition_outlier <- !qc$failed & !is.na(qc$sim_p) & qc$sim_p < SIM_P_FLAG
qc$low_depth <- qc$lib_clean < LOW_DEPTH & !qc$failed

cat(sprintf("\nlibrary size (clean): median %.0f, range %.0f-%.0f\n",
            median(lib), min(lib), max(lib)))
ok <- !qc$failed
kr <- kruskal.test(log10(qc$lib_clean[ok]) ~ factor(qc$row[ok]))$p.value
kc <- kruskal.test(log10(qc$lib_clean[ok]) ~ factor(qc$col_phys[ok]))$p.value
cat(sprintf("plate effect on depth (non-failed wells): row p = %.2g, column p = %.2g\n", kr, kc))
cat("median clean library size by row:\n"); print(tapply(qc$lib_clean[ok], qc$row[ok], median))
cat("\nFAILED wells (excluded downstream):\n")
print(qc[qc$failed, c("well", "lib_clean", "top_feature", "top_share", "similarity", "sim_p", "fail_reason")],
      row.names = FALSE)
cat("\ncomposition_outlier wells (KEPT; inspect, likely jackpot noise):\n")
print(qc[qc$composition_outlier, c("well", "lib_clean", "top_feature", "similarity", "sim_null_q05", "sim_p")],
      row.names = FALSE)

## ---- 4. Jackpot screen --------------------------------------------------------
live <- clean[, !qc$failed, drop = FALSE]
tot  <- rowSums(live); mx <- apply(live, 1, max)
jp   <- tot >= JACKPOT_MIN & mx / pmax(tot, 1) >= JACKPOT_SHARE
jackpots <- data.frame(feature = rownames(live)[jp], total_reads = tot[jp],
                       top_well = colnames(live)[apply(live[jp, , drop = FALSE], 1, which.max)],
                       top_share = round(mx[jp] / tot[jp], 3), row.names = NULL)
cat(sprintf("\njackpot features (>= %.0f%% of >= %d reads in one well): %d\n",
            100 * JACKPOT_SHARE, JACKPOT_MIN, nrow(jackpots)))
if (nrow(jackpots)) print(jackpots[order(-jackpots$total_reads), ], row.names = FALSE)
if (any(jp & grepl(TARGET_REGEX, rownames(live))))
  cat("WARNING: the target miRNA itself looks like a jackpot feature!\n")
write.table(jackpots, file.path(out_dir, "jackpot_features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nlow_depth wells (< %d reads; KEPT, check sensitivity): %s\n", LOW_DEPTH,
            paste(qc$well[qc$low_depth], collapse = ", ")))

## ---- Write ------------------------------------------------------------------
out_counts <- data.frame(feature = rownames(clean), clean, check.names = FALSE)
write.table(out_counts, file.path(out_dir, "counts_clean.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(map1[, c("feature", "final_feature")], file.path(out_dir, "feature_collapse_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc, file.path(out_dir, "cell_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(out_dir, "01_qc_plots.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
rows <- sort(unique(qc$row)); cols <- sort(unique(qc$col_phys))
M <- matrix(NA, length(rows), length(cols), dimnames = list(rows, cols))
M[cbind(match(qc$row, rows), match(qc$col_phys, cols))] <- log10(pmax(qc$lib_clean, 1))
image(seq_along(cols), seq_along(rows), t(M[rev(rows), , drop = FALSE]), axes = FALSE,
      col = hcl.colors(30, "viridis"), xlab = "physical column", ylab = "row",
      main = "log10 clean library size")
axis(1, seq_along(cols), cols); axis(2, seq_along(rows), rev(rows), las = 1)
text(rep(seq_along(cols), each = length(rows)), rep(seq_along(rows), length(cols)),
     round(t(M[rev(rows), , drop = FALSE])[cbind(rep(seq_along(cols), each = length(rows)),
                                                   rep(seq_along(rows), length(cols)))], 1),
     cex = 0.7, col = "white")
o <- order(qc$lib_clean)
plot(log10(qc$lib_clean), qc$similarity, pch = ifelse(qc$failed, 4, 19),
     col = ifelse(qc$failed, "red", ifelse(qc$composition_outlier, "purple",
                  ifelse(qc$low_depth, "orange", "grey30"))),
     xlab = "log10 clean library size", ylab = "Spearman similarity to consensus",
     main = "Composition vs depth-matched null")
lines(log10(qc$lib_clean[o]), qc$sim_null_q05[o], lty = 2)
legend("bottomright", c("kept", "low depth (kept)", "composition outlier (kept)", "failed", "null 5th pct"),
       pch = c(19, 19, 19, 4, NA), lty = c(NA, NA, NA, NA, 2),
       col = c("grey30", "orange", "purple", "red", "black"),
       bty = "n", cex = 0.8)
dev.off()
cat("\nwritten to", out_dir, ": counts_clean.tsv, feature_collapse_map.tsv, cell_qc.tsv, \n  jackpot_features.tsv, 01_qc_plots.pdf\n")
sink()
