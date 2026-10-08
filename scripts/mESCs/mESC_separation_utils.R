library(tidyverse)

# ---------------------------------------------------------------------------
# Shared separation-statistics + plotting logic, factored out of
# separate_KO_control.R so every method-specific script (ours, miReact,
# bayesReact, miTEA-HiRes) uses IDENTICAL math -- avoids the risk of subtly
# different formulas creeping in across scripts, the same duplication-risk
# principle applied throughout this whole project.
#
# Replicates the ORIGINAL script's specific conventions exactly, not a
# "corrected" or more standard version:
#   - Cohen's-d-style effect size using the KO group's OWN SD as the
#     denominator (not a pooled SD across both groups).
#   - Overlap % via 2*pnorm(-abs(d)/2)*100.
#
# PRIMARY METRIC (added): single-cutoff misclassification.
#   A cell is called "Control" if activity > cutoff, else "KO". The direction
#   is fixed a priori (controls = intact miRNAs = higher activity), so a score
#   that separates the wrong way gets ~50% (chance), unlike abs(d) above.
#   The cutoff minimises the BALANCED error (mean of the two per-group error
#   rates), so the cutoff cannot win by calling everything the larger group.
#   Rank-based, so it is comparable across methods with different score scales
#   (log2 ratio, -log10 p, posterior means), and needs no normality/equal-SD
#   assumption. overlap/2 is the same quantity under equal-variance normals.
#   bal_error_cv: cutoff picked on 4/5 of cells, applied to the held-out 1/5,
#   removing the optimism of choosing the best cutoff on the scored cells.
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"

sep_error <- function(score, is_ko) {
  stopifnot(length(score) == length(is_ko), !anyNA(score), !anyNA(is_ko))
  ctrl <- score[!is_ko]; ko <- score[is_ko]
  v   <- sort(unique(score))
  cut <- c(-Inf, (head(v, -1) + tail(v, -1)) / 2, Inf)       # cutoffs between distinct values
  ctrl_wrong <- round(ecdf(ctrl)(cut) * length(ctrl))         # controls <= cutoff -> called KO
  ko_wrong   <- round((1 - ecdf(ko)(cut)) * length(ko))       # KOs > cutoff -> called Control
  bal <- (ctrl_wrong / length(ctrl) + ko_wrong / length(ko)) / 2
  i <- which.min(bal)
  list(bal_error = bal[i], cutoff = cut[i], ctrl_called_ko = ctrl_wrong[i],
       ko_called_ctrl = ko_wrong[i], cells_wrong = ctrl_wrong[i] + ko_wrong[i])
}

# Note: calls set.seed(), so it resets the global RNG stream.
sep_error_cv <- function(score, is_ko, k = 5, seed = 1) {
  set.seed(seed)
  fold <- integer(length(score))
  fold[is_ko]  <- sample(rep_len(1:k, sum(is_ko)))            # stratified folds
  fold[!is_ko] <- sample(rep_len(1:k, sum(!is_ko)))
  wrong_ctrl <- 0; wrong_ko <- 0
  for (f in 1:k) {
    cut <- sep_error(score[fold != f], is_ko[fold != f])$cutoff
    wrong_ctrl <- wrong_ctrl + sum(score[fold == f & !is_ko] <= cut)
    wrong_ko   <- wrong_ko   + sum(score[fold == f &  is_ko] >  cut)
  }
  (wrong_ctrl / sum(!is_ko) + wrong_ko / sum(is_ko)) / 2
}

compute_separation <- function(activity_named_vec, pheno) {
  # activity_named_vec: named numeric vector, names = cell IDs matching pheno$cells
  activity_df <- data.frame(cells = names(activity_named_vec), activity = as.numeric(activity_named_vec))
  merged <- merge(activity_df, pheno, by = "cells")

  ko_activ <- merged[merged$Exp == "KO", ]
  control_activ <- merged[merged$Exp == "Control", ]

  d <- (mean(ko_activ$activity) - mean(control_activ$activity)) / sd(ko_activ$activity)
  overlap <- (2 * pnorm(-abs(d) / 2)) * 100
  diff <- mean(control_activ$activity) - mean(ko_activ$activity)

  is_ko <- merged$Exp == "KO"
  sep <- sep_error(merged$activity, is_ko)

  c(list(merged = merged, d = d, overlap = overlap, diff = diff,
         n_ko = nrow(ko_activ), n_control = nrow(control_activ)),
    sep, list(bal_error_cv = sep_error_cv(merged$activity, is_ko)))
}

plot_separation <- function(sep_result, title_prefix) {
  ggplot(sep_result$merged, aes(x = activity, fill = Exp)) +
    geom_density(color = "black", alpha = 0.8, linewidth = 1.3) +
    geom_vline(xintercept = sep_result$cutoff, linetype = "dashed", linewidth = 1) +
    theme(panel.border = element_rect(colour = "black", linewidth = 2, fill = NA),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "black", linetype = "dashed", linewidth = 0.2),
          panel.grid.major.x = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 16, hjust = 0.5),
          legend.title = element_blank()) +
    scale_fill_manual(values = c("royalblue4", "lightblue1")) +
    labs(title = sprintf("%s\nMisclassified = %d/%d cells (balanced error %.1f%%, CV %.1f%%)",
                          title_prefix, sep_result$cells_wrong,
                          sep_result$n_ko + sep_result$n_control,
                          100 * sep_result$bal_error, 100 * sep_result$bal_error_cv),
         subtitle = sprintf("Overlap (normal approx.) = %.3f%%   Mean difference = %.2f   Dashed line = cutoff",
                             sep_result$overlap, sep_result$diff),
         x = "Activity", y = "Density")
}

report_separation <- function(sep_result, label) {
  cat(sprintf("%-40s: misclassified = %3d (Ctrl->KO %3d, KO->Ctrl %3d), bal. error = %5.1f%%, CV = %5.1f%% | overlap = %6.3f%%, mean diff = %7.3f, d = %6.3f (n_KO=%d, n_Control=%d)\n",
              label, sep_result$cells_wrong, sep_result$ctrl_called_ko, sep_result$ko_called_ctrl,
              100 * sep_result$bal_error, 100 * sep_result$bal_error_cv,
              sep_result$overlap, sep_result$diff, sep_result$d,
              sep_result$n_ko, sep_result$n_control))
}
