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
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"

compute_separation <- function(activity_named_vec, pheno) {
  # activity_named_vec: named numeric vector, names = cell IDs matching pheno$cells
  activity_df <- data.frame(cells = names(activity_named_vec), activity = as.numeric(activity_named_vec))
  merged <- merge(activity_df, pheno, by = "cells")

  ko_activ <- merged[merged$Exp == "KO", ]
  control_activ <- merged[merged$Exp == "Control", ]

  d <- (mean(ko_activ$activity) - mean(control_activ$activity)) / sd(ko_activ$activity)
  overlap <- (2 * pnorm(-abs(d) / 2)) * 100
  diff <- mean(control_activ$activity) - mean(ko_activ$activity)

  list(merged = merged, d = d, overlap = overlap, diff = diff,
       n_ko = nrow(ko_activ), n_control = nrow(control_activ))
}

plot_separation <- function(sep_result, title_prefix) {
  ggplot(sep_result$merged, aes(x = activity, fill = Exp)) +
    geom_density(color = "black", alpha = 0.8, linewidth = 1.3) +
    theme(panel.border = element_rect(colour = "black", linewidth = 2, fill = NA),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "black", linetype = "dashed", linewidth = 0.2),
          panel.grid.major.x = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 18, hjust = 0.5),
          legend.title = element_blank()) +
    scale_fill_manual(values = c("royalblue4", "lightblue1")) +
    labs(title = sprintf("%s     Overlap = %.3f%%   Mean difference = %.2f",
                          title_prefix, sep_result$overlap, sep_result$diff),
         x = "Activity", y = "Density")
}

report_separation <- function(sep_result, label) {
  cat(sprintf("%-40s: overlap = %6.3f%%, mean diff = %7.3f, d = %6.3f (n_KO=%d, n_Control=%d)\n",
              label, sep_result$overlap, sep_result$diff, sep_result$d,
              sep_result$n_ko, sep_result$n_control))
}
