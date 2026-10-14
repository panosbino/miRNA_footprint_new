library(tidyverse)
library(readxl)

BASE_DIR <- "~/Desktop/Projects/miRNA_footprint_new"
DB_DIR   <- file.path(BASE_DIR, "resources/mouse")
OUT_DIR  <- file.path(BASE_DIR, "analysis/mESCs/targets")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Adjust file names to what you actually downloaded.
MIRDB_FILE      <- file.path(DB_DIR, "miRDB_v6.0_prediction_result.txt.gz")
MIRTARBASE_FILE <- file.path(DB_DIR, "mmu_MTI.csv")
TARBASE_FILE    <- file.path(DB_DIR, "Mus_musculus_TarBase-v9.tsv")
MATURE_FA       <- file.path(DB_DIR, "mature.fa")           # miRBase v22.1

# --- seed families -----------------------------------------------------------
# Family = all mouse mature miRNAs whose seed (nt SEED_START-SEED_END) equals
# the representative's. USER_LISTS is only compared against this.
SEED_START <- 2
SEED_END   <- 8   # TargetScan 7mer definition; 2-7 would merge miR291a and miR292a
FAMILY_REPS <- c(miR17   = "mmu-miR-17-5p",
                 miR291a = "mmu-miR-291a-3p",
                 miR292a = "mmu-miR-292a-3p",
                 miR7    = "mmu-miR-7a-5p",     # no mmu-miR-7-5p in mouse; 7a/7b
                 miR183  = "mmu-miR-183-5p")
CONTROL_FAMILIES <- c("miR7", "miR183")
USER_LISTS <- list(
  miR17   = c("mmu-miR-17-5p", "mmu-miR-20a-5p", "mmu-miR-93-5p",
              "mmu-miR-106a-5p", "mmu-miR-106b-5p"),
  miR291a = c("mmu-miR-291a-3p", "mmu-miR-294-3p", "mmu-miR-295-3p", "mmu-miR-302a-3p"),
  miR292a = c("mmu-miR-292a-3p", "mmu-miR-467a-5p")
)

fa   <- readLines(MATURE_FA)
hdr  <- grep("^>", fa)
ends <- c(hdr[-1] - 1, length(fa))
mature <- tibble(
  name = sub("^>(\\S+).*$", "\\1", fa[hdr]),
  seq  = map2_chr(hdr + 1, ends, \(s, e) paste(fa[s:e], collapse = ""))
) |>
  dplyr::filter(startsWith(name, "mmu-")) |>
  mutate(seq  = gsub("T", "U", toupper(seq)),
         seed = substr(seq, SEED_START, SEED_END))

families <- imap_dfr(FAMILY_REPS, \(rep, fam) {
  s <- mature$seed[mature$name == rep]
  if (length(s) != 1) stop(rep, " not found (or duplicated) in ", MATURE_FA)
  dplyr::filter(mature, seed == s) |> mutate(family = fam)
})
dup <- families |> dplyr::count(name) |> dplyr::filter(n > 1)
if (nrow(dup) > 0) stop("miRNA in more than one family (shared seed): ", paste(dup$name, collapse = ", "))

cat(sprintf("\nSeed families (nt %d-%d) from %s:\n", SEED_START, SEED_END, basename(MATURE_FA)))
print(dplyr::select(families, family, name, seed, seq), n = Inf)
for (fam in names(USER_LISTS)) {
  members <- families$name[families$family == fam]
  not_in_family <- setdiff(USER_LISTS[[fam]], members)
  not_in_list   <- setdiff(members, USER_LISTS[[fam]])
  if (length(not_in_family)) cat(sprintf("  %s: in your list but not in seed family (other seed, or name not in miRBase): %s\n",
                                         fam, paste(not_in_family, collapse = ", ")))
  if (length(not_in_list))   cat(sprintf("  %s: same seed but not in your list: %s\n",
                                         fam, paste(not_in_list, collapse = ", ")))
}

# query_mirna holds the family label; db_mirna the member it came from
alias_tbl <- families |>
  transmute(query_mirna = family, db_mirna = name, key = tolower(name))

# --- helpers -----------------------------------------------------------------
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))[1]
  if (is.na(hit)) stop(sprintf("None of [%s] found. Columns are: %s",
                               paste(candidates, collapse = ", "),
                               paste(names(df), collapse = ", ")))
  hit
}

match_mirnas <- function(df, mirna_col, db) {
  hits <- df |>
    mutate(key = tolower(.data[[mirna_col]])) |>
    inner_join(dplyr::select(alias_tbl, query_mirna, key), by = "key") |>
    dplyr::select(-key)
  cat(sprintf("\n[%s] database names matched:\n", db))
  print(dplyr::count(hits, query_mirna, db_name = .data[[mirna_col]]))
  for (fam in unique(alias_tbl$query_mirna)) {
    absent <- setdiff(alias_tbl$key[alias_tbl$query_mirna == fam], tolower(hits[[mirna_col]]))
    if (length(absent)) cat(sprintf("  %s: no entries in %s for %s\n", fam, db,
                                    paste(absent, collapse = ", ")))
  }
  hits
}

# --- miRDB (prediction) ------------------------------------------------------
# File has no header: miRNA, RefSeq accession, target score (50-100).
mirdb_all <- read_tsv(MIRDB_FILE, col_names = c("mirna", "refseq", "score"),
                      col_types = "ccd") |>
  dplyr::filter(startsWith(mirna, "mmu-")) |>
  mutate(refseq = sub("\\.\\d+$", "", refseq))

# Map all mouse RefSeq IDs (not just matched ones): needed for the overlap background
refseq_map <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                                    keys = unique(mirdb_all$refseq), keytype = "REFSEQ",
                                    columns = c("SYMBOL", "ENTREZID"))
mirdb_universe <- unique(na.omit(refseq_map$SYMBOL))

mirdb_out <- match_mirnas(mirdb_all, "mirna", "miRDB") |>
  left_join(refseq_map, by = c("refseq" = "REFSEQ"), relationship = "many-to-many") |>
  group_by(query_mirna, db_mirna = mirna, gene_symbol = SYMBOL, entrez_id = ENTREZID) |>
  summarise(score = max(score), refseq = paste(unique(refseq), collapse = ";"),
            .groups = "drop")
cat(sprintf("miRDB: %d RefSeq IDs could not be mapped to a gene (kept, symbol = NA)\n",
            sum(is.na(mirdb_out$gene_symbol))))
write_tsv(mirdb_out, file.path(OUT_DIR, "targets_mirdb.tsv"))

# --- miRTarBase (validated) --------------------------------------------------
mtb_all <- read_csv(MIRTARBASE_FILE, show_col_types = FALSE)
c_mir  <- pick_col(mtb_all, c("miRNA"))
c_gene <- pick_col(mtb_all, c("Target Gene"))
c_ent  <- pick_col(mtb_all, c("Target Gene (Entrez ID)", "Target Gene (Entrez Gene ID)"))
c_sup  <- pick_col(mtb_all, c("Support Type"))
c_exp  <- pick_col(mtb_all, c("Experiments"))
c_pmid <- pick_col(mtb_all, c("References (PMID)"))
c_sp   <- intersect("Species (Target Gene)", names(mtb_all))

if (length(c_sp) == 1) {
  n0 <- nrow(mtb_all); mtb_all <- dplyr::filter(mtb_all, .data[[c_sp]] == "mmu")
  cat(sprintf("miRTarBase: dropped %d rows with non-mouse target genes\n", n0 - nrow(mtb_all)))
}
mtb_universe <- unique(na.omit(mtb_all[[c_gene]]))

mtb_out <- match_mirnas(mtb_all, c_mir, "miRTarBase") |>
  group_by(query_mirna, db_mirna = .data[[c_mir]], gene_symbol = .data[[c_gene]],
           entrez_id = as.character(.data[[c_ent]])) |>
  summarise(support_type = paste(sort(unique(.data[[c_sup]])), collapse = "; "),
            experiments  = paste(sort(unique(unlist(strsplit(.data[[c_exp]], "//")))), collapse = "; "),
            n_pmids      = n_distinct(.data[[c_pmid]]),
            .groups = "drop")
write_tsv(mtb_out, file.path(OUT_DIR, "targets_mirtarbase.tsv"))

# --- TarBase v9 (validated, mostly CLIP) -------------------------------------
tb_all <- read_tsv(TARBASE_FILE, show_col_types = FALSE, guess_max = 1e5)
t_mir  <- pick_col(tb_all, c("mirna_name", "mirna", "miRNA"))
t_gene <- pick_col(tb_all, c("gene_name", "geneName"))
t_id   <- pick_col(tb_all, c("gene_id", "geneId", "ensembl_gene_id"))
t_meth <- pick_col(tb_all, c("experimental_method", "method"))
t_reg  <- intersect(c("regulation", "up_down", "positive_negative"), names(tb_all))[1]
t_sp   <- intersect(c("species"), names(tb_all))

if (length(t_sp) == 1) tb_all <- dplyr::filter(tb_all, grepl("mus musculus|^mmu$", .data[[t_sp]], ignore.case = TRUE))
tb_all <- mutate(tb_all, gene_symbol = sub("\\(mmu\\)$", "", .data[[t_gene]]))
tb_universe <- unique(na.omit(tb_all$gene_symbol))

tb <- match_mirnas(tb_all, t_mir, "TarBase")
if (!is.na(t_reg)) { cat("TarBase regulation values:\n"); print(table(tb[[t_reg]], useNA = "ifany")) }

tb_out <- tb |>
  group_by(query_mirna, db_mirna = .data[[t_mir]], gene_symbol,
           gene_id = .data[[t_id]]) |>
  summarise(methods    = paste(sort(unique(.data[[t_meth]])), collapse = "; "),
            regulation = if (is.na(t_reg)) NA_character_ else
                           paste(sort(unique(.data[[t_reg]])), collapse = "; "),
            n_records  = n(),
            .groups = "drop")
write_tsv(tb_out, file.path(OUT_DIR, "targets_tarbase.tsv"))

# --- summary -----------------------------------------------------------------
target_sets <- bind_rows(miRDB = mirdb_out, miRTarBase = mtb_out, TarBase = tb_out, .id = "db") |>
  dplyr::filter(!is.na(gene_symbol)) |>
  distinct(db, family = query_mirna, gene_symbol)

cat("\nUnique target genes per family and database:\n")
target_sets |>
  dplyr::count(family, db) |>
  pivot_wider(names_from = db, values_from = n, values_fill = 0) |>
  print()

# --- overlap between families, per database ----------------------------------
# Background = all genes with any mouse entry in that database. Fisher test is
# one-sided (more overlap than chance). Control pairs give an empirical
# baseline for what "unrelated families" overlap looks like.
universes <- list(miRDB = mirdb_universe, miRTarBase = mtb_universe, TarBase = tb_universe)

overlap <- target_sets |>
  group_by(db) |>
  group_modify(\(d, key) {
    sets <- split(d$gene_symbol, d$family)
    if (length(sets) < 2) return(tibble())
    U <- length(universes[[key$db]])
    map_dfr(combn(names(sets), 2, simplify = FALSE), \(p) {
      a <- sets[[p[1]]]; b <- sets[[p[2]]]
      k <- length(intersect(a, b))
      tab <- matrix(c(k, length(a) - k, length(b) - k, U - length(union(a, b))), nrow = 2)
      ft  <- fisher.test(tab, alternative = "greater")
      tibble(family_a = p[1], family_b = p[2],
             pair_type = c("focal-focal", "focal-control", "control-control")[
               1 + sum(p %in% CONTROL_FAMILIES)],
             n_a = length(a), n_b = length(b), n_shared = k,
             expected_shared = length(a) * length(b) / U,
             jaccard = k / length(union(a, b)),
             overlap_coef = k / min(length(a), length(b)),
             odds_ratio = unname(ft$estimate), p = ft$p.value)
    })
  }) |>
  ungroup() |>
  mutate(p_adj = p.adjust(p, method = "BH")) |>
  arrange(db, pair_type, family_a, family_b)

cat("\nPairwise target overlap (gene symbols):\n")
print(overlap, n = Inf, width = Inf)
write_tsv(overlap, file.path(OUT_DIR, "family_target_overlap.tsv"))
