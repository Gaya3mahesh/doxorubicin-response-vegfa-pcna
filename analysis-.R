# Doxorubicin sensitivity vs VEGFA and PCNA expression (PRISM 20Q2 + DepMap 25Q2)
# Run from the project folder. Needs these files in ./data:
#   prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv
#   OmicsExpressionProteinCodingGenesTPMLogp1Stranded.csv

library(tidyverse)

dr_file   <- "data/prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv"
expr_file <- "data/OmicsExpressionProteinCodingGenesTPMLogp1Stranded.csv"

dir.create("figures", showWarnings = FALSE)

# 1. Doxorubicin dose-response parameters from PRISM ---------------------------
dr <- read_csv(
  dr_file,
  col_select = c(broad_id, depmap_id, ccle_name, screen_id,
                 auc, ec50, ic50, r2, passed_str_profiling, name),
  show_col_types = FALSE
)

dox <- dr %>% filter(str_detect(str_to_lower(name), "doxorubicin"))

# A few doxorubicin rows have no cell line ID, so they are dropped.
# All remaining rows come from one screen (HTS002), one row per cell line.
dox_clean <- dox %>%
  filter(!is.na(depmap_id), !is.na(ccle_name), screen_id == "HTS002") %>%
  select(depmap_id, ccle_name, screen_id, auc, ec50, ic50, r2, passed_str_profiling)

stopifnot(nrow(count(dox_clean, depmap_id) %>% filter(n > 1)) == 0)

# 2. VEGFA and PCNA expression from DepMap -------------------------------------
get_genes <- function(file, genes) {
  hdr <- names(read_csv(file, n_max = 0, show_col_types = FALSE))
  pattern <- paste0("^(", paste(genes, collapse = "|"), ") \\(")
  idx <- c(1, grep(pattern, hdr))
  read_csv(file, col_select = all_of(idx), show_col_types = FALSE) %>%
    rename(depmap_id = 1) %>%
    rename_with(~ sub(" \\(.*", "", .x), .cols = -1)
}

genes <- get_genes(expr_file, c("VEGFA", "PCNA")) %>%
  filter(!is.na(depmap_id))

# 3. Join the two datasets -----------------------------------------------------
combined_all <- dox_clean %>%
  inner_join(genes, by = "depmap_id") %>%
  filter(!is.na(VEGFA), !is.na(PCNA), !is.na(ic50))

cat("Doxorubicin response records:", nrow(dox_clean), "\n")
cat("Matched cell lines:", n_distinct(combined_all$depmap_id), "\n")

# 4. The three cell lines I work with ------------------------------------------
my_dox <- combined_all %>%
  filter(str_detect(ccle_name, "^(MDAMB231|MCF7|NCIH358)_")) %>%
  select(ccle_name, VEGFA, PCNA, auc, ic50, r2)
print(my_dox)
write_csv(my_dox, "results_my_lines_doxorubicin.csv")

# 5. Spearman correlation: all lines -------------------------------------------
cor_one <- function(data, gene) {
  x <- data[[gene]]; y <- data$ic50
  keep <- complete.cases(x, y)
  test <- cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
  tibble(Gene = gene, N = sum(keep),
         Correlation = unname(test$estimate), P_value = test$p.value)
}

cor_all <- bind_rows(cor_one(combined_all, "VEGFA"),
                     cor_one(combined_all, "PCNA")) %>%
  mutate(Adjusted_P = p.adjust(P_value, method = "BH"))
print(cor_all)
write_csv(cor_all, "results_correlations_all_lines.csv")

# 6. Spearman correlation: breast and lung separately --------------------------
# Cancer type is taken from the CCLE name (for example MCF7_BREAST).
cancer_data <- combined_all %>%
  mutate(cancer_type = case_when(
    str_detect(str_to_upper(ccle_name), "BREAST") ~ "Breast",
    str_detect(str_to_upper(ccle_name), "LUNG")   ~ "Lung",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(cancer_type))

print(count(cancer_data, cancer_type))

cor_by_type <- cancer_data %>%
  group_by(cancer_type) %>%
  group_modify(~ bind_rows(cor_one(.x, "VEGFA"), cor_one(.x, "PCNA"))) %>%
  ungroup() %>%
  mutate(Adjusted_P = p.adjust(P_value, method = "BH"))
print(cor_by_type)
write_csv(cor_by_type, "results_correlations_by_cancer_type.csv")

# 7. Plots ---------------------------------------------------------------------
plot_dox <- function(df, gene, title, facet = FALSE) {
  p <- ggplot(df, aes(.data[[gene]], ic50)) +
    geom_point(alpha = 0.55) +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
    theme_minimal() +
    labs(title = title,
         x = paste(gene, "expression (DepMap 25Q2)"),
         y = "Doxorubicin IC50 (PRISM 20Q2)")
  if (facet) p <- p + facet_wrap(~ cancer_type, scales = "free")
  p
}

ggsave("figures/VEGFA_vs_doxorubicin_IC50.png",
       plot_dox(combined_all, "VEGFA", "VEGFA expression vs doxorubicin IC50"),
       width = 8, height = 6, dpi = 300)
ggsave("figures/PCNA_vs_doxorubicin_IC50.png",
       plot_dox(combined_all, "PCNA", "PCNA expression vs doxorubicin IC50"),
       width = 8, height = 6, dpi = 300)
ggsave("figures/VEGFA_doxorubicin_by_cancer_type.png",
       plot_dox(cancer_data, "VEGFA", "VEGFA vs doxorubicin IC50, breast and lung", TRUE),
       width = 10, height = 5, dpi = 300)
ggsave("figures/PCNA_doxorubicin_by_cancer_type.png",
       plot_dox(cancer_data, "PCNA", "PCNA vs doxorubicin IC50, breast and lung", TRUE),
       width = 10, height = 5, dpi = 300)
