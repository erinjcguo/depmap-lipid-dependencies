# ─────────────────────────────────────────────────────────────
# Lipid-metabolism gene dependencies in liver, colorectal and
# pancreatic cancer cell lines
#
# Data: DepMap Public 26Q3 CRISPR gene effect (Chronos)
# A gene effect score near 0 means knocking the gene out did
# little. A score near -1 means knocking it out hurt the cells
# about as much as removing a gene every cell needs.
# Convention: a score below -0.5 = the line "depends" on that gene.
# ─────────────────────────────────────────────────────────────

library(tidyverse)

# 1. Read the file downloaded from DepMap Custom Downloads
dat <- read_csv("CRISPR_DepMap_Public_26Q3Score_Chronos_subsetted.csv")

genes <- c("SCD", "FASN", "ACACA", "SREBF1", "SREBF2", "HMGCR")

# 2. Keep only the three cancer types, and reshape so there is
#    one row per cell line per gene (easier to plot and summarise)
tidy <- dat %>%
  filter(lineage_1 %in% c("Liver", "Bowel", "Pancreas")) %>%
  select(depmap_id, cell_line_display_name, lineage_1, all_of(genes)) %>%
  pivot_longer(all_of(genes), names_to = "gene", values_to = "gene_effect") %>%
  drop_na(gene_effect)

# 3. Summary table: how many lines in each group depend on each gene
summary_tbl <- tidy %>%
  group_by(gene, lineage_1) %>%
  summarise(
    n            = n(),
    median_score = round(median(gene_effect), 3),
    n_dependent  = sum(gene_effect < -0.5),
    pct_dependent = round(100 * mean(gene_effect < -0.5), 1),
    .groups = "drop"
  ) %>%
  arrange(gene, lineage_1)

print(summary_tbl)
write_csv(summary_tbl, "results_summary.csv")

# 4. Figure 1 — distribution of scores, by gene and cancer type
p1 <- ggplot(tidy, aes(x = lineage_1, y = gene_effect, fill = lineage_1)) +
  geom_hline(yintercept = -0.5, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  facet_wrap(~ gene) +
  labs(
    title    = "Lipid-metabolism gene dependencies across three cancer types",
    subtitle = "DepMap Public 26Q3 CRISPR gene effect (Chronos). Dashed line = dependency threshold (-0.5)",
    x = NULL, y = "Gene effect score"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

ggsave("fig1_boxplots.png", p1, width = 9, height = 6, dpi = 300)

# 5. Figure 2 — median score per gene and cancer type
p2 <- summary_tbl %>%
  ggplot(aes(x = lineage_1, y = gene, fill = median_score)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = median_score), size = 3.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Median\ngene effect") +
  labs(
    title = "Median gene effect by cancer type",
    subtitle = "More negative (blue) = stronger dependency",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11)

ggsave("fig2_heatmap.png", p2, width = 6, height = 4.5, dpi = 300)

# 6. Is the difference between groups bigger than chance?
#    Kruskal-Wallis is a rank-based test that does not assume
#    the scores are normally distributed.
tests <- tidy %>%
  group_by(gene) %>%
  summarise(p_value = kruskal.test(gene_effect ~ lineage_1)$p.value, .groups = "drop") %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))

print(tests)
write_csv(tests, "results_tests.csv")
