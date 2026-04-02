# ============================================================
# Figure 3C and 3D: genotype overlap of DE genes and
# subgenome-specific counts of drought-responsive genes
#
# Repository root: clover_drought
#
# Inputs are read from .data/RNAseq/
# Outputs are written to ./figures/Figure3/
setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/Drought_F3_paper/clover_drought")
# ============================================================

library(tidyverse)
library(DESeq2)
library(qvalue)
library(cowplot)
library(ggVennDiagram)

# ----------------------------
# Paths
# ----------------------------

rna_dir <- "./data/RNAseq"
counts_file <- file.path(rna_dir, "map_hap1.5_exon_20240228.txt")
output_dir <- file.path("figures", "Figure3")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(counts_file)) {
  stop("Missing required repository file: RNAseq/map_hap1.5_exon_20240228.txt")
}

# ----------------------------
# Read count matrix
# ----------------------------

cts <- read_tsv(counts_file, comment = "#", show_col_types = FALSE) %>%
  dplyr::select(-c(2:6))

sample_names <- basename(colnames(cts)[2:ncol(cts)]) %>%
  str_remove("_hap1.5_2passBasicAligned.sortedByCoord.out.bam") %>%
  str_remove("^map_")

cts_matrix <- as.matrix(cts[, 2:ncol(cts)])
rownames(cts_matrix) <- cts$Geneid
colnames(cts_matrix) <- sample_names

# ----------------------------
# Sample metadata
# ----------------------------

sample_info <- data.frame(
  genotype = rep(c("DMN010", "GFL007", "STL0701"), each = 6),
  treat = rep(rep(c("control", "drought"), each = 3), times = 3),
  row.names = colnames(cts_matrix)
) %>%
  mutate(
    genotype = factor(genotype, levels = c("DMN010", "GFL007", "STL0701")),
    treat = factor(treat, levels = c("control", "drought"))
  )

# ----------------------------
# Helper function:
# run DESeq2 for one genotype and return significant drought-response genes
# ----------------------------

run_treatment_de <- function(count_mat, coldata, genotype_label) {
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~ treat
  )
  
  dds <- DESeq(dds)
  
  res <- results(dds, name = "treat_drought_vs_control", alpha = 0.05)
  
  res_tbl <- as_tibble(res, rownames = "Gene") %>%
    mutate(
      subgenome = case_when(
        str_detect(Gene, "drTriRepe4Chr[1-8]g") ~ "occidentale",
        str_detect(Gene, "drTriRepe4Chr9g") | str_detect(Gene, "drTriRepe4Chr1[0-6]g") ~ "pallescens",
        TRUE ~ NA_character_
      ),
      regulation = case_when(
        log2FoldChange < 0 ~ "Negative",
        log2FoldChange >= 0 ~ "Positive",
        TRUE ~ NA_character_
      ),
      genotype = genotype_label
    )
  
  res_tbl %>%
    filter(padj < 0.05)
}

# ----------------------------
# Differential expression by genotype
# ----------------------------

res_DMN_treat_sig <- run_treatment_de(
  count_mat = cts_matrix[, 1:6],
  coldata = sample_info[1:6, , drop = FALSE],
  genotype_label = "DMN010"
)

res_GFL_treat_sig <- run_treatment_de(
  count_mat = cts_matrix[, 7:12],
  coldata = sample_info[7:12, , drop = FALSE],
  genotype_label = "GFL007"
)

res_STL_treat_sig <- run_treatment_de(
  count_mat = cts_matrix[, 13:18],
  coldata = sample_info[13:18, , drop = FALSE],
  genotype_label = "STL0701"
)

# ----------------------------
# Fig. 3D:
# Subgenome-specific counts of significant up/downregulated genes
# ----------------------------

res_subgenome <- bind_rows(
  res_DMN_treat_sig,
  res_GFL_treat_sig,
  res_STL_treat_sig
) %>%
  group_by(genotype, subgenome, regulation) %>%
  summarise(count = n(), .groups = "drop") %>%
  drop_na() %>%
  mutate(
    count = if_else(regulation == "Negative", -count, count),
    genotype = factor(
      genotype,
      levels = c("DMN010", "STL0701", "GFL007"),
      labels = c("DMN_010", "STL_0701", "GFL_007")
    )
  )

p_fig3D <- ggplot(
  data = res_subgenome,
  aes(x = genotype, y = count, fill = subgenome, group = subgenome)
) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "gray") +
  annotate(
    geom = "text",
    x = 0.5,
    y = 2850,
    hjust = 0,
    label = "Fisher’s Exact Test\nNot significant in all comparisons"
  ) +
  scale_fill_manual(values = c("darkolivegreen3", "darkslategray2")) +
  scale_y_continuous(
    breaks = c(-2000, -1000, 0, 1000, 2000),
    limits = c(-2600, 3000),
    labels = function(x) abs(x)
  ) +
  labs(y = "Downregulated genes      Upregulated genes", x = "") +
  theme_bw() +
  theme(legend.position = c(0.18, 0.2), legend.background = element_rect(fill = NA))

p_fig3D

ggsave(
  filename = file.path(output_dir, "Fig_3D_subgenome_dominance_treat.png"),
  plot = p_fig3D,
  width = 3.5,
  height = 3.5,
  dpi = 600
)

# ----------------------------
# Fisher's exact tests reported in original script
# ----------------------------

total_occidentale <- 48107
total_pallescens <- 48124

fisher_results <- res_subgenome %>%
  mutate(count = abs(count)) %>%
  group_by(genotype, regulation) %>%
  summarise(
    fisher_p_value = {
      sub_occ <- count[subgenome == "occidentale"]
      sub_pall <- count[subgenome == "pallescens"]
      
      non_expressed_occ <- total_occidentale - sub_occ
      non_expressed_pall <- total_pallescens - sub_pall
      
      contingency_table <- matrix(
        c(sub_occ, non_expressed_occ, sub_pall, non_expressed_pall),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
          Subgenome = c("occidentale", "pallescens"),
          Expression = c("Expressed", "Not Expressed")
        )
      )
      
      fisher.test(contingency_table)$p.value
    },
    .groups = "drop"
  )

print(fisher_results)

# ----------------------------
# Fig. 3C:
# Venn diagram of significant drought-responsive genes
# ----------------------------

venn_data <- list(
  DMN_010 = res_DMN_treat_sig$Gene,
  GFL_007 = res_GFL_treat_sig$Gene,
  STL_0701 = res_STL_treat_sig$Gene
)

p_fig3C <- ggVennDiagram(
  venn_data,
  label_alpha = 0,
  label_size = 3,
  set_size = 3,
  edge_size = 0.5
) +
  scale_fill_gradient(low = "grey90", high = "red") +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  scale_x_continuous(expand = expansion(mult = 0.1))

p_fig3C

ggsave(
  filename = file.path(output_dir, "Fig_3C_Venn_DE_genotypes.png"),
  plot = p_fig3C,
  width = 3.5,
  height = 3.5,
  dpi = 600
)
