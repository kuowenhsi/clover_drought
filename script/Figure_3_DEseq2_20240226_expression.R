# ============================================================
# Figure 3: subgenome PCA and selected gene expression plots
# Repository root: clover_drought
#
# Inputs are read from ./data/RNAseq/
# Outputs are written to ./figures/Figure3/
setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/Drought_F3_paper/clover_drought")
# ============================================================

library(tidyverse)
library(DESeq2)
library(ggh4x)

# ----------------------------
# Paths
# ----------------------------

rna_dir <- "./data/RNAseq"
orthogroup_file <- file.path(rna_dir, "Orthogroups", "Orthogroups.tsv")
counts_file <- file.path(rna_dir, "map_hap1.5_exon_20240228.txt")
output_dir <- file.path("figures", "Figure3")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_files <- c(orthogroup_file, counts_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required repository files:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

# ----------------------------
# Helper functions
# ----------------------------

read_orthogroups <- function(path) {
  ortho_data <- read_tsv(path, show_col_types = FALSE) %>%
    mutate(
      Arabidopsis = str_split(Arabidopsis, ", "),
      Hap1_v1.5_protein_chr16 = str_split(Hap1_v1.5_protein_chr16, ", "),
      Hap1_v1.5_protein_chr8 = str_split(Hap1_v1.5_protein_chr8, ", ")
    ) %>%
    filter(!is.na(Arabidopsis)) %>%
    dplyr::rename(
      occidentale = Hap1_v1.5_protein_chr8,
      pallescens = Hap1_v1.5_protein_chr16
    ) %>%
    mutate(
      Arabidopsis = lapply(Arabidopsis, function(x) str_split_i(x, "\\|", 1)),
      pallescens = lapply(pallescens, function(x) str_split_i(x, "\\.", 1)),
      occidentale = lapply(occidentale, function(x) str_split_i(x, "\\.", 1)),
      pallescens = lapply(pallescens, unique),
      occidentale = lapply(occidentale, unique)
    )
  
  ortho_data_unnested <- ortho_data %>%
    unnest(Arabidopsis) %>%
    unnest(occidentale) %>%
    unnest(pallescens) %>%
    filter(!(is.na(occidentale) & is.na(pallescens))) %>%
    mutate(Arabidopsis = str_split_i(Arabidopsis, "[.]", 1))
  
  list(
    ortho_data = ortho_data,
    ortho_data_unnested = ortho_data_unnested
  )
}

read_counts_matrix <- function(path) {
  cts <- read_tsv(path, comment = "#", show_col_types = FALSE) %>%
    dplyr::select(-c(2:6))
  
  sample_names <- basename(colnames(cts)[2:ncol(cts)]) %>%
    str_remove("_hap1.5_2passBasicAligned.sortedByCoord.out.bam") %>%
    str_remove("^map_")
  
  cts_matrix <- as.matrix(cts[, 2:ncol(cts)])
  rownames(cts_matrix) <- cts$Geneid
  colnames(cts_matrix) <- sample_names
  
  cts_matrix
}

make_subgenome_pca_plot <- function(expr_matrix, coldata, title_expr, show_legend = TRUE) {
  pca_data <- prcomp(expr_matrix)
  pca_scores <- as.data.frame(pca_data$x)
  
  variance_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2) * 100
  
  pca_scores <- pca_scores %>%
    mutate(
      genotype = factor(
        coldata$genotype,
        levels = c("DMN010", "STL0701", "GFL007"),
        labels = c("DMN_010", "STL_0701", "GFL_007")
      ),
      treatment = factor(coldata$treat, levels = c("control", "drought"))
    )
  
  p <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
    geom_point(
      aes(shape = genotype, fill = interaction(genotype, treatment)),
      size = 3,
      stroke = 0.2,
      show.legend = show_legend
    ) +
    scale_shape_manual(values = c(21, 22, 24)) +
    scale_fill_manual(values = c(
      "#EE0000CC", "#0000FFCC", "#008B00CC",
      "#EE000030", "#0000FF30", "#008B0030"
    )) +
    ggtitle(title_expr) +
    labs(
      x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
      y = paste0("PC2 (", round(variance_explained[2], 2), "%)")
    ) +
    theme_bw() +
    theme(
      legend.position = if (show_legend) c(0.85, 0.26) else "none",
      legend.title = element_blank(),
      legend.spacing = unit(0, "in"),
      legend.background = element_rect(fill = NA),
      plot.title = element_text(hjust = 0.5)
    )
  
  p
}

make_two_gene_plot <- function(dds_obj, gene_a, gene_b, title_expr, normalized = FALSE,show_legend = FALSE) {
  strip_colors <- c("Gene_A" = "darkolivegreen3", "Gene_B" = "darkslategray2")
  
  single_gene_data <- as_tibble(
    counts(dds_obj[c(gene_a, gene_b), ], normalized = normalized),
    rownames = "Gene"
  ) %>%
    pivot_longer(
      cols = 2:ncol(.),
      names_to = c("genotype", "treat", "rep"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    mutate(
      Gene = factor(
        Gene,
        levels = c(gene_a, gene_b),
        labels = c(
          str_replace(gene_a, "drTriRepe4", "occidentale - "),
          str_replace(gene_b, "drTriRepe4", "pallescens - ")
        )
      ),
      genotype = factor(
        genotype,
        levels = c("DMN010", "STL0701", "GFL007"),
        labels = c("DMN_010", "STL_0701", "GFL_007")
      )
    )
  
  ggplot(single_gene_data, aes(x = genotype, y = value)) +
    geom_point(
      aes(fill = interaction(genotype, treat), shape = genotype),
      stroke = 0.2,
      size = 3,
      position = position_dodge(width = 0.5),
      show.legend = show_legend
    ) +
    facet_wrap2(
      . ~ Gene,
      strip = strip_themed(
        background_x = list(
          element_rect(fill = strip_colors[["Gene_A"]]),
          element_rect(fill = strip_colors[["Gene_B"]])
        )
      )
    ) +
    scale_y_continuous("Normalized expression counts") +
    scale_x_discrete("") +
    scale_fill_manual(values = c(
      "#EE000090", "#0000FF90", "#008B0090",
      "#EE000030", "#0000FF30", "#008B0030"
    )) +
    scale_shape_manual(values = c(21, 22, 24)) +
    ggtitle(title_expr) +
    theme_bw() +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8),
      title = element_text(size = 9)
    )
}

# ----------------------------
# Orthogroups and counts
# ----------------------------

ortho_list <- read_orthogroups(orthogroup_file)
ortho_data <- ortho_list$ortho_data
ortho_data_unnested <- ortho_list$ortho_data_unnested

cts_matrix <- read_counts_matrix(counts_file)

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
# DESeq2
# ----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = cts_matrix,
  colData = sample_info,
  design = ~ genotype + treat + genotype:treat
)

dds_normalized <- DESeq(dds)
rld <- rlogTransformation(dds_normalized)

# ----------------------------
# PCA by subgenome
# ----------------------------

total_exp <- t(assay(rld))

column_to <- grep("^drTriRepe4Chr[1-8]g", colnames(total_exp))
column_tp <- grep("^drTriRepe4Chr(9|10|11|12|13|14|15|16)g", colnames(total_exp))

sub_to_exp <- total_exp[, column_to, drop = FALSE]
sub_tp_exp <- total_exp[, column_tp, drop = FALSE]

p_fig3A <- make_subgenome_pca_plot(
  expr_matrix = sub_to_exp,
  coldata = colData(rld),
  title_expr = expression(italic("T. occidentale") ~ "subgenome"),
  show_legend = TRUE
)

p_fig3B <- make_subgenome_pca_plot(
  expr_matrix = sub_tp_exp,
  coldata = colData(rld),
  title_expr = expression(italic("T. pallescens") ~ "subgenome"),
  show_legend = FALSE
)

p_fig3A
p_fig3B

ggsave(
  file.path(output_dir, "Fig_3A_T_occidentale_PCA.png"),
  plot = p_fig3A,
  width = 3.5,
  height = 3.5,
  dpi = 600
)

ggsave(
  file.path(output_dir, "Fig_3B_T_pallescens_PCA.png"),
  plot = p_fig3B,
  width = 3.5,
  height = 3.5,
  dpi = 600
)

# ----------------------------
# Fig. 3E: proline transporter ortholog pair
# Arabidopsis gene: AT2G36590
# ----------------------------

proline_pair <- ortho_data_unnested %>%
  filter(Arabidopsis == "AT2G36590") %>%
  distinct(Arabidopsis, occidentale, pallescens) %>%
  drop_na()

if (nrow(proline_pair) == 0) {
  stop("Could not find an ortholog pair for AT2G36590 in Orthogroups.tsv")
}

proline_gene_a <- proline_pair$occidentale[1]
proline_gene_b <- proline_pair$pallescens[1]

p_fig3E <- make_two_gene_plot(
  dds_obj = dds_normalized,
  gene_a = proline_gene_a,
  gene_b = proline_gene_b,
  title_expr = expression(italic("AT2G36590") ~ "(proline transporter)"),
  normalized = TRUE,
  show_legend = FALSE
)

p_fig3E

ggsave(
  file.path(output_dir, "Fig_3E_AT2G35690_proline_transporter.png"),
  plot = p_fig3E,
  width = 3.5,
  height = 3.5,
  dpi = 600
)

# ----------------------------
# Fig. 3F: beta-CAS ortholog pair
# ----------------------------

beta_cas_gene_a <- "drTriRepe4Chr7g205700"
beta_cas_gene_b <- "drTriRepe4Chr15g223500"

p_fig3F <- make_two_gene_plot(
  dds_obj = dds_normalized,
  gene_a = beta_cas_gene_a,
  gene_b = beta_cas_gene_b,
  title_expr = expression(beta ~ italic("CAS") ~ "genes"),
  normalized = TRUE,
  show_legend = FALSE
)

p_fig3F


ggsave(
  file.path(output_dir, "Fig_3F_beta_CAS.png"),
  plot = p_fig3F,
  width = 3.5,
  height = 3.5,
  dpi = 600
)

