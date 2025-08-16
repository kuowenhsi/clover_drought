# Load required libraries
library(tidyverse)
library(DESeq2)
library(qvalue)
library(cowplot)
library(biomaRt)

# Set working directory
setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/RNAseq")

#--------------------------------------------------
# Helper functions to identify orthologs
#--------------------------------------------------
get_Arabidopsis_gene <- function(x) sapply(ortho_data$Arabidopsis, \(y) any(grepl(x, y)))
get_pallescens_gene  <- function(x) sapply(ortho_data$pallescens,  \(y) any(grepl(x, y)))
get_occidentale_gene <- function(x) sapply(ortho_data$occidentale, \(y) any(grepl(x, y)))

analyze_contrast <- function(contrast_name) {
  # Run DESeq2 results
  res <- results(dds_normalized, name = contrast_name, alpha = 0.05)
  
  # Convert to tibble
  res_treat <- as_tibble(res, rownames = "Gene")
  
  # Diagnostic plots (optional)
  print(sum(res_treat$padj < 0.05, na.rm = TRUE) / nrow(res_treat))
  hist(res_treat$padj, main = paste("padj histogram for", contrast_name))
  hist(res_treat$pvalue, main = paste("pvalue histogram for", contrast_name))
  
  # Trifolium occidentale-specific genes
  res_treat_to <- res_treat %>%
    filter(padj < 0.05) %>%
    filter(str_detect(Gene, "drTriRepe4Chr[1-8]g")) %>%
    left_join(ortho_data_single, by = c("Gene" = "occidentale")) %>%
    drop_na()
  
  # Trifolium pallescens-specific genes
  res_treat_tp <- res_treat %>%
    filter(padj < 0.05) %>%
    filter(str_detect(Gene, "drTriRepe4Chr9g") | str_detect(Gene, "drTriRepe4Chr1[0-6]g")) %>%
    left_join(ortho_data_single, by = c("Gene" = "pallescens")) %>%
    drop_na()
  
  # List of genes of interest
  interedted_Gene_to <- c(res_treat_to$Gene, res_treat_to$pallescens)
  interedted_Gene_tp <- c(res_treat_tp$Gene, res_treat_tp$occidentale)
  
  # Final annotated and colored tables
  res_treat_to_final <- as_tibble(res, rownames = "Gene") %>%
    filter(Gene %in% interedted_Gene_to) %>%
    left_join(ortho_data_single_long, by = "Gene") %>%
    filter(!is.na(Orthogroup)) %>%
    mutate(point_color = case_when(
      (padj < 0.05) & (log2FoldChange > 0) ~ "red",
      (padj < 0.05) & (log2FoldChange < 0) ~ "blue",
      TRUE ~ "#11111111"
    )) %>%
    group_by(Orthogroup) %>%
    mutate(line_color = case_when(
      all(point_color != "#11111111") & (prod(log2FoldChange) > 0) ~ "green4",
      all(point_color != "#11111111") & (prod(log2FoldChange) < 0) ~ "purple2",
      TRUE ~ "#11111111"
    )) %>%
    ungroup() %>%
    mutate(
      point_color = factor(point_color, levels = c("#11111111", "red", "blue")),
      line_color = factor(line_color, levels = c("#11111111", "green4", "purple2"))
    ) %>%
    filter(!is.na(Orthogroup))%>%
    arrange(Orthogroup)
  
  res_treat_tp_final <- as_tibble(res, rownames = "Gene") %>%
    filter(Gene %in% interedted_Gene_tp) %>%
    left_join(ortho_data_single_long, by = "Gene") %>%
    filter(!is.na(Orthogroup)) %>%
    mutate(point_color = case_when(
      (padj < 0.05) & (log2FoldChange > 0) ~ "red",
      (padj < 0.05) & (log2FoldChange < 0) ~ "blue",
      TRUE ~ "#11111111"
    )) %>%
    group_by(Orthogroup) %>%
    mutate(line_color = case_when(
      all(point_color != "#11111111") & (prod(log2FoldChange) > 0) ~ "green4",
      all(point_color != "#11111111") & (prod(log2FoldChange) < 0) ~ "purple2",
      TRUE ~ "#11111111"
    )) %>%
    ungroup() %>%
    mutate(
      point_color = factor(point_color, levels = c("#11111111", "red", "blue")),
      line_color = factor(line_color, levels = c("#11111111", "green4", "purple2"))
    ) %>%
    filter(!is.na(Orthogroup))%>%
    arrange(Orthogroup)
  
  return(list(to = res_treat_to_final, tp = res_treat_tp_final))
}

plot_homeolog_expression <- function(data,
                                     y_annot_occ = 5,
                                     y_annot_pal = 4,
                                     y_annot_note = 8,
                                     note_text = "",
                                     y_label_text = "(Drought / Control)",
                                     row_filter = NULL) {
  # Subset the data if row_filter is provided
  plot_data <- if (!is.null(row_filter)) data[row_filter, ] else data
  
  ggplot(plot_data, aes(x = Subgenome, y = log2FoldChange)) +
    geom_line(aes(group = Orthogroup), color = plot_data$line_color) +
    geom_point(color = plot_data$point_color) +
    annotate("text", x = "occidentale", y = y_annot_occ,
             label = paste0(sum(data$Subgenome == "occidentale" & data$padj < 0.05, na.rm = TRUE), "/", sum(data$Subgenome == "occidentale"))) +
    annotate("text", x = "pallescens", y = y_annot_pal,
             label = paste0(sum(data$Subgenome == "pallescens" & data$padj < 0.05, na.rm = TRUE), "/", sum(data$Subgenome == "pallescens"))) +
    annotate("text", x = 0.5, y = y_annot_note, label = note_text, hjust = 0, vjust = 1, size = 3) +
    ylab(parse(text = sprintf("log[2]~%s", y_label_text)))+
    xlab("") +
    scale_x_discrete(labels = c("occidentale" = "Occid.", "pallescens" = "Palle.")) +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 10))
}


#--------------------------------------------------
# Load and process orthogroup data
#--------------------------------------------------
ortho_data <- read_tsv("./Orthogroups/Orthogroups.tsv") %>%
  mutate(across(c(Arabidopsis, Hap1_v1.5_protein_chr16, Hap1_v1.5_protein_chr8), str_split, ", ")) %>%
  filter(!is.na(Arabidopsis)) %>%
  dplyr::rename(occidentale = Hap1_v1.5_protein_chr8,
         pallescens  = Hap1_v1.5_protein_chr16) %>%
  mutate(
    Arabidopsis = map(Arabidopsis, ~ str_split_i(.x, "\\|", 1)),
    pallescens  = map(pallescens,  ~ unique(str_split_i(.x, "\\.", 1))),
    occidentale = map(occidentale, ~ unique(str_split_i(.x, "\\.", 1)))
  )

# Long orthogroup format for many-to-many mappings
ortho_data_unnested <- ortho_data %>%
  unnest(Arabidopsis)%>%
  unnest(occidentale)%>%
  unnest(pallescens)%>%
  filter(!(is.na(occidentale) & is.na(pallescens)))%>%
  mutate(Arabidopsis = str_split_i(Arabidopsis, "[.]", 1))

# One-to-one orthologs only
ortho_data_single <- ortho_data %>%
  dplyr::select(-Arabidopsis) %>%
  filter(map_int(occidentale, length) == 1,
         map_int(pallescens,  length) == 1) %>%
  mutate(across(c(occidentale, pallescens), unlist)) %>%
  drop_na()

ortho_data_single_long <- ortho_data_single %>%
  pivot_longer(cols = c(occidentale, pallescens), names_to = "Subgenome", values_to = "Gene")

# Many-to-many orthologs
ortho_data_multi <- ortho_data %>%
  dplyr::select(-Arabidopsis)%>%
  unnest(occidentale)%>%
  unnest(pallescens)%>%
  drop_na()

ortho_data_multi_long <- ortho_data_multi %>%
  pivot_longer(cols = c(occidentale, pallescens), names_to = "Subgenome", values_to = "Gene")

#--------------------------------------------------
# Load and process count matrix
#--------------------------------------------------
cts <- read_tsv("map_hap1.5_exon_20240228.txt", comment = "#") %>%
  dplyr::select(-c(2:6))

# Extract sample names
sample_names <- colnames(cts)[-1] %>%
  str_split_i("/", 9) %>%
  str_remove("_hap1.5_2passBasicAligned.sortedByCoord.out.bam") %>%
  str_remove("map_") %>%
  na.omit()

# Create expression matrix
cts_matrix <- as.matrix(cts[, -1])
rownames(cts_matrix) <- cts$Geneid
colnames(cts_matrix) <- sample_names

#--------------------------------------------------
# Create sample metadata and DESeq2 object
#--------------------------------------------------
sample_info <- data.frame(
  genotype = rep(c("DMN010", "GFL007", "STL0701"), each = 6),
  treat    = rep(rep(c("control", "drought"), each = 3), times = 3),
  row.names = sample_names
) %>%
  mutate(
    genotype = factor(genotype, levels = c("DMN010", "GFL007", "STL0701")),
    treat    = factor(treat, levels = c("control", "drought"))
  )

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = cts_matrix,
  colData   = sample_info,
  design    = ~ genotype + treat + genotype:treat
)

# Run DESeq2 normalization and model fitting
dds_normalized <- DESeq(dds)
resultsNames(dds_normalized)


#--------------------------------------------------
# Define the contrast names (excluding "Intercept")
#--------------------------------------------------
contrast_names <- resultsNames(dds_normalized)[-1]

# Create lists to store results
res_to_list <- list()
res_tp_list <- list()

# Loop through each contrast and apply the function
for (contrast in contrast_names) {
  cat("Analyzing contrast:", contrast, "\n")
  result <- analyze_contrast(contrast)
  res_to_list[[contrast]] <- result$to
  res_tp_list[[contrast]] <- result$tp
}


#-------------------------------------------------
# Prepare for plotting
#-------------------------------------------------
# === 1. Drought Treatment Main Effect ===
res_treat_to <- res_to_list[["treat_drought_vs_control"]]
res_treat_tp <- res_tp_list[["treat_drought_vs_control"]]

# Occidentale (To) perspective
p1 <- plot_homeolog_expression(
  data = res_treat_to,
  y_annot_occ = 5,
  y_annot_pal = 4,
  y_annot_note = 8,
  note_text = "Among their homeologous copies in the\nPalle. subgenome, 40 out of 85 also show DE."
)
p1
ggsave("Homeolog_expression_Control_Drought_from_To.png", plot = p1, width = 3.5, height = 3.5, dpi = 600)

# Pallescens (Tp) perspective
p2 <- plot_homeolog_expression(
  data = res_treat_tp,
  y_annot_occ = 5,
  y_annot_pal = 4.5,
  y_annot_note = 8,
  note_text = "There are 85 DE genes in Palle. subgenome.\nTheir homeologous copies in Occid. subgenome,\n40 out of 85, also show DE."
)
p2
ggsave("Homeolog_expression_Control_Drought_from_Tp.png", plot = p2, width = 3.5, height = 3.5, dpi = 600)

# === 2. Genotype: GFL007 vs DMN010 ===
res_gfl_to <- res_to_list[["genotype_GFL007_vs_DMN010"]]
res_gfl_tp <- res_tp_list[["genotype_GFL007_vs_DMN010"]]

p3 <- plot_homeolog_expression(
  data = res_gfl_to,
  y_annot_occ = 10.5,
  y_annot_pal = 7.5,
  y_annot_note = 8,
  note_text = "",
  y_label_text = "(GFL_007 / DMN_010)",
  row_filter = 1:300
)
p3
ggsave("Homeolog_expression_genotype_GFL007_vs_DMN010_from_To.png", plot = p3, width = 2.5, height = 2.5, dpi = 600)

p4 <- plot_homeolog_expression(
  data = res_gfl_tp,
  y_annot_occ = 6,
  y_annot_pal = 12.5,
  note_text = "",
  y_label_text = "(GFL_007 / DMN_010)",
  row_filter = 1:300
)
p4
ggsave("Homeolog_expression_genotype_GFL007_vs_DMN010_from_Tp.png", plot = p4, width = 2.5, height = 2.5, dpi = 600)

# === 3. Genotype: STL0701 vs DMN010 ===
res_stl_to <- res_to_list[["genotype_STL0701_vs_DMN010"]]
res_stl_tp <- res_tp_list[["genotype_STL0701_vs_DMN010"]]

p5 <- plot_homeolog_expression(
  data = res_stl_to,
  y_annot_occ = 9.5,
  y_annot_pal = 9,
  y_annot_note = 8,
  note_text = "",
  y_label_text = "(STL_0701 / DMN_010)",
  row_filter = 1:300
)
p5

ggsave("Homeolog_expression_genotype_STL0701_vs_DMN010_from_To.png", plot = p5, width = 2.5, height = 2.5, dpi = 600)

p6 <- plot_homeolog_expression(
  data = res_stl_tp,
  y_annot_occ = 6,
  y_annot_pal = 9,
  y_annot_note = 8,
  note_text = "",
  y_label_text = "(STL_0701 / DMN_010)",
  row_filter = 1:300
)
p6

ggsave("Homeolog_expression_genotype_STL0701_vs_DMN010_from_Tp.png", plot = p6, width = 2.5, height = 2.5, dpi = 600)

# === 4. Interaction: GFL007 × Drought ===
res_interact_gfl_to <- res_to_list[["genotypeGFL007.treatdrought"]]
res_interact_gfl_tp <- res_tp_list[["genotypeGFL007.treatdrought"]]

p7 <- plot_homeolog_expression(
  data = res_interact_gfl_to,
  y_annot_occ = 8,
  y_annot_pal = 5,
  y_annot_note = 8,
  note_text = "",
  y_label_text = "(GFL_007 / DMN_010)"
)
p7
ggsave("Homeolog_expression_genotypeGFL007.treatdrought_from_To.png", plot = p7, width = 2.5, height = 2.5, dpi = 600)

p8 <- plot_homeolog_expression(
  data = res_interact_gfl_tp,
  y_annot_occ = 8,
  y_annot_pal = 7,
  y_annot_note = 8,
  note_text = "",
  y_label_text = "(GFL_007 / DMN_010)"
)
p8
ggsave("Homeolog_expression_genotypeGFL007.treatdrought_from_Tp.png", plot = p8, width = 2.5, height = 2.5, dpi = 600)

# === 5. Interaction: STL0701 × Drought ===
res_interact_stl_to <- res_to_list[["genotypeSTL0701.treatdrought"]]
res_interact_stl_tp <- res_tp_list[["genotypeSTL0701.treatdrought"]]

p9 <- plot_homeolog_expression(
  data = res_interact_stl_to,
  y_annot_occ = 9,
  y_annot_pal = 8,
  y_annot_note = 8,
  note_text = "",
  y_label_text = "(STL_0701 / DMN_010)"
)
p9

ggsave("Homeolog_expression_genotypeSTL0701.treatdrought_from_To.png", plot = p9, width = 2.5, height = 2.5, dpi = 600)

p10 <- plot_homeolog_expression(
  data = res_interact_stl_tp,
  y_annot_occ = 9,
  y_annot_pal = 11,
  y_annot_note = 8,
  note_text = "",
  y_label_text = "(STL_0701 / DMN_010)"
)
p10
ggsave("Homeolog_expression_genotypeSTL0701.treatdrought_from_Tp.png", plot = p10, width = 2.5, height = 2.5, dpi = 600)

