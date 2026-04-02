# ============================================================
# Figure 6A-C: Chr4 QTL summary and allele-effect plots
#
# Repository root: clover_drought
#
# Inputs are read from ./data/
# Outputs are written to ./figures/Figure6/
setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/Drought_F3_paper/clover_drought")
# ============================================================

library(tidyverse)
library(qtl)
library(car)

# ----------------------------
# Paths
# ----------------------------

data_dir <- "data"
output_dir <- file.path("figures", "Figure6")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_files <- c(
  file.path(data_dir, "leaf_area_LOD_tb_20240303.csv"),
  file.path(data_dir, "flower_count_LOD_tb_20240303.csv"),
  file.path(data_dir, "dry_weight_LOD_tb_20240303.csv"),
  file.path(data_dir, "bwa_DG_F3_DP0_F3_m0.5_hardfilter_pre_imputed_homo_maf0.2_hww0.01_LM.csv"),
  file.path(data_dir, "DG_F3_physical_order_LM_pheno_20240219.csv"),
  file.path(data_dir, "dry_weight_log_BLUPs_20240826.csv"),
  file.path(data_dir, "flower_count_BLUPs_20240826.csv"),
  file.path(data_dir, "leaf_area_BLUPs_20240825.csv")
)

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
# Helper function
# ----------------------------
get_effect_plot <- function(cross_object,
                            marker,
                            Trt1,
                            Trt2 = NULL,
                            pos,
                            lod,
                            y_label = "pheno") {
  
  if (is.null(Trt2)) {
    if (str_detect(Trt1, "Control")) {
      Trt2 <- str_replace(Trt1, "Control", "Drought")
    }
    if (str_detect(Trt1, "Drought")) {
      Trt2 <- str_replace(Trt1, "Drought", "Control")
    }
  }
  
  Trt1 <- as.character(Trt1)
  Trt2 <- as.character(Trt2)
  
  if (!Trt1 %in% colnames(cross_object$pheno)) {
    stop("Phenotype column not found: ", Trt1)
  }
  if (!Trt2 %in% colnames(cross_object$pheno)) {
    stop("Phenotype column not found: ", Trt2)
  }
  
  allele_effect_Trt1 <- plotPXG(cross_object, marker, pheno.col = Trt1) %>%
    drop_na() %>%
    mutate(Trt = Trt1)
  
  allele_effect_Trt2 <- plotPXG(cross_object, marker, pheno.col = Trt2) %>%
    drop_na() %>%
    mutate(Trt = Trt2)
  
  allele_effect <- bind_rows(allele_effect_Trt1, allele_effect_Trt2) %>%
    mutate(marker_genotype = get(marker)) %>%
    mutate(
      marker_genotype = factor(marker_genotype, labels = c("AA", "AB", "BB")),
      Trt = factor(Trt, levels = c(Trt1, Trt2))
    ) %>%
    filter(!is.na(marker_genotype), !is.na(pheno))
  
  n_trt <- n_distinct(droplevels(allele_effect$Trt))
  n_geno <- n_distinct(droplevels(allele_effect$marker_genotype))
  
  if (n_trt < 2) {
    stop(
      "Only one treatment level remained after filtering: ",
      paste(unique(allele_effect$Trt), collapse = ", ")
    )
  }
  
  if (n_geno < 2) {
    stop(
      "Only one marker genotype level remained after filtering for marker ",
      marker
    )
  }
  
  model_effect <- lm(pheno ~ marker_genotype * Trt, data = allele_effect)
  ANOVA_test <- car::Anova(model_effect, type = "III")
  
  interaction_row <- grep("marker_genotype:Trt", rownames(ANOVA_test), value = TRUE)[1]
  interaction_p <- ANOVA_test[interaction_row, "Pr(>F)"]
  
  value_str <- format(sprintf("%.2e", interaction_p), scientific = TRUE)
  parts <- strsplit(value_str, "e")[[1]]
  coefficient <- parts[1]
  exponent <- as.integer(parts[2])
  
  expr <- paste(
    "Interaction:",
    deparse(bquote(italic(p) == .(coefficient) %*% 10^.(exponent)))
  )
  
  sample_size <- allele_effect %>%
    count(marker_genotype, name = "n") %>%
    arrange(marker_genotype)
  
  p <- ggplot(data = allele_effect, aes(x = marker_genotype, y = pheno)) +
    geom_point(
      aes(color = Trt),
      position = position_dodge(width = 0.2),
      alpha = 0.5,
      size = 1
    ) +
    stat_summary(
      geom = "line",
      fun = median,
      aes(color = Trt, group = Trt),
      position = position_dodge(width = 0.5),
      show.legend = FALSE
    ) +
    geom_boxplot(
      aes(color = Trt),
      position = position_dodge(width = 0.5),
      width = 0.1,
      show.legend = FALSE,
      outlier.shape = NA,
      linewidth = 0.2
    ) +
    stat_summary(
      geom = "point",
      fun = median,
      aes(color = Trt),
      shape = 18,
      size = 3,
      position = position_dodge(width = 0.5),
      show.legend = FALSE
    ) +
    geom_text(
      x = 0.5,
      y = max(allele_effect$pheno, na.rm = TRUE) * 0.95,
      label = expr,
      parse = TRUE,
      size = 3.5,
      fontface = "plain",
      check_overlap = TRUE,
      hjust = 0
    ) +
    scale_x_discrete(
      name = marker,
      labels = c(
        paste0("AA (", sample_size$n[[1]], ")"),
        paste0("AB (", sample_size$n[[2]], ")"),
        paste0("BB (", sample_size$n[[3]], ")")
      )
    ) +
    scale_y_continuous(oob = scales::oob_squish) +
    scale_color_manual(values = c("#00BFC4", "#F8766D")) +
    ylab(y_label) +
    ggtitle(paste("Pos =", round(pos, 2), "LOD =", round(lod, 2))) +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.background = element_rect(fill = NA)
    )
  
  return(p)
}

# ----------------------------
# Load QTL summary tables
# ----------------------------

trait_code_levels <- c(
  "Genetic distance",
  "flower_count",
  "leaf_area_w_2",
  "leaf_area_w_3",
  "leaf_area_w_4",
  "leaf_area_w_5",
  "shoot_weight",
  "root_weight"
)

trait_label_levels <- c(
  "Genetic distance",
  "Inflorescence",
  "Leaf area w2",
  "Leaf area w3",
  "Leaf area w4",
  "Leaf area w5",
  "Shoot mass",
  "Root mass"
)

total_LOD_tb <- bind_rows(
  read_csv(file.path(data_dir, "leaf_area_LOD_tb_20240303.csv"), show_col_types = FALSE) %>%
    mutate(trait_code = paste("leaf_area", trait, sep = "_")),
  read_csv(file.path(data_dir, "flower_count_LOD_tb_20240303.csv"), show_col_types = FALSE) %>%
    mutate(trait_code = paste("flower_count", trait, sep = "_")),
  read_csv(file.path(data_dir, "dry_weight_LOD_tb_20240303.csv"), show_col_types = FALSE) %>%
    mutate(trait_code = str_replace(trait, "w", "weight"))
) %>%
  mutate(
    trait_code = str_remove(str_remove(trait_code, "_Control"), "_Drought"),
    trait_label = factor(trait_code, levels = trait_code_levels, labels = trait_label_levels)
  )

# ----------------------------
# Load cross objects and phenotype data
# ----------------------------

DG_F3_data <- read.cross(
  format = "csv",
  dir = data_dir,
  file = "bwa_DG_F3_DP0_F3_m0.5_hardfilter_pre_imputed_homo_maf0.2_hww0.01_LM.csv",
  estimate.map = FALSE,
  genotypes = c("AA", "AB", "BB"),
  alleles = c("A", "B")
)

other_data <- read.cross(
  format = "csv",
  dir = data_dir,
  file = "DG_F3_physical_order_LM_pheno_20240219.csv",
  estimate.map = FALSE,
  genotypes = c("AA", "AB", "BB"),
  alleles = c("A", "B")
)

other_pheno <- other_data$pheno %>%
  select(1:4)

combined_blups <- read_csv(file.path(data_dir, "dry_weight_log_BLUPs_20240826.csv"), show_col_types = FALSE) %>%
  left_join(read_csv(file.path(data_dir, "flower_count_BLUPs_20240826.csv"), show_col_types = FALSE), by = "Genotype") %>%
  left_join(read_csv(file.path(data_dir, "leaf_area_BLUPs_20240825.csv"), show_col_types = FALSE), by = "Genotype")

DG_F3_pheno <- DG_F3_data$pheno %>%
  left_join(combined_blups, by = "Genotype") %>%
  left_join(other_pheno, by = "Genotype")

DG_F3_data$pheno <- DG_F3_pheno

# ----------------------------
# Build genetic map table
# ----------------------------

geneticMapList <- list()

for (chr in names(DG_F3_data$geno)) {
  markers <- names(DG_F3_data$geno[[chr]]$map)
  distances <- DG_F3_data$geno[[chr]]$map
  
  geneticMapList[[chr]] <- data.frame(
    Marker = markers,
    GeneticDistance = distances
  )
}

geneticMapList_tb <- bind_rows(geneticMapList, .id = "mid_chr") %>%
  transmute(
    mid_chr = mid_chr,
    mid_pos = as.double(GeneticDistance),
    trait_label = factor("Genetic distance", levels = trait_label_levels)
  )

# ----------------------------
# Restrict to chromosome 4
# ----------------------------

total_LOD_tb_Chr4 <- total_LOD_tb %>%
  filter(mid_chr == "drTriRepe4Chr4")

geneticMapList_tb_Chr4 <- geneticMapList_tb %>%
  filter(mid_chr == "drTriRepe4Chr4")

# ----------------------------
# Fig. 6A: Chr4 QTL summary plot
# ----------------------------

p_fig6A <- ggplot(data = total_LOD_tb_Chr4, aes(x = trait_label, y = mid_pos)) +
  geom_line(
    data = geneticMapList_tb_Chr4,
    aes(x = trait_label, y = mid_pos),
    inherit.aes = FALSE,
    color = "gray60"
  ) +
  geom_point(
    data = geneticMapList_tb_Chr4,
    aes(x = trait_label, y = mid_pos),
    inherit.aes = FALSE,
    shape = "-",
    size = 10,
    color = "gray40"
  ) +
  geom_rect(
    data = data.frame(
      trait_label = factor("Genetic distance", levels = trait_label_levels),
      ymin = 1355.5686,
      ymax = 1576.4467
    ),
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "blue",
    alpha = 0.2
  ) +
  geom_segment(
    data = filter(total_LOD_tb_Chr4, arrow_head == "first"),
    aes(xend = trait_label, y = low_pos, yend = high_pos, color = treat),
    arrow = arrow(type = "closed", length = grid::unit(0.15, "inches"), ends = "first"),
    lineend = "butt",
    linewidth = 1,
    show.legend = FALSE
  ) +
  geom_segment(
    data = filter(total_LOD_tb_Chr4, arrow_head == "last"),
    aes(xend = trait_label, y = low_pos, yend = high_pos, color = treat),
    arrow = arrow(type = "closed", length = grid::unit(0.15, "inches"), ends = "last"),
    lineend = "butt",
    linewidth = 1,
    show.legend = FALSE
  ) +
  geom_point(aes(color = treat)) +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  xlab("") +
  ylab("Genetic distance (cM)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = c(0.1, 0.9),
    legend.title = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_fig6A

ggsave(
  filename = file.path(output_dir, "Fig_6A_chr4_qtl_summary.png"),
  plot = p_fig6A,
  width = 3,
  height = 5,
  dpi = 600
)

# ----------------------------
# Fig. 6B: inflorescence allele-effect plot
# ----------------------------

p_fig6B <- get_effect_plot(
  cross_object = DG_F3_data,
  marker = total_LOD_tb_Chr4$mid_marker[[5]],
  Trt1 = "cum_W5_Control",
  pos = total_LOD_tb_Chr4$mid_pos[[5]],
  lod = total_LOD_tb_Chr4$mid_lod[[5]],
  y_label = "Inflorescence BLUPs"
) +
  theme(axis.title.x = element_text(margin = margin(t = 20)))

p_fig6B

ggsave(
  filename = file.path(output_dir, "Fig_6B_inflorescence_effect.png"),
  plot = p_fig6B,
  width = 3.5,
  height = 5,
  dpi = 600
)

# ----------------------------
# Fig. 6C: leaf area week 5 allele-effect plot
# ----------------------------

p_fig6C <- get_effect_plot(
  cross_object = DG_F3_data,
  marker = total_LOD_tb_Chr4$mid_marker[[4]],
  Trt1 = "Control_w_5",
  pos = total_LOD_tb_Chr4$mid_pos[[4]],
  lod = total_LOD_tb_Chr4$mid_lod[[4]],
  y_label = "Leaf area w5 BLUPs"
) +
  theme(axis.title.x = element_text(margin = margin(t = 20)))

p_fig6C

ggsave(
  filename = file.path(output_dir, "Fig_6C_leaf_area_w5_effect.png"),
  plot = p_fig6C,
  width = 3.5,
  height = 5,
  dpi = 600
)

