# ============================================================
# Figure 5 and related statistics
#
# Repository root: clover_drought
#
# Inputs are read from ./data/
# Outputs are written to ./figures/Figure5/ and ./results/Figure5/
setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/Drought_F3_paper/clover_drought")
# ============================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(Hmisc)
library(cocor)
library(qvalue)
library(readxl)

# ----------------------------
# Paths
# ----------------------------

data_dir <- "data"
output_fig_dir <- file.path("figures", "Figure5")
output_res_dir <- file.path("results", "Figure5")

dir.create(output_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_res_dir, showWarnings = FALSE, recursive = TRUE)

blup_flower_file <- file.path(data_dir, "flower_count_BLUPs_20240311.csv")
blup_leaf_file   <- file.path(data_dir, "leaf_area_BLUPs_20240311.csv")
blup_dw_file     <- file.path(data_dir, "dry_weight_BLUPs_20240311.csv")

flower_file <- file.path(data_dir, "total_flower_20220524.xlsx")
leaf_file   <- file.path(data_dir, "combined_imputed_leaf_area_20230518.tsv")
dw_file     <- file.path(data_dir, "total_DW_20220625.xlsx")

required_files <- c(
  blup_flower_file,
  blup_leaf_file,
  blup_dw_file,
  flower_file,
  leaf_file,
  dw_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required input files:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

# ----------------------------
# Helper functions
# ----------------------------

process_rcorr_result_triangle_factors <- function(rcorr_result, triangle = "lower") {
  if (!triangle %in% c("lower", "upper")) {
    stop("The 'triangle' parameter must be either 'lower' or 'upper'.")
  }
  
  cor_matrix <- rcorr_result$r
  p_matrix <- rcorr_result$P
  var_names <- rownames(cor_matrix)
  
  cor_df <- as.data.frame(as.table(cor_matrix))
  names(cor_df) <- c("row", "column", "cor")
  
  p_df <- as.data.frame(as.table(p_matrix))
  names(p_df) <- c("row", "column", "p")
  
  combined_df <- merge(cor_df, p_df, by = c("row", "column")) %>%
    mutate(
      row_index = match(row, var_names),
      column_index = match(column, var_names)
    )
  
  if (triangle == "lower") {
    combined_df <- combined_df %>% filter(row_index > column_index)
  } else {
    combined_df <- combined_df %>% filter(row_index < column_index)
  }
  
  combined_df %>%
    select(-row_index, -column_index) %>%
    mutate(
      significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE ~ ""
      ),
      row = factor(row, levels = var_names),
      column = factor(column, levels = var_names)
    )
}

fit_anova_table <- function(data, phenotype_cols, formula_rhs) {
  map_dfr(phenotype_cols, function(phenotype) {
    model_formula <- as.formula(paste(phenotype, formula_rhs))
    model <- lmer(model_formula, data = data)
    
    as.data.frame(anova(model)) %>%
      rownames_to_column("component") %>%
      mutate(phenotype = phenotype)
  })
}

estimate_heritability <- function(data, phenotype_cols, formula_rhs, genetic_term) {
  map_dfr(phenotype_cols, function(phenotype) {
    model_formula <- as.formula(paste(phenotype, formula_rhs))
    model <- lmer(model_formula, data = data)
    
    vc_df <- as.data.frame(VarCorr(model))
    
    genetic_variance <- vc_df %>%
      filter(grp == genetic_term) %>%
      pull(vcov)
    
    residual_variance <- attr(VarCorr(model), "sc")^2
    total_variance <- genetic_variance + residual_variance
    heritability <- genetic_variance / total_variance
    
    tibble(
      phenotype = phenotype,
      genetic_variance = genetic_variance,
      residual_variance = residual_variance,
      total_variance = total_variance,
      heritability = heritability
    )
  })
}

# ----------------------------
# Fig. 5A: bivariate correlation heatmap
# ----------------------------

BLUP_data <- read_csv(blup_flower_file, show_col_types = FALSE) %>%
  left_join(read_csv(blup_leaf_file, show_col_types = FALSE), by = "Genotype") %>%
  left_join(read_csv(blup_dw_file, show_col_types = FALSE), by = "Genotype")

BLUE_data_control <- BLUP_data %>%
  select(starts_with("Control")) %>%
  rename_with(~ str_remove(.x, "^Control_"))

BLUE_data_drought <- BLUP_data %>%
  select(starts_with("Drought")) %>%
  rename_with(~ str_remove(.x, "^Drought_"))

resH_control <- Hmisc::rcorr(as.matrix(BLUE_data_control))
resH_drought <- Hmisc::rcorr(as.matrix(BLUE_data_drought))

cor_control_result <- process_rcorr_result_triangle_factors(resH_control, triangle = "lower")
cor_drought_result <- process_rcorr_result_triangle_factors(resH_drought, triangle = "upper")

cor_merge_result <- bind_rows(cor_control_result, cor_drought_result) %>%
  complete(row, column)

p_fig5A <- ggplot(data = cor_merge_result, aes(x = row, y = column)) +
  geom_raster(aes(fill = cor)) +
  geom_text(aes(label = significance), color = "black", size = 8) +
  scale_fill_gradient2(
    name = NULL,
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  labs(fill = "Correlation Coefficient") +
  scale_x_discrete(
    "Control",
    labels = c(
      "flower" = "Inflorescence",
      "leaf_area_w2" = "Leaf area w2",
      "leaf_area_w3" = "Leaf area w3",
      "leaf_area_w4" = "Leaf area w4",
      "leaf_area_w5" = "Leaf area w5",
      "shoot_w" = "Shoot mass",
      "root_w" = "Root mass"
    )
  ) +
  scale_y_discrete(
    "Drought",
    labels = c(
      "flower" = "Inflorescence",
      "leaf_area_w2" = "Leaf area w2",
      "leaf_area_w3" = "Leaf area w3",
      "leaf_area_w4" = "Leaf area w4",
      "leaf_area_w5" = "Leaf area w5",
      "shoot_w" = "Shoot mass",
      "root_w" = "Root mass"
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_fixed()

ggsave(
  filename = file.path(output_fig_dir, "Fig_5A_bivariate_cor_plot.png"),
  plot = p_fig5A,
  width = 8,
  height = 6,
  dpi = 600
)

write_csv(
  cor_merge_result,
  file.path(output_res_dir, "Fig_5A_correlation_matrix_long.csv")
)

# ----------------------------
# Correlation-difference tests between control and drought
# ----------------------------

variable_names <- colnames(BLUE_data_control)

results_df <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (var1 in variable_names) {
  for (var2 in variable_names) {
    if (var1 != var2) {
      message(paste0("~", var1, " + ", var2, " | ", var1, " + ", var2))
      
      test_result <- cocor::cocor(
        formula = as.formula(paste0("~", var1, " + ", var2, " | ", var1, " + ", var2)),
        data = list(as.data.frame(BLUE_data_control), as.data.frame(BLUE_data_drought))
      )
      
      results_df <- bind_rows(
        results_df,
        data.frame(
          Variable1 = var1,
          Variable2 = var2,
          P_Value = test_result@fisher1925$p.value
        )
      )
    }
  }
}

results_df <- results_df %>%
  mutate(
    q_value = qvalue(P_Value, lambda = 0)$qvalue,
    significance_p = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    significance_q = case_when(
      q_value < 0.001 ~ "***",
      q_value < 0.01  ~ "**",
      q_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

write_csv(
  results_df,
  file.path(output_res_dir, "Fig_5_correlation_difference_tests.csv")
)

# ----------------------------
# rGE between control and drought
# ----------------------------

rGE_df <- data.frame(trait = NULL, rGE = NULL, p_value = NULL)

for (i in colnames(BLUE_data_control)) {
  rGE <- Hmisc::rcorr(BLUE_data_control[[i]], BLUE_data_drought[[i]])
  
  rGE_df <- bind_rows(
    rGE_df,
    data.frame(
      trait = i,
      rGE = rGE$r[1, 2],
      p_value = rGE$P[1, 2]
    )
  )
}

N <- 300

rGE_df <- rGE_df %>%
  mutate(
    FishersZ = 0.5 * log((1 + rGE) / (1 - rGE)),
    SE_Z = 1 / sqrt(N - 3),
    hypothetical_r = 0.9,
    Z_expected = 0.5 * log((1 + hypothetical_r) / (1 - hypothetical_r)),
    Z_score = (FishersZ - Z_expected) / SE_Z,
    p_value_new = 2 * (1 - pnorm(abs(Z_score)))
  ) %>%
  select(-Z_expected, -Z_score)

write_csv(
  rGE_df,
  file.path(output_res_dir, "Fig_5_rGE.csv")
)

# ----------------------------
# Mixed-model ANOVA summaries
# ----------------------------

fl_anova_data <- read_xlsx(flower_file) %>%
  mutate(
    plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300, length.out = n()),
    total_fl = W1 + W2 + W3 + W4 + W5
  ) %>%
  select(
    Genotype = unique_F3_id,
    Trt,
    Subgroup = plant_group,
    Replication = Rep,
    total_fl
  ) %>%
  mutate(
    Trt = factor(Trt),
    Subgroup = factor(Subgroup),
    Replication = factor(Replication)
  )

flower_anova <- fit_anova_table(
  data = fl_anova_data,
  phenotype_cols = "total_fl",
  formula_rhs = "~ Trt + (1|Replication) + (1|Genotype) + Trt:Genotype"
)

lf_raw <- read_tsv(leaf_file, show_col_types = FALSE)

lf_anova_data <- lf_raw %>%
  mutate(subgroup = rep(rep(c(1, 2), each = 1500), length.out = n())) %>%
  select(2, 5, 6, 9, 11, 12, 13) %>%
  left_join(
    lf_raw %>%
      mutate(subgroup = rep(rep(c(1, 2), each = 1500), length.out = n())) %>%
      filter(exp_week == 1) %>%
      select(
        unique_F3_id,
        Trt,
        subgroup,
        Rep,
        leaf_area_w1 = leaf_area
      ),
    by = c("unique_F3_id", "Trt", "subgroup", "Rep")
  ) %>%
  mutate(
    Relative_growth = leaf_area - leaf_area_w1,
    exp_week = paste0("w_", exp_week)
  ) %>%
  select(
    Genotype = unique_F3_id,
    Trt,
    Subgroup = subgroup,
    Replication = Rep,
    exp_week,
    Relative_growth
  ) %>%
  pivot_wider(names_from = exp_week, values_from = Relative_growth) %>%
  select(-w_1) %>%
  mutate(
    Trt = factor(Trt),
    Subgroup = factor(Subgroup),
    Replication = factor(Replication)
  )

leaf_anova_phenotypes <- names(lf_anova_data)[
  !(names(lf_anova_data) %in% c("Genotype", "Subgroup", "Replication", "Trt"))
]

leaf_anova <- fit_anova_table(
  data = lf_anova_data,
  phenotype_cols = leaf_anova_phenotypes,
  formula_rhs = "~ Trt + (1|Subgroup) + (1|Subgroup:Replication) + (1|Subgroup:Genotype) + Trt:Genotype"
)

dw_anova_data <- read_xlsx(dw_file) %>%
  mutate(
    plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300, length.out = n()),
    Trt_group_rep = paste(Trt, plant_group, Rep, sep = "_"),
    shoot_w = Shoot_gross - Shoot_bag,
    root_w = Root_gross - Root_bag
  ) %>%
  select(
    Genotype = unique_F3_id,
    Trt,
    Subgroup = plant_group,
    Replication = Rep,
    shoot_w,
    root_w
  ) %>%
  mutate(
    Trt = factor(Trt),
    Subgroup = factor(Subgroup),
    Replication = factor(Replication)
  )

dw_anova_phenotypes <- names(dw_anova_data)[
  !(names(dw_anova_data) %in% c("Genotype", "Subgroup", "Replication", "Trt"))
]

dw_anova <- fit_anova_table(
  data = dw_anova_data,
  phenotype_cols = dw_anova_phenotypes,
  formula_rhs = "~ Trt + (1|Subgroup) + (1|Subgroup:Replication) + (1|Subgroup:Genotype) + Trt:Genotype"
)

total_ANOVA_df <- bind_rows(flower_anova, leaf_anova, dw_anova)

write_csv(
  total_ANOVA_df,
  file.path(output_res_dir, "Fig_5_total_ANOVA_df.csv")
)

# ----------------------------
# Heritability: flower number
# ----------------------------

fl_h2_data <- read_xlsx(flower_file) %>%
  mutate(
    plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300, length.out = n()),
    total_fl = W1 + W2 + W3 + W4 + W5
  ) %>%
  select(
    Genotype = unique_F3_id,
    Trt,
    Subgroup = plant_group,
    Replication = Rep,
    total_fl
  ) %>%
  pivot_wider(names_from = Trt, values_from = total_fl) %>%
  mutate(
    Replication = rep(1:6, each = 150, length.out = n()),
    Subgroup = factor(Subgroup),
    Replication = factor(Replication)
  )

flower_h2_phenotypes <- names(fl_h2_data)[
  !(names(fl_h2_data) %in% c("Genotype", "Subgroup", "Replication"))
]

flower_h2 <- estimate_heritability(
  data = fl_h2_data,
  phenotype_cols = flower_h2_phenotypes,
  formula_rhs = "~ (1|Replication) + (1|Genotype)",
  genetic_term = "Genotype"
)

write_csv(
  flower_h2,
  file.path(output_res_dir, "Fig_5_flower_heritability.csv")
)

# ----------------------------
# Heritability: leaf area growth traits
# ----------------------------

lf_h2_data <- lf_raw %>%
  mutate(subgroup = rep(rep(c(1, 2), each = 1500), length.out = n())) %>%
  select(2, 5, 6, 9, 11, 12, 13) %>%
  left_join(
    lf_raw %>%
      mutate(subgroup = rep(rep(c(1, 2), each = 1500), length.out = n())) %>%
      filter(exp_week == 1) %>%
      select(
        unique_F3_id,
        Trt,
        subgroup,
        Rep,
        leaf_area_w1 = leaf_area
      ),
    by = c("unique_F3_id", "Trt", "subgroup", "Rep")
  ) %>%
  mutate(
    Relative_growth = leaf_area - leaf_area_w1,
    exp_week = paste0("w_", exp_week)
  ) %>%
  select(
    Genotype = unique_F3_id,
    Trt,
    Subgroup = subgroup,
    Replication = Rep,
    exp_week,
    Relative_growth
  ) %>%
  pivot_wider(names_from = c(Trt, exp_week), values_from = Relative_growth) %>%
  select(-Control_w_1, -Drought_w_1) %>%
  mutate(
    Subgroup = factor(Subgroup),
    Replication = factor(Replication)
  )

leaf_h2_phenotypes <- names(lf_h2_data)[
  !(names(lf_h2_data) %in% c("Genotype", "Subgroup", "Replication"))
]

leaf_h2 <- estimate_heritability(
  data = lf_h2_data,
  phenotype_cols = leaf_h2_phenotypes,
  formula_rhs = "~ (1|Subgroup) + (1|Subgroup:Replication) + (1|Subgroup:Genotype)",
  genetic_term = "Subgroup:Genotype"
)

write_csv(
  leaf_h2,
  file.path(output_res_dir, "Fig_5_leaf_area_heritability.csv")
)

# ----------------------------
# Heritability: dry weights
# ----------------------------

dw_h2_data <- read_xlsx(dw_file) %>%
  mutate(
    plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300, length.out = n()),
    Trt_group_rep = paste(Trt, plant_group, Rep, sep = "_"),
    shoot_w = Shoot_gross - Shoot_bag,
    root_w = Root_gross - Root_bag
  ) %>%
  select(
    Genotype = unique_F3_id,
    Trt,
    Subgroup = plant_group,
    Replication = Rep,
    shoot_w,
    root_w
  ) %>%
  pivot_wider(names_from = Trt, values_from = c(shoot_w, root_w)) %>%
  mutate(
    Subgroup = factor(Subgroup),
    Replication = factor(Replication)
  )

dw_h2_phenotypes <- names(dw_h2_data)[
  !(names(dw_h2_data) %in% c("Genotype", "Subgroup", "Replication"))
]

dw_h2 <- estimate_heritability(
  data = dw_h2_data,
  phenotype_cols = dw_h2_phenotypes,
  formula_rhs = "~ (1|Replication) + (1|Genotype)",
  genetic_term = "Genotype"
)

write_csv(
  dw_h2,
  file.path(output_res_dir, "Fig_5_dry_weight_heritability.csv")
)
