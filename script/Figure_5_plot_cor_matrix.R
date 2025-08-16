library(tidyverse)
library(lme4)
library(lmerTest)
library(Hmisc)
library(cocor)
library(qvalue)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper")


BLUP_data <- read_csv("flower_count_BLUPs_20240311.csv")%>%
  left_join(read_csv("leaf_area_BLUPs_20240311.csv"), by = "Genotype")%>%
  left_join(read_csv("dry_weight_BLUPs_20240311.csv"), by = "Genotype")

BLUE_data_control <- select(BLUP_data, starts_with("Control"))%>%
  rename_all(.funs = ~str_remove(., "Control_"))
BLUE_data_drought <- select(BLUP_data, starts_with("Drought"))%>%
  rename_all(.funs = ~str_remove(., "Drought_"))

resH_control <- Hmisc::rcorr(as.matrix(BLUE_data_control))
resH_drought <- rcorr(as.matrix(BLUE_data_drought))
resH_control



process_rcorr_result_triangle_factors <- function(rcorr_result, triangle = "lower") {
  # Ensure that the triangle parameter is valid
  if (!triangle %in% c("lower", "upper")) {
    stop("The 'triangle' parameter must be either 'lower' or 'upper'.")
  }
  
  # Extract the correlation and p-value matrices
  cor_matrix <- rcorr_result$r
  p_matrix <- rcorr_result$P
  
  # Get row and column names from the correlation matrix
  var_names <- rownames(cor_matrix)
  
  # Convert matrices to long format data frames
  cor_df <- as.data.frame(as.table(cor_matrix))
  names(cor_df) <- c("row", "column", "cor")
  p_df <- as.data.frame(as.table(p_matrix))
  names(p_df) <- c("row", "column", "p")
  
  # Combine the correlation and p-value data frames
  combined_df <- merge(cor_df, p_df, by = c("row", "column"))
  
  # Determine row and column indices based on the variable names
  combined_df$row_index <- match(combined_df$row, var_names)
  combined_df$column_index <- match(combined_df$column, var_names)
  
  # Filter for lower or upper triangle based on indices
  if (triangle == "lower") {
    combined_df <- combined_df %>% filter(row_index > column_index)
  } else { # upper triangle
    combined_df <- combined_df %>% filter(row_index < column_index)
  }
  
  # Remove the index columns
  combined_df <- combined_df %>% select(-row_index, -column_index)
  
  # Add significance levels based on p values
  combined_df <- combined_df %>%
    mutate(significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ))
  
  # Convert 'row' and 'column' back to factors with the original ordering
  combined_df$row <- factor(combined_df$row, levels = var_names)
  combined_df$column <- factor(combined_df$column, levels = var_names)
  
  return(combined_df)
}


cor_control_result <- process_rcorr_result_triangle_factors(resH_control, triangle = "lower")
cor_drought_result <- process_rcorr_result_triangle_factors(resH_drought, triangle = "upper")
cor_merge_result <- bind_rows(cor_control_result, cor_drought_result)%>%
  complete(row, column)


# Plot
p <- ggplot(data = cor_merge_result, aes(x = row, y = column)) +
  geom_raster(aes(fill = cor)) + # Color tiles by correlation coefficient
  geom_text(aes(label = significance), color = "black") + # Add significance levels
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-1, 1)) + # Adjust color gradient
  theme_minimal() + # Use a minimal theme
  labs(fill = "Correlation Coefficient") + # Label for the color scale
  scale_x_discrete("Control", label = c("Flower count", "Leaf area w2", "Leaf area w3", "Leaf area w4", "Leaf area w5", "Shoot mass", "Root mass"))+
  scale_y_discrete("Drought", label = c("Flower count", "Leaf area w2", "Leaf area w3", "Leaf area w4", "Leaf area w5", "Shoot mass", "Root mass"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.background = element_rect(fill = "white", color = NA))+
  coord_fixed()

p

ggsave("bivariate_cor_plot_20240312.png", width = 8, height = 6)

colnames(BLUE_data_control)
colnames(BLUE_data_drought)

# Perform the statistical test
test_result <- cocor::cocor(formula = ~w_2 + w_3 | w_2 + w_3,
                            data = list(as.data.frame(BLUE_data_control), as.data.frame(BLUE_data_drought)))
test_result
# Extract p-value from the test result
# The structure of test_result may vary based on the test used; you might need to adjust the code to correctly extract the p-value.
# Assuming the p-value is located in a list element named 'p.value'
p_value <- test_result@fisher1925$p.value

# Create a data frame to store the results
results_df <- data.frame(Variable1 = "flower",
                         Variable2 = "root_w",
                         P_Value = p_value)

# Assuming BLUE_data_control and BLUE_data_drought are already loaded in your workspace

# Get the column names (excluding 'flower' if you only want to compare it with other variables)
variable_names <- colnames(BLUE_data_control)

# Initialize an empty data frame for results
results_df <- data.frame(Variable1 = character(),
                         Variable2 = character(),
                         P_Value = numeric(),
                         stringsAsFactors = FALSE)

# Loop through all unique pairs of variables
for (var1 in variable_names) {
  for (var2 in variable_names) {
    if (var1 != var2) {  # Ensure not comparing a variable with itself
      # Perform the statistical test
      print(paste0("~", var1, " + ", var2, " | ", var1, " + ", var2))
      test_result <- cocor::cocor(formula = as.formula(paste0("~", var1, " + ", var2, " | ", var1, " + ", var2)),
                                  data = list(as.data.frame(BLUE_data_control), as.data.frame(BLUE_data_drought)))
      
      # Extract the p-value
      p_value <- test_result@fisher1925$p.value
      
      # Append the results to the results data frame
      results_df <- rbind(results_df, data.frame(Variable1 = var1,
                                                 Variable2 = var2,
                                                 P_Value = p_value))
    }
  }
}

# Inspect the first few rows of the results
results_df <- results_df %>%
  mutate(q_value = qvalue(P_Value, lambda=0)$qvalue)%>%
  mutate(significance_p = case_when(
    P_Value < 0.001 ~ "***",
    P_Value < 0.01  ~ "**",
    P_Value < 0.05  ~ "*",
    TRUE      ~ ""
  ))%>%
    mutate(significance_q = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE      ~ ""
    ))

###################  
rGE_df <- data.frame(trait = NULL, rGE = NULL, p_value = NULL)
for (i in colnames(BLUE_data_control)){
  rGE <- rcorr(BLUE_data_control[[i]], BLUE_data_drought[[i]])
  rGE_df <- bind_rows(rGE_df, data.frame(trait = i, rGE = rGE$r[1,2], p_value = rGE$P[1,2]))
}

N <- 300
rGE_df <- rGE_df  %>%
mutate(FishersZ = 0.5 * log((1 + rGE) / (1 - rGE)),
       SE_Z = 1 / sqrt(N - 3),
       # Hypothetical comparison against a Z value (not 1, as directly testing against 1 is impractical)
       # Calculate Z-score for comparison. Here, hypothetically comparing against 0.8 just for demonstration
       # In practice, specify the hypothetical_r you are interested in
       hypothetical_r = 0.9, # This is just for demonstration
       Z_expected = 0.5 * log((1 + hypothetical_r) / (1 - hypothetical_r)),
       Z_score = (FishersZ - Z_expected) / SE_Z,
       # Calculate p-value from Z-score, assuming a two-tailed test
       p_value_new = 2 * (1 - pnorm(abs(Z_score)))) %>%
  select(-Z_expected, -Z_score)

colnames(rGE_df)
write_csv(rGE_df, "rGE_20240312.csv")

plot(BLUE_data_control$flower, BLUE_data_drought$flower)
plot(BLUE_data_drought$flower, BLUE_data_drought$root_w)


#################
library(readxl)

fl_data <- read_xlsx("../GBS_mapping_population/Drought_exp/total_flower_20220524.xlsx")%>%
  mutate(plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300))%>%
  mutate(total_fl = W1 + W2 + W3 + W4 + W5)%>%
  select(Genotype = unique_F3_id, Trt, Subgroup = plant_group, Replication = Rep, total_fl)%>%
  mutate(Subgroup = as.factor(Subgroup), Replication = as.factor(Replication))

str(fl_data)

# List of phenotype columns
phenotype_cols <- names(fl_data)[!(names(fl_data) %in% c("Genotype", "Subgroup", "Replication", "Trt"))]
phenotype_cols
# Initialize an empty list to store BLUPs for each phenotype
ANOVA_list <- list()

# Loop over each phenotype column
for (phenotype in phenotype_cols) {
  phenotype = "total_fl"
  
  # Fit the mixed model
  model <- lmer(as.formula(paste(phenotype, "~ Trt + (1|Replication) + (1|Genotype) + Trt:Genotype")), data = fl_data)
  
  VarCorr(model)
  lmerTest::ranova(model)
  # Perform ANOVA on the model
  anova_result <- anova(model)
  
  # Convert the ANOVA result to a data frame
  anova_df <- as.data.frame(anova_result)%>%
    rownames_to_column("conponent")%>%
    mutate(phenotype = phenotype)
  
  ANOVA_list[[phenotype]] <- anova_df
}

# Combine all BLUPs data frames into a single data frame
# Assuming that the genotype order is consistent across all BLUPs calculations
total_ANOVA_df <- bind_rows(ANOVA_list)


lf_data_blup <- read_tsv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_mapping_population/Drought_exp/combined_imputed_leaf_area_20230518.tsv")%>%
  mutate(subgroup = rep(rep(c(1, 2), each = 1500), length = 9000))%>%
  select(c(2,5,6, 9, 11,12,13))%>%
  left_join(select(filter(., exp_week == 1), unique_F3_id, Trt, subgroup, leaf_area_w1 = leaf_area, Rep), by = c("unique_F3_id", "Trt", "subgroup", "Rep"))%>%
  mutate(Relative_growth = leaf_area - leaf_area_w1, exp_week = paste("w", exp_week, sep = "_"))%>%
  select(Genotype = unique_F3_id, Trt, Subgroup = subgroup, Replication = Rep, exp_week, Relative_growth)%>%
  pivot_wider(names_from = c("exp_week"), values_from = c("Relative_growth"))%>%
  select(-w_1, -w_1)



# List of phenotype columns
phenotype_cols <- names(lf_data_blup)[!(names(lf_data_blup) %in% c("Genotype", "Subgroup", "Replication", "Trt"))]
phenotype_cols
# Initialize an empty list to store BLUPs for each phenotype
ANOVA_list <- list()

# Loop over each phenotype column
for (phenotype in phenotype_cols) {
  # phenotype = "total_fl"
  
  # Fit the mixed model
  model <- lmer(as.formula(paste(phenotype, "~ Trt + (1|Subgroup) + (1|Subgroup:Replication) + (1|Subgroup:Genotype) + Trt:Genotype")), data = lf_data_blup)
  
  # VarCorr(model)
  # lmerTest::ranova(model)
  # Perform ANOVA on the model
  anova_result <- anova(model)
  
  # Convert the ANOVA result to a data frame
  anova_df <- as.data.frame(anova_result)%>%
    rownames_to_column("conponent")%>%
    mutate(phenotype = phenotype)
  
  ANOVA_list[[phenotype]] <- anova_df
}

# Combine all BLUPs data frames into a single data frame
# Assuming that the genotype order is consistent across all BLUPs calculations
total_ANOVA_df <- bind_rows(bind_rows(ANOVA_list), total_ANOVA_df)

dw_data <- read_xlsx("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_mapping_population/Drought_exp/total_DW_20220625.xlsx")%>%
  mutate(plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300))%>%
  mutate(Trt_group_rep = paste(Trt, plant_group, Rep, sep = "_"))%>%
  mutate(shoot_w = Shoot_gross - Shoot_bag, root_w = Root_gross - Root_bag)%>%
  select(Genotype = unique_F3_id, Trt, Subgroup = plant_group, Replication = Rep, shoot_w, root_w)%>%
  mutate(Subgroup = as.factor(Subgroup), Replication = as.factor(Replication))

# List of phenotype columns
phenotype_cols <- names(dw_data)[!(names(dw_data) %in% c("Genotype", "Subgroup", "Replication", "Trt"))]
phenotype_cols
# Initialize an empty list to store BLUPs for each phenotype
ANOVA_list <- list()

# Loop over each phenotype column
for (phenotype in phenotype_cols) {
  # phenotype = "total_fl"
  
  # Fit the mixed model
  model <- lmer(as.formula(paste(phenotype, "~ Trt + (1|Subgroup) + (1|Subgroup:Replication) + (1|Subgroup:Genotype) + Trt:Genotype")), data = dw_data)
  
  # VarCorr(model)
  # lmerTest::ranova(model)
  # Perform ANOVA on the model
  anova_result <- anova(model)
  
  # Convert the ANOVA result to a data frame
  anova_df <- as.data.frame(anova_result)%>%
    rownames_to_column("conponent")%>%
    mutate(phenotype = phenotype)
  
  ANOVA_list[[phenotype]] <- anova_df
}

# Combine all BLUPs data frames into a single data frame
# Assuming that the genotype order is consistent across all BLUPs calculations
total_ANOVA_df <- bind_rows(bind_rows(ANOVA_list), total_ANOVA_df)

write_csv(total_ANOVA_df, "total_ANOVA_df_20240312.csv")

var_out <- VarCorr(model)
str(var_out)
##################

fl_data <- read_xlsx("../GBS_mapping_population/Drought_exp/total_flower_20220524.xlsx")%>%
  mutate(plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300))%>%
  mutate(total_fl = W1 + W2 + W3 + W4 + W5)%>%
  select(Genotype = unique_F3_id, Trt, Subgroup = plant_group, Replication = Rep, total_fl)%>%
  pivot_wider(names_from = c("Trt"), values_from = c("total_fl"))%>%
  mutate(Replication = rep(1:6, each = 150))%>%
  mutate(Subgroup = as.factor(Subgroup), Replication = as.factor(Replication))

str(fl_data)

# List of phenotype columns
phenotype_cols <- names(fl_data)[!(names(fl_data) %in% c("Genotype", "Subgroup", "Replication"))]
phenotype_cols
# Initialize an empty list to store BLUPs for each phenotype
heritability_list <- list()

# Loop over each phenotype column
for (phenotype in phenotype_cols) {
  # phenotype = "Control"
  
  
  # Fit the mixed model
  model <- lmer(as.formula(paste(phenotype, "~ (1|Replication) + (1|Genotype)")), data = fl_data)
  
  # Extract variance components (these are actually standard deviations)
  var_components_sd <- VarCorr(model)
  print(var_components_sd)

  genetic_variance <- var_components_sd$`Genotype`[1,1]
  
  # Calculate total variance, which includes genetic variance, variance from other random effects, and residual variance
  # Don't forget to square the residual standard deviation to get the residual variance
  residual_variance <- attr(var_components_sd, "sc")^2
  total_variance <- genetic_variance + residual_variance
  
  # Calculate heritability
  heritability <- genetic_variance / total_variance
  
  
  # Store the BLUPs data frame in the list, named by phenotype
  heritability_list[[phenotype]] <- heritability
}

heritability_list


lf_data_blup <- read_tsv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_mapping_population/Drought_exp/combined_imputed_leaf_area_20230518.tsv")%>%
  mutate(subgroup = rep(rep(c(1, 2), each = 1500), length = 9000))%>%
  select(c(2,5,6, 9, 11,12,13))%>%
  left_join(select(filter(., exp_week == 1), unique_F3_id, Trt, subgroup, leaf_area_w1 = leaf_area, Rep), by = c("unique_F3_id", "Trt", "subgroup", "Rep"))%>%
  mutate(Relative_growth = leaf_area - leaf_area_w1, exp_week = paste("w", exp_week, sep = "_"))%>%
  select(Genotype = unique_F3_id, Trt, Subgroup = subgroup, Replication = Rep, exp_week, Relative_growth)%>%
  pivot_wider(names_from = c("Trt", "exp_week"), values_from = c("Relative_growth"))%>%
  select(-Control_w_1, -Drought_w_1)


# List of phenotype columns
phenotype_cols <- names(lf_data_blup)[!(names(lf_data_blup) %in% c("Genotype", "Subgroup", "Replication"))]
phenotype_cols
# Initialize an empty list to store BLUPs for each phenotype
heritability_list <- list()

# Loop over each phenotype column
for (phenotype in phenotype_cols) {
  # phenotype = "Control"
  
  
  # Fit the mixed model
  model <- lmer(as.formula(paste(phenotype, "~ (1|Subgroup) + (1|Subgroup:Replication) + (1|Subgroup:Genotype)")), data = lf_data_blup)
  
  # Extract variance components (these are actually standard deviations)
  var_components_sd <- VarCorr(model)
  print(var_components_sd)
  
  genetic_variance <- var_components_sd$`Subgroup:Genotype`[1,1]
  
  # Calculate total variance, which includes genetic variance, variance from other random effects, and residual variance
  # Don't forget to square the residual standard deviation to get the residual variance
  residual_variance <- attr(var_components_sd, "sc")^2
  total_variance <- genetic_variance + residual_variance
  
  # Calculate heritability
  heritability <- genetic_variance / total_variance
  
  
  # Store the BLUPs data frame in the list, named by phenotype
  heritability_list[[phenotype]] <- heritability
}

heritability_list

dw_data <- read_xlsx("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_mapping_population/Drought_exp/total_DW_20220625.xlsx")%>%
  mutate(plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300))%>%
  mutate(Trt_group_rep = paste(Trt, plant_group, Rep, sep = "_"))%>%
  mutate(shoot_w = Shoot_gross - Shoot_bag, root_w = Root_gross - Root_bag)%>%
  select(Genotype = unique_F3_id, Trt, Subgroup = plant_group, Replication = Rep, shoot_w, root_w)%>%
  pivot_wider(names_from = c("Trt"), values_from = c("shoot_w", "root_w"))%>%
  mutate(Subgroup = as.factor(Subgroup), Replication = as.factor(Replication))


str(dw_data)
# List of phenotype columns
phenotype_cols <- names(dw_data)[!(names(dw_data) %in% c("Genotype", "Subgroup", "Replication"))]
phenotype_cols
# Initialize an empty list to store BLUPs for each phenotype
heritability_list <- list()

# Loop over each phenotype column
for (phenotype in phenotype_cols) {
  # phenotype = "Control"
  
  
  # Fit the mixed model
  # model <- lmer(as.formula(paste(phenotype, "~ (1|Subgroup) + (1|Subgroup:Replication) + (1|Subgroup:Genotype)")), data = dw_data)
  
  model <- lmer(as.formula(paste(phenotype, "~ (1|Replication) + (1|Genotype)")), data = dw_data)
  
  # Extract variance components (these are actually standard deviations)
  var_components_sd <- VarCorr(model)
  print(var_components_sd)
  
  # genetic_variance <- var_components_sd$`Subgroup:Genotype`[1,1]
  genetic_variance <- var_components_sd$Genotype[1,1]
  
  # Calculate total variance, which includes genetic variance, variance from other random effects, and residual variance
  # Don't forget to square the residual standard deviation to get the residual variance
  residual_variance <- attr(var_components_sd, "sc")^2
  total_variance <- genetic_variance + residual_variance
  
  # Calculate heritability
  heritability <- genetic_variance / total_variance
  
  
  # Store the BLUPs data frame in the list, named by phenotype
  heritability_list[[phenotype]] <- heritability
}

heritability_list
