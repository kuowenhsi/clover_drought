library(tidyverse)
library(readxl)
library(cowplot)
library(ungeviz)
library(see)
library(car)
library(agricolae)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_mapping_population/Drought_exp/parental_lines")


get_relative_growth_rate <- function(x){
  x <- x %>%
    group_by(ID) %>%
    mutate(Growth_rate = (Green_area - lead(Green_area))/as.numeric(Date - lead(Date)))%>%
    mutate(Relative_growth = Green_area - first(Green_area, order_by = Date))%>%
    ungroup()
  
  return(x)
  
}


# make a note here

plot_data <- read_xlsx("Parental_drought_exp_20210912.xlsx", na = c("", "NA")) %>%
  filter(Date > as.Date("2021-08-01"))%>%
  filter(Accession != "no_plant")%>%
  arrange(Date)%>%
  group_by(ID)%>%
  mutate(cum_inflo = cumsum(Inflo))%>%
  mutate(Day_after_Trt = as.numeric(difftime(Date, min(Date), units = "days")))

water_data <- read_xlsx("Lam_waterpotential.xlsx", na = c("", "NA"))%>%
  pivot_longer(cols = c(W1, W2, T1, T2)) %>%
  mutate(leaf_id = rep(1:2, 96), name = str_sub(name, 1, 1)) %>%
  pivot_wider(id_cols = c(ID, Accession, Trt, Rep, leaf_id), names_from = name, values_from = value)%>%
  mutate(Water_potential = -W/10)

dryW_data <- read_xlsx("Dry_weight.xlsx", na = c("", "NA")) %>%
  mutate(DW_shoot_root = Shoot/Root,
         DW_root_shoot = Root/Shoot,
         DW_total = Shoot+Root)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper")

lm_model <- lm(Shoot ~ Accession + Trt + Accession:Trt, data = dryW_data)
ANOVA_test <- (Anova(lm_model, type = "II"))
value_str <- format(sprintf("%.2e", ANOVA_test$`Pr(>F)`[[3]]), scientific = TRUE)

# Extract the coefficient and exponent parts
parts <- strsplit(value_str, "e")[[1]]
coefficient <- parts[1]
exponent <- as.integer(parts[2])
expr <- paste("Interaction:", deparse(bquote(italic(p) == .(coefficient) %*% 10^.(exponent))))

post_result <- (HSD.test(lm_model, trt = c("Accession", "Trt"), alpha = 0.05, group = TRUE))$groups%>%
  rownames_to_column("items")%>%
  separate(items, into = c("Accession", "Trt"), sep = ":")%>%
  left_join(dryW_data %>% group_by(Accession, Trt) %>% summarise(y_max = max(Shoot)), by = c("Accession", "Trt"))

p1 <- ggplot(data = dryW_data, aes(x = Accession, y = Shoot))+
  geom_violin(aes(group = interaction(Accession, Trt), fill = Trt), color = NA, scale = "width", alpha = 0.7)+
  geom_boxplot(aes(group =interaction(Accession, Trt)), fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.9))+
  geom_point(aes(group =interaction(Accession, Trt)),position = position_dodge(width = 0.9), size = 0.5)+
  geom_text(x = 0.5, y = max(dryW_data$Shoot)*1, label = expr, parse = TRUE, size = 3.5, fontface = "plain",check_overlap = TRUE, hjust = 0)+
  geom_text(data = post_result, aes(y = y_max + max(dryW_data$Shoot)*0.04, label = groups, group = interaction(Accession, Trt)), position = position_dodge(width = 0.9))+
  scale_x_discrete("", labels = c("DMN_010", "GFL_007", "STL_0701"))+
  ylab("Shoot dry mass (g)")+
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank())
p1


# ggsave("parental_shoot_dw.png", width = 3.5, height = 3.5, dpi = 900)
################################################################

perform_analysis <- function(data, response_variable, Y_label, x_position = 0.1, y_position = 0.9, y_pad) {
  
  # data = water_data; response_variable = "Water_potential"

  formula_str <- as.formula(paste(response_variable, "~ Accession + Trt + Accession:Trt"))
  lm_model <- lm(formula_str, data = data)
  ANOVA_test <- Anova(lm_model, type = "II")
  
  value_str <- format(sprintf("%.2e", ANOVA_test$`Pr(>F)`[[3]]), scientific = TRUE)
  parts <- strsplit(value_str, "e")[[1]]
  coefficient <- parts[1]
  exponent <- as.integer(parts[2])
  
  
  if (ANOVA_test$`Pr(>F)`[[3]] < 0.05){
    post_result <- HSD.test(lm_model, trt = c("Accession", "Trt"), alpha = 0.05, group = TRUE)$groups %>%
      rownames_to_column("items") %>%
      separate(items, into = c("Accession", "Trt"), sep = ":") %>%
      left_join(data %>% group_by(Accession, Trt) %>% summarise(y_max = max(!!as.symbol(response_variable), na.rm = TRUE)), by = c("Accession", "Trt"))
    
    expr <- paste(deparse(bquote(italic(p) == .(coefficient) %*% 10^.(exponent))))
  } else {
    data_control <- filter(data, Trt == "control")
    lm_model <- lm(as.formula(paste(response_variable, "~ Accession")), data = data_control)
    post_result_1 <- HSD.test(lm_model, trt = c("Accession"), alpha = 0.05, group = TRUE)$groups %>%
      rownames_to_column("Accession") %>%
      mutate(Trt = "control") %>%
      left_join(data %>% group_by(Accession, Trt) %>% summarise(y_max = max(!!as.symbol(response_variable), na.rm = TRUE)), by = c("Accession", "Trt"))
    data_drought <- filter(data, Trt == "drought")
    lm_model <- lm(as.formula(paste(response_variable, "~ Accession")), data = data_drought)
    post_result_2 <- HSD.test(lm_model, trt = c("Accession"), alpha = 0.05, group = TRUE)$groups %>%
      rownames_to_column("Accession") %>%
      mutate(Trt = "drought") %>%
      left_join(data %>% group_by(Accession, Trt) %>% summarise(y_max = max(!!as.symbol(response_variable))), by = c("Accession", "Trt"))
    post_result <- bind_rows(post_result_1, post_result_2)
    
    expr <- paste(paste(deparse(bquote(italic(p) == .(coefficient) %*% 10^.(exponent)))), "~ '(n.s.)'")
  }
  
  y_npc <- min(data[[response_variable]], na.rm = TRUE) + y_position*(max(data[[response_variable]], na.rm = TRUE) - min(data[[response_variable]], na.rm = TRUE))
  x_npc <- 3.5*x_position
  
  p <- ggplot(data = data, aes(x = Accession, y = !!as.symbol(response_variable))) +
    geom_violin(aes(group = interaction(Accession, Trt), fill = Trt), color = NA, scale = "width", alpha = 0.7) +
    geom_boxplot(aes(group =interaction(Accession, Trt)), fill = NA, width = 0.2, outlier.shape = NA, position = position_dodge(width = 0.9)) +
    geom_point(aes(group =interaction(Accession, Trt)), position = position_dodge(width = 0.9), size = 0.5) +
    annotate(geom = "text", x = x_npc, y = y_npc, label = expr, parse = TRUE, size = 3.5, fontface = "plain", check_overlap = TRUE, hjust = 0) +
    geom_text(data = post_result, aes(y = y_max + abs(max(data[[response_variable]], na.rm = TRUE))*0.04 + y_pad, label = groups, group = interaction(Accession, Trt)), position = position_dodge(width = 0.9)) +
    scale_x_discrete("", labels = c("DMN_010", "GFL_007", "STL_0701")) +
    ylab(Y_label) +
    scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank())
  
  return(p)
}

# Example of calling the function
# Replace 'your_dataset' with your actual dataset variable
# perform_analysis(your_dataset, "Shoot")



###############################################################

p2 <- perform_analysis(dryW_data, "Root", "Root dry mass (g)", 0.15, 1, 0)
p2


# ggsave("parental_root_dw.png", width = 3.5, height = 3.5, dpi = 900)

p3 <- perform_analysis(dryW_data, "DW_total", "Shoot + Root dry mass (g)", 0.15, 1, 0)
p3

# ggsave("parental_DW_total.png", width = 3.5, height = 3.5, dpi = 900)


p4 <- perform_analysis(dryW_data, "DW_shoot_root", "Shoot/Root ratio", 0.6, 1, 0)
p4


# ggsave("parental_DW_shoot_root_total.png", width = 3.5, height = 3.5, dpi = 900)



p5 <- perform_analysis(dryW_data, "DW_root_shoot", "Root/Shoot ratio", 0.15, 1, 0)
p5


# ggsave("parental_DW_root_shoot.png", width = 3.5, height = 3.5, dpi = 900)


p6 <- perform_analysis(water_data, "Water_potential", expression(psi[w]~(Mpa)), 0.15, 0.05, 0.1)
p6



# ggsave("parental_water_potential.png", width = 3.5, height = 3.5, dpi = 900)



p7 <- perform_analysis(water_data, "T", expression(Leaf~thickness~(mu*m)), 0.6, 0.05, 0)
p7


# ggsave("parental_leaf_thickness.png", width = 3.5, height = 3.5, dpi = 900)



growth_data <- plot_data %>%
  filter(Green_area > 0)%>%
  get_relative_growth_rate()



p <- ggplot(data = plot_data, aes(x = Day_after_Trt, y = VWC)) +
  stat_summary(fun = mean, geom = "line", aes(color = Trt), size = 1)+
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(color = Trt))+
  labs(x = "Days after treatment (treatment started 7 days after propagation)", y = "Volumetric water content (%)")+
  scale_color_manual("",values = c("skyblue2", "brown2"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), legend.position = c(0.05,0.93), legend.title = element_blank(), legend.background = element_rect(fill = NA))+
  facet_wrap(.~Accession, nrow = 1)
p
# ggsave("parental_VWC_plot.png", width = 8, height = 3)

##########

mean_growth_data <- growth_data %>%
  group_by(ID, Accession, Trt, Rep) %>%
  summarise(mean_growth_rate = mean(Growth_rate, na.rm = TRUE))


p8 <- perform_analysis(mean_growth_data, "mean_growth_rate", expression(Average~growth~rate~(mm**2/Day)), 0.15, 1, 0)
p8


# ggsave("parental_mean_growth_rate.png", width = 3.5, height = 3.5, dpi = 900)


total_IF_data <- plot_data %>%
  filter(Date == max(Date))


p9 <- perform_analysis(total_IF_data, "cum_inflo", "Total inflorescence (count)", 0.6, 1, 0)
p9


# ggsave("parental_Total_inflorescence.png", width = 3.5, height = 3.5, dpi = 900)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_mapping_population/Drought_exp/parental_lines/LICOR/clean")



data_0916 <- read_xlsx("20210916_.xlsx")
data_0917 <- read_xlsx("20210917_.xlsx")
data_0918 <- read_xlsx("20210918_.xlsx")
data_0919 <- read_xlsx("20210919_.xlsx")
data_0922 <- read_xlsx("20210922_.xlsx")

LICOR_data <- bind_rows(data_0916, data_0917, data_0918, data_0919, data_0922) %>% 
  mutate(WUE = Photo/Trmmol) %>% # umol CO2/mol H2O
  group_by(ID, Accession, Trt, Bio_rep) %>%
  summarise_all(.funs = mean) %>%
  select(-Tec_rep, -Obs)

head(LICOR_data)

p10 <- perform_analysis(data = LICOR_data, response_variable = "Photo", Y_label = expression(P[n]~(mu*mol~m^{-2}~s^{-2})), 0.15, 1, 0)
p10

# ggsave("parental_Pn_20240326.png", width = 3.5, height = 3.5, dpi = 900)


p11 <- perform_analysis(data = LICOR_data, response_variable = "Cond", Y_label = expression(G[s]~(mol~m^{-2}~s^{-1})), 0.15, 1, 0)
p11

# ggsave("parental_Gs_20240326.png", width = 3.5, height = 3.5, dpi = 900)

p12 <- perform_analysis(data = LICOR_data, response_variable = "WUE", Y_label = expression(WUE~(mu*mol~CO[2]~mol^{-1}~H[2]*O)), 0.15, 1, 0)
p12

# ggsave("parental_WUE_20240326.png", width = 3.5, height = 3.5, dpi = 900)

setwd("/Users/kuowenhsi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/Drought_F3_paper")

p_comb <- plot_grid(plotlist = lapply(paste0("p", 1:12), get), align = "hv", nrow = 4, ncol = 3, labels = paste0("(", LETTERS[1:12], ")"), label_y = 1.05)

p_comb <- ggdraw(plot_grid(
  plotlist = lapply(paste0("p", 1:12), get),
  align = "hv",
  nrow = 4,
  ncol = 3,
  labels = paste0("(", LETTERS[1:12], ")"),
  label_y = 1.05
)) + theme(plot.margin = margin(t = 20, unit = "pt"), plot.background = element_rect(fill = "white", color = NA))
ggsave("Figure_1_comb_20250511.png", width = 10.5, height = 13.4, dpi = 900)
