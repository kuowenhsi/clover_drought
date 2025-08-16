library(tidyverse)
library(grid)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper/clover_drought")

### Leaf area and cyanogenesis

lf_data_cyano <- read_tsv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_mapping_population/Drought_exp/combined_imputed_leaf_area_20230518.tsv")%>%
  mutate(subgroup = rep(rep(c(1, 2), each = 1500), length = 9000))%>%
  select(c(2,5,6, 9, 11,12,13))%>%
  left_join(select(filter(., exp_week == 1), unique_F3_id, Trt, subgroup, leaf_area_w1 = leaf_area, Rep), by = c("unique_F3_id", "Trt", "subgroup", "Rep"))%>%
  mutate(Relative_growth = leaf_area - leaf_area_w1, exp_week = paste("w", exp_week, sep = "_"))%>%
  select(Genotype = unique_F3_id, Trt, Subgroup = subgroup, Replication = Rep, exp_week, Relative_growth)%>%
  filter(exp_week != "w_1")%>%
  left_join(other_pheno, by = "Genotype")%>%
  mutate(Cyanotype = case_when((Ac == 1)&(Li == 1) ~ "AcLi",
                               (Ac == 1)&(Li == 0) ~ "Acli",
                               (Ac == 0)&(Li == 1) ~ "acLi",
                               (Ac == 0)&(Li == 0) ~ "acli"))%>%
  mutate(Replication = rep(1:6, each = 1200))%>%
  mutate(Subgroup = as.factor(Subgroup), Replication = as.factor(Replication))



cyano_model_w2 <- lmerTest::lmer(Relative_growth ~ 1 + Cyanotype + Trt + Cyanotype:Trt + (1|Replication), data = filter(lf_data_cyano, exp_week == "w_2"))
summary(cyano_model_w2)
anova_result <- anova(cyano_model_w2, type = "II")
anova_result
write.csv(anova_result, "anova_result_table.csv", row.names = FALSE)

cyano_model_w3 <- lmerTest::lmer(Relative_growth ~ 1 + Cyanotype + Trt + Cyanotype:Trt + (1|Subgroup) + (1|Subgroup:Replication), data = filter(lf_data_cyano, exp_week == "w_3"))
summary(cyano_model_w3)
anova_result <- anova(cyano_model_w3, type = "II")
write.csv(anova_result, "anova_result_table.csv", row.names = FALSE)

cyano_model_w4 <- lmerTest::lmer(Relative_growth ~ 1 + Cyanotype + Trt + Cyanotype:Trt + (1|Subgroup) + (1|Subgroup:Replication), data = filter(lf_data_cyano, exp_week == "w_4"))
summary(cyano_model_w4)
anova_result <- anova(cyano_model_w4, type = "II")
write.csv(anova_result, "anova_result_table.csv", row.names = FALSE)

cyano_model_w5 <- lmerTest::lmer(Relative_growth ~ 1 + Cyanotype + Trt + Cyanotype:Trt + (1|Subgroup) + (1|Subgroup:Replication), data = filter(lf_data_cyano, exp_week == "w_5"))
summary(cyano_model_w5)
anova_result <- anova(cyano_model_w5, type = "II")
write.csv(anova_result, "anova_result_table.csv", row.names = FALSE)


library(ggh4x)

unique(fl_cyano_data$Subgroup)
unique(fl_cyano_data$Replication)

p <- ggplot(data = filter(lf_data_cyano, exp_week == "w_4"), aes(x = interaction(Trt, Cyanotype), y = Relative_growth))+
  geom_point(aes(group = Cyanotype),size = 0.3, color = "gray65")+
  geom_violin(aes(group = interaction(Trt, Cyanotype), fill = Cyanotype), color = NA, alpha = 0.7, scale = "width")+
  geom_boxplot(aes(group = interaction(Trt, Cyanotype)),width = 0.1, fill = NA, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_brewer(palette = "Set2")+
  facet_wrap2(.~ Replication, nrow = 2, strip = strip_themed(background_x = elem_list_rect(fill = c("#ABEBC6", "#F9E79F"), color =c("#ABEBC6", "#F9E79F"))),
              labeller = as_labeller(c(`1` = 'S1R1 - 2021-09-30', `2` = 'S2R1 - 2021-11-11', `3` = 'S1R2 - 2021-12-24', `4` = 'S2R2 - 2022-02-18', `5` = 'S1R3 - 2022-04-01', `6` = 'S2R3 - 2022-05-12')))+
  theme_bw()+
  xlab("")+
  ylab(expression("Leaf area ("*"mm"^2*")"))+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.8), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))

p

ggsave("leaf_area_cyanotype_rep_w4.png", width = 10, height = 8)


# Convert the plot to a grob
g <- ggplotGrob(p)

g$grobs[[2]]$children$panel.border..rect.10780$gp$col <- "#ABEBC6"
g$grobs[[3]]$children$panel.border..rect.10793$gp$col <- "#F9E79F"
g$grobs[[4]]$children$panel.border..rect.10806$gp$col <- "#ABEBC6"
g$grobs[[5]]$children$panel.border..rect.10819$gp$col <- "#F9E79F"
g$grobs[[6]]$children$panel.border..rect.10832$gp$col <- "#ABEBC6"
g$grobs[[7]]$children$panel.border..rect.10845$gp$col <- "#F9E79F"


# Draw the plot with customized panel borders
png("leaf_area_cyanotype_rep.png", width = 10, height = 8, units = "in", res = 600)
grid.draw(g)
dev.off()


##########
# Cyanogenesis and flower count

##################
fl_cyano_data <- left_join(fl_data, other_pheno, by = "Genotype")%>%
  mutate(Cyanotype = case_when((Ac == 1)&(Li == 1) ~ "AcLi",
                               (Ac == 1)&(Li == 0) ~ "Acli",
                               (Ac == 0)&(Li == 1) ~ "acLi",
                               (Ac == 0)&(Li == 0) ~ "acli"))%>%
  select(c(1:5,6), Cyanotype)%>%
  mutate(Replication = rep(1:6, each = 150))%>%
  pivot_longer(names_to = "Trt", cols = c("Control", "Drought"), values_to = "flower_count")%>%
  mutate(Subgroup = as.factor(Subgroup), Replication = as.factor(Replication))


cyano_model <- lmerTest::lmer(flower_count ~ 1 + Cyanotype + Trt + Cyanotype:Trt + (1|Replication), data = fl_cyano_data)
summary(cyano_model)
anova_result <- anova(cyano_model, type = "II")
anova_result
write.csv(anova_result, "anova_result_table.csv", row.names = FALSE)

library(ggh4x)

unique(fl_cyano_data$Subgroup)
unique(fl_cyano_data$Replication)

p <- ggplot(data = fl_cyano_data, aes(x = interaction(Trt, Cyanotype), y = flower_count))+
  geom_point(aes(group = Cyanotype),size = 0.3, color = "gray65")+
  geom_violin(aes(group = interaction(Trt, Cyanotype), fill = Cyanotype), color = NA, alpha = 0.7, scale = "width")+
  geom_boxplot(aes(group = interaction(Trt, Cyanotype)),width = 0.1, fill = NA, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_brewer(palette = "Set2")+
  facet_wrap2(.~ Replication, nrow = 2, strip = strip_themed(background_x = elem_list_rect(fill = c("#ABEBC6", "#F9E79F"), color =c("#ABEBC6", "#F9E79F"))),
              labeller = as_labeller(c(`1` = 'S1R1 - 2021-09-30', `2` = 'S2R1 - 2021-11-11', `3` = 'S1R2 - 2021-12-24', `4` = 'S2R2 - 2022-02-18', `5` = 'S1R3 - 2022-04-01', `6` = 'S2R3 - 2022-05-12')))+
  theme_bw()+
  xlab("")+
  ylab("Inflorescence count")+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.8), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))

p

ggsave("flower_count_cyanotype_rep.png", width = 10, height = 8)


library(grid)
# Convert the plot to a grob
g <- ggplotGrob(p)

g$grobs[[2]]$children$panel.border..rect.6475$gp$col <- "#ABEBC6"
g$grobs[[3]]$children$panel.border..rect.6488$gp$col <- "#F9E79F"
g$grobs[[4]]$children$panel.border..rect.6501$gp$col <- "#ABEBC6"
g$grobs[[5]]$children$panel.border..rect.6514$gp$col <- "#F9E79F"
g$grobs[[6]]$children$panel.border..rect.6527$gp$col <- "#ABEBC6"
g$grobs[[7]]$children$panel.border..rect.6540$gp$col <- "#F9E79F"


# Draw the plot with customized panel borders
png("flower_count_cyanotype_rep.png", width = 10, height = 8, units = "in", res = 600)
grid.draw(g)
dev.off()


####################
# Cyanogenesis and dry weight


dw_cyano_data <- read_xlsx("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_mapping_population/Drought_exp/total_DW_20220625.xlsx")%>%
  mutate(plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300))%>%
  mutate(Trt_group_rep = paste(Trt, plant_group, Rep, sep = "_"))%>%
  mutate(shoot_w = Shoot_gross - Shoot_bag, root_w = Root_gross - Root_bag)%>%
  mutate(root_w = case_when(root_w < -0.1 ~ as.numeric(NA), TRUE ~ root_w))%>%
  select(Genotype = unique_F3_id, Trt, Subgroup = plant_group, Replication = Rep, shoot_w, root_w)%>%
  mutate(Subgroup = as.factor(Subgroup), Replication = as.factor(Replication))%>%
  left_join(other_pheno, by = "Genotype")%>%
  mutate(Cyanotype = case_when((Ac == 1)&(Li == 1) ~ "AcLi",
                               (Ac == 1)&(Li == 0) ~ "Acli",
                               (Ac == 0)&(Li == 1) ~ "acLi",
                               (Ac == 0)&(Li == 0) ~ "acli"))

str(dw_cyano_data)
table(filter(dw_cyano_data, Replication == 1, Trt == "Control")$Cyanotype) + table(filter(dw_cyano_data, Replication == 2, Trt == "Control")$Cyanotype)

# Full model
cyano_model <- lmer(root_w ~ 1 + Cyanotype + (1|Replication), data = dw_cyano_data)

# Reduced model (without the fixed effect of interest)
cyano_model_reduced <- lmer(shoot_w ~ 1 + (1|Subgroup) + (1|Subgroup:Replication), data = dw_cyano_data)

# Perform the likelihood ratio test
anova(cyano_model_reduced, cyano_model)



cyano_model <- lmerTest::lmer(shoot_w ~ 1 + Cyanotype + Trt + Cyanotype:Trt + (1|Replication), data = dw_cyano_data)
summary(cyano_model)
VarCorr(cyano_model)
anova_result <- anova(cyano_model, type = "II")
anova_result
write.csv(anova_result, "anova_result_table.csv", row.names = FALSE)

p <- ggplot(data = dw_cyano_data, aes(x = interaction(Trt, Cyanotype), y = shoot_w))+
  geom_point(aes(group = Cyanotype),size = 0.3, color = "gray65")+
  geom_violin(aes(group = interaction(Trt, Cyanotype), fill = Cyanotype), color = NA, alpha = 0.7, scale = "width")+
  geom_boxplot(aes(group = interaction(Trt, Cyanotype)),width = 0.1, fill = NA, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_brewer(palette = "Set2")+
  facet_wrap2(.~ Replication, nrow = 2, strip = strip_themed(background_x = elem_list_rect(fill = c("#ABEBC6", "#F9E79F"), color =c("#ABEBC6", "#F9E79F"))),
              labeller = as_labeller(c(`1` = 'S1R1 - 2021-09-30', `2` = 'S2R1 - 2021-11-11', `3` = 'S1R2 - 2021-12-24', `4` = 'S2R2 - 2022-02-18', `5` = 'S1R3 - 2022-04-01', `6` = 'S2R3 - 2022-05-12')))+
  theme_bw()+
  xlab("")+
  ylab("Shoot dry mass (g)")+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.8), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))

p

ggsave("shoot_dw_cyanotype_rep.png", width = 10, height = 8)

library(grid)
# Convert the plot to a grob
g <- ggplotGrob(p)

g$grobs[[2]]$children$panel.border..rect.60047$gp$col <- "#ABEBC6"
g$grobs[[3]]$children$panel.border..rect.60060$gp$col <- "#F9E79F"
g$grobs[[4]]$children$panel.border..rect.60073$gp$col <- "#ABEBC6"
g$grobs[[5]]$children$panel.border..rect.60086$gp$col <- "#F9E79F"
g$grobs[[6]]$children$panel.border..rect.60099$gp$col <- "#ABEBC6"
g$grobs[[7]]$children$panel.border..rect.60112$gp$col <- "#F9E79F"


# Draw the plot with customized panel borders
png("shoot_dw_cyanotype_rep.png", width = 10, height = 8, units = "in", res = 600)
grid.draw(g)
dev.off()


####

p <- ggplot(data = dw_cyano_data, aes(x = interaction(Trt, Cyanotype), y = root_w))+
  geom_point(aes(group = Cyanotype),size = 0.3, color = "gray65")+
  geom_violin(aes(group = interaction(Trt, Cyanotype), fill = Cyanotype), color = NA, alpha = 0.7, scale = "width")+
  geom_boxplot(aes(group = interaction(Trt, Cyanotype)),width = 0.1, fill = NA, outlier.shape = NA, show.legend = FALSE)+
  scale_fill_brewer(palette = "Set2")+
  facet_wrap2(.~ Replication, nrow = 2, strip = strip_themed(background_x = elem_list_rect(fill = c("#ABEBC6", "#F9E79F"), color =c("#ABEBC6", "#F9E79F"))),
              labeller = as_labeller(c(`1` = 'S1R1 - 2021-09-30', `2` = 'S2R1 - 2021-11-11', `3` = 'S1R2 - 2021-12-24', `4` = 'S2R2 - 2022-02-18', `5` = 'S1R3 - 2022-04-01', `6` = 'S2R3 - 2022-05-12')))+
  theme_bw()+
  xlab("")+
  ylab("Root dry mass (g)")+
  theme(panel.border = element_rect(colour = "black", linewidth = 1), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))

p

ggsave("root_dw_cyanotype_rep.png", width = 10, height = 8)

library(grid)
# Convert the plot to a grob
g <- ggplotGrob(p)

g$grobs[[2]]$children$panel.border..rect.65353$gp$col <- "#ABEBC6"
g$grobs[[3]]$children$panel.border..rect.65366$gp$col <- "#F9E79F"
g$grobs[[4]]$children$panel.border..rect.65379$gp$col <- "#ABEBC6"
g$grobs[[5]]$children$panel.border..rect.65392$gp$col <- "#F9E79F"
g$grobs[[6]]$children$panel.border..rect.65405$gp$col <- "#ABEBC6"
g$grobs[[7]]$children$panel.border..rect.65418$gp$col <- "#F9E79F"


# Draw the plot with customized panel borders
png("root_dw_cyanotype_rep.png", width = 10, height = 8, units = "in", res = 600)
grid.draw(g)
dev.off()
