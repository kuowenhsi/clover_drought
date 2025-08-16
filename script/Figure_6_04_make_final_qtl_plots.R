library(tidyverse)
library(lme4)



# ploting the final combined plots

total_LOD_tb <- bind_rows(
  read_csv("./data/leaf_area_LOD_tb_20240303.csv") %>% mutate(trait = (paste("leaf_area", trait, sep = "_"))),
  read_csv("./data/flower_count_LOD_tb_20240303.csv") %>% mutate(trait = paste("flower_count", trait, sep = "_")),
  read_csv("./data/dry_weight_LOD_tb_20240303.csv")%>% mutate(trait = str_replace(trait, "w", "weight"))
)%>%
  mutate(trait = str_remove(str_remove(trait, "_Control"), "_Drought"))%>%
  mutate(trait = factor(trait, levels = c("Genetic distance", "flower_count", "leaf_area_w_2", "leaf_area_w_3", "leaf_area_w_4", "leaf_area_w_5", "shoot_weight",  "root_weight")))


geneticMapList_tb_s <- filter(geneticMapList_tb, mid_chr %in% total_LOD_tb$mid_chr)

p <- ggplot(data = total_LOD_tb, aes(x = trait, y = mid_pos))+
  geom_line(data = geneticMapList_tb_s, color = "gray60")+
  geom_point(data = geneticMapList_tb_s, shape = "-", size = 10, color = "gray40")+
  geom_rect(data = data.frame(mid_chr = "drTriRepe4Chr4", trait = "Genetic distance", mid_pos = 1250), xmin = -Inf , xmax = Inf, ymin = 1355.5686, ymax = 1576.4467, fill = "blue", alpha = .2)+
  geom_segment(data = filter(total_LOD_tb, arrow_head == "first"), aes(xend = trait, y = low_pos, yend = high_pos, color = treat),
               arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "first"), lineend = "butt", linewidth = 1,  show.legend = FALSE)+
  geom_segment(data = filter(total_LOD_tb, arrow_head == "last"), aes(xend = trait, y = low_pos, yend = high_pos, color = treat),
               arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last"), lineend = "butt", linewidth = 1, show.legend = FALSE)+
  geom_point(aes(color = treat))+
  scale_color_manual(values = c("#00BFC4", "#F8766D"))+
  facet_wrap(.~mid_chr, nrow = 1)+
  xlab("")+
  ylab("Genetic distance (cM)")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), legend.position = c(0.9, 0.9), legend.title = element_blank(), plot.background = element_rect(fill = "white", color = NA))

p  
ggsave("./figures/Total_BLUP_p0.1_20240304.png", width = 8, height = 6)

plotPXG(DG_F3_data, "drTriRepe4Chr4:42683586", pheno.col=6)

for (i in 1:nrow(LOD_tb)){
  p <- get_effect_plot(DG_F3_data, LOD_tb$mid_marker[[i]], Trt1 = LOD_tb$trait[[i]], pos = LOD_tb$mid_pos[[i]], lod = LOD_tb$mid_lod[[i]], y_label = "Leaf area BLUPs")
  # ggsave(paste("leaf_area_BLUP_p0.1_", LOD_tb$trait[[i]], "_",str_replace(LOD_tb$mid_marker[[i]], ":", "."), ".png", sep = ""), width = 5, height = 5)
}






#########################

get_effect_plot(DG_F3_data, LOD_tb$mid_marker[[1]], Trt1 = LOD_tb$trait[[1]], pos = LOD_tb$mid_pos[[1]], lod = LOD_tb$mid_lod[[1]])

library(ggplot2)
library(scales) # For show_col()

# Generate a factor with a specified number of levels
num_levels <- 2
factor_levels <- factor(1:num_levels)

# Use the default ggplot2 color palette to get colors for these levels
colors <- hue_pal()(num_levels)

# Display the colors
show_col(colors)
