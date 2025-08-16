library(qtl)
library(tidyverse)
library(qtltools)
library(readxl)
library(future.apply)
library(progressr)
library(ASMap)
library(lme4)
library(ggh4x)

# Set up future to use multicore parallelism
plan(multisession)

# Create a progress handler
handlers(global = TRUE)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper/clover_drought")


DG_F3_data <- read.cross("csv", "/Users/kuowenhsi/OneDrive - Washington University in St. Louis/leaf_mark_paper/VCF_bwa/", "bwa_DG_F3_DP0_F3_m0.5_hardfilter_pre_imputed_homo_maf0.2_hww0.01_LM.csv", estimate.map=FALSE,
                         genotypes=c("AA", "AB", "BB"), alleles = c("A", "B"))

###############################

other_data <- read.cross(file = "DG_F3_physical_order_LM_pheno_20240219.csv", dir = "/Users/kuowenhsi/OneDrive - Washington University in St. Louis/leaf_mark_paper/",format = "csv", estimate.map=FALSE, genotypes=c("AA", "AB", "BB"), alleles = c("A", "B"))

other_pheno <- other_data$pheno%>%select(1:4)

###############################
summaryMap(DG_F3_data)

plotMap(DG_F3_data)

DG_F3_data <- jittermap(DG_F3_data)

DG_F3_data <- DG_F3_data %>%
  calc.genoprob(step=0, error.prob=0.001,map.function=c("kosambi")) %>%
  sim.geno(step=0, n.draws=16, error.prob=0.001)

dw_data <- read_xlsx("./data/total_DW_20220625.xlsx")%>%
  mutate(plant_group = rep(c(1, 2, 1, 2, 1, 2), each = 300))%>%
  mutate(Trt_group_rep = paste(Trt, plant_group, Rep, sep = "_"))%>%
  mutate(shoot_w = Shoot_gross - Shoot_bag, root_w = Root_gross - Root_bag)%>%
  select(Genotype = unique_F3_id, Trt, Replication = Rep, drymass_shoot = shoot_w, drymass_root = root_w)%>%
  mutate(drymass_root = drymass_root + abs(min(drymass_root, na.rm = TRUE)))%>%
  # mutate(shoot_w = log(shoot_w + 1), root_w = log(root_w + 1))%>%
  pivot_wider(names_from = c("Trt"), values_from = c("drymass_shoot", "drymass_root"))%>%
  mutate(Replication = as.factor(Replication))

write_csv(dw_data, "./data/drymass_shootroot_raw_20240826.csv")

str(dw_data)

# List of phenotype columns
phenotype_cols <- names(dw_data)[!(names(dw_data) %in% c("Genotype", "Subgroup", "Replication"))]

# Initialize an empty list to store BLUPs for each phenotype
blups_list <- list()

# Loop over each phenotype column
for (phenotype in phenotype_cols) {
  # phenotype = "Control_w_2"
  
  # Fit the mixed model
  model <- lmer(as.formula(paste(phenotype, "~ (1|Replication) + (1|Genotype)")), data = dw_data)
  
  # Calculate BLUPs for genotypes
  blups <- ranef(model)[["Genotype"]]
  
  # Convert the BLUPs to a data frame and reset row names
  blups_df <- as.data.frame(blups)
  # rownames(blups_df) <- lf_data_blup$Genotype
  names(blups_df) <- phenotype
  
  # Store the BLUPs data frame in the list, named by phenotype
  blups_list[[phenotype]] <- blups_df
}

# Combine all BLUPs data frames into a single data frame
# Assuming that the genotype order is consistent across all BLUPs calculations
combined_blups <- bind_cols(blups_list)%>%
  rownames_to_column(var = "Genotype")

write_csv(combined_blups, "./data/dry_weight_log_BLUPs_20240826.csv")
#######################


DG_F3_pheno <- DG_F3_data$pheno %>%
  left_join(combined_blups, by = c("Genotype"))%>%
  left_join(other_pheno, by = "Genotype")


DG_F3_data$pheno <- DG_F3_pheno

summary(DG_F3_data)
names(DG_F3_pheno)[2:5]

qtl_results_list <- list()
for (i in names(DG_F3_pheno)[2:5]){
  qtl_results_list[[i]] <- scanone(DG_F3_data, pheno.col=i, model="normal", method = "hk")
  plot(qtl_results_list[[i]], ylab="LOD score", alternate.chrid=TRUE, main = paste0("Biomass ", i))
}



length(qtl_results_list)

qtl_sig_list <- list()
data_list <- names(DG_F3_pheno)[2:5]
# Set up the progress bar
with_progress({
  p <- progressor(along = data_list)
  
  # Use future_lapply to apply the function in parallel with progress updates
  qtl_sig_list <- future_lapply(data_list, function(x) {
    p()  # Update the progress bar
    scanone(DG_F3_data, pheno.col=x, model="normal", method = "hk", n.perm = 1000)
  })
})

names(qtl_sig_list) <- names(DG_F3_pheno)[2:5]

length(qtl_sig_list)
lapply(qtl_results_list, summary)
lapply(qtl_sig_list, summary)

qtl_aboveSig_list <- list()
for (i in 1:4){
  qtl_aboveSig_list[[i]] <- summary(qtl_results_list[[i]], perms=qtl_sig_list[[i]], alpha=0.1, pvalues=TRUE)
}

names(qtl_aboveSig_list) <- names(DG_F3_pheno)[2:5]
length(qtl_aboveSig_list)
qtl_aboveSig_list

#############################


LOD_list <- vector(mode = "list", length = 4)
names(LOD_list) <- names(DG_F3_pheno)[2:5]
for (i in names(LOD_list)){
  
  significant_peaks <- qtl_aboveSig_list[[i]]
  print(i)
  if (nrow(significant_peaks) == 0) {
    LOD_list[[i]] = NULL
    next}
  # Loop through each significant peak
  ci_results <- tibble(low_marker=NULL, mid_marker=NULL, high_marker=NULL, low_chr=NULL, mid_chr=NULL, high_chr=NULL, low_pos=NULL, mid_pos=NULL, high_pos=NULL, low_lod=NULL, mid_lod=NULL, high_lod=NULL, additive_effect = NULL)
  for (h in 1:nrow(significant_peaks)) {
    
    # Extract chromosome and position for the current peak
    current_chr <- significant_peaks$chr[h]
    current_pos <- significant_peaks$pos[h]
    
    # Calculate the 1-LOD drop confidence interval for the current peak
    ci_result <- lodint(qtl_results_list[[i]], chr=current_chr, drop=1)
    
    qtl_object <- makeqtl(DG_F3_data, chr=current_chr, pos=current_pos, what = "prob")
    additive_effect <- (fitqtl(DG_F3_data, pheno.col = i, qtl = qtl_object, method = "hk", get.ests = TRUE))$ests$ests[[2]]
    
    
    # Store the result
    ci_results <- bind_rows(ci_results, tibble(low_marker = rownames(ci_result)[[1]], mid_marker = rownames(ci_result)[[2]], high_marker=rownames(ci_result)[[3]], low_chr=ci_result$chr[[1]], mid_chr=ci_result$chr[[2]], high_chr=ci_result$chr[[3]], low_pos=ci_result$pos[[1]], mid_pos=ci_result$pos[[2]], high_pos=ci_result$pos[[3]], low_lod=ci_result$lod[[1]], mid_lod=ci_result$lod[[2]], high_lod=ci_result$lod[[3]], additive_effect = additive_effect))
    
  }
  
  LOD_list[[i]] <- ci_results
  
}
length(LOD_list)
LOD_list 

geneticMapList <- list()

# Extract the genetic map for each chromosome
for (chr in names(DG_F3_data$geno)) {
  # Extract marker names and genetic distances for the current chromosome
  markers <- names(DG_F3_data$geno[[chr]]$map)
  distances <- DG_F3_data$geno[[chr]]$map
  
  # Create a data frame for the current chromosome
  df <- data.frame(Marker = markers, GeneticDistance = distances)
  
  # Add the data frame to the list
  geneticMapList[[chr]] <- df
}
geneticMapList[[1]]

geneticMapList_tb <- bind_rows(geneticMapList, .id = "chr")%>%
  select(c(1,3))%>%
  rename(mid_pos = GeneticDistance, mid_chr = chr)%>%
  mutate(mid_pos = as.double(mid_pos), trait = "Genetic distance")%>%
  remove_rownames()

###############

LOD_tb <- bind_rows(LOD_list, .id = "trait")%>%
  filter(trait != "leaf_mark_q") %>%
  mutate(arrow_head = case_when(additive_effect < 0 ~ "first", additive_effect > 0 ~ "last"),
         treat = str_split_i(trait, "_", 3))%>%
  # mutate(trait = str_remove(str_remove(str_remove(trait, "_w_Control"), "_w_Drought"), "area_"))%>%
  # mutate(trait = paste("S", str_split_i(trait, "_", 2), "_R", str_split_i(trait, "_", 3), "_", str_split_i(trait, "_", 1), sep = ""))%>%
  mutate(trait = factor(trait))

x_labels <- LOD_tb %>%
  mutate(trait = str_remove(str_remove(str_remove(trait, "_w_Control"), "_w_Drought"), "area_"))%>%
  mutate(trait = paste("S", str_split_i(trait, "_", 2), "_R", str_split_i(trait, "_", 3), "_", str_split_i(trait, "_", 1), sep = ""))%>%
  pull(trait)

geneticMapList_tb_s <- filter(geneticMapList_tb, mid_chr %in% LOD_tb$mid_chr)%>%
  mutate(mid_chr = factor(mid_chr, levels = paste0("drTriRepe4Chr", 1:16)))

write_csv(LOD_tb, "./data/dry_weight_log_LOD_tb_20240826.csv")

p <- ggplot(data = LOD_tb, aes(x = trait, y = mid_pos))+
  geom_line(data = geneticMapList_tb_s, color = "gray60")+
  geom_point(data = geneticMapList_tb_s, shape = "-", size = 10, color = "gray40")+
  geom_segment(data = filter(LOD_tb, arrow_head == "first"), aes(xend = trait, y = low_pos, yend = high_pos, color = treat),
               arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "first"), lineend = "butt", linewidth = 1,  show.legend = FALSE)+
  geom_segment(data = filter(LOD_tb, arrow_head == "last"), aes(xend = trait, y = low_pos, yend = high_pos, color = treat),
               arrow = arrow(type = "closed", length = unit(0.15, "inches"), ends = "last"), lineend = "butt", linewidth = 1, show.legend = FALSE)+
  geom_point(aes(color = treat))+
  scale_color_manual(values = c("#00BFC4", "#F8766D"))+
  # scale_x_discrete(labels = c(x_labels[1:2], "", x_labels[3], "", x_labels[4],"",x_labels[5],"",x_labels[6]))+
  facet_grid(.~mid_chr, scales = "free_x", space = "free_x")+
  xlab("")+
  ylab("Genetic distance (cM)")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), legend.position = c(0.9, 0.9), legend.title = element_blank(), plot.background = element_rect(fill = "white", color = NA))

p 
ggsave("./figures/dry_weight_log_BLUP_p0.1_20240826.png", width = 6, height = 6)



get_effect_plot <- function(cross_object = DG_F3_data, marker, Trt1, Trt2 = NULL, pos = mid_pos, lod = mid_lod, y_label = "pheno") {
  
  # cross_object = DG_F3_data
  # marker=LOD_tb$mid_marker[[1]]
  # Trt1 = LOD_tb$trait[[1]]
  # Trt2 = NULL
  # pos = LOD_tb$mid_pos[[1]]
  # lod = LOD_tb$mid_lod[[1]]
  
  if (is.null(Trt2)){
    if (str_detect(Trt1, "Control")){
      Trt2 = str_replace(Trt1, "Control", "Drought")
    }
    if (str_detect(Trt1, "Drought")){
      Trt2 = str_replace(Trt1, "Drought", "Control")
    }
  }
  
  Trt1 <- as.character(Trt1)
  Trt2 <- as.character(Trt2)
  
  allele_effect_Trt1 <- plotPXG(cross_object,marker, pheno.col=Trt1) %>% drop_na() %>% mutate(Trt = Trt1)
  allele_effect_Trt2 <- plotPXG(cross_object,marker, pheno.col=Trt2) %>% drop_na() %>% mutate(Trt = Trt2)
  allele_effect <- bind_rows(allele_effect_Trt1, allele_effect_Trt2)
  
  model_effect <- lm(pheno ~ get(marker) + Trt + get(marker):Trt, data = allele_effect)
  ANOVA_test <- car::Anova(model_effect, type = "III")
  value_str <- format(sprintf("%.2e", ANOVA_test$`Pr(>F)`[[4]]), scientific = TRUE)
  
  # Extract the coefficient and exponent parts
  parts <- strsplit(value_str, "e")[[1]]
  coefficient <- parts[1]
  exponent <- as.integer(parts[2])
  expr <- paste("Interaction:", deparse(bquote(italic(p) == .(coefficient) %*% 10^.(exponent))))
  
  sample_size <- allele_effect %>%
    group_by(get(marker))%>%
    count()%>%
    rename(marker = "get(marker)")%>%
    ungroup()
  
  p <- ggplot(data = allele_effect, aes(x = get(marker), y = pheno))+
    geom_point(aes(color = Trt), position = position_dodge(width = 0.2), alpha = 0.8)+
    stat_summary(geom = "line", fun = "mean", aes(color = Trt, group = Trt),position = position_dodge(width = 0.2), show.legend = FALSE)+
    geom_boxplot(aes(color = Trt), position = position_dodge(width = 0.5), width = 0.1, show.legend = FALSE, outlier.shape = NA)+
    stat_summary(geom = "point", fun = "mean", aes(color = Trt), shape = 18, size = 6,position = position_dodge(width = 0.2), show.legend = FALSE)+ 
    geom_text(x = 0.5, y = max(allele_effect$pheno)*0.95, label = expr, parse = TRUE, size = 3.5, fontface = "plain",check_overlap = TRUE, hjust = 0)+
    scale_x_discrete(name = marker, labels = c(paste0("AA (", sample_size$n[[1]], ")"), paste0("AB (", sample_size$n[[2]], ")"), paste0("BB (", sample_size$n[[3]], ")")))+
    scale_y_continuous(oob = scales::oob_squish)+
    scale_color_manual(values = c("#00BFC4", "#F8766D"))+
    ylab(y_label)+
    ggtitle(paste(Trt1,"Pos =", round(pos, digits = 2), "LOD =", round(lod, digits = 2), sep = " "))+
    theme_bw()+
    theme(legend.position = c(0.8,0.95), legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 12), legend.background = element_rect(fill = NA))
  
  return(p)
}


for (i in 1:nrow(LOD_tb)){
  p <- get_effect_plot(DG_F3_data, LOD_tb$mid_marker[[i]], Trt1 = LOD_tb$trait[[i]], pos = LOD_tb$mid_pos[[i]], lod = LOD_tb$mid_lod[[i]], y_label = "Dry weight BLUPs")
  ggsave(paste("./figures/dry_weight_log_BLUP_p0.1_", LOD_tb$trait[[i]], "_",str_replace(LOD_tb$mid_marker[[i]], ":", "."), ".png", sep = ""), width = 5, height = 5)
}
