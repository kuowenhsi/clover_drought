library(tidyverse)
library(data.table)
library(qvalue)
library(cowplot)
library(ggrepel)
library(IRanges)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper")

################

chr_len_temp <- read_tsv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_mapping_population/F2_QTL_mapping/output/HiFi_HiC_LM_combined_v_1.5.fasta.fai", col_names = c("chr", "end")) %>%
  select(1:2)%>%
  mutate(start = 1)%>%
  filter(!str_detect(chr, "drTriRepe4Chr0c"))%>%
  mutate(chr = as.integer(str_remove(chr, "drTriRepe4Chr")))%>%
  arrange(chr) %>%
  mutate(lag_pos = lag(end, default = 0))%>%
  mutate(pos_pad = cumsum(lag_pos))%>%
  mutate(padded_start = start + pos_pad - 1, padded_end = end + pos_pad)%>%
  mutate(padded_chr_pos = (padded_start + padded_end)/2)

#################


one_drop_data <- list()

one_drop_data[[1]] <- read.csv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper/flower_count_LOD_tb_20240303.csv")%>%
  mutate(dataset = "flower_count")

one_drop_data[[2]] <- read.csv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper/dry_weight_LOD_tb_20240303.csv")%>%
  mutate(dataset = "dry_weight")

one_drop_data[[3]] <- read.csv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/Drought_F3_paper/leaf_area_LOD_tb_20240303.csv")%>%
  mutate(dataset = "leaf_area")


one_drop_data <- bind_rows(one_drop_data)
colnames(one_drop_data)

one_drop_data_pRange <- one_drop_data %>%
  mutate(lowmarker_pos = as.integer(str_remove(low_marker, "drTriRepe4Chr.:")), highmarker_pos = as.integer(str_remove(high_marker, "drTriRepe4Chr.:")), chr = as.integer(str_remove(mid_chr, "drTriRepe4Chr")), mid_marker_pos = as.integer(str_remove(mid_marker, "drTriRepe4Chr.:")))%>%
  left_join(chr_len_temp, by = "chr")%>%
  select(dataset, trait, chr, maxLod = mid_lod, lowmarker_pos, highmarker_pos, pos_pad, mid_marker_pos, treat)%>%
  mutate(padded_lowmarker_pos = lowmarker_pos + pos_pad, padded_highmarker_pos = highmarker_pos + pos_pad, ID = seq(nrow(.)),
         padded_mid_marker_pos = mid_marker_pos + pos_pad)%>%
  mutate(chr = factor(chr, levels = 1:16), name = "QTL mapping")
  

p_qtl <- ggplot(data = filter(one_drop_data_pRange, chr == 4), aes(xmin = padded_lowmarker_pos, xmax = padded_highmarker_pos, y = trait))+
  geom_rect(data = filter(chr_len_temp, chr == 4), aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_linerange(aes(color = treat), linewidth = 3)+
  geom_point(aes(x = padded_mid_marker_pos), size = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous(expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  # annotate("segment", x=124910454, y=4.5, xend=124910454, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  # annotate("segment", x=249995270, y=4.5, xend=249995270, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  # annotate("segment", x=582870000, y=4.5, xend=582870000, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  # annotate("segment", x=713443783, y=4.5, xend=713443783, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  # annotate("segment", x=874139073, y=4.5, xend=874139073, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  ylab("")+
  xlab("")+
  scale_color_manual(values = c("#00BFC4", "#F8766D"))+
  scale_y_discrete(labels = c("Leaf area w2", "Leaf area w3", "Leaf area w4", "Leaf area w5", "Flower count", "Root weight", "Shoot weight"))+
  # coord_cartesian(ylim = c(1, 3), clip="off") +
  theme_bw()+
  theme(plot.margin = unit(c(2,0.2,0.2,0.2), "lines"), legend.position = c(0.1, 0.8), legend.title = element_blank(), legend.background = element_rect(fill = NA, color = NA))+
  # theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
  #       axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())+
  guides(fill = "none")

p_qtl

qtl_range <- IRanges(start = one_drop_data_pRange$padded_lowmarker_pos, end = one_drop_data_pRange$padded_highmarker_pos, names = one_drop_data_pRange$pheno_group)
qtl_range


##################

Ac_snps <- read_tsv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/Figure_5/PC10_ADD/bwa_415_20230527_GWAS_PC10.Ac_CNV.glm.linear")%>%
  filter(LOG10_P > -log10(5e-8), `#CHROM` == "chr_02")


Li_snps <- read_tsv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/Figure_5/PC10_ADD/bwa_415_20230527_GWAS_PC10.Li_CNV.glm.linear")%>%
  filter(LOG10_P > -log10(5e-8), `#CHROM` == "chr_12")

AcLi_snp_ID <- c(Ac_snps$ID, Li_snps$ID)


###################

LFMM_PCA_data <- fread("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/LFMM_output_20230606.txt", header = FALSE)

unique(LFMM_PCA_data$V5)
unique(LFMM_PCA_data$V6)

# aridity
LFMM_data_sub <- LFMM_PCA_data[(LFMM_PCA_data$V5 %in% unique(LFMM_PCA_data$V5)[2]) & LFMM_PCA_data$V6 == 3,] %>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V5)%>% ## group by ENV
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))

nrow(filter(LFMM_data_sub, qvalue < 0.05))
LFMM_data_sub_sig <- filter(LFMM_data_sub, qvalue < 0.05)
LFMM_range <- IRanges(start = LFMM_data_sub_sig$padded_pos, end = LFMM_data_sub_sig$padded_pos, names = rep("LFMM", 442))
LFMM_range

lfmm_sig_Aridity_3 <- filter(LFMM_data_sub, qvalue < 0.05) %>% ungroup %>% select(V3) %>% separate(V3, into = c("chr", "position"), sep = ":")
write_tsv(lfmm_sig_Aridity_3, "lfmm_sig_Aridity_3_ID.txt", col_names = FALSE)  
  
  
ppop <- ggplot(data = filter(LFMM_data_sub, V1 == 4), aes(x = padded_pos, y = -log10(qvalue)))+
  geom_rect(data = filter(chr_len_temp, chr == 4), aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_point(color = filter(LFMM_data_sub, V1 == 4)$dot_color, size = 0.5)+
  # geom_point(data = filter(LFMM_data_sub, V3 %in% AcLi_snp_ID), color = "red", size = 0.5)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("Chr. 4 (Mbp)", expand = c(0,0), breaks = 5e6*(1:12) + 174244438, labels = 5*c(1:12))+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))+
  guides(fill = "none")

ppop
###################

PCadapt_data <- read_tsv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/Figure_4/PCadapt_output_20230613.txt")%>%
  left_join(chr_len_temp, by = c("CHR" = "chr"))%>%
  mutate(padded_pos = POSITION + pos_pad, qvalue = qvalue(pvalue)$qvalues)%>%
  mutate(dot_color = case_when(pvalue < 5e-8 ~ "black", TRUE ~ "gray80"))%>%
  mutate(name = "PCadapt")

p_QQ_PCadapt <- ggplot(data = PCadapt_data, aes(sample = -log(pvalue)))+
  stat_qq(aes(x = after_stat(theoretical)/log(10), y = after_stat(sample)/log(10)), distribution = qexp, size = 0.5)+
  stat_qq_line(aes(x = after_stat(x)/log(10), y = after_stat(y)/log(10)), distribution = qexp)+
  geom_hline(yintercept = -log10(5e-8), color = "red")+
  xlab(expression("Theoretical"~"-"*log[10]*"(p)"))+
  ylab(expression("Observed"~"-"*log[10]*"(p)"))+
  labs(title = "Q-Q plot for the PCadapt result")+
  theme_bw()+
  theme(panel.grid = element_blank())

p_QQ_PCadapt
ggsave(filename = "p_QQ_PCadapt.png", width = 5, height = 5, path = "./Figure_4/" )


nrow(filter(PCadapt_data, pvalue < 5e-8))
PCadapt_sig <- filter(PCadapt_data, pvalue < 5e-8)%>%
  select(CHR, POSITION) %>%
  mutate(CHR = paste("chr", str_pad(CHR, 2, "left", "0"), sep = "_"))

write_tsv(PCadapt_sig, "PCadapt_sig_pvalue.txt", col_names = FALSE)

tail(PCadapt_data)

PCadapt_data_sig <- PCadapt_data %>%
  filter(pvalue < 5e-8)

PCadapt_range <- IRanges(start = PCadapt_data_sig$padded_pos, end = PCadapt_data_sig$padded_pos, names = rep("PCadapt", 213))

# write_tsv(PCadapt_data_sig, "PCadapt_data_sig_SNPs_20230622.txt")

pcap <- ggplot(data = filter(PCadapt_data, chr == 4), aes(x = padded_pos, y = -log10(pvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_point(color = PCadapt_data$dot_color, size = 0.5)+
  geom_hline(yintercept = -log10(5e-8), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(p value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))+
  guides(fill = "none")


p_combined <- plot_grid(p_qtl, ppop, align = "v", ncol = 1, rel_heights = c(3, 3), axis = 'l', labels = c("(A)", "(B)"))

ggsave("GEA_aridity_chr4_20240326.png", width = 8, , height = 4.5, path = "./Figure_4/")
########################################
qtl_range
qtl_range_r <- reduce(qtl_range)
QTL_LFMM_intersect <- findOverlaps(qtl_range, LFMM_range)
reduce(qtl_range[unique(QTL_LFMM_intersect@from)])

########### annotated genes #############

# cand_genes <- read_csv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/05_candidate_gene/20230617_lfmm_degreeDay5_result_transfered.csv")%>%
#   filter(!is.na(gene_name))%>%
#   mutate(gene_pos = (hit_start + hit_end)/2, chr = as.integer(str_remove(hit_chr, "drTriRepe4Chr")))%>%
#   group_by(chr, gene_name, description, url)%>%
#   summarize(gene_pos = round(mean(gene_pos)))%>%
#   arrange(chr, gene_pos)
# 
# write_excel_csv(cand_genes, "cand_genes_for_plot.csv")

cand_genes_pad <- read_csv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/05_candidate_gene/cand_genes_lfmm_for_plot_20240109.csv")%>%
  left_join(chr_len_temp, by = "chr")%>%
  mutate(padded_pos = gene_pos + pos_pad, nudge_dist = case_when(HJUST == 0 ~ 50e5, HJUST == 1 ~ -50e5))


# cand_genes_PCadapt <- read_tsv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/05_candidate_gene/20230622_PCadapt_result_transfered.txt")%>%
#   filter(!is.na(gene_name))%>%
#   mutate(gene_pos = (hit_start + hit_end)/2, chr = as.integer(str_remove(hit_chr, "drTriRepe4Chr")))%>%
#   group_by(chr, gene_name, description, url)%>%
#   summarize(gene_pos = round(mean(gene_pos)))%>%
#   arrange(chr, gene_pos)
# 
# write_excel_csv(cand_genes_PCadapt, "cand_genes_PCadapt_for_plot.csv")

cand_genes_PCadapt_pad <- read_csv("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/05_candidate_gene/cand_genes_PCadapt_for_plot_20240109.csv")%>%
  left_join(chr_len_temp, by = "chr")%>%
  mutate(padded_pos = gene_pos + pos_pad, nudge_dist = case_when(HJUST == 0 ~ 50e5, HJUST == 1 ~ -50e5))

cand_genes_pad <- cand_genes_pad %>%
  left_join(cand_genes_PCadapt_pad %>% select(chr, gene_name) %>% mutate(also_PCadapt = "YES"), by = c("chr", "gene_name")) %>%
  mutate(font_color = case_when(also_PCadapt == "YES" ~ "red", TRUE ~ "blue3"))

cand_genes_PCadapt_pad <- cand_genes_PCadapt_pad %>%
  left_join(cand_genes_pad %>% select(chr, gene_name) %>% mutate(also_lfmm = "YES"), by = c("chr", "gene_name")) %>%
  mutate(font_color = case_when(also_lfmm == "YES" ~ "red", TRUE ~ "blue3"))






##########

ppop <- ggplot(data = LFMM_data_sub, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_segment(data = cand_genes_pad, aes(x = padded_pos, xend = padded_pos, yend = Y), y = 0, linewidth = 1, alpha = 0.5, color = "skyblue")+
  geom_text(data = cand_genes_pad, aes(x = padded_pos + nudge_dist, y = Y, label = gene_name, hjust = HJUST), color = cand_genes_pad$font_color, size = 3)+
  geom_point(color = LFMM_data_sub$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_sub, V3 %in% AcLi_snp_ID), color = "red", size = 0.5)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank(),plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))+
  guides(fill = "none")

ppop

ggsave("LFMM_k3_growingDegDays5_20240109.png", width = 10, height = 4)

pcap <- ggplot(data = PCadapt_data, aes(x = padded_pos, y = -log10(pvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_segment(data = cand_genes_PCadapt_pad, aes(x = padded_pos, xend = padded_pos, yend = Y), y = 0, linewidth = 1, alpha = 0.5, color = "skyblue")+
  geom_text(data = cand_genes_PCadapt_pad, aes(x = padded_pos + nudge_dist, y = Y, label = gene_name, hjust = HJUST), color = cand_genes_PCadapt_pad$font_color, size = 3)+
  geom_point(color = PCadapt_data$dot_color, size = 0.5)+
  geom_hline(yintercept = -log10(5e-8), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(p value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))+
  guides(fill = "none")

pcap

p_combined <- plot_grid(p_qtl, ppop,  pcap, align = "v", ncol = 1, rel_heights = c(2, 4, 4), axis = 'l', labels = c("(A)", "(B)", "(C)"))

ggsave("LFMM_test_20240311.png", width = 8, , height = 8, path = "./Figure_4/")

############ per chromosome ##############

for (i in 1:16){

  one_drop_data_pRange_chr <- one_drop_data_pRange %>% filter(chr == i)
  LFMM_data_sub_chr <- LFMM_data_sub %>% filter(V1 == i)
  PCadapt_data_chr <- PCadapt_data %>%filter(CHR == i)
  cand_genes_chr <- cand_genes_pad %>% filter(chr == i)
  cand_genes_chr_pca <- cand_genes_PCadapt_pad %>% filter(chr == i)
  
  print("subset finished")
  print(i)
  
  p_qtl <- ggplot(data = one_drop_data_pRange_chr, aes(xmin = lowmarker_pos, xmax = highmarker_pos, y = pheno_group))+
    geom_linerange(aes(color = pheno_group), linewidth = 3)+
    scale_x_continuous(expand = c(0.1,0.1), limits = c(1, filter(chr_len_temp, chr == i)$end))+
    ylab("")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())+
    guides(fill = "none", color = "none")

  ppop <- ggplot(data = LFMM_data_sub_chr, aes(x = V2, y = -log10(qvalue)))+
    geom_segment(data = cand_genes_chr, aes(x = gene_pos, xend = gene_pos, yend = Y), y = 0, linewidth = 1, alpha = 0.5, color = "skyblue")+
    geom_text(data = cand_genes_chr, aes(x = gene_pos + 5e5, y = Y, label = gene_name, hjust = 0), color = cand_genes_chr$font_color, size = 3)+
    geom_point(color = LFMM_data_sub_chr$dot_color, size = 1)+
    geom_point(data = filter(LFMM_data_sub_chr, V3 %in% AcLi_snp_ID), color = "red", size = 1)+
    # geom_text_repel(data = cand_genes_chr, aes(x = gene_pos, label = gene_name), y = max(-log10(LFMM_data_sub_chr$qvalue), na.rm = TRUE) - 0.3, size = 3)+
    geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
    scale_x_continuous(paste("chr", i), expand = c(0.1,0.1))+
    scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
    theme_bw()+
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())+
    guides(fill = "none")
  
  pcap <- ggplot(data = PCadapt_data_chr, aes(x = POSITION, y = -log10(pvalue)))+
    geom_segment(data = cand_genes_chr_pca, aes(x = gene_pos, xend = gene_pos, yend = Y), y = 0, linewidth = 1, alpha = 0.5, color = "skyblue")+
    geom_text(data = cand_genes_chr_pca, aes(x = gene_pos + 5e5, y = Y, label = gene_name, hjust = 0), color = cand_genes_chr_pca$font_color, size = 3)+
    geom_point(color = PCadapt_data_chr$dot_color, size = 1)+
    geom_hline(yintercept = -log10(5e-8), color = "red", linewidth = 0.5, alpha = 1)+
    scale_x_continuous(paste("chr", i, "(Mbp)"), expand = c(0.1,0.1), breaks = scales::breaks_extended(7), labels = function(x)format(x/1e6, nsmall = 2))+
    scale_y_continuous(expression("-"*log[10]*"(p value)"), expand = c(0,0,0.1,0.1))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    guides(fill = "none")
  
  p_combined <- plot_grid(p_qtl, ppop,  pcap, align = "v", ncol = 1, rel_heights = c(1, 4, 4), axis = 'l', labels = c("(A)", "(B)", "(C)"), hjust = 0)
  
  ggsave(paste0("GEA_20240109_chr", i,".png"), width = 8, , height = 9)
    
}




######################################
####### 16 ENVIRON variables #########
######################################

for (i in 1:3){

LFMM_data_sub <- LFMM_PCA_data[(LFMM_PCA_data$V5 %in% unique(LFMM_PCA_data$V5)[1:16]) & LFMM_PCA_data$V6 == i,] %>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V5)%>% ## group by ENV
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))

tail(LFMM_data_sub)


ppop <- ggplot(data = LFMM_data_sub, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_point(color = LFMM_data_sub$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_sub, V3 %in% AcLi_snp_ID), color = "red", size = 0.5)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none")+
  facet_wrap(.~V5, ncol = 1, scales = "free_x")

ggsave(paste0("LFMM_noBWAnoVBC_k", i, "_16_ENVIRON_20230608.png"), width = 10, height = 30)

}

########## degree 5 days #############

unique(LFMM_PCA_data$V5)

LFMM_data_sub <- LFMM_PCA_data[(LFMM_PCA_data$V5 == "current_30arcsec_growingDegDays5"),] %>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V6)%>% ## group by k
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))

tail(LFMM_data_sub)


ppop <- ggplot(data = LFMM_data_sub, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_point(color = LFMM_data_sub$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_sub, V3 %in% AcLi_snp_ID), color = "red", size = 0.5)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none")+
  facet_wrap(.~V6, ncol = 1, scales = "free_x")

ggsave(paste0("LFMM_noBWAnoVBC_k123_growingDegDays5_20230608.png"), width = 10, height = 12)

LFMM_data_sub_sig <- LFMM_data_sub %>%
  filter(qvalue < 0.05, V6 == 3)

write_csv(LFMM_data_sub_sig, "LFMM_data_sub_sig_k3_growingDegDays5_20230617.csv")


### chr_02 ####

LFMM_data_sub <- LFMM_PCA_data[(LFMM_PCA_data$V5 == "current_30arcsec_growingDegDays5"),] %>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V6)%>% ## group by k
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))%>%
  filter(V1 == 2)

tail(LFMM_data_sub)


ppop <- ggplot(data = LFMM_data_sub, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_point(color = LFMM_data_sub$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_sub, V3 %in% AcLi_snp_ID), color = "red", size = 0.5)+
  geom_text(data = filter(LFMM_data_sub, qvalue < 0.05), aes(label = V3), size = 2.5, nudge_y = 0.3)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_x_continuous("chr 2", expand = c(0,0))+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none")+
  facet_wrap(.~V6, ncol = 1, scales = "free_x")

ggsave(paste0("LFMM_noBWAnoVBC_k123_chr02_growingDegDays5_20230608.png"), width = 10, height = 12)


### chr_12 ####

LFMM_data_sub <- LFMM_PCA_data[(LFMM_PCA_data$V5 == "current_30arcsec_growingDegDays5"),] %>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V6)%>% ## group by k
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))%>%
  filter(V1 == 12)

tail(LFMM_data_sub)


ppop <- ggplot(data = LFMM_data_sub, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_point(color = LFMM_data_sub$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_sub, V3 %in% AcLi_snp_ID), color = "red", size = 0.5)+
  geom_text(data = filter(LFMM_data_sub, qvalue < 0.05), aes(label = V3), size = 2.5, nudge_y = 0.3)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_x_continuous("chr 12", expand = c(0,0))+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none")+
  facet_wrap(.~V6, ncol = 1, scales = "free_x")

ggsave(paste0("LFMM_noBWAnoVBC_k123_chr12_growingDegDays5_20230608.png"), width = 10, height = 12)

######################################
####### 19 BIOCLIM variables #########
######################################

for (i in 1:3){

LFMM_data_sub <- LFMM_PCA_data[(LFMM_PCA_data$V5 %in% unique(LFMM_PCA_data$V5)[17:35]) & LFMM_PCA_data$V6 == i,] %>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V5)%>% ## group by ENV
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))

tail(LFMM_data_sub)


ppop <- ggplot(data = LFMM_data_sub, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_point(color = LFMM_data_sub$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_sub, V3 %in% AcLi_snp_ID), color = "red", size = 0.5)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none")+
  facet_wrap(.~V5, ncol = 1, scales = "free_x")

ggsave(paste0("LFMM_noBWAnoVBC_k", i, "_19_BIOCLIM_20230607.png"), width = 10, height = 30)

}



###################

LFMM_data <- fread("LFMM_output_noBWA_noVBC_20230522.txt", header = FALSE)

LFMM_data_k1 <- LFMM_data[LFMM_data$V6 == 3,]  ## filter for k value
unique(LFMM_data_k1$V5)

#####################

LFMM_data_k1 <- LFMM_data_k1%>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V5)%>%
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))%>%
  mutate(AcLi_color = case_when(V3 == "Ac" ~ "yellow", V3 == "Li" ~ "blue", TRUE ~ as.character(NA)),
         AcLi_shape = case_when(V3 == "Ac" ~ 15, V3 == "Li" ~ 17, TRUE ~ as.integer(NA)))%>%
  left_join(AcLi_snps, by = c("V3" = "snp"))

tail(LFMM_data_k1)

LFMM_data_k1_AcLi <- filter(LFMM_data_k1, V3 == "Ac" | V3 == "Li")


####

ppop <- ggplot(data = LFMM_data_k1, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_point(color = LFMM_data_k1$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_k1, !is.na(trait)), color = "red", size = 0.5)+
  geom_point(data = LFMM_data_k1_AcLi, color = LFMM_data_k1_AcLi$AcLi_color, size = 2, shape = LFMM_data_k1_AcLi$AcLi_shape)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none")+
  facet_wrap(.~V5, ncol = 1, scales = "free")



p_combined <- plot_grid(p_qtl, ppop, align = "v", ncol = 1, rel_heights = c(1, 45), axis = 'l')

ggsave("LFMM_k3_noBWA_noVBC_20230523.png", width = 10, , height = 50, limitsize = FALSE)



######################


LFMM_data_k1_sub <- LFMM_data_k1[LFMM_data_k1$V5 == "BIO1" | LFMM_data_k1$V5 == "BIO12",]%>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V5)%>%
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))%>%
  mutate(AcLi_color = case_when(V3 == "Ac" ~ "yellow", V3 == "Li" ~ "blue", TRUE ~ as.character(NA)),
         AcLi_shape = case_when(V3 == "Ac" ~ 15, V3 == "Li" ~ 17, TRUE ~ as.integer(NA)))%>%
  left_join(AcLi_snps, by = c("V3" = "snp"))
  
tail(LFMM_data_k1_sub)

LFMM_data_k1_sub_AcLi <- filter(LFMM_data_k1_sub, V3 == "Ac" | V3 == "Li")





 ##################### 

ppop <- ggplot(data = LFMM_data_k1_sub, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_point(color = LFMM_data_k1_sub$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_k1_sub, !is.na(trait)), color = "red", size = 0.5)+
  geom_point(data = LFMM_data_k1_sub_AcLi, color = LFMM_data_k1_sub_AcLi$AcLi_color, size = 2, shape = LFMM_data_k1_sub_AcLi$AcLi_shape)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none")+
  facet_wrap(.~V5, ncol = 1, scales = "free_x")



p_combined <- plot_grid(p_qtl, ppop, align = "v", ncol = 1, rel_heights = c(1, 3), axis = 'l')

ggsave("LFMM_test_20230520.png", width = 10, , height = 4)


###################### QQ plot ####################

p <- ggplot(data = LFMM_data_k1_sub, aes(sample = -log(V4)))+
  stat_qq(aes(x = after_stat(theoretical)/log(10), y = after_stat(sample)/log(10)), distribution = qexp)+
  stat_qq_line(aes(x = after_stat(x)/log(10), y = after_stat(y)/log(10)), distribution = qexp)+
  xlab(expression("Theoretical"~"-"*log[10]*"(p)"))+
  ylab(expression("Observed"~"-"*log[10]*"(p)"))+
  theme_bw()+
  theme(panel.grid = element_blank())
  

p


################################################


### PC1 from k = 1-5

LFMM_data_PC1 <- LFMM_PCA_data[LFMM_PCA_data$V5 == "PC1",] %>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V6)%>% ## group by k
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))%>%
  mutate(AcLi_color = case_when(V3 == "Ac" ~ "yellow", V3 == "Li" ~ "blue", TRUE ~ as.character(NA)),
         AcLi_shape = case_when(V3 == "Ac" ~ 15, V3 == "Li" ~ 17, TRUE ~ as.integer(NA)))%>%
  left_join(AcLi_snps, by = c("V3" = "snp"))

tail(LFMM_data_PC1)

LFMM_data_PC1_AcLi <- filter(LFMM_data_PC1, V3 == "Ac" | V3 == "Li")


ppop <- ggplot(data = LFMM_data_PC1, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_point(color = LFMM_data_PC1$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_PC1, !is.na(trait)), color = "red", size = 0.5)+
  geom_point(data = LFMM_data_PC1_AcLi, color = LFMM_data_PC1_AcLi$AcLi_color, size = 2, shape = LFMM_data_PC1_AcLi$AcLi_shape)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none")+
  facet_wrap(.~V6, ncol = 1, scales = "free_x")

ggsave("LFMM_PC1_k12345_noBWA_noVBC_20230523.png", width = 10, height = 10)


##################

### PC2 from k = 1-5

LFMM_data_PC2 <- LFMM_PCA_data[LFMM_PCA_data$V5 == "PC2",] %>%
  left_join(chr_len_temp, by = c("V1" = "chr"))%>%
  group_by(V6)%>% ## group by k
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues)%>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))%>%
  left_join(AcLi_snps, by = c("V3" = "snp"))

tail(LFMM_data_PC2)

LFMM_data_PC2_AcLi <- filter(LFMM_data_PC2, V3 == "Ac" | V3 == "Li")


ppop <- ggplot(data = LFMM_data_PC2, aes(x = padded_pos, y = -log10(qvalue)))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE)+
  geom_point(color = LFMM_data_PC2$dot_color, size = 0.5)+
  geom_point(data = filter(LFMM_data_PC2, !is.na(trait)), color = "red", size = 0.5)+
  geom_point(data = LFMM_data_PC2_AcLi, color = LFMM_data_PC2_AcLi$AcLi_color, size = 2, shape = LFMM_data_PC2_AcLi$AcLi_shape)+
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5)+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(q value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none")+
  facet_wrap(.~V6, ncol = 1, scales = "free_x")

ggsave("LFMM_PC2_k12345_noBWA_noVBC_20230523.png", width = 10, height = 10)
