library(tidyverse)
library(DESeq2)
library(qvalue)
library(cowplot)
library(biomaRt)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/RNAseq")

get_Arabidopsis_gene <- function(x){
  sapply((sapply(ortho_data$Arabidopsis, function(y) grepl(x, y))), any)
}

get_pallescens_gene <- function(x){
  sapply((sapply(ortho_data$pallescens, function(y) grepl(x, y))), any)
}

get_occidentale_gene <- function(x){
  sapply((sapply(ortho_data$occidentale, function(y) grepl(x, y))), any)
}


ortho_data <- read_tsv("./Orthogroups/Orthogroups.tsv")%>%
  mutate(Arabidopsis = str_split(Arabidopsis, ", "), Hap1_v1.5_protein_chr16 = str_split(Hap1_v1.5_protein_chr16, ", "), Hap1_v1.5_protein_chr8 = str_split(Hap1_v1.5_protein_chr8, ", "))%>%
  filter(!is.na(Arabidopsis))%>%
  dplyr::rename(occidentale = Hap1_v1.5_protein_chr8, pallescens = Hap1_v1.5_protein_chr16)%>%
  mutate(Arabidopsis = lapply(Arabidopsis, function(x)str_split_i(x,"\\|",1)))%>%
  mutate(pallescens = lapply(pallescens, function(x)str_split_i(x,"\\.",1)))%>% ## remove alternative splicing
  mutate(occidentale = lapply(occidentale, function(x)str_split_i(x,"\\.",1)))%>% ## remove alternative splicing
  mutate(pallescens = lapply(pallescens, unique))%>% ## remove alternative splicing
  mutate(occidentale = lapply(occidentale, unique)) ## remove alternative splicing


ortho_data_unnested <- ortho_data %>%
  unnest(Arabidopsis)%>%
  unnest(occidentale)%>%
  unnest(pallescens)%>%
  filter(!(is.na(occidentale) & is.na(pallescens)))%>%
  mutate(Arabidopsis = str_split_i(Arabidopsis, "[.]", 1))

ortho_data_single <- ortho_data %>%
  dplyr::select(-Arabidopsis)%>%
  filter(sapply(occidentale, length) == 1)%>%
  filter(sapply(pallescens, length) == 1)%>%
  mutate_at(2:3, unlist)%>%
  drop_na()

ortho_data_single_long <- ortho_data_single %>%
  pivot_longer(cols = 2:3, names_to = "Subgenome", values_to = "Gene")

ortho_data_multi <- ortho_data %>%
  dplyr::select(-Arabidopsis)%>%
  unnest(occidentale)%>%
  unnest(pallescens)%>%
  drop_na()

ortho_data_multi_long <- ortho_data_multi %>%
  pivot_longer(cols = 2:3, names_to = "Subgenome", values_to = "Gene")


cts <- read_tsv("map_hap1.5_exon_20240228.txt", comment = "#")%>%
  dplyr::select(-c(2:6))

colnames(cts)
str_split_i(colnames(cts), pattern = "/", i = 9) %>% str_remove(pattern = "_hap1.5_2passBasicAligned.sortedByCoord.out.bam") %>% str_remove(pattern = "map_") %>% na.omit()

cts_matrix <- as.matrix(cts[,2:19])
rownames(cts_matrix) <- cts$Geneid
colnames(cts_matrix) <- str_split_i(colnames(cts), pattern = "/", i = 9) %>% str_remove(pattern = "_hap1.5_2passBasicAligned.sortedByCoord.out.bam") %>% str_remove(pattern = "map_") %>% na.omit()
cts_matrix[1:10, 1:18]

sample_info <- data.frame(genotype = rep(c("DMN010", "GFL007", "STL0701"), each = 6), treat = rep(rep(c("control", "drought"), each = 3), times = 3), row.names = colnames(cts_matrix))%>%
  mutate(genotype = factor(genotype, levels = c("DMN010", "GFL007", "STL0701")), treat = factor(treat, levels = c("control", "drought")))



dds <- DESeqDataSetFromMatrix(countData = cts_matrix,
                              colData = sample_info,
                              design = ~ genotype + treat + genotype:treat)


dds_normalized <- DESeq(dds)
dds_normalized
resultsNames(dds_normalized)

####################"treat_drought_vs_control"

res <- results(dds_normalized, name = "treat_drought_vs_control", alpha = 0.05)

res_treat <- as_tibble(res, rownames = "Gene") 

sum(res_treat$padj < 0.05, na.rm = TRUE)/nrow(res_treat)
hist(res_treat$padj)
hist(res_treat$pvalue)

At_ortholog_both_sig <- ortho_data_unnested %>%
  filter((occidentale %in% res_treat_to$Gene) | (pallescens %in% res_treat_tp$Gene))

tair_mart <- useMart(biomart = 'plants_mart',
                     host = 'https://plants.ensembl.org', dataset = 'athaliana_eg_gene')

listAttributes(tair_mart)

annot <- getBM(
  values = ortho_data_unnested$Arabidopsis,
  mart = tair_mart,
  attributes = c('ensembl_gene_id', 'entrezgene_id',
                 'description', 'external_gene_name', "go_id", "name_1006", "definition_1006"),
  filters = 'ensembl_gene_id')

proline_gene <- annot %>%
  filter(str_detect(description, "proline transporter"))

proline_gene <- ortho_data_unnested %>%
  filter(Arabidopsis %in% proline_gene$ensembl_gene_id)%>%
  distinct()

#####

proline_gene
library(ggh4x)

for (i in 1:nrow(proline_gene)){

  Gene_A <- proline_gene[i, 4] %>% as.character()
  Gene_B <- proline_gene[i, 3] %>% as.character()
  strip_colors <- c("Gene_A" = "darkolivegreen3", "Gene_B" = "darkslategray2")
  
  single_gene_data <- as_tibble(counts(dds_normalized[c(Gene_A, Gene_B),]), rownames = "Gene")%>%
    pivot_longer(cols = 2:19, names_to = c("genotype", "treat", "rep"), names_sep = "_", values_to = "value") %>%
    mutate(Gene = factor(Gene, levels = c(Gene_A, Gene_B)), genotype = factor(genotype, labels = c("DMN_010", "GFL_007", "STL_0701")))
  
  p <- ggplot(data = single_gene_data, aes(x = genotype, y = value))+
    geom_point(aes(color = treat, shape = genotype))+
    facet_wrap2(.~Gene, strip = strip_themed(background_x = list(element_rect(fill = strip_colors[["Gene_A"]]), element_rect(fill = strip_colors[["Gene_B"]]))))+
    scale_y_continuous("")+
    scale_x_discrete("")+
    scale_color_manual(values = c("#00BFC4", "#F8766D"))+
    scale_shape_manual(values = c(16, 17, 18))+
    ggtitle(paste0(proline_gene$Arabidopsis[[i]], " (proline transporter)"))+
    theme_bw()+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 8), title = element_text(size = 9))
    
  
  ggsave(paste0("./res_treat/PROT_",i , ".png"), width = 3.5, height = 3.5)
}
################################################################
# beta-CAS


Gene_A <- "drTriRepe4Chr7g205700"
Gene_B <- "drTriRepe4Chr15g223500"
strip_colors <- c("Gene_A" = "darkolivegreen3", "Gene_B" = "darkslategray2")

single_gene_data <- as_tibble(counts(dds_normalized[c(Gene_A, Gene_B),]), rownames = "Gene")%>%
  pivot_longer(cols = 2:19, names_to = c("genotype", "treat", "rep"), names_sep = "_", values_to = "value") %>%
  mutate(Gene = factor(Gene, levels = c(Gene_A, Gene_B)), genotype = factor(genotype, labels = c("DMN_010", "GFL_007", "STL_0701")))

p <- ggplot(data = single_gene_data, aes(x = genotype, y = value))+
  geom_point(aes(color = treat, shape = genotype))+
  facet_wrap2(.~Gene, strip = strip_themed(background_x = list(element_rect(fill = strip_colors[["Gene_A"]]), element_rect(fill = strip_colors[["Gene_B"]]))))+
  scale_y_continuous("")+
  scale_x_discrete("")+
  scale_color_manual(values = c("#00BFC4", "#F8766D"))+
  scale_shape_manual(values = c(16, 17, 18))+
  ggtitle(expression(beta~italic("CAS")~"genes"))+
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 8), title = element_text(size = 9))

p

ggsave(paste0("./res_treat/beta_CAS" , ".png"), width = 3.5, height = 3.5)



############ Cyanogenesis gene expression ######################

res_trt <- as_tibble(counts(dds_normalized[c("drTriRepe4Chr2g112600", "drTriRepe4Chr2g111200", "drTriRepe4Chr2g110400", "drTriRepe4Chr12g338400"),]), rownames = "Gene")%>%
  pivot_longer(2:19, names_to = c("genotype", "treat", "rep"), names_sep = "_") %>%
  mutate(Gene = factor(Gene, levels = c("drTriRepe4Chr2g112600", "drTriRepe4Chr2g111200", "drTriRepe4Chr2g110400", "drTriRepe4Chr12g338400"), labels = c("CYP79D15", "CYP736A187", "UGT85K17", "Li")))

p <- ggplot(data = res_trt, aes(x = genotype, y = value))+
  geom_point(aes(color = treat, group = treat), position = position_dodge(width = 0.8))+
  scale_x_discrete("")+
  scale_y_log10("",)+
  scale_color_manual(values = c("#00BFC4", "#F8766D"))+
  facet_wrap(Gene ~ ., nrow = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), strip.text = element_text(face = "italic"), legend.position = c(0.92, 0.9), legend.title = element_blank())
p


ggsave("Ac_Li_gene_expression_geno_trt.png", width = 8, height = 5, dpi = 600)



################################################################
############ Cyanogenesis gene expression II ######################

res_trt <- as_tibble(counts(dds_normalized[c("drTriRepe4Chr2g112600", "drTriRepe4Chr2g111200", "drTriRepe4Chr2g110400"),]), rownames = "Gene")%>%
  pivot_longer(2:19, names_to = c("genotype", "treat", "rep"), names_sep = "_") %>%
  mutate(Gene = factor(Gene, levels = c("drTriRepe4Chr2g112600", "drTriRepe4Chr2g111200", "drTriRepe4Chr2g110400"), labels = c("CYP79D15", "CYP736A187", "UGT85K17")))

p <- ggplot(data = res_trt, aes(x = genotype, y = value))+
  geom_point(aes(color = treat, group = treat), position = position_dodge(width = 0.8))+
  scale_x_discrete("")+
  scale_y_log10("",)+
  scale_color_manual(values = c("#00BFC4", "#F8766D"))+
  facet_wrap(Gene ~ ., nrow = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), strip.text = element_text(face = "italic"), legend.position = c(0.09, 0.9), legend.title = element_blank())
p


ggsave("Ac_Li_gene_expression_geno_trt2.png", width = 8, height = 5, dpi = 600)



################################################################

intersect(res_treat$Gene, res_interaction_GFL$Gene)

ortho_data_drTriRepe4Chr6g039500 <- ortho_data %>%
  filter(get_occidentale_gene("drTriRepe4Chr6g039500"))


res_treat_sig <- res_treat %>%
  filter(padj < 0.05)

single_gene_data <- as_tibble(counts(dds_normalized[c(res_interaction_GFL_treat[[1]], "drTriRepe4Chr3g551400"),]), rownames = "Gene")%>%
  pivot_longer(2:19, names_to = c("genotype", "treat", "rep"), names_sep = "_")

p <- ggplot(data = single_gene_data, aes(x = genotype, y = value))+
  geom_point(aes(color = treat))+
  facet_wrap(.~Gene)
p

plotCounts(dds_normalized, gene = "drTriRepe4Chr3g461900", intgroup = "genotype:treat")
plotCounts(dds_normalized, gene = res_interaction_GFL_treat[[1]], intgroup = "genotype:treat")

####################

for (i in resultsNames(dds_normalized)[-1]){
  
  # Assuming 'dds_normalized' is your DESeqDataSet object and the analysis is already done.
  # change this for each comparison
  res <- results(dds_normalized, name = i)
  
  # Convert results to data frame
  res_df <- as.data.frame(res)
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id")%>%
    filter(!is.na(pvalue) & !is.na(log2FoldChange)) %>%
    mutate(significant = ifelse(pvalue < 0.01, "yes", "no"))%>%
    mutate(subgenome = case_when(str_detect(gene_id, "drTriRepe4Chr[1-8]g") ~ "To",
                                 str_detect(gene_id, "drTriRepe4Chr9g|drTriRepe4Chr10g|drTriRepe4Chr1[1-6]g") ~ "Tp",
                                 TRUE ~ "scaffolds"))%>%
    mutate(color_s = paste(significant, subgenome, sep = "_"))
  
  # Volcano plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = color_s)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("yes_To" = "red", "yes_Tp" = "blue"),
                       breaks = c("yes_To", "yes_Tp"),
                       labels = c("T. occidentale subgenome", "T. pallescens subgenome"))+
    labs(title = i,
         x = "Log2 Fold Change",
         y = "-Log10 p-value") +
  scale_x_continuous(limits = c(-15, 15), oob = scales::oob_squish)+
    theme_bw()+
    theme(legend.title = element_blank(), legend.position = c(0.8, 0.9))
  
  ggsave(paste0("volcano_", i, ".png"), width = 5, height = 5)
}




# Generate MA plot
plotMA(res, main = "MA Plot")

####################
rld <- rlogTransformation(dds_normalized)
rld$genotype
rld$treat
t(assay(rld))[1:10, 1:10]
plotPCA(rld, intgroup=c("genotype", "treat"))


# res <- results(dds_normalized, name = "treat_drought_vs_control")

# sig_genes <- res[which(res$padj < 0.05), ]
rld_sig <- rld[rownames(rld) %in% rownames(sig_genes), ]
# plotPCA(rld_sig, intgroup=c("genotype", "treat"))



pca_data <- prcomp(t(assay(rld))) # Ensure scaling is applied if needed

# Extract PCA scores
pca_scores <- as.data.frame(pca_data$x)

# Add genotype and treatment information
pca_scores$genotype <- colData(rld)$genotype
pca_scores$treatment <- colData(rld)$treat

# Optionally, name the PCA dimensions for easier reference
names(pca_scores)[1:2] <- c("PC1", "PC2")

# Calculate the percentage of variance explained by each PC
variance_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2) * 100
variance_explained_pc1 <- round(variance_explained[1], 2) # for PC1
variance_explained_pc2 <- round(variance_explained[2], 2) # for PC2

pca_scores <- pca_scores %>%
  mutate(genotype = factor(genotype, labels = c("DMN_010", "GFL_007", "STL_0701")))


p <- ggplot(pca_scores, aes(x = PC1, y = PC2, shape = genotype, color = treatment)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(16, 17, 18)) + # Adjust as needed for the number of genotypes
  scale_color_manual(values = c("#00BFC4", "#F8766D"))+ # Adjust colors as needed for treatments
  theme_minimal() +
  labs(x = paste("PC1 (", variance_explained_pc1, "%)", sep=""),
       y = paste("PC2 (", variance_explained_pc2, "%)", sep="")) +
  theme_bw()+
  theme(legend.position = c(0.85,0.25))
p

ggsave("PCA_all_expression_20250511.png", width = 5, height = 5)

###########################

res <- results(dds_normalized, name = "treat_drought_vs_control")


sig_genes <- res[which(res$padj < 0.05), ]
rld_sig <- rld[rownames(rld) %in% rownames(sig_genes), ]

############################

total_exp <- t(assay(rld))
# Assuming 'my_matrix' is your matrix
column_names <- colnames(total_exp)
length(unique(column_names))

# Finding column indices with names matching "drTriRepe4Chr[1-8]"
column_to <- grep("^drTriRepe4Chr[1-8]g", column_names)
length(column_to)

column_tp <- grep("^drTriRepe4Chr(9|10|11|12|13|14|15|16)g", column_names)
length(column_tp)

# Selecting the matched columns
sub_to_exp <- total_exp[, column_to]
sub_tp_exp <- total_exp[, column_tp]

unique(sub("g....00", "",colnames(sub_to_exp)))
unique(sub("g....00", "",colnames(sub_tp_exp)))
# Now 'selected_columns' contains the columns of 'my_matrix' that match the pattern

######################
pca_data <- prcomp(sub_to_exp) # Ensure scaling is applied if needed

# Extract PCA scores
pca_scores <- as.data.frame(pca_data$x)

# Add genotype and treatment information
pca_scores$genotype <- colData(rld)$genotype
pca_scores$treatment <- colData(rld)$treat

# Optionally, name the PCA dimensions for easier reference
names(pca_scores)[1:2] <- c("PC1", "PC2")

# Calculate the percentage of variance explained by each PC
variance_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2) * 100
variance_explained_pc1 <- round(variance_explained[1], 2) # for PC1
variance_explained_pc2 <- round(variance_explained[2], 2) # for PC2

pca_scores <- pca_scores %>%
  mutate(genotype = factor(genotype, labels = c("DMN_010", "GFL_007", "STL_0701")))

p <- ggplot(pca_scores, aes(x = PC1, y = PC2, shape = genotype, color = treatment)) +
  geom_point(size = 3) +
  scale_shape_manual("", values = c(16, 17, 18)) + # Adjust as needed for the number of genotypes
  scale_color_manual("", values = c("#00BFC4", "#F8766D"), labels = c("Control", "Drought"))+ # Adjust colors as needed for treatments
  theme_minimal() +
  labs(title = expression(italic("T. occidentale")~"subgenome"),
    x = paste("PC1 (", variance_explained_pc1, "%)", sep=""),
       y = paste("PC2 (", variance_explained_pc2, "%)", sep="")) +
  theme_bw()+
  # theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(0.85,0.26), plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.spacing = unit(0, "in"), legend.background = element_rect(fill = NA))
p

ggsave("PCA_To_subgenome_expression_20250511.png", width = 3.5, height = 3.5, dpi = 600)

#######################




########################
######################
pca_data <- prcomp(sub_tp_exp) # Ensure scaling is applied if needed

# Extract PCA scores
pca_scores <- as.data.frame(pca_data$x)

# Add genotype and treatment information
pca_scores$genotype <- colData(rld)$genotype
pca_scores$treatment <- colData(rld)$treat

# Optionally, name the PCA dimensions for easier reference
names(pca_scores)[1:2] <- c("PC1", "PC2")

# Calculate the percentage of variance explained by each PC
variance_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2) * 100
variance_explained_pc1 <- round(variance_explained[1], 2) # for PC1
variance_explained_pc2 <- round(variance_explained[2], 2) # for PC2

pca_scores <- pca_scores %>%
  mutate(genotype = factor(genotype, labels = c("DMN_010", "GFL_007", "STL_0701")))

p <- ggplot(pca_scores, aes(x = PC1, y = PC2, shape = genotype, color = treatment)) +
  geom_point(size = 3) +
  scale_shape_manual("", values = c(16, 17, 18)) + # Adjust as needed for the number of genotypes
  scale_color_manual("", values = c("#00BFC4", "#F8766D"), labels = c("Control", "Drought"))+ # Adjust colors as needed for treatments
  theme_minimal() +
  labs(title = expression(italic("T. pallescens")~"subgenome"),
       x = paste("PC1 (", variance_explained_pc1, "%)", sep=""),
       y = paste("PC2 (", variance_explained_pc2, "%)", sep="")) +
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  # theme(legend.position = c(0.85,0.26), plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.spacing = unit(0, "in"), legend.background = element_rect(fill = NA))
p

ggsave("PCA_Tp_subgenome_expression_20250511.png", width = 3.5, height = 3.5)



sub_to_exp <- 


plotPCA(rld, intgroup=c("treat"))

plotCounts(dds_normalized, gene = "drTriRepe4Chr1g000200", intgroup = "genotype")
plotCounts(dds_normalized, gene = "gene-QL285_007457", intgroup = "dev_stage")
plotCounts(dds_normalized, gene = "gene-QL285_069672", intgroup = "dev_stage")
plotCounts(dds_normalized, gene = "gene-QL285_069674", intgroup = "dev_stage")
gene-QL285_069672
