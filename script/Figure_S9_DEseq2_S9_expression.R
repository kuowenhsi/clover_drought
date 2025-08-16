library(tidyverse)
library(DESeq2)
setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/RNAseq")

cts <- read_tsv("S10_count_hap1_1.5_NCBI_20231024.txt", comment = "#")%>%
  select(-c(2:7))

str_split_i(colnames(cts), pattern = "/", i = 8) %>% str_split_i(pattern = "_", i = 2) %>% na.omit()

cts_matrix <- as.matrix(cts[,2:13])
rownames(cts_matrix) <- cts$Geneid
colnames(cts_matrix) <- str_split_i(colnames(cts), pattern = "/", i = 8) %>% str_split_i(pattern = "_", i = 2) %>% na.omit()
cts_matrix[1:10, 1:12]

sample_info <- data.frame(dev_stage = c("stolon", "root", "stolon", "root", "youngLeaf", "matureLeaf", "stolon", "root", "matureLeaf", "youngLeaf", "matureLeaf", "youngLeaf"), row.names = colnames(cts_matrix))

colnames(cts_matrix)

dds <- DESeqDataSetFromMatrix(countData = cts_matrix,
                              colData = sample_info,
                              design = ~ dev_stage)

dds_normalized <- DESeq(dds)
dds_normalized
res <- results(dds_normalized)
res


res_dev <- as_tibble(counts(dds_normalized[c("gene-QL285_007426", "gene-QL285_007412", "gene-QL285_007404", "gene-QL285_069672"),]), rownames = "Gene")%>%
  pivot_longer(2:13, names_to = "accession")%>%
  mutate(Tissue = rep(c("stolon", "root", "stolon", "root", "youngLeaf", "matureLeaf", "stolon", "root", "matureLeaf", "youngLeaf", "matureLeaf", "youngLeaf"), 4)) %>%
  mutate(Gene = factor(Gene, levels = c("gene-QL285_007426", "gene-QL285_007412", "gene-QL285_007404", "gene-QL285_069672"), labels = c("CYP79D15", "CYP736A187", "UGT85K17", "Li")))%>%
  filter(Gene != "Li")

p <- ggplot(data = res_dev, aes(x = Tissue, y = value))+
  geom_point()+
  scale_x_discrete("")+
  scale_y_log10("")+
  facet_wrap(Gene ~ ., nrow = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1), strip.text = element_text(face = "italic"))
p


ggsave("Ac_Li_gene_expression_dif_dev2.png", width = 8, height = 5, dpi = 600)

## CYP79D15 drTriRepe4Chr2g112600.1 
plotCounts(dds_normalized, gene = "gene-QL285_007426", intgroup = "dev_stage")

## CYP79D15 drTriRepe4Chr2g115700.1
plotCounts(dds_normalized, gene = "gene-QL285_007457", intgroup = "dev_stage")

## CYP736A187 drTriRepe4Chr2g111200.1 

plotCounts(dds_normalized, gene = "gene-QL285_007412", intgroup = "dev_stage")

## CYP736A187 drTriRepe4Chr2g114400.1

plotCounts(dds_normalized, gene = "gene-QL285_007444", intgroup = "dev_stage")

## UGT85K17 drTriRepe4Chr2g110400.1

plotCounts(dds_normalized, gene = "gene-QL285_007404", intgroup = "dev_stage")

## UGT85K17 drTriRepe4Chr2g113800.1

plotCounts(dds_normalized, gene = "gene-QL285_007438", intgroup = "dev_stage")


## Li drTriRepe4Chr12g338400.1
plotCounts(dds_normalized, gene = "gene-QL285_069672", intgroup = "dev_stage")


## Li drTriRepe4Chr12g338600.1
plotCounts(dds_normalized, gene = "gene-QL285_069674", intgroup = "dev_stage")

