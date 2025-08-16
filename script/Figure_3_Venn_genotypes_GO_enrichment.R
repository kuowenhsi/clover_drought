library(tidyverse)
library(DESeq2)
library(qvalue)
library(cowplot)
library(ggVennDiagram)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/RNAseq")


cts <- read_tsv("map_hap1.5_exon_20240228.txt", comment = "#")%>%
  select(-c(2:6))

colnames(cts)
str_split_i(colnames(cts), pattern = "/", i = 9) %>% str_remove(pattern = "_hap1.5_2passBasicAligned.sortedByCoord.out.bam") %>% str_remove(pattern = "map_") %>% na.omit()

cts_matrix <- as.matrix(cts[,2:19])
rownames(cts_matrix) <- cts$Geneid
colnames(cts_matrix) <- str_split_i(colnames(cts), pattern = "/", i = 9) %>% str_remove(pattern = "_hap1.5_2passBasicAligned.sortedByCoord.out.bam") %>% str_remove(pattern = "map_") %>% na.omit()
cts_matrix[1:10, 1:18]

sample_info <- data.frame(genotype = rep(c("DMN010", "GFL007", "STL0701"), each = 6), treat = rep(rep(c("control", "drought"), each = 3), times = 3), row.names = colnames(cts_matrix))%>%
  mutate(genotype = factor(genotype, levels = c("DMN_010", "GFL_007", "STL_0701")), treat = factor(treat, levels = c("control", "drought")))

#############################

dds_DMN <- DESeqDataSetFromMatrix(countData = cts_matrix[,1:6],
                              colData = sample_info[1:6,],
                              design = ~ treat)

dds_DMN_normalized <- DESeq(dds_DMN)
dds_DMN_normalized
resultsNames(dds_DMN_normalized)

res_DMN <- results(dds_DMN_normalized, name = "treat_drought_vs_control", alpha = 0.05)

res_DMN_treat <- as_tibble(res_DMN, rownames = "Gene")%>%
  mutate(subgenome = case_when(str_detect(Gene, "drTriRepe4Chr[1-8]g") ~ "occidentale", str_detect(Gene, "drTriRepe4Chr9g") | str_detect(Gene, "drTriRepe4Chr1[0-6]g") ~ "pallescens"))%>%
  mutate(regulation = case_when(log2FoldChange < 0 ~ "Negative", log2FoldChange >= 0 ~ "Positive"))%>%
  mutate(genotype = "DMN")
res_DMN_treat_sig <- res_DMN_treat %>%
  filter(padj < 0.05)

#############################################


dds_GFL <- DESeqDataSetFromMatrix(countData = cts_matrix[,7:12],
                                  colData = sample_info[7:12,],
                                  design = ~ treat)

dds_GFL_normalized <- DESeq(dds_GFL)
dds_GFL_normalized
resultsNames(dds_GFL_normalized)

res_GFL <- results(dds_GFL_normalized, name = "treat_drought_vs_control", alpha = 0.05)

res_GFL_treat <- as_tibble(res_GFL, rownames = "Gene") %>%
  mutate(subgenome = case_when(str_detect(Gene, "drTriRepe4Chr[1-8]g") ~ "occidentale", str_detect(Gene, "drTriRepe4Chr9g") | str_detect(Gene, "drTriRepe4Chr1[0-6]g") ~ "pallescens"))%>%
  mutate(regulation = case_when(log2FoldChange < 0 ~ "Negative", log2FoldChange >= 0 ~ "Positive"))%>%
  mutate(genotype = "GFL")
res_GFL_treat_sig <- res_GFL_treat %>%
  filter(padj < 0.05)

############################################


dds_STL <- DESeqDataSetFromMatrix(countData = cts_matrix[,13:18],
                                  colData = sample_info[13:18,],
                                  design = ~ treat)

dds_STL_normalized <- DESeq(dds_STL)
dds_STL_normalized
resultsNames(dds_STL_normalized)

res_STL <- results(dds_STL_normalized, name = "treat_drought_vs_control", alpha = 0.05)

res_STL_treat <- as_tibble(res_STL, rownames = "Gene") %>%
  mutate(subgenome = case_when(str_detect(Gene, "drTriRepe4Chr[1-8]g") ~ "occidentale", str_detect(Gene, "drTriRepe4Chr9g") | str_detect(Gene, "drTriRepe4Chr1[0-6]g") ~ "pallescens"))%>%
  mutate(regulation = case_when(log2FoldChange < 0 ~ "Negative", log2FoldChange >= 0 ~ "Positive"))%>%
  mutate(genotype = "STL")
res_STL_treat_sig <- res_STL_treat %>%
  filter(padj < 0.05)
res_STL_treat_sig
#########################################

res_subgenome <- bind_rows(res_DMN_treat_sig, res_GFL_treat_sig, res_STL_treat_sig)%>%
  group_by(genotype, subgenome, regulation)%>%
  summarize(count = n())%>%
  drop_na()%>%
  ungroup()%>%
  mutate(count = case_when(regulation == "Negative" ~ -count, TRUE ~ count), genotype = factor(genotype, labels = c("DMN_010", "GFL_007", "STL_0701")))

p <- ggplot(data = res_subgenome, aes(x = genotype, y = count, fill = subgenome, group = subgenome))+
  geom_col(position = position_dodge())+
  geom_hline(yintercept = 0, linewidth = 0.5, color = "gray")+
  annotate(geom = "text", x = 0.5, y = 2000, hjust = 0, label = "Fisher’s Exact Test\nNot significant in all comparisons")+
  scale_fill_manual(values = c("darkolivegreen3", "darkslategray2"))+
  scale_y_continuous(breaks = c(-2000, -1000, 0, 1000, 2000), labels = function(x){abs(x)})+
  labs(y = "Downregulated genes      Upregulated genes", x = "")+
  theme_bw()+
  theme(legend.position = c(0.2, 0.2))
p

ggsave("subgenome_dominance_treat_20250509.png", width = 3.5, height = 3.5, dpi = 600)

# Define the total number of genes in each subgenome
total_occidentale <- 48107
total_pallescens <- 48124

# Perform Fisher’s Exact Test within each genotype and regulation category
results <- res_subgenome %>%
  mutate(count = abs(count))%>%
  group_by(genotype, regulation) %>%
  summarise(
    fisher_p_value = {
      # Create contingency table
      sub_occ <- count[subgenome == "occidentale"]
      sub_pall <- count[subgenome == "pallescens"]
      
      # Compute non-expressed genes
      non_expressed_occ <- total_occidentale - sub_occ
      non_expressed_pall <- total_pallescens - sub_pall
      
      # Construct contingency table
      contingency_table <- matrix(c(sub_occ, non_expressed_occ,
                                    sub_pall, non_expressed_pall),
                                  nrow = 2, byrow = TRUE,
                                  dimnames = list(
                                    Subgenome = c("occidentale", "pallescens"),
                                    Expression = c("Expressed", "Not Expressed")
                                  ))
      
      # Perform Fisher's exact test
      test_result <- fisher.test(contingency_table)
      
      # Return p-value
      test_result$p.value
    },
    .groups = "drop"
  )

# Print the results
print(results)
#########################################

Venn_data <- list(DMN_010 = res_DMN_treat_sig$Gene, GFL_007 = res_GFL_treat_sig$Gene, STL_0701 = res_STL_treat_sig$Gene)

ggVennDiagram(Venn_data, label_alpha = 0, label_size =3, set_size = 3, edge_size = 0.5)+
  scale_fill_gradient(low="grey90",high = "red")+
  theme(plot.background = element_rect(fill = "white", color = NA))+
  scale_x_continuous(expand = expansion(mult = .1))

ggsave("Venn_DE_genotypes_20250509.png", width = 3.5, height = 3.5, dpi = 600)

ggVennDiagram(Venn_data, force_upset = TRUE)

shared_gene <- Reduce(intersect, Venn_data)
shared_gene

only_DMN <- setdiff(Venn_data$DMN, union(Venn_data$GFL, Venn_data$STL))
only_GFL <- setdiff(Venn_data$GFL, union(Venn_data$DMN, Venn_data$STL))
only_STL <- setdiff(Venn_data$STL, union(Venn_data$DMN, Venn_data$GFL))

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
  filter(!is.na(Arabidopsis))%>%
  mutate(Arabidopsis = str_split_i(Arabidopsis, "[.]", 1))

ortho_data_unnested_shared_gene <- ortho_data_unnested %>%
  pivot_longer(col = 3:4, names_to = "subgenome", values_to = "Gene")%>%
  filter(Gene %in% shared_gene)%>%
  distinct()

ortho_data_unnested_only_DMN <- ortho_data_unnested %>%
  pivot_longer(col = 3:4, names_to = "subgenome", values_to = "Gene")%>%
  filter(Gene %in% only_DMN)%>%
  distinct()

ortho_data_unnested_only_GFL <- ortho_data_unnested %>%
  pivot_longer(col = 3:4, names_to = "subgenome", values_to = "Gene")%>%
  filter(Gene %in% only_GFL)%>%
  distinct()

ortho_data_unnested_only_STL <- ortho_data_unnested %>%
  pivot_longer(col = 3:4, names_to = "subgenome", values_to = "Gene")%>%
  filter(Gene %in% only_STL)%>%
  distinct()


unique(ortho_data_unnested_shared_gene$Gene)
ortho_data_unnested_shared_gene

shared_GENES <- ortho_data_unnested_shared_gene$Arabidopsis %>% unique()
only_DMN_GENES <- ortho_data_unnested_only_DMN$Arabidopsis %>% unique()
only_GFL_GENES <- ortho_data_unnested_only_GFL$Arabidopsis %>% unique()
only_STL_GENES <- ortho_data_unnested_only_STL$Arabidopsis %>% unique()
##########################

# library(biomaRt)
# tair_mart <- useMart(biomart = 'plants_mart',
#                      host = 'https://plants.ensembl.org', dataset = 'athaliana_eg_gene')
# 
# annot <- getBM(
#   values = GENES,
#   mart = tair_mart,
#   attributes = c('ensembl_gene_id', 'entrezgene_id',
#                  'description', 'external_gene_name'),
#   filters = 'ensembl_gene_id')
# 
# annot
# 
# ortho_data_unnested_shared_gene_annot <- ortho_data_unnested_shared_gene %>%
#   mutate(Arabidopsis = str_split_i(Arabidopsis, "[.]", 1))%>%
#   left_join(annot, by = c("Arabidopsis" = "ensembl_gene_id"))
# 
# library(writexl)
# 
# write_xlsx(ortho_data_unnested_shared_gene_annot, "DE_genes_intersect_DMN_GFL_STL.xlsx")


library(topGO)

background_genes <- ortho_data_unnested$Arabidopsis %>% unique()

# Define ALL_genes and Enriched_genes
ALL_genes <- background_genes
Enriched_genes <- shared_GENES

# Connect to the database
mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "athaliana_eg_gene",
                         host = 'plants.ensembl.org')

# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id",
                                        "go_id"), mart = mart)
# Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]

# Convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))

all.genes <- sort(unique(ALL_genes))

Enriched_genes <- factor(as.numeric(all.genes %in% Enriched_genes))
names(Enriched_genes) = all.genes

# Create topGO objects for each ontology
go.objBP = new("topGOdata", ontology='BP',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

go.objCC = new("topGOdata", ontology='CC',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

go.objMF = new("topGOdata", ontology='MF',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

# Run enrichment test
resultsBP <- runTest(go.objBP, algorithm = "elim", statistic = "fisher")
resultsCC <- runTest(go.objCC, algorithm = "elim", statistic = "fisher")
resultsMF <- runTest(go.objMF, algorithm = "elim", statistic = "fisher")

# Get top significant GO terms (adjust `numTerms` as needed)
numTerms <- 10  # Number of terms to display

# Generate tables
results.tabBP <- GenTable(object = go.objBP, elimFisher = resultsBP, topNodes = numTerms)
results.tabCC <- GenTable(object = go.objCC, elimFisher = resultsCC, topNodes = numTerms)
results.tabMF <- GenTable(object = go.objMF, elimFisher = resultsMF, topNodes = numTerms)

results.tabBP$Ontology <- "Biological process"
results.tabCC$Ontology <- "Cellular component"
results.tabMF$Ontology <- "Molecular function"
GO_results <- rbind(results.tabBP, results.tabCC, results.tabMF)%>%
  mutate(elimFisher = -log10(as.numeric(elimFisher)))

p <- ggplot(GO_results, aes(y = reorder(Term, elimFisher), x = Significant)) +
  geom_bar(stat = "identity", aes(fill = elimFisher)) +
  labs(title = "DMN GFL STL shared differentially expressed genes", x = "Gene count", y = "") +
  scale_fill_fermenter("-log(p)", type = "seq", palette = "OrRd", direction = 1)+
  facet_wrap(.~Ontology, scales = "free", ncol = 1)+
  theme_bw()
p

ggsave("GO_shared.png", width = 8, height = 10, dpi = 600)

###############

# Define ALL_genes and Enriched_genes
ALL_genes <- background_genes
Enriched_genes <- only_DMN_GENES

all.genes <- sort(unique(ALL_genes))

Enriched_genes <- factor(as.numeric(all.genes %in% Enriched_genes))
names(Enriched_genes) = all.genes

# Create topGO objects for each ontology
go.objBP = new("topGOdata", ontology='BP',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

go.objCC = new("topGOdata", ontology='CC',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

go.objMF = new("topGOdata", ontology='MF',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

# Run enrichment test
resultsBP <- runTest(go.objBP, algorithm = "elim", statistic = "fisher")
resultsCC <- runTest(go.objCC, algorithm = "elim", statistic = "fisher")
resultsMF <- runTest(go.objMF, algorithm = "elim", statistic = "fisher")

# Get top significant GO terms (adjust `numTerms` as needed)
numTerms <- 10  # Number of terms to display

# Generate tables
results.tabBP <- GenTable(object = go.objBP, elimFisher = resultsBP, topNodes = numTerms)
results.tabCC <- GenTable(object = go.objCC, elimFisher = resultsCC, topNodes = numTerms)
results.tabMF <- GenTable(object = go.objMF, elimFisher = resultsMF, topNodes = numTerms)

results.tabBP$Ontology <- "Biological process"
results.tabCC$Ontology <- "Cellular component"
results.tabMF$Ontology <- "Molecular function"
GO_results <- rbind(results.tabBP, results.tabCC, results.tabMF)%>%
  mutate(elimFisher = -log10(as.numeric(elimFisher)))

p <- ggplot(GO_results, aes(y = reorder(Term, elimFisher), x = Significant)) +
  geom_bar(stat = "identity", aes(fill = elimFisher)) +
  labs(title = "DMN only differentially expressed genes", x = "Gene count", y = "") +
  scale_fill_fermenter("-log(p)",type = "seq", palette = "OrRd", direction = 1)+
  facet_wrap(.~Ontology, scales = "free", ncol = 1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1))

p

ggsave("GO_DMN_only.png", width = 8, height = 10, dpi = 600)

##################################

# Define ALL_genes and Enriched_genes
ALL_genes <- background_genes
Enriched_genes <- only_GFL_GENES

all.genes <- sort(unique(ALL_genes))

Enriched_genes <- factor(as.numeric(all.genes %in% Enriched_genes))
names(Enriched_genes) = all.genes

# Create topGO objects for each ontology
go.objBP = new("topGOdata", ontology='BP',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

go.objCC = new("topGOdata", ontology='CC',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

go.objMF = new("topGOdata", ontology='MF',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

# Run enrichment test
resultsBP <- runTest(go.objBP, algorithm = "elim", statistic = "fisher")
resultsCC <- runTest(go.objCC, algorithm = "elim", statistic = "fisher")
resultsMF <- runTest(go.objMF, algorithm = "elim", statistic = "fisher")

# Get top significant GO terms (adjust `numTerms` as needed)
numTerms <- 10  # Number of terms to display

# Generate tables
results.tabBP <- GenTable(object = go.objBP, elimFisher = resultsBP, topNodes = numTerms)
results.tabCC <- GenTable(object = go.objCC, elimFisher = resultsCC, topNodes = numTerms)
results.tabMF <- GenTable(object = go.objMF, elimFisher = resultsMF, topNodes = numTerms)

results.tabBP$Ontology <- "Biological process"
results.tabCC$Ontology <- "Cellular component"
results.tabMF$Ontology <- "Molecular function"
GO_results <- rbind(results.tabBP, results.tabCC, results.tabMF)%>%
  mutate(elimFisher = -log10(as.numeric(elimFisher)))

p <- ggplot(GO_results, aes(y = reorder(Term, elimFisher), x = Significant)) +
  geom_bar(stat = "identity", aes(fill = elimFisher)) +
  labs(title = "GFL only differentially expressed genes", x = "Gene count", y = "") +
  scale_fill_fermenter("-log(p)",type = "seq", palette = "OrRd", direction = 1)+
  facet_wrap(.~Ontology, scales = "free", ncol = 1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1))

p

ggsave("GO_GFL_only.png", width = 8, height = 10, dpi = 600)

##################################

# Define ALL_genes and Enriched_genes
ALL_genes <- background_genes
Enriched_genes <- only_STL_GENES

all.genes <- sort(unique(ALL_genes))

Enriched_genes <- factor(as.numeric(all.genes %in% Enriched_genes))
names(Enriched_genes) = all.genes

# Create topGO objects for each ontology
go.objBP = new("topGOdata", ontology='BP',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

go.objCC = new("topGOdata", ontology='CC',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

go.objMF = new("topGOdata", ontology='MF',
               allGenes = Enriched_genes,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)

# Run enrichment test
resultsBP <- runTest(go.objBP, algorithm = "elim", statistic = "fisher")
resultsCC <- runTest(go.objCC, algorithm = "elim", statistic = "fisher")
resultsMF <- runTest(go.objMF, algorithm = "elim", statistic = "fisher")

# Get top significant GO terms (adjust `numTerms` as needed)
numTerms <- 10  # Number of terms to display

# Generate tables
results.tabBP <- GenTable(object = go.objBP, elimFisher = resultsBP, topNodes = numTerms)
results.tabCC <- GenTable(object = go.objCC, elimFisher = resultsCC, topNodes = numTerms)
results.tabMF <- GenTable(object = go.objMF, elimFisher = resultsMF, topNodes = numTerms)

results.tabBP$Ontology <- "Biological process"
results.tabCC$Ontology <- "Cellular component"
results.tabMF$Ontology <- "Molecular function"
GO_results <- rbind(results.tabBP, results.tabCC, results.tabMF)%>%
  mutate(elimFisher = -log10(as.numeric(elimFisher)))

p <- ggplot(GO_results, aes(y = reorder(Term, elimFisher), x = Significant)) +
  geom_bar(stat = "identity", aes(fill = elimFisher)) +
  labs(title = "STL only differentially expressed genes", x = "Gene count", y = "") +
  scale_fill_fermenter("-log(p)",type = "seq", palette = "OrRd", direction = 1)+
  facet_wrap(.~Ontology, scales = "free", ncol = 1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 1))

p

ggsave("GO_STL_only.png", width = 8, height = 10, dpi = 600)
