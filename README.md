American Journal of Botany
Manuscript Number: AJB-D-25-00324R2

Full Title: Phenotypic and genetic bases of variable drought stress response in a widely adapted
allotetraploid species

Article Type: Special Issue Article
Corresponding Author: Kenneth M. Olsen
Washington University In St Louis: Washington University in St Louis
St. Louis, MO UNITED STATES OF AMERICA
Corresponding Author's Institution: Washington University In St Louis: Washington University in St Louis
Corresponding Author E-Mail: kolsen@wustl.edu

This repository contains figure-generation scripts, processed input data, and selected result tables used for phenotypic, QTL, and RNA-seq analyses in the drought study.

Repository structure
clover_drought/
├── data/                  # processed input tables and supporting files
│   ├── RNAseq/            # RNA-seq count matrix and orthogroup files
│   ├── LICOR/             # cleaned and raw LICOR files
│   └── ...
├── figures/               # generated figure outputs
│   ├── Figure3/
│   ├── Figure4/
│   ├── Figure5/
│   ├── Figure6/
│   └── ...
├── results/               # exported statistical summaries
│   └── Figure5/
├── script/                # figure-generation scripts
└── README.md
Directory overview
data/

Processed input files used by the analysis scripts.

R scripts used to generate manuscript figures.

Current scripts:

Figure_1_make_parents_plot.R
Figure_2_making_plot_20250511.R
Figure_3_DEseq2_20240226_expression.R
Figure_3_Venn_genotypes_GO_enrichment.R
Figure_4_generate_hoemologous_plots_20250505.R
Figure_5_plot_cor_matrix.R
Figure_6_make_final_qtl_plots.R
