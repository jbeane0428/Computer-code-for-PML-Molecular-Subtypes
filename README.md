# Computer-code-for-PML-Molecular-Subtypes
R scripts and input data to produce results in the manuscript:

Molecular subtyping reveals immune alterations associated with progression of bronchial premalignant lesions.

Authors:  Jennifer Beane*, Sarah A. Mazzilli, Joshua D. Campbell, Grant Duclos, Kostyantyn Krysan, Christopher Moy, Catalina Perdomo, Michael Schaffer, Gang Liu, Sherry Zhang, Hangqio Liu, Jessica Vick, Samjot S. Dhillon, Suso J. Platero, Steven M. Dubinett, Christopher Stevenson, Mary E. Reid, Marc E. Lenburg, Avrum E. Spira

Nature Communications 2019

Input data:
Endobronchial biopsy and brush data:  GSE109743_exp_count.txt
Endobronchial biopsy and brush phenotypic data:  GSE109743_SAMPLES.txt
NTCU mouse model data:  GSE111091_exp_count.txt
NTCU mouse model phenotypic data:  GSE111091_SAMPLES.txt
TCGA SCC data:  TCGA_SCC_cpm.txt (Campbell et al. Nature Genetics 2016, Log2 plus 1 TPM)

R scripts:
processData_and_runWGCNA.R (processes data and compute WGCNA modules for each dataset)
build_modules.R (Constructs consensus modules in Figure 1)
score_samples.R
biopsy_prediction_validaton_set.R
brush_prediction.R

Output data from running the scripts:

Residual matrices for the datasets used
bx.discovery.resid.rds
br.discovery.resid.rds
bx.validation.resid.rds
br.validation.resid.rds
mouse.resid.rds
tcga.scc.resid.rds

WGCNA output

Consensus Module Construction
correlated_modules_bx_discovery.rds
correlated_modules_br_discovery.rds
correlated_modules_mouse.rds
correlated_modules_tcga_scc.rds
wgcna_gene_list.rds (List of genes in each of the 9 consensus gene modules)


