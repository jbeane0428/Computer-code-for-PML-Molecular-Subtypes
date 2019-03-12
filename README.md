# Computer-code-for-PML-Molecular-Subtypes
R scripts and input data to produce results in the manuscript:

Molecular subtyping reveals immune alterations associated with progression of bronchial premalignant lesions.

Authors:  Jennifer Beane*, Sarah A. Mazzilli, Joshua D. Campbell, Grant Duclos, Kostyantyn Krysan, Christopher Moy, Catalina Perdomo, Michael Schaffer, Gang Liu, Sherry Zhang, Hangqio Liu, Jessica Vick, Samjot S. Dhillon, Suso J. Platero, Steven M. Dubinett, Christopher Stevenson, Mary E. Reid, Marc E. Lenburg, Avrum E. Spira

Nature Communications 2019

Input data:
Endobronchial biopsy and brush data:  GSE109743_exp_count.txt\n
Endobronchial biopsy and brush phenotypic data:  GSE109743_SAMPLES.txt\n
NTCU mouse model data:  GSE111091_exp_count.txt
NTCU mouse model phenotypic data:  GSE111091_SAMPLES.txt
TCGA SCC data:  TCGA_SCC_cpm.txt (Campbell et al. Nature Genetics 2016, Log2 plus 1 TPM)
Genes used to predict genomic smoking status:  smoking_signature_genes.txt (derived from GSE7895)

R scripts:
processData_and_runWGCNA.R (processes data and compute WGCNA modules for each dataset)
build_modules.R (Constructs consensus modules in Figure 1)
conduct_consensus_clustering.R (Performs consensus clustering on consensus modules)
predict_smoking_status.R (Predicts smoking status for each patient/timepoint pair)
select_msubtype_genes.R (Selects genes for molecular subtype classifier)
msubtype_classifier.R (Predict molecular subtype in validation set and Proliferative or not in brushes)
WGCNA_wrapper.R (used in processData_and_runWGCNA.R)
wv_alg.R (used in predict_smoking_status.R)

Output data from running the scripts:

Residual matrices for the datasets used
bx.discovery.resid.rds
br.discovery.resid.rds
bx.validation.resid.rds
br.validation.resid.rds
mouse.resid.rds
tcga.scc.resid.rds

WGCNA output
bx.discovery.wgcna.rds
br.discovery.wgcna.rds
mouse.wgcna.rds
tcga.scc.wgcna.rds

Consensus Module Construction
correlated_modules_bx_discovery.rds
correlated_modules_br_discovery.rds
correlated_modules_mouse.rds
correlated_modules_tcga_scc.rds
wgcna_gene_list.rds (List of genes in wgcna modules for each dataset above)
combined_gene_set.rds (List of genes in each of the 9 consensus gene modules - Ensembl IDs)
combined_gene_set_symbols.rds (List of genes in each of the 9 consensus gene modules - HGNC symbols)

Consensus Clustering
bx_discovery_resid_ConsensusCluster.rda

Molecular Subtype Classifer
molecular.subtype.classifier.genes.rds
bx_discovery_msubtype_classifier_results.rds
bx_validation_msubtype_classifier_results.rds
bx_discovery_msubtype_classifier_results_Prolif.or.Not.rds
br_discovery_msubtype_classifier_results_Prolif.or.Not.rds
bx_validation_msubtype_classifier_results_Prolif.or.Not.rds
br_validation_msubtype_classifier_results_Prolif.or.Not.rds

