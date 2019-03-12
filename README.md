# Computer-code-for-PML-Molecular-Subtypes
R scripts and input data to produce results in the manuscript:

Molecular subtyping reveals immune alterations associated with progression of bronchial premalignant lesions.

Authors:  Jennifer Beane*, Sarah A. Mazzilli, Joshua D. Campbell, Grant Duclos, Kostyantyn Krysan, Christopher Moy, Catalina Perdomo, Michael Schaffer, Gang Liu, Sherry Zhang, Hangqio Liu, Jessica Vick, Samjot S. Dhillon, Suso J. Platero, Steven M. Dubinett, Christopher Stevenson, Mary E. Reid, Marc E. Lenburg, Avrum E. Spira

Nature Communications 2019

Input data:<br />
Endobronchial biopsy and brush data:  GSE109743_exp_count.txt<br />
Endobronchial biopsy and brush phenotypic data:  GSE109743_SAMPLES.txt<br />
NTCU mouse model data:  GSE111091_exp_count.txt<br />
NTCU mouse model phenotypic data:  GSE111091_SAMPLES.txt<br />
TCGA SCC data:  TCGA_SCC_cpm.txt (Campbell et al. Nature Genetics 2016, Log2 plus 1 TPM)<br />
Genes used to predict genomic smoking status:  smoking_signature_genes.txt and smoking_signature_data.txt (derived from GSE7895)<br />

R scripts:<br />
processData_and_runWGCNA.R (processes data and compute WGCNA modules for each dataset)<br />
build_modules.R (Constructs consensus modules in Figure 1)<br />
conduct_consensus_clustering.R (Performs consensus clustering on consensus modules)<br />
predict_smoking_status.R (Predicts smoking status for each patient/timepoint pair)<br />
select_msubtype_genes.R (Selects genes for molecular subtype classifier)<br />
msubtype_classifier.R (Predict molecular subtype in validation set and Proliferative or not in brushes)<br />
WGCNA_wrapper.R (used in processData_and_runWGCNA.R)<br />
wv_alg.R (used in predict_smoking_status.R)<br />

Output data from running the scripts:<br />

Residual matrices for the datasets used<br />
bx.discovery.resid.rds<br />
br.discovery.resid.rds<br />
bx.validation.resid.rds<br />
br.validation.resid.rds<br />
mouse.resid.rds<br />
tcga.scc.resid.rds<br />

WGCNA output<br />
bx.discovery.wgcna.rds<br />
br.discovery.wgcna.rds<br />
mouse.wgcna.rds<br />
tcga.scc.wgcna.rds<br />

Consensus Module Construction<br />
correlated_modules_bx_discovery.rds<br />
correlated_modules_br_discovery.rds<br />
correlated_modules_mouse.rds<br />
correlated_modules_tcga_scc.rds<br />
wgcna_gene_list.rds (List of genes in wgcna modules for each dataset above)<br />
combined_gene_set.rds (List of genes in each of the 9 consensus gene modules - Ensembl IDs)<br />
combined_gene_set_symbols.rds (List of genes in each of the 9 consensus gene modules - HGNC symbols)<br />

Consensus Clustering<br />
bx_discovery_resid_ConsensusCluster.rda<br />

Molecular Subtype Classifer<br />
molecular.subtype.classifier.genes.rds<br />
bx_discovery_msubtype_classifier_results.rds<br />
bx_validation_msubtype_classifier_results.rds<br />
bx_discovery_msubtype_classifier_results_Prolif.or.Not.rds<br />
br_discovery_msubtype_classifier_results_Prolif.or.Not.rds<br />
bx_validation_msubtype_classifier_results_Prolif.or.Not.rds<br />
br_validation_msubtype_classifier_results_Prolif.or.Not.rds<br />

