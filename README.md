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
TCGA SCC data:  TCGA_SCC_cpm.txt

R scripts:
process_data.R
run_wgcna.R
create_modules.R
score_samples.R
biopsy_prediction_validaton_set.R
brush_prediction.R


