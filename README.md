# em_paper
All the scripts below rely on a data frame "epiGenesDF", which corresponds to Supplemental Table 2 in the paper.

Scripts for the variation-tolerance analysis:
1. EM_paper_variation_tolerance_objects.R
This script depends on a data frame "exactab", which is downloaded from the exac browser.
2. EM_paper_variation_tolerance_figures.R

Scripts for co-expression analyses:

1. EM_paper_gtex_read.R
This script reads the gtex raw gene-level quantifications from GTEx into R, and transforms them into log2(RPM+1).
2. EM_paper_28_raw_expression_matrices.R
This script obtains the expression matrices subsequently used in the co-expression analyses. It depends on the log2(RPM+1) expression data obtained using the script "EM_paper_gtex_read.R", but also on the log2(RPKM+1) expression data (for filtering on expression level), where the RPKM expression matrix is downloaded from the GTEx portal. The "GTEx_Data_V6_Annotations_SampleAttributesDS.txt" annotation file is also downloaded from the GTEx portal.

3. EM_paper_coexpression_functions.R
4. EM_paper_main_coexpression_analysis.R
5. EM_paper_coexpression_neuro_enrichment.R
The above 3 scripts form the core of the co-expression analysis, where the log2(RPM+1) expression matrices obtained via the script "EM_paper_28_raw_expression_matrices.R" are preprocessed, co-expression networks are estimated for each tissue, and subsequently the results are integrated across tissues. Finally, the enrichment of genes associated with neurological dysfunction in the highly co-expressed group is tested.

6. EM_paper_coexpression_cormat_thresholding.R
This script estimates the co-expression network for each tissue by thresholding the corresponding correlation matrix.
7. EM_paper_coexpression_subsampling.R
This script examines the robustness of the co-expression analysis to sample outliers.





