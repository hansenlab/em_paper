# em_paper
All the scripts below rely on a data frame "epiGenesDF", which corresponds to Supplemental Table 2 in the paper.

* Scripts for the variation-tolerance analysis:
1. EM_paper_variation_tolerance_objects.R.  
This script depends on a data frame "exactab", which is downloaded from the ExAC browser. The "isEM" column is defined based on the epiGenesDF data frame, and the "isTF" based on the list of transcription factors downloaded from Barrera et al., 2016.
2. EM_paper_variation_tolerance_figures.R. 
3. EM_paper_accessory_EM_subunit_analysis_and_figures.R
4. EM_paper_substrate_specificities_analysis_and_figures.R.  

* Scripts for the local domain constraint analysis:
1. EM_paper_domain_constraint_functions.R. 
2. EM_paper_domain_constraint_binary.R. 
3. EM_paper_domain_constraint_analysis_quantitative_version.R.  
4. EM_paper_domain_constraint_quantitative_version_figures.R.   
The above 4 scripts depend on a bed file "EMdomainsvCCR_2.bed", which is provided as supplemental table 12 in our paper.
5. EM_paper_map_protein_domain_to_genomic_coordinates.R
This script was used to define the protein domain genomic coordinates of the EM genes used in the local mutational constraint analysis.

* Scripts for the co-expression analysis:

1. EM_paper_gtex_read.R.   
This script reads the gtex raw gene-level quantifications from GTEx into R, and transforms them into log2(RPM+1).
2. EM_paper_28_raw_expression_matrices.R.   
This script obtains the expression matrices subsequently used in the co-expression analyses. It depends on the log2(RPM+1) expression data obtained using the script "EM_paper_gtex_read.R", but also on the log2(RPKM+1) expression data (for filtering on expression level), where the RPKM expression matrix is downloaded from the GTEx portal. The "GTEx_Data_V6_Annotations_SampleAttributesDS.txt" annotation file is also downloaded from the GTEx portal.

3. EM_paper_coexpression_functions.R. 
4. EM_paper_main_coexpression_analysis.R. 
5. EM_paper_coexpression_neuro_enrichment.R.   
The above 3 scripts form the core of the co-expression analysis, where the log2(RPM+1) expression matrices obtained via the script "EM_paper_28_raw_expression_matrices.R" are preprocessed, co-expression networks are estimated for each tissue, and subsequently the results are integrated across tissues. Finally, the enrichment of genes associated with neurological dysfunction in the highly co-expressed group is tested.

6. EM_paper_coexpression_cormat_thresholding.R.   
This script estimates the co-expression network for each tissue by thresholding the corresponding correlation matrix.
7. EM_paper_coexpression_subsampling.R.   
This script examines the robustness of the co-expression analysis to sample outliers.


* Scripts for the tissue specificity and expression level analyses:
1. EM_paper_tissue_specificity_analysis.R
2. EM_paper_tissue_specificity_analysis_batch_correction.R
3. EM_paper_expression_levels_tissue_specificity_figures.R.  
These scripts examine the tissue specificity and expression levels of EM genes, within the 28 tissues also used in the co-expression analysis. They use the log2(RPKM+1) expression data, where the RPKM expression matrix is downloaded from the GTEx portal. They also use the "GTEx_Data_V6_Annotations_SampleAttributesDS.txt" annotation file, which is also downloaded from the GTEx portal.

* Scripts for the trans-acting factor binding analysis:   
1. EM_paper_TF_chipseq_read.R
2. EM_paper_k562_analysis_and_figures.R






