load(file = "epigenetic_list/objects/epiGenesDF.rda")
load(file = "epigenetic_list/objects/exactab.rda")
library(matrixStats)
library(cummeRbund)

##define the following function
getMedianExprLevels <- function(gene_names, tissue_name){
  selected_tissue <- which(tissue_mat %in% c(tissue_name))
  gtex_selected_tissue <- gtex.cds.mat[, selected_tissue, drop=FALSE]
  median_expr_levels <- rowMedians(gtex_selected_tissue)
  
  gtex_input_genes <- which(gtex.cds.df$gene_names %in% gene_names)
  ###now remove duplicates
  genes <- gtex.cds.df$gene_names[which(gtex.cds.df$gene_names %in% gene_names)]
  gtex_input_genes <- gtex_input_genes[-which(duplicated(genes))]
  ###
  median_expr_levels[gtex_input_genes]
}

###import the expression data
load(file = "data/gtex.cds.rpkm.rda")
load(file = "data/gtex.cds.df.rda")
gtex.cds.mat <- gtex.cds.rpkm #this is necessary given the way the above functions are defined

#also import the meta data
options(stringsAsFactors = FALSE)
gtex_samples <- read.delim("data/GTEx_Data_V6_Annotations_SampleAttributesDS.txt")
#for some reason, the annotation file contains sample IDs that are not in the gtex matrix column names. remove those sample ids
notinmat <- which(!(gtex_samples$SAMPID %in% colnames(gtex.cds.rpkm)))
sampids <- gtex_samples$SAMPID[-notinmat]             
tissueids <- gtex_samples$SMTSD[-notinmat]
tissue_mat <- vector()
for (i in 1:8555){tissue_mat[i] <- tissueids[which(sampids %in% c(colnames(gtex.cds.rpkm)[i]))]} #this is basically instantaneous, but probably still not optimal since it's a for loop
gtex_tissues <- unique(tissue_mat)
distinct_gtex_tissues <- gtex_tissues[c(1,3,6,12,13,21,27,29,33,34,35,36,37,38,39,40,
                                        41,42,43,44,46,47,48,49,50,51,52,53)]
###define classes of genes that will be used
epiGenes <- epiGenesDF$Gene_name_gtex
library(gdata)
tflist <- read.xls('data/Barrera_2016_Science_TableS2.xls')
names(tflist) <- c("ensembl", "Gene.Symbol")
tflist[,1] <- as.character(tflist[,1])
tflist[,2] <- as.character(tflist[,2])
tf_genes <- tflist$Gene.Symbol
tf_genes <- tf_genes[-c(1,2)]

#and now keep only those that are in gtex
input_epig <- gtex.cds.df$gene_names[which(gtex.cds.df$gene_names %in% epiGenes)]
input_tf <- gtex.cds.df$gene_names[which(gtex.cds.df$gene_names %in% tf_genes)]
input_allother <- gtex.cds.df$gene_names[-which(gtex.cds.df$gene_names %in% c(epiGenes, tf_genes))]



###get data frames with rows being genes and columns tissues
###each entry is median expr level (across individuals) of a given gene in a given tissue
###uncorrected here just means that in this tissue specificity analysis we will not correct for batch effects etc
epi_medians_uncorrected <- sapply(distinct_gtex_tissues, getMedianExprLevels, gene_names = input_epig)
epi_medians_uncorrected <- as.data.frame(epi_medians_uncorrected)
tf_medians_uncorrected <- sapply(distinct_gtex_tissues, getMedianExprLevels, gene_names = input_tf)
tf_medians_uncorrected <- as.data.frame(tf_medians_uncorrected)
all_other_genes_medians_uncorrected <- sapply(distinct_gtex_tissues, getMedianExprLevels, gene_names = input_allother)
all_other_genes_medians_uncorrected <- as.data.frame(all_other_genes_medians_uncorrected)


#####name the rows (corresponding to genes) of each dataframe
#epi
gtex_input_genes <- which(gtex.cds.df$gene_names %in% input_epig)
###now remove duplicates
genes <- gtex.cds.df$gene_names[which(gtex.cds.df$gene_names %in% input_epig)]
gtex_input_genes <- gtex_input_genes[-which(duplicated(genes))]

genenames <- gtex.cds.df$gene_names[gtex_input_genes]
rownames(epi_medians_uncorrected) <- genenames

#tf
gtex_input_genes <- which(gtex.cds.df$gene_names %in% input_tf)
###now remove duplicates
genes <- gtex.cds.df$gene_names[which(gtex.cds.df$gene_names %in% input_tf)]
gtex_input_genes <- gtex_input_genes[-which(duplicated(genes))]

genenames <- gtex.cds.df$gene_names[gtex_input_genes]
rownames(tf_medians_uncorrected) <- genenames

#all other
gtex_input_genes <- which(gtex.cds.df$gene_names %in% input_allother)
###now remove duplicates
genes <- gtex.cds.df$gene_names[which(gtex.cds.df$gene_names %in% input_allother)]
gtex_input_genes <- gtex_input_genes[-which(duplicated(genes))]

genenames <- gtex.cds.df$gene_names[gtex_input_genes]
rownames(all_other_genes_medians_uncorrected) <- genenames

citric_genes <- read.delim(file = "data/protein_class_Citric.tsv") #remember strings as factors
citric_genes <- citric_genes[[1]]

citric_medians_uncorrected <- all_other_genes_medians_uncorrected[which(rownames(all_other_genes_medians_uncorrected) 
                                 %in% citric_genes), , drop = FALSE]


epi_prob_vectors_uncorrected <- makeprobs(t(as.matrix(epi_medians_uncorrected)))
tf_prob_vectors_uncorrected <- makeprobs(t(as.matrix(tf_medians_uncorrected)))
allother_prob_vectors_uncorrected <- makeprobs(t(as.matrix(all_other_genes_medians_uncorrected)))
citric_prob_vectors_uncorrected <- makeprobs(t(as.matrix(citric_medians_uncorrected)))

EM_spec_scores_uncorrected <- sapply(distinct_gtex_tissues, function (xx) 
  1 - JSdistFromP(t(epi_prob_vectors_uncorrected), diag(1, nrow = 28)[which(distinct_gtex_tissues == xx), ]))
EM_spec_scores_uncorrected <- as.matrix(EM_spec_scores_uncorrected)



#######
em_distance_from_uniform_uncorrected <- JSdistFromP(t(epi_prob_vectors_uncorrected), rep(1/28, 28))
em_distance_from_uniform_uncorrected <- em_distance_from_uniform_uncorrected[-which(
  rownames(epi_medians_uncorrected) 
  %in% epiGenesDF$Gene_name_gtex_rpkm[which(epiGenesDF$Gene_name %in% em_and_tf)])]
em_new_rownames <- rownames(epi_medians_uncorrected)[-which(
  rownames(epi_medians_uncorrected) 
  %in% epiGenesDF$Gene_name_gtex_rpkm[which(epiGenesDF$Gene_name %in% em_and_tf)])]

tf_distance_from_uniform_uncorrected <- JSdistFromP(t(tf_prob_vectors_uncorrected), rep(1/28, 28))
#tf_distance_from_uniform[which(is.na(tf_distance_from_uniform))] <- rep(1, 
#                                                                        length(which(is.na(tf_distance_from_uniform))))
tf_distance_from_uniform_uncorrected <- tf_distance_from_uniform_uncorrected[-which(
  rownames(tf_medians_uncorrected) %in% epiGenesDF$Gene_name_gtex_rpkm[which(epiGenesDF$Gene_name 
                                                                             %in% em_and_tf)])]
tf_new_rownames <- rownames(tf_medians_uncorrected)[-which(
  rownames(tf_medians_uncorrected) 
  %in% epiGenesDF$Gene_name_gtex_rpkm[which(epiGenesDF$Gene_name %in% em_and_tf)])]

allother_distance_from_uniform_uncorrected <- JSdistFromP(t(allother_prob_vectors_uncorrected), rep(1/28, 28))

citric_distance_from_uniform_uncorrected <- JSdistFromP(t(citric_prob_vectors_uncorrected), rep(1/28, 28))


###the following will be used when selecting genes expressed at a certain cutoff
epi_medians_of_medians <- rowMedians(as.matrix(epi_medians_uncorrected))
tf_medians_of_medians <- rowMedians(as.matrix(tf_medians_uncorrected))
all_other_medians_of_medians <- rowMedians(as.matrix(all_other_genes_medians_uncorrected))


