######get the raw expression matrix for the 28 gtex tissues used
options(stringsAsFactors = FALSE)
load(file = "data/gtex.cds.rpkm.rda")
gtex.cds.mat <- gtex.cds.rpkm
load(file = "data/gtex.cds.df.rda")

gtex_samples <- read.delim("data/GTEx_Data_V6_Annotations_SampleAttributesDS.txt")
#for some reason, the annotation file contains sample IDs that are not in the gtex matrix column names. remove those sample ids
notinmat <- which(!(gtex_samples$SAMPID %in% colnames(gtex.cds.rpkm)))
sampids <- gtex_samples$SAMPID[-notinmat]             
tissueids <- gtex_samples$SMTSD[-notinmat]
tissue_mat <- vector()
for (i in 1:8555){tissue_mat[i] <- tissueids[which(sampids %in% c(colnames(gtex.cds.rpkm)[i]))]} #this is pretty fast, but probably still not optimal since it's a for loop
gtex_tissues <- unique(tissue_mat)

#define function that will be used in the following line
whichExpressedAboveZero <- function(tissue_name, expr_mat_rpkm){
  gtex.cds.mat <- expr_mat_rpkm #note that this gtex.cds.mat is for the function's local environment, not for the global enrivonment
  selected_tissue <- which(tissue_mat %in% c(tissue_name))
  gtex_selected_tissue <- gtex.cds.mat[, selected_tissue, drop=FALSE]
  median_expr_levels <- rowMedians(gtex_selected_tissue)
  genes_to_keep <- which(median_expr_levels > 0)
  genes <- gtex.cds.df$gene_names[genes_to_keep]
  if (!(length(which(duplicated(genes))) %in% c(0))){
    genes <- genes[-which(duplicated(genes))]
  } 
  genes
}


which.expressed <- sapply(gtex_tissues, function(xx) whichExpressedAboveZero(xx, gtex.cds.mat)) ##this will be used later because the coexpression networks will be constructed with only the expressed genes in that tissue


###load the RPM data (because we have already removed NAs and RPM vs RPKM doesn't make a difference for coexpression)
load(file="data/gtex.cds.rda") #the expression matrix loaded is again named gtex.cds.mat, for convenience

#I now need to run the following lines again
options(stringsAsFactors = FALSE)
gtex_samples <- read.delim("data/GTEx_Data_V6_Annotations_SampleAttributesDS.txt")
#for some reason, the annotation file contains sample IDs that are not in the gtex matrix column names. remove those sample ids
notinmat <- which(!(gtex_samples$SAMPID%in%colnames(gtex.cds.mat)))
sampids <- gtex_samples$SAMPID[-notinmat]             
tissueids <- gtex_samples$SMTSD[-notinmat]
tissue_mat <- vector()
for (i in 1:8550){tissue_mat[i] <- tissueids[which(sampids==colnames(gtex.cds.mat)[i])]} #this is pretty fast, but probably still not optimal since it's a for loop
gtex_tissues <- unique(tissue_mat)

distinct_gtex_tissues <- gtex_tissues[c(1,3,6,12,13,21,27,29,33,34,35,36,37,38,39,40,
                                        41,42,43,44,46,47,48,49,50,51,52,53)]



#for a specific tissue coexpression network the following function should be used:
getExprMat <- function(tissue_name, expr_mat_rpm){
  gtex.cds.mat <- expr_mat_rpm #again note this gtex.cds.mat is only in the local environment of the function, not in the global environment
  gtex_selected_tissue <- gtex.cds.mat[which(gtex.cds.df$gene_names %in%  
                                               which.expressed[[which(gtex_tissues %in% c(tissue_name))]]), 
                                       which(tissue_mat %in% c(tissue_name)), drop = FALSE]
  rownames(gtex_selected_tissue) <- gtex.cds.df$gene_names[which(gtex.cds.df$gene_names 
                                                                 %in% which.expressed[[which(gtex_tissues %in% c(tissue_name))]])]
  gtex_selected_tissue <- gtex_selected_tissue[rownames(gtex_selected_tissue)[-which(duplicated(rownames(gtex_selected_tissue)))], 
                                               , drop = FALSE]
}

raw_expr.mat_list <- lapply(distinct_gtex_tissues, function(xx) getExprMat(xx, gtex.cds.mat))


