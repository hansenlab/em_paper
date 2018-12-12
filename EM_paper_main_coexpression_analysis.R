load(file = "data/gtex_28_tissues_raw_expr_mat_list.rda")

mean_centered_raw.expr_mat_list <- lapply(raw_expr.mat_list, function(xx) {
  exprmat <- xx
  rowmeans <- rowMeans(exprmat)
  exprmat_centered <- sweep(exprmat , 1, FUN = "-", rowmeans)
  rowstdeviations <- rowSds(exprmat)
  exprmat_standardized <- sweep(exprmat_centered, 1, FUN="/", rowstdeviations)
  exprmat_standardized
}) 

clean_expr_mat_list <- lapply(mean_centered_raw.expr_mat_list, removePCsFromExprMat, number_of_pcs = 4)

saveSoftPowerPlot <- function(expr_mat, tissue, file_name){
  quartz(file = paste0("coexpression_softpowers/softpowers_pcremoval/", file_name, ".pdf"), 
         width = 7, height = 7, type = "pdf")
  determineSoftPower(expr_mat, tissue)
  dev.off()
}

filenames <- c("adsub", "adrenal", "arttibial", "cereb", "cortex", "mammarytissue", "colon_transv", 
               "esoph_mucosa", "heartleftventr", "kidneycortex", "liver", "lung", "minorsalivarygland", 
               "muscle_skel", "nervetibial", "ovary", "pancreas", "pituitary", "prostate", 
               "skin_notsunexposed", "small_intestine_terminal_ileum", "spleen", "stomach", 
               "testis", "thyroid", "uterus", "vagina", "wholeblood")


determineSoftPower <- function(expr_mat, tissue_name){
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(t(expr_mat), powerVector = powers, verbose = 0)
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste0("Scale independence_", tissue_name));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.90,col="red")
}


saveSoftPowerPlot <- function(expr_mat, tissue, file_name){
  quartz(file = paste0("coexpression_softpowers/softpowers2/", file_name, ".pdf"), 
         width = 7, height = 7, type = "pdf")
  determineSoftPower(expr_mat, tissue)
  dev.off()
}

for (i in 1:length(distinct_gtex_tissues)){
  saveSoftPowerPlot(clean_expr_mat_list[[i]], distinct_gtex_tissues[[i]], filenames[[i]])
}

###after looking at the softpower plots
#pcremoval_softpowers <- c(7, 6, 7, 8, 14, 8, 10, 6, 8, 7, 7, 8, 7, 6, 9, 8, 
#                          8, 9, 9, 9, 9, 7, 8, 10, 8, 8, 5) those are wrong I had forgotten to center the data before PC removal
clean_expr_mat_list_pcremoval <- list()
for (i in 1:8){
  clean_expr_mat_list_pcremoval[[i]] <- clean_expr_mat_list[[i]]
}
clean_expr_mat_list_pcremoval[[9]] <- clean_expr_mat_list[[23]]
for (i in 10:23){
  clean_expr_mat_list_pcremoval[[i]] <- clean_expr_mat_list[[i-1]]
}
for (i in 24:28){
  clean_expr_mat_list_pcremoval[[i]] <- clean_expr_mat_list[[i]]
}
pcremoval_softpowers <- c(6, 6, 9, 6, 7, 7, 5, 5, 6, 7, 9, 7, 4, 9, 7, 10, 7, 6, 9, 7, 8, 6, 6, 7, 7, 6, 6, 10) #note I have moved stomach after esophagus - mucosa in the ordering of the tissues and thus in the ordering of softpowers as well

softpowers <- pcremoval_softpowers


#nowcreate a list with the clean matrices by keeping only the common genes
#load EM genes
load(file = "epigenetic_list/objects/epiGenesDF.rda")
epiGenes <- epiGenesDF$Gene_name_gtex

list.of.expr.mats.to.use <- createExprMatListWithCommonGenes(clean_expr_mat_list_pcremoval, epiGenes)


tissue_specific_modules_list <- list()
for (i in 1:length(distinct_gtex_tissues)){
  tissue_specific_modules_list[[i]] <- getModules(epiGenes, list.of.expr.mats.to.use[[i]], 
                                                  removeCrossMappable = FALSE,
                                                  "TOM", "modules", 0.20, softpowers[i])
}

names(tissue_specific_modules_list) <- distinct_gtex_tissues[c(1:8, 23, 9:22, 24:28)]
module_membership_matrix <- as.matrix(as.data.frame(tissue_specific_modules_list))
rownames(module_membership_matrix) <- rownames(list.of.expr.mats.to.use[[1]]) #doesn't matter which matrix from the list I use all have the same names

module_colors <- unique(as.vector(module_membership_matrix))
rest_of_the_colors <- module_colors[-which(module_colors %in% c("turquoise", "blue", "brown", "green", "yellow",
                                                                "red", "grey"))]

module_membership_matrix[which(module_membership_matrix %in% c("turquoise"))] <- 10
module_membership_matrix[which(module_membership_matrix %in% c("blue"))] <- 40
module_membership_matrix[which(module_membership_matrix %in% c("brown"))] <- 70
module_membership_matrix[which(module_membership_matrix %in% c("green"))] <- 100
module_membership_matrix[which(module_membership_matrix %in% c("yellow"))] <- 130
module_membership_matrix[which(module_membership_matrix %in% c("red"))] <- 160
module_membership_matrix[which(module_membership_matrix %in% c("grey"))] <- 190

if (!(length(rest_of_the_colors) %in% c(0))){
  for (i in 1:length(rest_of_the_colors)){
    module_membership_matrix[which(module_membership_matrix %in% c(rest_of_the_colors[i]))] <- (190+i)}
}

class(module_membership_matrix) <- "numeric" 

###now perform the module membership analysis across tissues
library(stringdist)
#load(file = "figures/coexpression/ruv_module_membership_matrix.rda")
#in order to get accurate distances between genes, replace the number corresponding to the grey module (which represents
#singletons) with a negative number specific to each row. Otherwise, the hamming distance considers the singletons
#as being coexpressed in the same module
for (i in 1:nrow(module_membership_matrix)){
  module_membership_matrix[i,][which(module_membership_matrix[i,] 
                                     %in% c(190))] <- (-i)}

module_membership_matrix_rowlist <- as.list(data.frame(t(module_membership_matrix)))
distancelist <- lapply(module_membership_matrix_rowlist, 
                       function(xx) sapply(1:nrow(module_membership_matrix), 
                                           function(x) 
                                             seq_dist(module_membership_matrix_rowlist[[x]], xx, 
                                                      method = "hamming")))


em_distancevector <- as.numeric(lapply(distancelist, 
                                       function(xx) length(which(xx < 19))))

inverse_distance_list <- lapply(distancelist, function(xx) rep(28, length(distancelist)) - xx)
flattened <- unlist(inverse_distance_list)
adjmatrix <- t(matrix(flattened, nrow = length(distancelist)))


adjmatrix[which(adjmatrix > 14)] <- 30
adjmatrix[which(adjmatrix < 14 & adjmatrix > 9)] <- 31
adjmatrix[which(!(adjmatrix %in% c(30, 31)))] <- 32

adjmatrix <- adjmatrix[order(em_distancevector, decreasing = TRUE), order(em_distancevector, decreasing = TRUE), drop = FALSE]
rownames(adjmatrix) <- rownames(module_membership_matrix)[order(em_distancevector, decreasing = TRUE)]


##########plot partnership heatmap
colfunc <- colorRampPalette(c("orange", "lightskyblue3", "gray25"))
heatmap(adjmatrix, colfunc(3), scale = "none", Colv = NA, Rowv = NA, 
        labRow = "", labCol = "", add.expr = c(segments(76, 0, 76, 76, lwd = 5), 
                                               segments(0, 76, 76, 76, lwd = 5), 
                                               segments(159, 0, 159, 159, lwd = 5), 
                                               segments(0, 159, 159, 159, lwd = 5)))


