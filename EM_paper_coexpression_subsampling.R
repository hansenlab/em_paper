####construct the network for EM genes by subsampling to ensure robustness to outliers etc


####the following is calculated with the module_membership_matrix that was built using
###the coexpression network which was constructed using all of the samples in each tissue
original_highly_coexpressed <- rownames(module_membership_matrix)[which(em_distancevector > 74)]

###now see what happens if we subsample
getEMnetworkAfterSubsampling <- function(expr.mat.list, genes){
  clean_expr_mat_list_pcremoval_subsampled <- lapply(expr.mat.list, function(xx) {
    mat <- xx
    if (ncol(mat) < 100 & ncol(mat) > 50){
      mat <- mat[, sample(ncol(mat), ncol(mat) - 20), drop = FALSE]
    }
    else if (ncol(mat) < 30){
      mat <- mat[, sample(ncol(mat), ncol(mat) - 5), drop = FALSE]
    }
    else {
      mat <- mat[, sample(ncol(mat), ncol(mat)/2), drop = FALSE]
    }
    mat
  })
  list.of.expr.mats.to.use <- createExprMatListWithCommonGenes(clean_expr_mat_list_pcremoval_subsampled, genes)
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
  
  list(module_membership_matrix, length(which(em_distancevector > 74)), length(which(original_highly_coexpressed 
                                                                                     %in% rownames(module_membership_matrix)[which(em_distancevector > 74)])))
}



getEMnetworkAfterSubsampling(clean_expr_mat_list_pcremoval, epiGenes)








####
em_network_after_subsampling <- get(load(file = "data/em_network_after_subsampling.rda"))

subsampling_mod_memb_matrix_list <- lapply(seq(1, 900, by = 3), function(xx) em_network_after_subsampling[[xx]])

subsampling_coexpressed_neighbors_list <- lapply(subsampling_mod_memb_matrix_list, function(xx) 
  getRandomCoexpressionMeasure(xx, "get_number_of_coexpressed_genes", 19, 74))

subsampling_overlap_with_original_highly_coexpressed <- lapply(subsampling_mod_memb_matrix_list, function(xx){
  module_membership_matrix_rowlist <- as.list(data.frame(t(xx)))
  distancelist <- lapply(module_membership_matrix_rowlist, 
                         function(yy) sapply(1:nrow(xx), 
                                             function(x) 
                                               seq_dist(module_membership_matrix_rowlist[[x]], yy, 
                                                        method = "hamming")))
  
  em_distancevector <- as.numeric(lapply(distancelist, 
                                         function(yy) length(which(yy < 19))))
  length(which(rownames(xx)[which(em_distancevector > 74)] %in% rownames(adjmatrix)[1:74]))
})

