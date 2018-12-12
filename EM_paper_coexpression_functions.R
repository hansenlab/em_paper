removePCsFromExprMat <- function(expr_mat, number_of_pcs){
  t(removePrincipalComponents(t(expr_mat), number_of_pcs)) ###unclear how the number of PCs to remove should be determined
}

###function for network construction and module detection
getModules <- function(gene_names, expr_matrix, removeCrossMappable = c(TRUE, FALSE), dissimilarity = c("bicor", "TOM", "SubsetTOM", "Full TOM"), 
                       output = c("module summary", "modules"), modulemerge_threshold, soft_power){
  expressed <- rownames(expr_matrix)[which(rownames(expr_matrix) %in% gene_names)]
  expr.mat <- expr_matrix[which(rownames(expr_matrix) %in% expressed), , drop = FALSE]
  softPower <- soft_power; ###this needs to be determined above, will probably be different depending on the samples
  adjacency <- adjacency(t(expr.mat), power = softPower, type = "unsigned");
  if (removeCrossMappable == TRUE){
    ##this requires a data frame where each row is a pair of genes with crossmapping potential
    ##and the first column is the name of one gene and the second column the name of the other gene
    for (i in 1:length(crossmappable_genes$name1)){
      rowindex <- which(rownames(adjacency) %in% crossmappable_genes$name1[i])
      colindex <- which(colnames(adjacency) %in% crossmappable_genes$name2[i])
      adjacency[rowindex, colindex] <- 0
      adjacency[colindex, rowindex] <- 0
      #adjacency <- make.symmetric(adjacency)
    }
  }
  #  if (removeCrossMappable == TRUE){
  ##this requires a data frame where each row is a pair of genes with crossmapping potential
  ##and the first column is the name of one gene and the second column the name of the other gene
  
  #    adjacency[1:240, 1:240] <- 0
  #adjacency[200:250, 1:150] <- 0
  #  }
  if (dissimilarity %in% c("TOM"))
  {TOM <- TOMsimilarity(adjacency);
  diss <- 1-TOM} else if (dissimilarity %in% "SubsetTOM"){
    TOM <- subsetTOM(t(expr_matrix), subset = which(rownames(expr_matrix) %in% gene_names), power = softPower, networkType = "unsigned")
    diss <- 1-TOM
  } else if (dissimilarity %in% "Full TOM"){
    adjacency <- adjacency(t(expr_matrix), power = softPower, type = "unsigned")
    TOM <- TOMsimilarity(adjacency)
    TOM <- TOM[which(rownames(expr_matrix) %in% gene_names), 
               which(rownames(expr_matrix) %in% gene_names), drop = FALSE]
    diss <- 1-TOM
  } else {
    s <- abs(bicor(t(expr.mat)))
    a <- s^softPower
    diss <- 1-a
  }
  
  geneTree <- hclust(as.dist(diss), method = "average");
  minModuleSize <- 15 ###this has to be determined based on the ribo genes, I choose the minimum module size that still leads to almost all the ribo genes being in one module
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = diss,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = minModuleSize);
  dynamicColors <- labels2colors(dynamicMods)
  if (length(unique(as.vector(dynamicColors))) %in% c(1)){
    if (output == "module summary"){table(dynamicColors)}
    else if (output == "modules"){dynamicColors}
  } else if (!(length(unique(as.vector(dynamicColors))) %in% c(1))){
    MEList <- moduleEigengenes(t(expr.mat), colors = dynamicColors)
    MEs <- MEList$eigengenes
    MEDiss <- 1-cor(MEs);
    METree <- hclust(as.dist(MEDiss), method = "average");
    merge <- mergeCloseModules(t(expr.mat), dynamicColors, 
                               cutHeight = modulemerge_threshold, verbose = 3)
    if (output == "module summary"){table(merge$colors)}
    else if (output == "modules"){merge$colors}
  }
}

createExprMatListWithCommonGenes <- function(list_of_expr.mats, gene_names){
  rownames_list <- lapply(list_of_expr.mats, function(xx) rownames(xx))
  rownames_intersection <- Reduce(intersect, rownames_list)
  new_list_of_expr.mats <- lapply(list_of_expr.mats, function(xx) xx[which(rownames(xx) %in% rownames_intersection), , drop = FALSE])
  new_list_of_expr.mats <- lapply(new_list_of_expr.mats, function(xx) xx[which(rownames(xx) %in% gene_names), , drop = FALSE])
  new_list_of_expr.mats
}


getRandModuleMembershipMatrix <- function(number_of_genes, onlyHighlyExpressed = c(TRUE, FALSE), genes_above_cutoff){
  if (onlyHighlyExpressed %in% c(FALSE)){
    randgeneindices <- sample(12000, number_of_genes)
    randgenes <- rownames(clean_expr_mat_list[[1]])[randgeneindices]
  } else{
    #highly_expressed <- rownames(all_other_genes_medians)[which(all_other_medians_of_medians > expression_cutoff)] 
    rand_samp_highly_expressed <- sample(genes_above_cutoff, number_of_genes)
    #the vector with names of genes expressed above the cutoff is defined outside the function
    ###need to have calculated this vector across those 28 tissues by first getting a matrix where each entry is median 
    ###(across individuals) expression of a given gene in a given tissue, and then getting a vector where each 
    ###coordinate corresponds to a gene and is the row median calculated from that matrix
    randgenes <- rownames(clean_expr_mat_list_pcremoval[[1]])[which(rownames(clean_expr_mat_list_pcremoval[[1]]) 
                                                                    %in% rand_samp_highly_expressed)]
    #crossmappable_genes <- getCrossmappablePairs(randgenes)
  }
  rand.expr.mats.to.use <- createExprMatListWithCommonGenes(clean_expr_mat_list_pcremoval, randgenes)
  rand_tissue_specific_modules_list <- lapply(1:length(distinct_gtex_tissues), function(xx) getModules(randgenes, rand.expr.mats.to.use[[xx]],
                                                                                                       removeCrossMappable = FALSE,
                                                                                                       "TOM", "modules", 0.20, softpowers[xx]))
  #  rand_tissue_specific_modules_list <- list()
  #  for (i in 1:length(distinct_gtex_tissues)){
  #    rand_tissue_specific_modules_list[[i]] <- getModules(randgenes, rand.expr.mats.to.use[[i]], 
  #                                                         "TOM", "modules", 0.20, softpowers[i])
  #    }
  names(rand_tissue_specific_modules_list) <- distinct_gtex_tissues
  rand_module_membership_matrix <- as.matrix(as.data.frame(rand_tissue_specific_modules_list))
  rownames(rand_module_membership_matrix) <- rownames(rand.expr.mats.to.use[[1]]) #doesn't matter which matrix from the list I use all have the same names
  
  #howmanysingletons_inrand <- apply(rand_module_membership_matrix, 1, function(xx) length(which(xx == "grey")))
  #median(howmanysingletons_inrand)
  #hist(howmanysingletons_inrand)
  #sometimes there's also other colors so make sure those are handled as well, otherwise there will be NAs
  #in the matrix and those lead to inaccurate distances
  module_colors <- unique(as.vector(rand_module_membership_matrix))
  rest_of_the_colors <- module_colors[-which(module_colors %in% c("turquoise", "blue", "brown", "green", "yellow",
                                                                  "red", "grey"))]
  
  rand_module_membership_matrix[which(rand_module_membership_matrix %in% c("turquoise"))] <- 10
  rand_module_membership_matrix[which(rand_module_membership_matrix %in% c("blue"))] <- 40
  rand_module_membership_matrix[which(rand_module_membership_matrix %in% c("brown"))] <- 70
  rand_module_membership_matrix[which(rand_module_membership_matrix %in% c("green"))] <- 100
  rand_module_membership_matrix[which(rand_module_membership_matrix %in% c("yellow"))] <- 130
  rand_module_membership_matrix[which(rand_module_membership_matrix %in% c("red"))] <- 160
  rand_module_membership_matrix[which(rand_module_membership_matrix %in% c("grey"))] <- 190
  
  if (!(length(rest_of_the_colors) %in% c(0))){
    for (i in 1:length(rest_of_the_colors)){
      rand_module_membership_matrix[which(rand_module_membership_matrix %in% c(rest_of_the_colors[i]))] <- (190+i)}
  }
  
  class(rand_module_membership_matrix) <- "numeric" 
  for (i in 1:nrow(rand_module_membership_matrix)){
    rand_module_membership_matrix[i,][which(rand_module_membership_matrix[i,] 
                                            %in% c(190))] <- (-i)}
  rand_module_membership_matrix
}


getRandomCoexpressionMeasure <- function(rand_mod_memb_matrix, output = c("plot_density", "insert_density",
                                                                          "get_number_of_coexpressed_genes", 
                                                                          "percentiles"), 
                                         distance_threshold, coexpression_threshold){  
  rand_rowlist <- as.list(data.frame(t(rand_mod_memb_matrix)))
  rand_distancelist <- mclapply(rand_rowlist, 
                                function(xx) sapply(1:nrow(rand_mod_memb_matrix), 
                                                    function(x) seq_dist(rand_rowlist[[x]], xx, 
                                                                         method = "hamming")), mc.cores = 4)
  rand_distancevector <- as.numeric(lapply(rand_distancelist, 
                                           function(xx) length(which(xx < distance_threshold))))
  
  if (output %in% c("plot_density")){plot(density(rand_distancevector, from = 0), 
                                          xlab = "number of stable coexpressed neighbors in more than half of the tissues", 
                                          main = "random coexpression across tissues", 
                                          xlim = c(0, 50), ylim = c(0, 0.25))}
  else if (output %in% c("insert_density")){lines(density(rand_distancevector, from = 0))}
  else if (output %in% c("get_number_of_coexpressed_genes")){
    length(which(rand_distancevector > coexpression_threshold))
  } else if (output %in% c("percentiles")){
    c(quantile(rand_distancevector, 0.95), quantile(rand_distancevector, 0.85), 
      quantile(rand_distancevector, 0.75))
  }
}


getDistanceVector <- function(mod_memb_matrix, tissue_cutoff){
  for (i in 1:nrow(mod_memb_matrix)){
    mod_memb_matrix[i,][which(mod_memb_matrix[i,] 
                              %in% c(190))] <- (-i)}
  
  module_membership_matrix_rowlist <- as.list(data.frame(t(mod_memb_matrix)))
  distancelist <- lapply(module_membership_matrix_rowlist, 
                         function(xx) sapply(1:nrow(mod_memb_matrix), 
                                             function(x) 
                                               seq_dist(module_membership_matrix_rowlist[[x]], xx, 
                                                        method = "hamming")))
  
  
  distance_vector <- as.numeric(lapply(distancelist, 
                                       function(xx) length(which(xx < tissue_cutoff))))
  distance_vector
  
}
