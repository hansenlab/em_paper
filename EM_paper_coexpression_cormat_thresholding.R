###coexpression network construction with correlation matrix thresholding
library(matrixStats)
library(Rgraphviz)

getCorVectorFromCorMat <- function(cor_mat){
  m <- cor_mat
  m[lower.tri(m)]<-NA
  diag(m)<-NA
  cors <- as.vector(m)
  cors <- cors[-which(is.na(cors))]
  cors
}


getGenesInMaxConnComp <- function(expr_mat, threshold){
  corMat <- cor(t(expr_mat))
  round(corMat,2)
  diag(corMat) <- 0
  #round(rowMaxs(corMat), 2)
  #round(rowMins(corMat), 2)
  #corMat <- corMat[order(rownames(corMat)), order(colnames(corMat))]
  adj <- abs(corMat)
  adj[adj > threshold] <- 1
  adj[adj < threshold] <- 0
  diag(adj) <- 0
  gA <- graphAM(adjMat = adj, edgemode = "undirected")
  gNEL <- as(gA, "graphNEL")
  complengths <- lengths(connComp(gNEL))
  maxcompsize <- max(complengths)
  connComp(gNEL)[which(lengths(connComp(gNEL)) %in% maxcompsize)][[1]]}

getCompSize <- function(expr_mat, threshold, output = c("max", "all")){
  corMat <- cor(t(expr_mat))
  round(corMat,2)
  diag(corMat) <- 0
  #round(rowMaxs(corMat), 2)
  #round(rowMins(corMat), 2)
  #corMat <- corMat[order(rownames(corMat)), order(colnames(corMat))]
  adj <- abs(corMat)
  adj[adj > threshold] <- 1
  adj[adj < threshold] <- 0
  diag(adj) <- 0
  gA <- graphAM(adjMat = adj, edgemode = "undirected")
  gNEL <- as(gA, "graphNEL")
  complengths <- lengths(connComp(gNEL))
  if (output %in% c("max")){
    max(complengths)
  } else {
    complengths
  }
}

clean_expr_mat_list_onlyabove05RPKM <- lapply(clean_expr_mat_list_pcremoval, function(xx) 
  xx[which(rownames(xx) %in% rownames(all_other_genes_medians)[
    which(all_other_medians_of_medians > 0.5)]), , drop = FALSE])

getRandCompSize <- function(expr_mat, number_of_genes, threshold, output = c("max", "all")){
  randsample <- sample(1:nrow(expr_mat), number_of_genes)
  getCompSize(expr_mat[randsample, , drop = FALSE], threshold, output)
}


#can determine thresholds in an unbiased way by looking at the distribution of correlations between
#random genes (and the distribution of ribosomal genes)
#
thresholds <- sapply(1:27, function(xx) {
  randsample <- sample(1:nrow(clean_expr_mat_list_onlyabove05RPKM[[xx]]), 2000)
  randcor <- cor(t(clean_expr_mat_list_onlyabove05RPKM[[xx]][randsample, , drop = FALSE]))
  randcorvector <- getCorVectorFromCorMat(randcor)
  
  round(quantile(randcorvector, 0.998), 2) #this procedure is sensitive to the exact cutoff
})



###the two entries in the above vector that are equal to 1 indicate that the corresponding columns of the
###maxconncomp_membership_matrix constructed below should consist of all 0s
maxconncomp_membership_vector_list <- lapply(1:length(distinct_gtex_tissues[c(1:8, 23, 9:22, 24:27)]), function(xx) {
  epimat <- clean_expr_mat_list_pcremoval[[xx]][which(rownames(clean_expr_mat_list_pcremoval[[xx]]) 
                                                      %in% rownames(module_membership_matrix)), , 
                                                drop = FALSE]
  maxconncompgenes <- getGenesInMaxConnComp(epimat, thresholds[xx])
  output_vector <- sapply(rownames(epimat), function(x) 
    ifelse(x %in% maxconncompgenes, 1, 0))
  output_vector
})

maxconncomp_membership_matrix <- t(do.call(rbind, maxconncomp_membership_vector_list))

maxconncomp_membership_matrix[, 5] <- 0
maxconncomp_membership_matrix[, 11] <- 0

for (i in 1:nrow(maxconncomp_membership_matrix)){
  maxconncomp_membership_matrix[i,][which(maxconncomp_membership_matrix[i,] 
                                          %in% c(0))] <- (-i)}

maxconncomp_membership_matrix_rowlist <- as.list(data.frame(t(maxconncomp_membership_matrix)))
distancelist_after_cormat_thresholding <- lapply(maxconncomp_membership_matrix_rowlist, 
                                                 function(xx) sapply(1:nrow(maxconncomp_membership_matrix), 
                                                                     function(x) 
                                                                       seq_dist(maxconncomp_membership_matrix_rowlist[[x]], xx, 
                                                                                method = "hamming")))

em_distancevector_after_cormat_thresholding <- as.numeric(lapply(distancelist_after_cormat_thresholding, 
                                                                 function(xx) length(which(xx < 19))))

rownames(maxconncomp_membership_matrix) <- rownames(module_membership_matrix)

length(which(em_distancevector_after_cormat_thresholding > 74))

length(which(rownames(maxconncomp_membership_matrix)[which(em_distancevector_after_cormat_thresholding > 74)]
             %in% rownames(adjmatrix)[1:74]))





###also do the following if visual inspection of the random graphs is required
plotRandomGraph <- function(expr_mat, number_of_genes, threshold){
  randsample <- sample(1:nrow(expr_mat), number_of_genes)
  corMat <- cor(t(expr_mat[randsample, , drop = FALSE]))
  round(corMat,2)
  diag(corMat) <- 0
  #round(rowMaxs(corMat), 2)
  #round(rowMins(corMat), 2)
  #corMat <- corMat[order(rownames(corMat)), order(colnames(corMat))]
  rownames(corMat) <- colnames(corMat) <- rownames(expr_mat)[randsample]
  adj <- abs(corMat)
  adj[adj > threshold] <- 1
  adj[adj < threshold] <- 0
  diag(adj) <- 0
  gA <- graphAM(adjMat = adj, edgemode = "undirected")
  gNEL <- as(gA, "graphNEL")
  maintitle <- paste0(threshold, "; ", "maximum_component_size = ", 
                      getCompSize(expr_mat[randsample, , drop = FALSE], threshold, "max"))
  plot(gA, "neato", main = maintitle)
}

getConnectivity <- function(expr_mat, threshold){
  corMat <- cor(t(expr_mat))
  round(corMat,2)
  diag(corMat) <- 0
  #round(rowMaxs(corMat), 2)
  #round(rowMins(corMat), 2)
  #corMat <- corMat[order(rownames(corMat)), order(colnames(corMat))]
  rownames(corMat) <- colnames(corMat) <- rownames(expr_mat)
  adj <- abs(corMat)
  adj[adj > threshold] <- 1
  adj[adj < threshold] <- 0
  diag(adj) <- 0
  gA <- graphAM(adjMat = adj, edgemode = "undirected")
  gNEL <- as(gA, "graphNEL")
  complengths <- lengths(connComp(gNEL))
  maxcompsize <- max(complengths)
  maxconncompgenes <- connComp(gNEL)[which(lengths(connComp(gNEL)) %in% maxcompsize)][[1]]
  mean(graph::degree(gNEL)[which(nodes(gNEL) %in% maxconncompgenes)])
}

getRandConnectivity <- function(expr_mat, number_of_genes, threshold){
  randsample <- sample(1:nrow(expr_mat), number_of_genes)
  getConnectivity(expr_mat[randsample, , drop = FALSE], threshold)
}



####plot the connectivity properties of the random graphs arising from the above choice of thresholds
#first size of maximum connected component
quartz(file = "em_vs_random_maxconncomponent_sizes.pdf", width = 8, height = 8, type = "pdf")
op <- par(mfrow = c(3, 3))
sapply(c(1:8), function(xx) {
  em_maxconncompsize <- getCompSize(clean_expr_mat_list_pcremoval[[xx]][which(rownames(clean_expr_mat_list_pcremoval[[xx]]) 
                                                                              %in% rownames(module_membership_matrix)), 
                                                                        , drop = FALSE], thresholds[xx], "max")
  plot(density(replicate(300, getRandCompSize(clean_expr_mat_list_onlyabove05RPKM[[xx]], 270, 
                                              thresholds[xx], "max"))), 
       xlab = "maximum connected\ncomponent size", ylab = "Density", cex.lab = 1.1, 
       lwd = 1.7, xlim = c(0, 180), ylim = c(0, 0.055), col = 4, xaxt = 'n', yaxt = 'n',
       main = distinct_gtex_tissues[xx], bty = 'l')
  abline(v = em_maxconncompsize, col = "orange", lwd = 1.7)
  axis(1, at = c(0, 70, 140), cex.axis = 1.1)
  axis(2, at = c(0, 0.050), cex.axis = 1.1)
  #legend <- legend("top", legend = c("EM genes", "Random genes\n(expression level = EM genes)"), 
  #                 col = c(2,4), lty = "solid", lwd = 2, cex = 0.45)
})
#plot(0.00011, 1, xlab = "", xaxt = 'n', yaxt = 'n', ylab = '', bty = 'n', pch = 19, col = rgb(1,1,1))
em_maxconncompsize <- getCompSize(clean_expr_mat_list_pcremoval[[10]][which(rownames(clean_expr_mat_list_pcremoval[[10]]) 
                                                                            %in% rownames(module_membership_matrix)), 
                                                                      , drop = FALSE], thresholds[10], "max")
plot(density(replicate(300, getRandCompSize(clean_expr_mat_list_onlyabove05RPKM[[10]], 270, 
                                            thresholds[10], "max"))), 
     xlab = "maximum connected\ncomponent size", ylab = "Density", cex.lab = 1.1, 
     lwd = 1.7, xlim = c(0, 180), ylim = c(0, 0.055), col = 4, xaxt = 'n', yaxt = 'n',
     main = distinct_gtex_tissues[9], bty = 'l')
abline(v = em_maxconncompsize, col = "orange", lwd = 1.7)
axis(1, at = c(0, 70, 140), cex.axis = 1.1)
axis(2, at = c(0, 0.050), cex.axis = 1.1)
legend <- legend("top", legend = c("EM genes", "Random genes\n(expression level = EM genes)"), 
                 col = c("orange", 4), lty = "solid", lwd = 2, bty = 'n')
par(op)
dev.off()



###now connectivity (i.e. average degree within the maximum connected component)
quartz(file = "em_and_random_maxconncompconnectivity.pdf", width = 8, height = 8, type = "pdf")
op <- par(mfrow = c(3, 3))
sapply(c(1:8), function(xx) {
  em_maxconncompconnectivity <- getConnectivity(clean_expr_mat_list_pcremoval[[xx]][which(rownames(clean_expr_mat_list_pcremoval[[xx]]) 
                                                                                          %in% rownames(module_membership_matrix)), 
                                                                                    , drop = FALSE], thresholds[xx])
  plot(density(replicate(300, getRandConnectivity(clean_expr_mat_list_onlyabove05RPKM[[xx]], 270, thresholds[xx]))), 
       xlab = "Average node degree", ylab = "Density", cex.lab = 1.1, lwd = 1.7, 
       xlim = c(0, 8), col = 4, ylim = c(0, 1.2), xaxt = 'n', yaxt = 'n',
       main = distinct_gtex_tissues[xx], bty = 'l')
  abline(v = em_maxconncompconnectivity, col = "orange", lwd = 1.7)
  axis(1, at = c(0, 4, 8), cex.axis = 1.1)
  axis(2, at = c(0, 1), cex.axis = 1.1)
  #legend <- legend("top", legend = c("EM genes", "Random genes\n(expression level = EM genes)"), 
  #                 col = c(2,4), lty = "solid", lwd = 2, cex = 0.45)
})
#plot(0.00011, 1, xlab = "", xaxt = 'n', yaxt = 'n', ylab = '', bty = 'n', pch = 19, col = rgb(1,1,1))
em_maxconncompconnectivity <- getConnectivity(clean_expr_mat_list_pcremoval[[10]][which(rownames(clean_expr_mat_list_pcremoval[[10]]) 
                                                                                        %in% rownames(module_membership_matrix)), 
                                                                                  , drop = FALSE], thresholds[10])
plot(density(replicate(300, getRandConnectivity(clean_expr_mat_list_onlyabove05RPKM[[10]], 270, thresholds[10]))), 
     xlab = "Average node degree", ylab = "Density", cex.lab = 1.1, lwd = 1.7, 
     xlim = c(0, 8), ylim = c(0, 1.2), col = 4, xaxt = 'n', yaxt = 'n',
     main = distinct_gtex_tissues[9], bty = 'l')
abline(v = em_maxconncompconnectivity, col = "orange", lwd = 1.7)
axis(1, at = c(0, 4, 8), cex.axis = 1.1)
axis(2, at = c(0, 1), cex.axis = 1.1)
legend <- legend("top", legend = c("EM genes", "Random genes\n(expression level = EM genes)"), 
                 col = c("orange", 4), lty = "solid", lwd = 2, bty = 'n')
par(op)
dev.off()




