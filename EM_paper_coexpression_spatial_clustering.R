####assess chromosomal clustering of EM genes
getEpiChromDistances <- function(chromosome, genes, gene_relationship){
  ###gene_relationship has to do with whether we are dealing with immediate neighbors (gene_relationship = 2)
  ###or genes separated by one gene (gene_relationship = 3) etc
  chrom <- chromosome
  if (length(which(epiGenesDF$Chr %in% chrom & epiGenesDF$Gene_name_gtex_rpm %in% genes)) < gene_relationship){
    paste0("less than ", gene_relationship, " genes in that chromosome") 
  } else {
    start_positions <- epiGenesDF$Start_Pos[which(epiGenesDF$Chr %in% chrom & epiGenesDF$Gene_name_gtex_rpm %in% genes)]
    end_positions <- epiGenesDF$End_Pos[which(epiGenesDF$Chr %in% chrom & epiGenesDF$Gene_name_gtex_rpm %in% genes)]
    start_positions <- sort(start_positions)
    end_positions <- sort(end_positions)
    distances <- sapply(gene_relationship:length(start_positions), function(xx) 
      start_positions[xx] - end_positions[xx-(gene_relationship-1)])
    sorted_distances <- sort(distances)
    sorted_distances
  }
}

all_chr_distances <- sapply(sapply(c(1:22, "X", "Y"), function(x) paste0("chr", x)), function(xx) 
  getEpiChromDistances(xx, rownames(adjmatrix)[1:74], gene_relationship = 2))

coexpressed_chrs <- sapply(sapply(c(1:22, "X", "Y"), function(x) paste0("chr", x)), function(xx) 
  length(which(epiGenesDF$Chr %in% xx & epiGenesDF$Gene_name_gtex_rpm %in% rownames(adjmatrix)[1:74])))
coexpressed_chrs <- as.numeric(coexpressed_chrs)
chisq.test(coexpressed_chrs/74, p = rep(1/24, 24), simulate.p.value = TRUE)

quartz(file = "em_chr_distances_less_than_1Mb_orange.pdf", width = 6.5, height = 3, type = "pdf")
plot(rep(1, length(all_chr_distances[[1]])), all_chr_distances[[1]], xlab = "Chromosome", 
     ylab = "Coexpressed EM distances (Mb)", cex.lab = 1.2, yaxt = 'n', xaxt = 'n', main = "", pch = 19, 
     ylim = c(0, 10^6), xlim = c(0.8, 48.5), cex = 0.35, col = "orange")
for (xx in 2:24){
  if (!(typeof(all_chr_distances[[xx]]) %in% c("character"))){
    points(rep(xx+(xx-1), length(all_chr_distances[[xx]])), all_chr_distances[[xx]], pch = 19, 
           cex = 0.35, col = "orange")
  }
  axis(1, at = seq(1, 48, by = 2), labels = c(1:22, "X", "Y"), cex.axis = 0.62)
  axis(2, at = c(10^5, 5*10^5, 9*10^5), labels = c(0.1, 0.5, 0.9), cex.axis = 0.8)
  #axis(2, at = c(0.1*10^7, 8*10^7, 16.5*10^7), labels = c(1000, 80000, 165000))
}
dev.off()

quartz(file = "em_chr_distances_more_than_1Mb_orange.pdf", width = 6.5, height = 3, type = "pdf")
plot(rep(1, length(all_chr_distances[[1]])), all_chr_distances[[1]], xlab = "Chromosome", 
     ylab = "Coexpressed EM distances (Mb)", cex.lab = 1.2, yaxt = 'n', xaxt = 'n', main = "", pch = 19, 
     ylim = c(10^6, 8*10^7), xlim = c(0.8, 48.5), cex = 0.35, col = "orange")
for (xx in 2:24){
  if (!(typeof(all_chr_distances[[xx]]) %in% c("character"))){
    points(rep(xx+(xx-1), length(all_chr_distances[[xx]])), all_chr_distances[[xx]], pch = 19, 
           cex = 0.35, col = "orange")
  }
  axis(1, at = seq(1, 48, by = 2), labels = c(1:22, "X", "Y"), cex.axis = 0.62)
  axis(2, at = c(0.1*10^7, 4*10^7, 8*10^7), labels = c(1, 40, 80))
}
dev.off()

all_chr_distances <- sapply(sapply(c(1:22, "X", "Y"), function(x) paste0("chr", x)), function(xx) 
  getEpiChromDistances(xx, rownames(adjmatrix)[75:157], gene_relationship = 2))

coexpressed_chrs <- sapply(sapply(c(1:22, "X", "Y"), function(x) paste0("chr", x)), function(xx) 
  length(which(epiGenesDF$Chr %in% xx & epiGenesDF$Gene_name_gtex_rpm %in% rownames(adjmatrix)[75:157])))
coexpressed_chrs <- as.numeric(coexpressed_chrs)
chisq.test(coexpressed_chrs/83, p = rep(1/24, 24), simulate.p.value = TRUE)

quartz(file = "em_chr_distances_less_than_1Mb_lightskyblue.pdf", width = 6.5, height = 3, type = "pdf")
plot(rep(1, length(all_chr_distances[[1]])), all_chr_distances[[1]], xlab = "Chromosome", 
     ylab = "Coexpressed EM distances (Mb)", cex.lab = 1.2, yaxt = 'n', xaxt = 'n', main = "", pch = 19, 
     ylim = c(0, 10^6), xlim = c(0.8, 48.5), cex = 0.35, col = "lightskyblue3")
for (xx in 2:24){
  if (!(typeof(all_chr_distances[[xx]]) %in% c("character"))){
    points(rep(xx+(xx-1), length(all_chr_distances[[xx]])), all_chr_distances[[xx]], pch = 19, 
           cex = 0.35, col = "lightskyblue3")
  }
  axis(1, at = seq(1, 48, by = 2), labels = c(1:22, "X", "Y"), cex.axis = 0.62)
  axis(2, at = c(10^5, 5*10^5, 9*10^5), labels = c(0.1, 0.5, 0.9), cex.axis = 0.8)
  #axis(2, at = c(0.1*10^7, 8*10^7, 16.5*10^7), labels = c(1000, 80000, 165000))
}
dev.off()

quartz(file = "em_chr_distances_more_than_1Mb_lightskyblue.pdf", width = 6.5, height = 3, type = "pdf")
plot(rep(1, length(all_chr_distances[[1]])), all_chr_distances[[1]], xlab = "Chromosome", 
     ylab = "Coexpressed EM distances (Mb)", cex.lab = 1.2, yaxt = 'n', xaxt = 'n', main = "", pch = 19, 
     ylim = c(10^6, 16.5*10^7), xlim = c(0.8, 48.5), cex = 0.35, col = "lightskyblue3")
for (xx in 2:24){
  if (!(typeof(all_chr_distances[[xx]]) %in% c("character"))){
    points(rep(xx+(xx-1), length(all_chr_distances[[xx]])), all_chr_distances[[xx]], pch = 19, 
           cex = 0.35, col = "lightskyblue3")
  }
  axis(1, at = seq(1, 48, by = 2), labels = c(1:22, "X", "Y"), cex.axis = 0.62)
  axis(2, at = c(0.1*10^7, 8*10^7, 16.5*10^7), labels = c(1, 80, 165))
}
dev.off()