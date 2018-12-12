######First, use recount to get a list of genes expressed in this cell line:
abstract_search('K562')


download_study("SRP010061")

load(file.path('SRP010061/rse_gene.Rdata'))

expr_matrix <- getRPKM(read_counts(rse_gene))

median_expr <- rowMedians(expr_matrix)

rownames(expr_matrix) <- gsub("[.].*$", "", rownames(expr_matrix))

#symbols <- sapply(rownames(expr_matrix), function(xx) gene_info$SYMBOL[which(gene_info$ENSEMBL == xx)])

gencode <- gsub('\\..*', '', names(recount_genes))

gene_info <- select(org.Hs.eg.db, gencode, c('ENTREZID', 'GENENAME', 'SYMBOL',
                                             'ENSEMBL'), 'ENSEMBL')

expressed_genes <- gene_info$SYMBOL[which(gene_info$ENSEMBL 
                                          %in% rownames(expr_matrix)[which(median_expr > 0)])]

expressed_genes <- expressed_genes[-which(is.na(expressed_genes))]

expressed_genes <- unique(expressed_genes)


#######Now get granges objects for the TSS of EM genes and all other (minus ribo) genes
em_hg19_coords <- as.data.frame(read.table("em_hg19_coords.bed",header = FALSE, 
                                           sep="\t", stringsAsFactors=FALSE, quote=""))


names_of_converted_in_order <- epiGenesDF$Gene_name[-c(121, 188)]
strand_of_converted_in_order <- epiGenesDF$Strand[-c(121, 188)]

colnames(em_hg19_coords) <- c("chr", "start", "end")

#add 1 to the start because ucsc works with 0 based coordinates
em_hg19_coords$start <- em_hg19_coords$start + 1

em_hg19_coords$gene_name <- names_of_converted_in_order
em_hg19_coords$strand <- strand_of_converted_in_order

em_plusminus5kb_tss <- em_hg19_coords


em_plusminus5kb_tss$start <- sapply(1:293, function(xx) {
  if (em_plusminus5kb_tss$strand[xx] %in% c("+")){
    minus <- em_hg19_coords$start[xx] - 5000
  } else if(em_hg19_coords$strand[xx] %in% c("-")) {
    minus <- em_hg19_coords$end[xx] - 5000
  }
  minus
})

em_plusminus5kb_tss$end <- sapply(1:293, function(xx) {
  if (em_plusminus5kb_tss$strand[xx] %in% c("+")){
    plus <- em_hg19_coords$start[xx] + 5000
  } else if (em_plusminus5kb_tss$strand[xx] %in% c("-")){
    plus <- em_hg19_coords$end[xx] + 5000
  }
  plus
})


em_tss_granges <- GRanges(seqnames = em_plusminus5kb_tss$chr, 
                          IRanges(start = em_plusminus5kb_tss$start, 
                                  width = (em_plusminus5kb_tss$end - em_plusminus5kb_tss$start+1)), 
                          gene_name = em_plusminus5kb_tss$gene_name)

#ribo genes serve as positive controls because they are also strongly co-expressed
ribo_coords <- read_csv("ribo_coords.csv")
ribo_coords <- ribo_coords[-which(duplicated(ribo_coords$geneSymbol)), , drop = FALSE]
ribo_coords <- ribo_coords[-which(!(ribo_coords$geneSymbol %in% ribo_genes)), , drop = FALSE]

colnames(ribo_coords) <- c("chr", "strand", "start", "end", "gene_name")

ribo_hg19_coords <- ribo_coords

ribo_hg19_coords$start <- ribo_hg19_coords$start + 1

ribo_plusminus5kb_tss <- ribo_hg19_coords


ribo_plusminus5kb_tss$start <- sapply(1:76, function(xx) {
  if (ribo_plusminus5kb_tss$strand[xx] %in% c("+")){
    minus <- ribo_hg19_coords$start[xx] - 5000
  } else if(ribo_hg19_coords$strand[xx] %in% c("-")) {
    minus <- ribo_hg19_coords$end[xx] - 5000
  }
  minus
})

ribo_plusminus5kb_tss$end <- sapply(1:76, function(xx) {
  if (ribo_plusminus5kb_tss$strand[xx] %in% c("+")){
    plus <- ribo_hg19_coords$start[xx] + 5000
  } else if (ribo_plusminus5kb_tss$strand[xx] %in% c("-")){
    plus <- ribo_hg19_coords$end[xx] + 5000
  }
  plus
})

ribo_tss_granges <- GRanges(seqnames = ribo_plusminus5kb_tss$chr, 
                            IRanges(start = ribo_plusminus5kb_tss$start, 
                                    width = (ribo_plusminus5kb_tss$end - ribo_plusminus5kb_tss$start+1)), 
                            gene_name = ribo_plusminus5kb_tss$gene_name)


all_coords <- read_csv("all_coords.csv")
all_coords <- all_coords[-which(duplicated(all_coords$geneSymbol)), , drop = FALSE]

colnames(all_coords) <- c("chr", "strand", "start", "end", "gene_name")

all_hg19_coords <- all_coords

all_hg19_coords$start <- all_hg19_coords$start + 1

all_plusminus5kb_tss <- all_hg19_coords


all_plusminus5kb_tss$start <- sapply(1:length(all_coords$gene_name), function(xx) {
  if (all_plusminus5kb_tss$strand[xx] %in% c("+")){
    minus <- all_hg19_coords$start[xx] - 5000
  } else if(all_hg19_coords$strand[xx] %in% c("-")) {
    minus <- all_hg19_coords$end[xx] - 5000
  }
  minus
})

all_plusminus5kb_tss$end <- sapply(1:length(all_coords$gene_name), function(xx) {
  if (all_plusminus5kb_tss$strand[xx] %in% c("+")){
    plus <- all_hg19_coords$start[xx] + 5000
  } else if (all_plusminus5kb_tss$strand[xx] %in% c("-")){
    plus <- all_hg19_coords$end[xx] + 5000
  }
  plus
})

all_tss_granges <- GRanges(seqnames = all_plusminus5kb_tss$chr, 
                           IRanges(start = all_plusminus5kb_tss$start, 
                                   width = (all_plusminus5kb_tss$end - all_plusminus5kb_tss$start+1)), 
                           gene_name = all_plusminus5kb_tss$gene_name)


all_tss_granges <- all_tss_granges[-which(all_tss_granges$gene_name %in% c(em_tss_granges$gene_name, 
                                                                           ribo_tss_granges$gene_name))]




getGeneTFPairsMatrix <- function(tss_granges, chip_granges_list, removeZeros = c(TRUE, FALSE)){
  gene_tf_pairs <- lapply(chip_granges_list, function(xx) {
    overlaps <- countOverlaps(tss_granges, xx)
    em_tf_vector <- rep(0, length(tss_granges$gene_name))
    em_tf_vector[which(overlaps > 0)] <- 1
    em_tf_vector
  })
  genetfpairs_matrix <- do.call(cbind, gene_tf_pairs)
  rownames(genetfpairs_matrix) <- tss_granges$gene_name
  if (removeZeros == TRUE){
    for (i in 1:nrow(genetfpairs_matrix)){
      genetfpairs_matrix[i,][which(genetfpairs_matrix[i,] ==0 )] <- (-i)}
  }
  genetfpairs_matrix
}


expr_em_tss_granges <- em_tss_granges[which(em_tss_granges$gene_name %in% expressed_genes)]
expr_all_tss_granges <- all_tss_granges[which(all_tss_granges$gene_name %in% expressed_genes)]


em_tfs_k562 <- getGeneTFPairsMatrix(expr_em_tss_granges, final_k562_granges_list, removeZeros = FALSE)
colnames(em_tfs_k562) <- tfs_k562
all_tfs_k562 <- getGeneTFPairsMatrix(expr_all_tss_granges, final_k562_granges_list, removeZeros = FALSE)
colnames(all_tfs_k562) <- tfs_k562




#########function that calculates the enrichment of a given TF at the promoters of one category of genes
#########versus another
getTFinGeneCategoryEnrichment <- function(gene_tf_pairs_matrix, category_1, category_2){
  category1_matrix <- gene_tf_pairs_matrix[which(rownames(gene_tf_pairs_matrix) 
                                                 %in% category_1), , drop = FALSE]
  category2_matrix <- gene_tf_pairs_matrix[which(rownames(gene_tf_pairs_matrix) 
                                                 %in% category_2), , drop = FALSE]
  enrichment_list <- lapply(1:dim(gene_tf_pairs_matrix)[2], function(xx) {
    n11 <-  apply(category1_matrix, 2, function(x) length(which(x == 1)))[xx]
    n12 <-  apply(category2_matrix, 2, function(x) length(which(x == 1)))[xx]
    n21 <- apply(category1_matrix, 2, function(x) length(which(x == 0)))[xx]
    n22 <- apply(category2_matrix, 2, function(x) length(which(x == 0)))[xx]
    ft <- fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))
    c(ft$p.value, ft$estimate, colnames(gene_tf_pairs_matrix)[xx], 
      round(n11/dim(category1_matrix)[1], 2), round(n12/dim(category2_matrix)[1], 2))
  })
  df <- do.call(rbind.data.frame, enrichment_list)
  colnames(df) <- c("p_val", "OR", "TF", "percentage_1", "percentage_2")
  df$p_val <- as.numeric(df$p_val)
  #df$adjusted_p_val <- p.adjust(df$p_val, method = "bonferroni")
  df$OR <- as.numeric(df$OR)
  df$percentage_1 <- as.numeric(df$percentage_1)
  df$percentage_2 <- as.numeric(df$percentage_2)
  df
}

em_column_zeros_k562 <- apply(em_tfs_k562, 2, function(xx) length(which(xx == 0))) #exclude TFs that don't bind anything
em_tfs_k562_no_column_zeros <- em_tfs_k562[, -which(em_column_zeros_k562 > (dim(em_tfs_k562)[1]-11)), 
                                           drop = FALSE]

em_tfs_k562_no_column_zeros <- em_tfs_k562_no_column_zeros[which(rownames(em_tfs_k562_no_column_zeros) %in% rownames(adjmatrix)), , 
                                                           drop = FALSE] #restrict to the rownames that are in the original adjmatrix


highly_coexp_vs_not_coexp <- getTFinGeneCategoryEnrichment(em_tfs_k562_no_column_zeros, 
                                                           rownames(adjmatrix)[1:74], rownames(adjmatrix)[158:270])


######get two kinds of null distributions. One by shuffling the labels of highly co-expressed and 
######not co-expressed EM genes, and one by comparing to random genes

#label shuffling
shuffled_enrichments <- replicate(2000, {
  sample1 <- sample(1:242, 166)
  sample2 <- sample1[1:72]
  sample3 <- sample1[73:166]
  
  getTFinGeneCategoryEnrichment(em_tfs_k562_no_column_zeros, 
                                k562_expr_rownames_adjmat[sample2], 
                                k562_expr_rownames_adjmat[sample3])
})

shuffled_pvals <- lapply(seq(1, 10000, by = 5), function(xx) shuffled_enrichments[[xx]])
shuffled_ors <- lapply(seq(2, 10000, by = 5), function(xx) shuffled_enrichments[[xx]])


#random genes
all_column_zeros_k562 <- apply(all_tfs_k562, 2, function(xx) length(which(xx == 0))) #exclude TFs that don't bind anything
all_tfs_k562_no_column_zeros <- all_tfs_k562[, -which(all_column_zeros_k562 > (dim(all_tfs_k562)[1]-51)), 
                                             drop = FALSE]

indices_to_sample1 <- which(rownames(all_tfs_k562_no_column_zeros) 
                            %in% symbols[which(median_expr > 1.2)])
indices_to_sample2 <- which(rownames(all_tfs_k562_no_column_zeros) 
                            %in% symbols[which(median_expr > 0.4)])

random_enrichments <- replicate(2000, {
  sample1 <- sample(indices_to_sample1, 72)
  sample2 <- sample(indices_to_sample2, 94)
  
  getTFinGeneCategoryEnrichment(all_tfs_k562_no_column_zeros, 
                                rownames(all_tfs_k562_no_column_zeros)[sample1], 
                                rownames(all_tfs_k562_no_column_zeros)[sample2])
})


random_pvals <- lapply(seq(1, 10000, by = 5), function(xx) random_enrichments[[xx]])
random_ors <- lapply(seq(2, 10000, by = 5), function(xx) random_enrichments[[xx]])



quartz(file = "EM_k562_random_distribution.pdf", width = 3, height = 3, type = "pdf")
plot(density(sapply(1:2000, function(xx) 
  length(which(random_pvals[[xx]] < 0.05 & log2(random_ors[[xx]]) > 1))), from = 0), 
  main = "K562", xlab = "OR > 2 & pval < 0.05", bty = 'l', lwd = 2, 
  yaxt = 'n', xaxt = 'n', cex.lab = 1)
axis(1, at = c(10, 50, 90), cex.axis = 1)
axis(2, at = c(0, 0.07), cex.axis = 1)

abline(v = length(which(highly_coexp_vs_not_coexp$p_val < 0.05 & log2(highly_coexp_vs_not_coexp$OR) > 1)), 
       col = "dark orange", lwd = 2)
legend <- legend("topright", legend = c("random (2000 samples)", "EM (observed)"), 
                 col = c(1, "dark orange"), bty = 'n', lty = "solid", lwd = 1.5, cex = 0.4)


dev.off()

quartz(file = "EM_k562_shuffled_distribution.pdf", width = 3, height = 3, type = "pdf")
plot(density(sapply(1:2000, function(xx) 
  length(which(shuffled_pvals[[xx]] < 0.05 & log2(shuffled_ors[[xx]]) > 1))), from = 0), 
  main = "K562", xlab = "OR > 2 & pval < 0.05", bty = 'l', lwd = 2, 
  yaxt = 'n', xaxt = 'n', cex.lab = 1)
axis(1, at = c(10, 50), cex.axis = 1)
axis(2, at = c(0, 0.10), cex.axis = 1)

abline(v = length(which(highly_coexp_vs_not_coexp$p_val < 0.05 & log2(highly_coexp_vs_not_coexp$OR) > 1)), 
       col = "dark orange", lwd = 2)
legend <- legend("top", legend = c("EM shuffled (2000 times)", "EM (observed)"), 
                 col = c(1, "dark orange"), bty = 'n', lty = "solid", lwd = 1.5, cex = 0.4)


dev.off()


