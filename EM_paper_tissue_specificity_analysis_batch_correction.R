load(file = "data/gtex.cds.rpkm.rda")
load(file = "data/gtex.cds.df.rda")
gtex.cds.mat <- gtex.cds.rpkm



options(stringsAsFactors = FALSE)
gtex_samples <- read.delim("data/GTEx_Data_V6_Annotations_SampleAttributesDS.txt")
#for some reason, the annotation file contains sample IDs that are not in the gtex matrix column names. remove those sample ids
notinmat <- which(!(gtex_samples$SAMPID %in% colnames(gtex.cds.rpkm)))
sampids <- gtex_samples$SAMPID[-notinmat]             


##get rins
rins <- gtex_samples$SMRIN[-notinmat]
#rin_mat <- vector()
#for (i in 1:8555){rin_mat[i] <- rins[which(sampids %in% c(colnames(gtex.cds.rpkm)[i]))]} #this is basically instantaneous, but probably still not optimal since it's a for loop
rin_mat <- sapply(1:8555, function(xx) 
  rins[which(sampids %in% c(colnames(gtex.cds.rpkm)[xx]))])

##get tissue ids
tissueids <- gtex_samples$SMTSD[-notinmat]
tissue_mat <- vector()
#for (i in 1:8555){tissue_mat[i] <- tissueids[which(sampids %in% c(colnames(gtex.cds.rpkm)[i]))]} #this is basically instantaneous, but probably still not optimal since it's a for loop
tissue_mat <- sapply(1:8555, function(xx) 
  tissueids[which(sampids %in% c(colnames(gtex.cds.rpkm)[xx]))])

#keep only the 28 tissues I'll use
tissue_ids_to_use <- tissue_mat[which(tissue_mat %in% distinct_gtex_tissues)]
rins_to_use <- rin_mat[which(tissue_mat %in% distinct_gtex_tissues)]

gtex.cds.rpkm.to.use <- gtex.cds.rpkm[, which(tissue_mat %in% distinct_gtex_tissues), drop = FALSE]

tissue_id_rin_df <- data.frame(tissue = tissue_ids_to_use, rin = rins_to_use)

mod <- model.matrix(~as.factor(tissue) + rin, 
                    data = tissue_id_rin_df)
mod0 <- model.matrix(~rin, 
                     data = tissue_id_rin_df)

number.of.zeroes <- apply(gtex.cds.rpkm.to.use, 1, function(xx) length(which(xx == 0)))

gtex.cds.rpkm.to.use.no.all.zeroes <- gtex.cds.rpkm.to.use[-which(number_of_zeroes == 5421), 
                                                           , drop = FALSE]


sva.tissue.spec.2 <- sva(gtex.cds.rpkm.to.use.no.all.zeroes, mod, mod0)

save(sva.tissue.spec.2, file = "/users/lboukas/sva.tissue.spec.2.rda")

mod.with.svs <- cbind(mod0, sva.tissue.spec.2$sv)

corrected.gtex.cds.rpkm <- removeBatchEffect(gtex.cds.rpkm.to.use.no.all.zeroes, 
                                             covariates = mod.with.svs[, -1, drop = FALSE], 
                                             design = model.matrix(~as.factor(tissue), 
                                                                   data = tissue_id_rin_df))



uncorrected.gtex.cds.mat <- gtex.cds.rpkm.to.use.no.all.zeroes

gtex.cds.df <- gtex.cds.df[-which(number.of.zeroes == 5421), , drop = FALSE]

neg_values <- which(gtex.cds.mat < 0)
original_value_in_neg_values <- uncorrected.gtex.cds.mat[neg_values]
summary(original_value_in_neg_values)
plot(density(original_value_in_neg_values), xlim = c(0,1))

medians_of_rows_with_neg_values <- rowMedians(gtex.cds.mat[rows_with_neg_values, , drop = FALSE])


###Set the negative values to 0

gtex.cds.mat[which(gtex.cds.mat < 0)] <- 0

###Remove all uniformly low counts
uniformly_low_counts <- apply(gtex.cds.mat, 1, function(xx) length(which(xx > 0.01)))
uniformly_low_counts <- which(uniformly_low_counts == 0)
gtex.cds.mat <- gtex.cds.mat[-uniformly_low_counts, ,drop = FALSE]

gtex.cds.df <- gtex.cds.df[-uniformly_low_counts, , drop = FALSE]

tissue_mat <- tissue_mat[which(tissue_mat %in% distinct_gtex_tissues)]

############and now proceed as with the uncorrected data
epi_medians <- sapply(distinct_gtex_tissues, getMedianExprLevels, gene_names = input_epig)
epi_medians <- as.data.frame(epi_medians)
tf_medians <- sapply(distinct_gtex_tissues, getMedianExprLevels, gene_names = input_tf)
tf_medians <- as.data.frame(tf_medians)
all_other_genes_medians <- sapply(distinct_gtex_tissues, getMedianExprLevels, gene_names = input_allother)
all_other_genes_medians <- as.data.frame(all_other_genes_medians)

####make the figures using the corrected version of expression levels
quartz(file = "EM_specificities_all_tissues_batch_corrected.pdf", width = 4.5, height = 4, type = "pdf")
scores.vec <- as.vector(EM_spec_scores)
ind.vec <- as.vector(sapply(1:28, function(xx) rep(xx, length(em_distance_from_uniform))))
#op <- par(mar = c(15,4,4,2) + 0.1)
plot(ind.vec, scores.vec, pch = 19, col = rgb(1, 0, 0, 1/2), 
     xlab = "Other Tissues", ylab = "Tissue Specificity Score", yaxt = 'n', xaxt = 'n', ylim = c(0,1), 
     main = "", cex = 0.5, cex.lab = 1.4, bty = 'l')
axis(1, at = 28, labels = "Testis", las = 2, cex.axis = 1.4)
axis(2, at = c(0.2, 0.8), cex.axis = 1.4)
#polygon(c(27.8, 27.8, 28.2, 28.2), c(0, 1, 1, 0), col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
#par(op)
dev.off()

#
quartz(file = "EMvsTCAcyclegenes_specificity_batch_corrected.pdf", width = 3.25, height = 3, type = "pdf")
plot(density(em_distance_from_uniform, from = 0, to = 1), col = 2, lwd = 1.7, ylab = "Density", 
     xlab = "Tissue Specificity Score", 
     xaxt ='n', yaxt = 'n', cex.lab = 1.4, main = "", ylim = c(0, 14), bty = 'l')
lines(density(citric_distance_from_uniform, from = 0, to = 1), col = 1, lwd = 1.7)
axis(2, at = c(2, 8, 14), cex.axis = 1.4)
axis(1, at = c(0.1, 0.5, 0.9), cex.axis = 1.4)
legend <- legend("top", legend = c("TCA cycle genes", "EM genes"), 
                 col = c(1, 2), lty = "solid", cex = 0.7, bty = 'n', lwd = 1.4)
#polygon(c(0.9, 0.9, 1.05, 1.05), c(0, 1.9, 1.9, 0), col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
dev.off()




#
quartz(file = "EMvsTFvsAllother_specificity_batch_corrected.pdf", width = 3.25, height = 3, type = "pdf")
plot(density(em_distance_from_uniform, from = 0, to = 1), col = 2, lwd = 1.7, ylab = "Density", 
     xlab = "Tissue Specificity Score", 
     xaxt ='n', yaxt = 'n', cex.lab = 1.4, main = "", ylim = c(0, 14), bty = 'l')
lines(density(tf_distance_from_uniform, from = 0, to = 1), col = 3, lwd = 1.7)
lines(density(allother_distance_from_uniform, from = 0, to = 1), col = 4, lwd = 1.7)
axis(2, at = c(2, 8, 14), cex.axis = 1.4)
axis(1, at = c(0.1, 0.5, 0.9), cex.axis = 1.4)
legend <- legend("top", legend = c("All other genes", "TF genes", "EM genes"), 
                 col = c(4, 3, 2), lty = "solid", cex = 0.7, bty = 'n', lwd = 1.4)
#polygon(c(0.9, 0.9, 1.05, 1.05), c(0, 1.9, 1.9, 0), col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
dev.off()


#####pLI EMvsTFvsAllOther (only genes with specificity score less than 0.1 after correction for batch effects)
em_plis_low_spec <- getpLIs(rownames(epi_medians)[which(em_distance_from_uniform < 0.1)], EMinput = TRUE, name_to_use = "gtex_rpkm")
tf_plis_low_spec <- getpLIs(rownames(tf_medians)[which(tf_distance_from_uniform < 0.1)], EMinput = FALSE)
allother_plis_low_spec <- getpLIs(rownames(all_other_genes_medians)[which(allother_distance_from_uniform < 0.1)], EMinput = FALSE)

quartz(file = "EMvsTFvsAllother_plis_onlylowspec_batch_corrected.pdf", width = 3.25, height = 3, type = "pdf")
plot(density(allother_plis_low_spec, from = 0, to = 1), 
     col = 4, lwd = 1.7, ylab = "Density", 
     xlab = "pLI", 
     xaxt ='n', yaxt = 'n', cex.lab = 1.4, main = "", bty = 'l')
lines(density(em_plis_low_spec, from = 0, to = 1), 
      col = 2, lwd = 1.7)
lines(density(tf_plis_low_spec, from = 0, to = 1), col = 3, lwd = 1.7)
axis(2, at = c(1, 2), cex.axis = 1.4)
axis(1, at = c(0.1, 0.5, 0.9), cex.axis = 1.4)
legend <- legend("top", legend = c("All other genes", "TF genes", "EM genes"), 
                 col = c(4, 3, 2), lty = "solid", cex = 0.7, bty = 'n', lwd = 1.4)
polygon(c(0.9, 0.9, 1.05, 1.05), c(0, 2.2, 2.2, 0), col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
dev.off()





