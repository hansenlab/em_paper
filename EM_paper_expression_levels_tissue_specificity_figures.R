#####expression levels in different tissues figure
quartz(file = "EMvsTFvsAllother_expression_boxplot.pdf", width = 8, height = 3.7, type = "pdf")
#op <- par(mar = c(4,4,4,2) + 0.1)
op <- par(mar = c(5.5, 5.5, 2, 2) + 0.1)
boxplot(all_other_genes_medians_uncorrected[, order(colMedians(as.matrix(epi_medians_uncorrected)), decreasing = TRUE)], 
        las = 2, ylab = "Median RPKM", yaxt = 'n', 
        xlim=c(0, 86), ylim = c(0, 15), outline = FALSE, boxfill=rgb(1, 1, 1, alpha=1), 
        border=rgb(1, 1, 1, alpha=1), at = seq(1, 4*28, by = 4), xaxt ='n', notch = FALSE, 
        medlwd = 0, xlab = "Tissues 1-28", cex.axis = 1.5, cex.lab = 1.8, frame = FALSE)
axis(2, at = c(2, 4), cex.axis = 1.4)
#axis(1, at = seq(1, 4*28, by = 4), labels = colnames(epi_medians), las = 2)
#axis(1, at = seq(1, 4*28, by = 4), labels = rep("", 28))

boxplot(tf_medians_uncorrected[, order(colMedians(as.matrix(epi_medians_uncorrected)), decreasing = TRUE)], xlab = "", 
        outline = FALSE, add = TRUE, at = seq(1, 28, by = 1), 
        yaxt = 'n', xaxt = 'n', col = "3", notch = FALSE, whisklty = 0, staplelty = 0, medlty = 0, boxlty = 0)
boxplot(all_other_genes_medians_uncorrected[, order(colMedians(as.matrix(epi_medians_uncorrected)), decreasing = TRUE)], xlab = "", 
        outline = FALSE, add = TRUE, at = seq(30, 57, by = 1), 
        yaxt = 'n', xaxt = 'n', col = "4", notch = FALSE, whisklty = 0, staplelty = 0, medlty = 0, boxlty = 0)
boxplot(epi_medians_uncorrected[, order(colMedians(as.matrix(epi_medians_uncorrected)), decreasing = TRUE)], xlab = "", 
        outline = FALSE, add = TRUE, at = seq(59, 86, by = 1), 
        yaxt = 'n', xaxt = 'n', col = "2", notch = FALSE, whisklty = 0, staplelty = 0, medlty = 0, boxlty = 0)
legend <- legend("topleft", c("EM", "TF", "All other"), 
                 fill = c(2, 3, 4), horiz = FALSE, cex = 0.8)
par(op)
dev.off()



###make the figures for tissue specificity
EM_spec_scores_uncorrected <- EM_spec_scores_uncorrected[, c(1:23, 25:28, 24), drop = FALSE]

quartz(file = "EM_specificities_all_tissues.pdf", width = 4.5, height = 4, type = "pdf")
scores.vec <- as.vector(EM_spec_scores_uncorrected)
ind.vec <- as.vector(sapply(1:28, function(xx) rep(xx, length(em_distance_from_uniform_uncorrected))))
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
quartz(file = "EMvsTFvsAllother_specificity.pdf", width = 3.25, height = 3, type = "pdf")
plot(density(em_distance_from_uniform_uncorrected, from = 0, to = 1), col = 2, lwd = 1.7, ylab = "Density", 
     xlab = "Tissue Specificity Score", 
     xaxt ='n', yaxt = 'n', cex.lab = 1.4, main = "", ylim = c(0, 14), bty = 'l')
lines(density(tf_distance_from_uniform_uncorrected, from = 0, to = 1), col = 3, lwd = 1.7)
lines(density(allother_distance_from_uniform_uncorrected, from = 0, to = 1), col = 4, lwd = 1.7)
axis(2, at = c(2, 8, 14), cex.axis = 1.4)
axis(1, at = c(0.1, 0.5, 0.9), cex.axis = 1.4)
legend <- legend("top", legend = c("All other genes", "TF genes", "EM genes"), 
                 col = c(4, 3, 2), lty = "solid", cex = 0.7, bty = 'n', lwd = 1.4)
#polygon(c(0.9, 0.9, 1.05, 1.05), c(0, 1.9, 1.9, 0), col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
dev.off()

#
quartz(file = "EMvsTCAcyclegenes_specificity.pdf", width = 3.25, height = 3, type = "pdf")
plot(density(em_distance_from_uniform_uncorrected, from = 0, to = 1), col = 2, lwd = 1.7, ylab = "Density", 
     xlab = "Tissue Specificity Score", 
     xaxt ='n', yaxt = 'n', cex.lab = 1.4, main = "", ylim = c(0, 14), bty = 'l')
lines(density(citric_distance_from_uniform_uncorrected, from = 0, to = 1), col = 1, lwd = 1.7)
axis(2, at = c(2, 8, 14), cex.axis = 1.4)
axis(1, at = c(0.1, 0.5, 0.9), cex.axis = 1.4)
legend <- legend("top", legend = c("TCA cycle genes", "EM genes"), 
                 col = c(1, 2), lty = "solid", cex = 0.7, bty = 'n', lwd = 1.4)
#polygon(c(0.9, 0.9, 1.05, 1.05), c(0, 1.9, 1.9, 0), col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
dev.off()


#######pLIs of EMvsTFvsAllother only for genes with low tissue specificity
em_plis_low_spec <- getpLIs(em_new_rownames[which(em_distance_from_uniform_uncorrected < 0.1)], EMinput = TRUE, name_to_use = "gtex_rpkm")
tf_plis_low_spec <- getpLIs(tf_new_rownames[which(tf_distance_from_uniform_uncorrected < 0.1)], EMinput = FALSE)
allother_plis_low_spec <- getpLIs(rownames(all_other_genes_medians_uncorrected)[which(allother_distance_from_uniform_uncorrected < 0.1)], EMinput = FALSE)


quartz(file = "EMvsTFvsAllother_plis_onlylowspec_allothergenesplotted.pdf", width = 3.25, height = 3, type = "pdf")
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






#########plot expression levels of high vs low pLI EM genes
#first define this function
compareExprAccordingtopLI <- function(gene_names, main_label){
  input_genes <- gtex.cds.df$gene_names[which(gtex.cds.df$gene_names %in% gene_names)]
  input_genes <- input_genes[which(input_genes %in% exactab$gene)]
  medians <- as.matrix(sapply(gtex_tissues, getMedianExprLevels, gene_names = input_genes))
  rownames(medians) <- gtex.cds.df$gene_names[which(gtex.cds.df$gene_names %in% input_genes)]
  pLIs <- getpLIs(rownames(medians), EMinput = EM, name_to_use = "gtex_rpkm")
  high_pli_medians <- medians[which(pLIs > 0.9), , drop=FALSE]
  high_pli_medians_of_medians <- rowMedians(high_pli_medians)
  plot(density(high_pli_medians_of_medians), col = rgb(1,0,0,1/2), 
       yaxt='n', xlab = "Expression", ylab ="", main = main_label)
  polygon(density(high_pli_medians_of_medians), col = rgb(1,0,0,1/2))
  
  low_pli_medians <- medians[-which(pLIs > 0.9), , drop=FALSE]
  low_pli_medians_of_medians <- rowMedians(low_pli_medians)
  lines(density(low_pli_medians_of_medians), col = 2)
  polygon(density(low_pli_medians_of_medians), col = rgb(0,0,1,1/2))
  
  legend <- legend("topright", legend = c("pLI > 0.9", "pLI < 0.9"), col=c(2,4), lty="solid")
}

#make the figure
compareExprAccordingtopLI(epiGenes, "EM genes")

