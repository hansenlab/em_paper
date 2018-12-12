load(file = "em_adjmatrix.rda")
load(file = "epiGenesDF.rda")

getEnrichmentInHighlyCoexpressed <- function(em_category, coexpressed, total){
  n11 <- length(which(em_category %in% coexpressed))
  n12 <- length(which(em_category %in% total)) - n11
  
  n21 <- length(coexpressed) - n11
  n22 <- length(total) - length(coexpressed) - n12
  
  fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))
}


getOR <- function(cutoff, category){
  n11 <- length(which(epiGenesDF$Gene_name_gtex[which(category 
                                                          %in% TRUE)] %in% rownames(adjmatrix)[1:cutoff]))
  n12 <- length(which(epiGenesDF$Gene_name_gtex[which(category
                                                          %in% TRUE)] %in% rownames(adjmatrix)[158:270]))
  
  n21 <- length(rownames(adjmatrix)[1:cutoff]) - n11
  n22 <- length(rownames(adjmatrix)[158:270]) - n12
  
  fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))
}


###
n_and_c <- epiGenesDF$Gene_name_gtex[which(epiGenesDF$is.neuro.associated == TRUE 
                                               & epiGenesDF$is.ca.associated == TRUE)]
c_only <- epiGenesDF$Gene_name_gtex[which(epiGenesDF$is.neuro.associated == FALSE & 
                                                epiGenesDF$is.ca.associated == TRUE)]
all_no_n_and_c <- rownames(adjmatrix)[-which(rownames(adjmatrix) %in% n_and_c)]
n_only <- epiGenesDF$Gene_name_gtex[which(epiGenesDF$is.neuro.associated == TRUE & 
                                                epiGenesDF$is.ca.associated == FALSE)]
coexpressed_no_n_and_c <- rownames(adjmatrix)[1:74][-which(rownames(adjmatrix)[1:74] %in% n_and_c)]

any_disease <- epiGenesDF$Gene_name_gtex[which(epiGenesDF$is.neuro.associated == TRUE 
                             | epiGenesDF$is.ca.associated == TRUE | epiGenesDF$is.MDEM == TRUE)]


all_plis_in_network <- sapply(rownames(adjmatrix), function(xx) 
  epiGenesDF$exac.pLI[which(epiGenesDF$Gene_name_gtex_rpm == xx)])

high_plis_in_network <- rownames(adjmatrix)[which(all_plis_in_network > 0.9)]

highly_coexpressed_high_plis <- high_plis_in_network[which(high_plis_in_network 
                                                           %in% rownames(adjmatrix)[1:74])]

neuro_associated_high_plis <- high_plis_in_network[which(high_plis_in_network 
                               %in% epiGenesDF$Gene_name_gtex[which(epiGenesDF$is.neuro.associated == TRUE)])]


ci1 <- getEnrichmentInHighlyCoexpressed(n_only, coexpressed_no_n_and_c, all_no_n_and_c)
ci2 <- getEnrichmentInHighlyCoexpressed(c_only, coexpressed_no_n_and_c, all_no_n_and_c)
ci3 <- getEnrichmentInHighlyCoexpressed(neuro_associated_high_plis, 
                                        highly_coexpressed_high_plis, high_plis_in_network) ###this is only considering enrichment of neuro genes with high pLI
ci4 <- getEnrichmentInHighlyCoexpressed(any_disease, rownames(adjmatrix)[1:74], rownames(adjmatrix)) 
ci7 <- getEnrichmentInHighlyCoexpressed(c(n_and_c, n_only), rownames(adjmatrix)[1:74], rownames(adjmatrix))###this is considering enrichment of all neuro without excluding ca
ci8 <- getEnrichmentInHighlyCoexpressed(c(n_and_c, c_only), rownames(adjmatrix)[1:74], rownames(adjmatrix))###this is considering enrichment of all ca without excluding ca
ci9 <- getEnrichmentInHighlyCoexpressed(n_and_c, rownames(adjmatrix)[1:74], rownames(adjmatrix))###this is considering enrichment of neuro and ca

quartz(file = "or_coexpression/ORs_neuro_drivers_coexpresion.pdf", width = 3.25, height = 3, type = "pdf")
op <- par(mar = c(6.5, 4.8, 1, 2) + 0.1)

plot(3, log2(ci4$estimate), pch = 19, cex = 1.2, xaxt = 'n', yaxt = 'n', xlab = "",
     ylab = "log2(OR) for\nenrichment", ylim = c(-1, 4.7), xlim = c(3, 15), bty = 'l', main = "")
segments(3, log2(ci4$conf.int[1]), 3, log2(ci4$conf.int[2]), lwd = 2.5)

points(5, log2(ci7$estimate), pch = 19, cex = 1.2)
segments(5, log2(ci7$conf.int[1]), 5, log2(ci7$conf.int[2]), lwd = 2.5)

points(7, log2(ci8$estimate), pch = 19, cex = 1.2)
segments(7, log2(ci8$conf.int[1]), 7, log2(ci8$conf.int[2]), lwd = 2.5)

points(9, log2(ci1$estimate), pch = 19, cex = 1.2)
segments(9, log2(ci1$conf.int[1]), 9, log2(ci1$conf.int[2]), lwd = 2.5)

points(11, log2(ci2$estimate), pch = 19, cex = 1.2)
segments(11, log2(ci2$conf.int[1]), 11, log2(ci2$conf.int[2]), lwd = 2.5)

points(13, log2(ci9$estimate), pch = 19, cex = 1.2)
segments(13, log2(ci9$conf.int[1]), 13, log2(ci9$conf.int[2]), lwd = 2.5)

points(15, log2(ci3$estimate), pch = 19, cex = 1.2)
segments(15, log2(ci3$conf.int[1]), 15, log2(ci3$conf.int[2]), lwd = 2.5)

abline(h = 0, lty = "longdash", lwd = 1.7)

abline(h = 1, lty = "longdash", lwd = 1.7)

axis(1, at = c(3, 5, 7, 9, 11, 13, 15), labels = c("Any disease", "Neuro", "Ca", 
                                                   "Neuro (no ca)", "Ca (no neuro)", "Neuro and Ca",
                                                   "high pLI Neuro\nvs other high pLI"), cex.axis = 0.8, las = 2)
axis(2, at = c(0, 1, 2, 3), cex.axis = 0.8, las = 1)
par(op)
dev.off()



######now perform a sensitivity analysis where we relax the cutoff for the definition of the
######highly co-expressed group and we look at what happens with the strength of the enrichment


df_neuro_vs_gray <- data.frame(OR = rep("", length(20:157)), ci_upper = rep("", length(20:157)), 
                               ci_lower = rep("", length(20:157)))
df_neuro_vs_gray$OR <- sapply(20:157, function(xx) 
  ors <- getOR(xx, epiGenesDF$is.neuro.associated)$estimate)
df_neuro_vs_gray$ci_lower <- sapply(20:157, function(xx) 
  ors <- getOR(xx, epiGenesDF$is.neuro.associated)$conf.int[1])
df_neuro_vs_gray$ci_upper <- sapply(20:157, function(xx) 
  ors <- getOR(xx, epiGenesDF$is.neuro.associated)$conf.int[2])


quartz(file = "or_coexpression/highly_coexpressed_vs_gray_various_cutoffs_log.pdf", width = 2.75, height = 3.5, type = "pdf")
op <- par(mar = c(5.5, 5.5, 1, 2) + 0.1)
plot(20:157, log2(df_neuro_vs_gray$OR), pch = 19, cex = 0.3, col = 1, xlab = "Highly Coexpressed\nSize", 
     ylab = "log2(OR) for Neuro\nenrichment", xaxt = 'n', main = "", cex.lab = 1, yaxt = 'n',
     xlim = c(20, 157), type = 'l', lwd = 2, bty = 'l', ylim = c(0, 5.5))
polygon(c(20:157, rev(20:157)), 
        c(log2(df_neuro_vs_gray$ci_lower), rev(log2(df_neuro_vs_gray$ci_upper))), 
        col = rgb(1,0,0,4/9))
axis(1, at = c(20, 140), cex.axis = 1)
axis(2, at = c(0, 1, 3, 5), cex.axis = 1)
abline(h = 0, lty = "longdash", lwd = 1.7)

abline(h = 1, lty = "longdash", lwd = 1.7)
par(op)
dev.off()





###compute enrichment of dual function in those associated with both neuro disease and ca
bothassoc <- epiGenesDF$Gene_name[which(epiGenesDF$is.ca.associated == TRUE
                                        & epiGenesDF$is.neuro.associated == TRUE)]

dual_func <- epiGenesDF$Gene_name[which(epiGenesDF$Epi_function %in% 
                                          c("Writer;Reader", "Remodeler;Reader", "Eraser;Reader"))]



n11 <- length(which(dual_func %in% bothassoc))
n12 <- length(dual_func) - n11

n21 <- length(bothassoc) - n11
n22 <- 295 - length(bothassoc) - n12

fisher.test(matrix(c(n11, n12, n21, n22), nrow = 2))

