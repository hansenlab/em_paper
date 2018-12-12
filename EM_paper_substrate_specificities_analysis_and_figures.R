#####EM genes with complementary activity


######Histone Methylation
H3K4_writers_hm <- c("KMT2A", "KMT2B", "KMT2C", "KMT2D", "SETD1A", "SETD1B", 
                     "SETD7", "SMYD1", "SMYD2", "ASH1L", "PRDM9")
H3K27_writers_hm <- c("EZH1", "EZH2")
H3K36_writers_hm <- c("NSD1", "WHSC1", "WHSC1L1", "SETD2", "SMYD2", "ASH1L", "SETD3", "SETMAR")
H3K9_writers_hm <- c("PRDM2", "EHMT1", "EHMT2", "SETDB1", "SUV39H1")

H3K4_erasers_hm <- c("KDM1A", "KDM1B", "KDM5A", "KDM5B", "KDM5C", "KDM5D", "NO66")
H3K27_erasers_hm <- c("KDM6A", "UTY", "KDM6B", "KDM7A", "PHF8")
H3K36_erasers_hm <- c("KDM2A", "KDM2B", "KDM4A", "KDM4B", "KDM4C", "KDM4D")
H3K9_erasers_hm <- c("KDM3A", "KDM3B", "JMJD1C", "KDM4A", "KDM4B", "KDM4C", "KDM4D", "PHF8", "PHF2")


######Histone Acetylation
H3K27_writers_hac <- c("EP300", "CREBBP")

H3K9_erasers_hac <- c("SIRT1", "SIRT2")


########
EM_specificities <- data.frame(gene_name = c(H3K4_writers_hm, H3K27_writers_hm, 
                                             H3K36_writers_hm, H3K9_writers_hm, H3K4_erasers_hm, 
                                             H3K27_erasers_hm, H3K36_erasers_hm, H3K9_erasers_hm, 
                                             H3K27_writers_hac, H3K9_erasers_hac), 
                               specificity = c(rep("H3K4_methylation_writer", length(H3K4_writers_hm)), 
                                               rep("H3K27_methylation_writer", length(H3K27_writers_hm)), 
                                               rep("H3K36_methylation_writer", length(H3K36_writers_hm)), 
                                               rep("H3K9_methylation_writer", length(H3K9_writers_hm)), 
                                               rep("H3K4_methylation_eraser", length(H3K4_erasers_hm)), 
                                               rep("H3K27_methylation_eraser", length(H3K27_erasers_hm)), 
                                               rep("H3K36_methylation_eraser", length(H3K36_erasers_hm)), 
                                               rep("H3K9_methylation_eraser", length(H3K9_erasers_hm)), 
                                               rep("H3K27_acetylation_writer", length(H3K27_writers_hac)), 
                                               rep("H3K9_acetylation_eraser", length(H3K9_erasers_hac))))

write_csv(EM_specificities, "EM_substrate_specificities.csv")
########


####
quartz(file = "em_pli_same_mark_enzymes.pdf", width = 4.5, height = 4.5, type = "pdf")
op <- par(mar = c(11,5,2,1) + 0.1)
plot(1, col = rgb(1,1,1), xlim = c(2, 20), ylab = "pLI", xlab = "", xaxt = 'n',
     yaxt = 'n', main = "", bty = 'u', ylim = c(-0.1,1.1))
polygon(c(0, 0, 22, 22), c(0.9, 1.2, 1.2, 0.9), 
        col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
plotpLIasPoints(H3K4_writers_hm, insert = TRUE, 2, rgb(1,0,0,0.8), EM = TRUE)
plotpLIasPoints(H3K4_erasers_hm, insert = TRUE, 4, rgb(1,0,0,0.8), EM = TRUE)
plotpLIasPoints(H3K27_writers_hm, insert = TRUE, 6, rgb(1,0,0,0.8), EM = TRUE)
plotpLIasPoints(H3K27_erasers_hm, insert = TRUE, 8, rgb(1,0,0,0.8), EM = TRUE)
plotpLIasPoints(H3K36_writers_hm, insert = TRUE, 10, rgb(1,0,0,0.8), EM = TRUE)
plotpLIasPoints(H3K36_erasers_hm, insert = TRUE, 12, rgb(1,0,0,0.8), EM = TRUE)
plotpLIasPoints(H3K9_writers_hm, insert = TRUE, 14, rgb(1,0,0,0.8), EM = TRUE)
plotpLIasPoints(H3K9_erasers_hm, insert = TRUE, 16, rgb(1,0,0,0.8), EM = TRUE)
plotpLIasPoints(H3K27_writers_hac, insert = TRUE, 18, rgb(1,0,0,0.8), EM = TRUE)
plotpLIasPoints(H3K9_erasers_hac, insert = TRUE, 20, rgb(1,0,0,0.8), EM = TRUE)
axis(1, at = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20), 
     labels = c("H3K4 meth writers", "H3K4 meth erasers", "H3K27 meth writers", "H3K27 meth erasers", 
                "H3K36 meth writers", "H3K36 meth erasers",
                "H3K9 meth writers", "H3K9 meth erasers",
                "H3K27 ac writers", "H3K9 ac erasers"),
     las = 2, cex.axis = 1.2)
axis(2, at = c(0.1, 0.9), cex.axis = 1.2)
abline(v=3, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=5, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=7, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=9, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=11, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=13, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=15, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=17, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
abline(v=19, lwd = 1.7, lty = "longdash", col = rgb(0, 0, 0, 0.4))
par(op)
dev.off()


###
same_mark_enzymes_list <- list(H3K4_writers_hm, H3K4_erasers_hm, H3K27_writers_hm, H3K27_erasers_hm, 
                               H3K36_writers_hm, H3K36_erasers_hm, H3K9_writers_hm, H3K9_erasers_hm, 
                               H3K27_writers_hac, H3K9_erasers_hac)

same_mark_enzymes_plis <- lapply(same_mark_enzymes_list, function(xx) 
  getpLIs(xx, EMinput = TRUE, name_to_use = "exac"))

same_mark_enzymes_plis_greaterThan0.9 <- lapply(same_mark_enzymes_plis, function(xx) 
  round(length(which(xx > 0.9))/length(xx), 2))




