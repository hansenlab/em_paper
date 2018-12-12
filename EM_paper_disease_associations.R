####number of EM genes with different disease associations
only_neuro_associated <- epiGenesDF$Gene_name[which(epiGenesDF$is.neuro.associated == TRUE &
                                                      epiGenesDF$is.ca.associated == FALSE)]
only_drivers <- epiGenesDF$Gene_name[which(epiGenesDF$is.ca.associated == TRUE & 
                                             epiGenesDF$is.MDEM == FALSE &
                                             epiGenesDF$is.neuro.associated == FALSE)]
neuro_and_drivers <- epiGenesDF$Gene_name[which(epiGenesDF$is.ca.associated == TRUE & 
                                                  epiGenesDF$is.neuro.associated == TRUE)]

mdem_only <- epiGenesDF$Gene_name[which(epiGenesDF$is.MDEM == TRUE 
                                        & epiGenesDF$is.neuro.associated == FALSE
                                        & epiGenesDF$is.ca.associated == FALSE)]

no_disease <- epiGenesDF$Gene_name[which(epiGenesDF$is.neuro.associated == FALSE
                                         & epiGenesDF$is.ca.associated == FALSE
                                         & epiGenesDF$is.MDEM == FALSE)]



quartz(file = "disease_association_mdem_new.pdf",width = 3, height = 4, type = "pdf")
op <- par(mar = c(9.3, 4, 1, 0) + 0.1)
barplot(c(100*round(length(mdem_only)/295, 2), 100*round(length(only_neuro_associated)/295, 2), 
          100*round(length(neuro_and_drivers)/295, 2), 100*round(length(only_drivers)/295, 2), 
          100*round(length(no_disease)/295, 2)), 
        names.arg = c("MDEM w/o Neuro", "Neuro", "Neuro & Ca", "Ca only", 
                      "No Disease"), 
        ylim = c(0, 80), yaxt = 'n', ylab = "Percentage",
        col = rep(rgb(1,0,0,7/9), 5), las = 2, 
        cex.names = 1.25, xlim = c(0, 1.5), width = c(0.2, 0.2, 0.2, 0.2, 0.2), cex.lab = 1.25)
axis(2, at = c(0, 25, 50, 75), cex.axis = 1.1)
par(op)
dev.off()




######disease-associated vs non disease-associated variation-tolerance
disease_associated <- epiGenesDF$Gene_name_exac[which(epiGenesDF$Gene_name %in% 
                                                        unique(c(mdem_only, only_neuro_associated, 
                                                                 neuro_and_drivers, only_drivers, ajhg_neuro)))]
non_disease_associated <- epiGenesDF$Gene_name_exac[which(epiGenesDF$Gene_name %in% 
                                                            unique(c(no_disease)))]


disease_associated_plis <- getpLIs(disease_associated, EMinput = TRUE, name_to_use = "exac")
non_disease_associated_plis <- getpLIs(non_disease_associated, EMinput = TRUE, name_to_use = "exac")

quartz(file = "em_disease_vs_nondisease_pli.pdf",width = 3.25, height = 3, type = "pdf")
plot(density(disease_associated_plis, from = 0, to = 1), col = "grey50", lwd = 1.7, ylab = "Density", xlab = "pLI", 
     xaxt ='n', yaxt = 'n', cex.lab = 1.4, main = "", bty = 'l')
lines(density(non_disease_associated_plis, from = 0, to = 1), col = 1, lwd = 1.7)
axis(2, at = c(1.5, 3), labels = c(1.5, 3), cex.axis = 1.4, las = 1)
axis(1, at = c(0.1, 0.5, 0.9), cex.axis = 1.4)
legend <- legend("top", legend = c("disease associated EM", "not disease associated EM"), 
                 col = c("grey50", 1), lty = "solid", cex = 0.55, bty = 'n', lwd = 1.4)
polygon(c(0.9, 0.9, 1.05, 1.05), c(0, 3.5, 3.5, 0), col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
dev.off()




