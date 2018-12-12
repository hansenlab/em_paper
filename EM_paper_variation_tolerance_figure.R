quartz(file = "EMvsTFvsAllother_plis_allothergenes_plotted.pdf", width = 3.25, height = 3, type = "pdf")
plot(density(allother_plis, from = 0, to = 1), 
     col = 4, lwd = 1.7, ylab = "Density", 
     xlab = "pLI", 
     xaxt ='n', yaxt = 'n', cex.lab = 1.4, main = "", bty = 'l')
lines(density(em_plis, 
              from = 0, to = 1), 
      col = 2, lwd = 1.7)
lines(density(tf_plis, from = 0, to = 1), col = 3, lwd = 1.7)
axis(2, at = c(2, 4), cex.axis = 1.4)
axis(1, at = c(0.1, 0.5, 0.9), cex.axis = 1.4)
legend <- legend("top", legend = c("All other genes", "TF genes", "EM genes"), 
                 col = c(4, 3, 2), lty = "solid", cex = 0.7, bty = 'n', lwd = 1.4)
polygon(c(0.9, 0.9, 1.05, 1.05), c(0, 2.6, 2.6, 0), col = adjustcolor("gray80", alpha.f = 0.5), lty = 0)
dev.off()



quartz(file = "pli_categories_violin_plots_updated.pdf", width = 3.2, height = 3.2, type = "pdf")
op <- par(mar = c(7,4,0.75,0.25) + 0.1)
col_for_points <- 1
col_for_density <- rgb(1, 0, 0, 0.5)
plot(1, col = rgb(1,1,1), xlim = c(1.7, 12.3), ylab = "pLI", xlab = "", xaxt = 'n',
     yaxt = 'n', main = "", bty = 'l', ylim = c(0,1))
polygon(c(1.2, 1.2, 12.8, 12.8), c(0.9, 1.05, 1.05, 0.9), col = adjustcolor("gray70", alpha.f = 0.5), lty = 0)
vioplot(getpLIs(writers, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 2, wex = 1.5)
plotpLIasPoints(writers, insert = TRUE, 2, col_for_points, c(1.4, 12.5), EM = TRUE)

vioplot(getpLIs(erasers, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 4, wex = 1.5)
plotpLIasPoints(erasers, insert = TRUE, 4, col_for_points, c(1.4, 12.5), EM = TRUE)

plotpLIasPoints(remodelers, insert = TRUE, 6, col_for_points, c(1.4, 12.5), EM = TRUE)

vioplot(getpLIs(readers, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 8, wex = 1.5)
plotpLIasPoints(readers, insert = TRUE, 8, col_for_points, c(1.4, 12.5), EM = TRUE)

vioplot(getpLIs(dualEMfunction, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 10, wex = 1.5)
plotpLIasPoints(dualEMfunction, insert = TRUE, 10, col_for_points, c(1.4, 12.5), EM = TRUE)

vioplot(getpLIs(EMandTF, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 12, wex = 1.5)
plotpLIasPoints(EMandTF, insert = TRUE, 12, col_for_points, c(1.4, 12.5), EM = TRUE)

axis(1, at = c(2, 4, 6, 8, 10, 12), 
     labels = c("Writers", "Erasers", "Remodelers", "Readers", "Dual Function", "EM and TF"),
     las = 2)
axis(2, at = c(0, 0.5, 1))
abline(v=3, lwd = 1.7, lty = "longdash")
abline(v=5, lwd = 1.7, lty = "longdash")
abline(v=7, lwd = 1.7, lty = "longdash")
abline(v=9, lwd = 1.7, lty = "longdash")
abline(v=11, lwd = 1.7, lty = "longdash")
par(op)
dev.off()





quartz(file = "pli_subcategories_violin_plots_updated.pdf", width = 4.5, height = 3.75, type = "pdf")
op <- par(mar = c(11,5,2,1) + 0.1)
col_for_points <- 1
col_for_density <- rgb(1, 0, 0, 0.5)
plot(1, col = rgb(1,1,1), xlim = c(2, 22), ylab = "pLI", xlab = "", xaxt = 'n',
     yaxt = 'n', main = "", bty = 'l', ylim = c(0,1))
polygon(c(1.2, 1.2, 22.8, 22.8), c(0.9, 1.05, 1.05, 0.9), col = adjustcolor("gray70", alpha.f = 0.5), lty = 0)
vioplot(getpLIs(hmts, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 2, wex = 1.5)
plotpLIasPoints(hmts, insert = TRUE, 2, col_for_points, c(1.4, 22.5), EM = TRUE)


vioplot(getpLIs(hats, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 4, wex = 1.5)
plotpLIasPoints(hats, insert = TRUE, 4, col_for_points, EM = TRUE)


plotpLIasPoints(dnmts, insert = TRUE, 6, col_for_points, EM = TRUE)

vioplot(getpLIs(hdms, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 8, wex = 1.5)
plotpLIasPoints(hdms, insert = TRUE, 8, col_for_points, EM = TRUE)

vioplot(getpLIs(hdacs, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at =10, wex = 1.5)
plotpLIasPoints(hdacs, insert = TRUE, 10, col_for_points, EM = TRUE)


plotpLIasPoints(dnmerasers, insert = TRUE, 12, col_for_points, EM = TRUE)

vioplot(getpLIs(hmrsonly, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 14, wex = 1.5)
plotpLIasPoints(hmrsonly, insert = TRUE, 14, col_for_points, EM = TRUE)

vioplot(getpLIs(harsonly, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 16, wex = 1.5)
plotpLIasPoints(harsonly, insert = TRUE, 16, col_for_points, EM = TRUE)

plotpLIasPoints(dnmrsonly, insert = TRUE, 18, col_for_points, EM = TRUE)

plotpLIasPoints(dnumrsonly, insert = TRUE, 20, col_for_points, EM = TRUE)

vioplot(getpLIs(dnmrsonly, EMinput = TRUE, name_to_use = "exac"), 
        col = col_for_density, drawRect = FALSE, border = rgb(1,1,1), add = TRUE, at = 22, wex = 1.5)
plotpLIasPoints(dualreaders, insert = TRUE, 22, col_for_points, EM = TRUE)


axis(1, at = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22), 
     labels = c("Hist meth writers", "Hist ac writers", "DNAm writers", "Hist meth erasers", 
                "Hist ac erasers", "DNAm erasers",
                "Hist meth readers", "Hist ac readers",
                "DNAm readers", "Unmeth CpG readers", "Dual Readers"),
     las = 2, cex.axis = 1.2)
axis(2, at = c(0.1, 0.9), cex.axis = 1.2)
abline(v=3, lwd = 1.7, lty = "longdash")
abline(v=5, lwd = 1.7, lty = "longdash")
abline(v=7, lwd = 1.7, lty = "longdash")
abline(v=9, lwd = 1.7, lty = "longdash")
abline(v=11, lwd = 1.7, lty = "longdash")
abline(v=13, lwd = 1.7, lty = "longdash")
abline(v=15, lwd = 1.7, lty = "longdash")
abline(v=17, lwd = 1.7, lty = "longdash")
abline(v=19, lwd = 1.7, lty = "longdash")
abline(v=21, lwd = 1.7, lty = "longdash")
par(op)
dev.off()

