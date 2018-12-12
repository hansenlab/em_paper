quartz(file = "collective_domain_constraint_pli.pdf", width = 3.75, height = 4.25, type = "pdf")
op <- par(mar = c(7, 6.2, 1, 1) + 0.1)
plot(rep(2, length(em_domain_constraint_vector_high_pli)) + rnorm(length(em_domain_constraint_vector_high_pli), 0, 0.2), 
     em_domain_constraint_vector_high_pli + rnorm(length(em_domain_constraint_vector_high_pli), 0, 0.05), xlim = c(1.4, 3.8), 
     pch = 19, col = "grey50", ylim = c(-1, 100), cex = 0.7, xaxt = 'n', xlab = "", 
     ylab = "EM-specific\nDomain Constraint", yaxt = 'n', cex.lab = 1.8)
points(rep(3.2, length(em_domain_constraint_vector_low_pli)) + rnorm(length(em_domain_constraint_vector_low_pli), 0, 0.2), 
       em_domain_constraint_vector_low_pli + rnorm(length(em_domain_constraint_vector_low_pli), 0, 0.05), 
       pch = 19, col = 1, cex = 0.7)
axis(2, at = c(10, 50, 90), labels = c(10, 50, 90), cex.axis = 1.8)
axis(1, at = c(2, 3.2), labels = c("High pLI\ngenes", "Low pLI\ngenes"), cex.axis = 1.8, las = 2)
par(op)

dev.off()


quartz(file = "em_vs_non_em_domain_constraint.pdf", width = 3, height = 4.25, type = "pdf")
op <- par(mar = c(9, 7, 1, 1) + 0.1)
boxplot(diffs_high_pli, at = 2, xlim = c(1.8, 2.7), col = "grey50", boxlty = 0, whisklty = 1, staplelwd = 0, outpch = 19, outcex = 0.4, 
        ylab = "Mean domain constraint\ndifference (EM-specific\nvs non-EM-specific)", yaxt = 'n', xaxt = 'n', cex.lab = 1.4, frame = TRUE)
boxplot(diffs_low_pli, at = 2.5, add = TRUE, col = 1, boxlty = 0, whisklty = 1, staplelwd = 0, outpch = 19, outcex = 0.4, 
        yaxt = 'n', frame = TRUE)
#axis(2, at = c(-50, -25, 0, 25, 50), labels = c(-50, -25, 0, 25, 50), cex.axis = 1.4, las = 1)
axis(2, at = c(-30, 0, 30), labels = c(-30, 0, 30), cex.axis = 1.4, las = 1)
axis(1, at = c(2, 2.5), labels = c("High pLI genes", "Low pLI genes"), cex.axis = 1.4, las = 2)
par(op)

dev.off()


####run this before the next figure
constraint_range_for_em_domains_appearing_multiple_times <- lapply(em_domains_appearing_multiple_times, function(xx) {
  mat <- em_domain_constraint_matrix_high_pli[grep(xx, rownames(em_domain_constraint_matrix_high_pli)), ]
  output_for_domain <- sapply(1:ncol(mat), function(x) {
    if (length(which(!(mat[, x] %in% c(190)))) > 1){
      max(mat[, x][which(!(mat[, x] %in% c(190)))]) - min(mat[, x][which(!(mat[, x] %in% c(190)))])
    }
  })
  unlist(output_for_domain)
})


quartz(file = "domain_constraint_range_when_multiple_copies.pdf", width = 4.5, height = 4, type = "pdf")
op <- par(mar = c(12, 5.9, 0.5, 0.5) + 0.1)
plot(rep(1, length(constraint_range_for_em_domains_appearing_multiple_times[[1]])), 
     constraint_range_for_em_domains_appearing_multiple_times[[1]], pch = 19, cex = 0.3, 
     ylim = c(0, 100), xlim = c(1, 9), 
     xlab = "", 
     ylab = "Constraint\nRange", xaxt = 'n', yaxt = 'n', col = 1, main = "", cex.lab = 1.4)
sapply(c(2:4, 6:10), function(xx) {
  if (xx < 5){
    points(rep(xx, length(constraint_range_for_em_domains_appearing_multiple_times[[xx]])), 
           constraint_range_for_em_domains_appearing_multiple_times[[xx]], pch = 19, cex = 0.3)
  } else {
    points(rep(xx-1, length(constraint_range_for_em_domains_appearing_multiple_times[[xx]])), 
           constraint_range_for_em_domains_appearing_multiple_times[[xx]], pch = 19, cex = 0.3)
  }
}) 
axis(1, at = 1:9, labels = c(em_domains_appearing_multiple_times[-5]), las = 2, cex.axis = 1.2)
axis(2, at = c(10, 50, 90), cex.axis = 1.4)
par(op)

dev.off()


