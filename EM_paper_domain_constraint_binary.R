#####binary interpretation of CCR domain constraint. Specifically, if a domain has greater than 10% of bases above
#####the 90th CCR percentile, we call that domain constrained. If not, we call it "unconstrained"
getDomainPercentilesforGene <- function(gene_name, domains_vector,
                                        domains_appearing_multiple_times, domain_max_multiplicities){
  df <- domain_constraint[which(domain_constraint$gene_name %in% c(gene_name)), ]
  all_domains_constraint_vector_for_gene <- sapply(domains_vector, function(xx) {
    if (length(which(df$domain_name %in% c(xx))) %in% c(0)){
      if(xx %in% domains_appearing_multiple_times){
        rep(190, domain_max_multiplicities[which(domains_vector %in% c(xx))])
      } else {190}
    } else {
      df <- df[which(df$domain_name %in% c(xx)), ]
      method_1 <- grep("domain_1", df$domain_number)
      method_2 <- grep("domain_2", df$domain_number)
      method_3 <- grep("domain_3", df$domain_number)
      method_4 <- grep("domain_4", df$domain_number)
      method_0 <- grep("domain_0", df$domain_number)
      
      if (!(length(method_1) %in% c(0))){
        output_vector_for_domain <- sapply(unique(df$domain_number[method_1]), function(x) {
          df <- df[which(df$domain_number %in% c(x)), ]
          ###here I decide on a way to represent the mutational constraint of a domain. One option is the weighted sum of weighted percentiles, another option
          ###are quantiles. The rationale behind quantiles is that even if a small portion of a domain is highly constrained, that argues in favor of the 
          ###domain being functionally important
          weighted_constraint <- sum(df$weighted_pct*df$bp_overlap)/sum(df$bp_overlap)
          quantile_based_constraint <- quantile(df$weighted_pct, 0.95)
          #quantile_based_constraint
          #max(df$weighted_pct)
          df_above_90th_pct <- df[which(df$weighted_pct > 90), ]
          exact_percentile <- 100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
          ifelse(exact_percentile>10, 1, 0)
          #weighted_constraint
        })
        if (xx %in% domains_appearing_multiple_times){
          output_vector_for_domain <- c(output_vector_for_domain, 
                                        rep(190, domain_max_multiplicities[which(domains_vector %in% c(xx))] - length(output_vector_for_domain)))
        }
        output_vector_for_domain
      } else {
        if (!(length(method_2) %in% c(0))){
          output_vector_for_domain <- sapply(unique(df$domain_number[method_2]), function(x) {
            df <- df[which(df$domain_number %in% c(x)), ]
            ###here I decide on a way to represent the mutational constraint of a domain. One option is the weighted sum of weighted percentiles, another option
            ###are quantiles. The rationale behind quantiles is that even if a small portion of a domain is highly constrained, that argues in favor of the 
            ###domain being functionally important
            weighted_constraint <- sum(df$weighted_pct*df$bp_overlap)/sum(df$bp_overlap)
            quantile_based_constraint <- quantile(df$weighted_pct, 0.95)
            #quantile_based_constraint
            #max(df$weighted_pct)
            df_above_90th_pct <- df[which(df$weighted_pct > 90), ]
            exact_percentile <- 100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
            #weighted_constraint
            ifelse(exact_percentile>10, 1, 0)
          })
          if (xx %in% domains_appearing_multiple_times){
            output_vector_for_domain <- c(output_vector_for_domain, 
                                          rep(190, domain_max_multiplicities[which(domains_vector %in% c(xx))] - length(output_vector_for_domain)))
          }
          output_vector_for_domain
        } else {
          if (!(length(method_3) %in% c(0))){
            output_vector_for_domain <- sapply(unique(df$domain_number[method_3]), function(x) {
              df <- df[which(df$domain_number %in% c(x)), ]
              ###here I decide on a way to represent the mutational constraint of a domain. One option is the weighted sum of weighted percentiles, another option
              ###are quantiles. The rationale behind quantiles is that even if a small portion of a domain is highly constrained, that argues in favor of the 
              ###domain being functionally important. Another option is percent of bases having constraint above the 90th percentile.
              weighted_constraint <- sum(df$weighted_pct*df$bp_overlap)/sum(df$bp_overlap)
              quantile_based_constraint <- quantile(df$weighted_pct, 0.95)
              #quantile_based_constraint
              #max(df$weighted_pct)
              df_above_90th_pct <- df[which(df$weighted_pct > 90), ]
              exact_percentile <- 100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
              #weighted_constraint
              ifelse(exact_percentile>10, 1, 0)
            })
            if (xx %in% domains_appearing_multiple_times){
              output_vector_for_domain <- c(output_vector_for_domain, 
                                            rep(190, domain_max_multiplicities[which(domains_vector %in% c(xx))] - length(output_vector_for_domain)))
            }
            output_vector_for_domain
          } else {
            if (!(length(method_4) %in% c(0))){
              output_vector_for_domain <- sapply(unique(df$domain_number[method_4]), function(x) {
                df <- df[which(df$domain_number %in% c(x)), ]
                ###here I decide on a way to represent the mutational constraint of a domain. One option is the weighted sum of weighted percentiles, another option
                ###are quantiles. The rationale behind quantiles is that even if a small portion of a domain is highly constrained, that argues in favor of the 
                ###domain being functionally important
                weighted_constraint <- sum(df$weighted_pct*df$bp_overlap)/sum(df$bp_overlap)
                quantile_based_constraint <- quantile(df$weighted_pct, 0.95)
                #quantile_based_constraint
                #max(df$weighted_pct)
                df_above_90th_pct <- df[which(df$weighted_pct > 90), ]
                exact_percentile <- 100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
                #weighted_constraint
                ifelse(exact_percentile>10, 1, 0)
              })
              if (xx %in% domains_appearing_multiple_times){
                output_vector_for_domain <- c(output_vector_for_domain, 
                                              rep(190, domain_max_multiplicities[which(domains_vector %in% c(xx))] - length(output_vector_for_domain)))
              }
              output_vector_for_domain
            } else {
              output_vector_for_domain <- sapply(unique(df$domain_number[method_0]), function(x) {
                df <- df[which(df$domain_number %in% c(x)), ]
                ###here I decide on a way to represent the mutational constraint of a domain. One option is the weighted sum of weighted percentiles, another option
                ###are quantiles. The rationale behind quantiles is that even if a small portion of a domain is highly constrained, that argues in favor of the 
                ###domain being functionally important
                weighted_constraint <- sum(df$weighted_pct*df$bp_overlap)/sum(df$bp_overlap)
                quantile_based_constraint <- quantile(df$weighted_pct, 0.95)
                #quantile_based_constraint
                #max(df$weighted_pct)
                df_above_90th_pct <- df[which(df$weighted_pct > 90), ]
                exact_percentile <- 100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
                #weighted_constraint
                ifelse(exact_percentile>10, 1, 0)
              })
              if (xx %in% domains_appearing_multiple_times){
                output_vector_for_domain <- c(output_vector_for_domain, 
                                              rep(190, domain_max_multiplicities[which(domains_vector %in% c(xx))] - length(output_vector_for_domain)))
              }
              output_vector_for_domain
            } 
          }
        }
      }
    }
  })
  unlist(all_domains_constraint_vector_for_gene)
}

###########Figures
quartz(file = "collective_domain_constraint_binary.pdf", width = 3, height = 3, type = "pdf")
op <- par(mar = c(4, 4.8, 1, 0) + 0.1)
barplot(c(length(which(em_domain_constraint_vector_high_pli == 1)), 
          length(which(em_domain_constraint_vector_high_pli == 0)),
          0), 
        names.arg = c("High pLI\ngenes", "Low pLI\ngenes", ""), 
        ylim = c(0, 300), yaxt = 'n', ylab = "Total Number of\nEM-specific Domains",
        col = c(rep("white", 3)), las = 1, 
        cex.names = 0.75, xlim = c(0, 1.5), width = c(0.33,0.691,0.4), cex.lab = 0.8, 
        border = "white")
barplot(c(length(which(em_domain_constraint_vector_high_pli == 1)), 
          length(which(em_domain_constraint_vector_high_pli == 0)),
          0,
          length(which(em_domain_constraint_vector_low_pli == 1)),
          length(which(em_domain_constraint_vector_low_pli == 0))), 
        names.arg = c("", "", "", 
                      "", ""), 
        ylim = c(0, 300), yaxt = 'n', ylab = "",
        col = c(rgb(1,0,0,0.7), "royalblue1", "white", rgb(1,0,0,0.7), "royalblue1"), las = 1, 
        cex.names = 1.25, xlim = c(0, 1.5), width = c(0.2, 0.2, 0.1, 0.2, 0.2), cex.lab = 2, 
        border = "white", add = TRUE)
axis(2, at = c(0, 150, 300), cex.axis = 0.7)
legend <- legend("topright", legend = c("contrained", "unconstrained"), 
                 fill = c(rgb(1,0,0,0.7), "royalblue1"), border = "white", bty = 'n', cex = 0.7)
par(op)

dev.off()



############
number_of_constrained_em_domains <- apply(em_domain_constraint_list_high_pli, 2, 
                                          function(xx) length(which(xx == 1)))

number_of_constrained_other_domains <- apply(other_domain_constraint_list_high_pli, 2, 
                                             function(xx) length(which(xx == 1)))

diff_in_number_of_constrained_domains <- number_of_constrained_em_domains - number_of_constrained_other_domains

sorted_diffs <- sort(diff_in_number_of_constrained_domains)



quartz(file = "em_vs_other_domain_constraint_binary.pdf", width = 3.35, height = 3, type = "pdf")
op <- par(mar = c(4, 4.8, 0.5, 0.5) + 0.1)
plot(1:max(which(sorted_diffs < 0)), 
     sorted_diffs[which(sorted_diffs < 0)], 
     col = "royalblue1", pch = 19, cex = 0.3, xlim = c(0,145), ylim = c(-6,10), bty = 'l',
     xaxt = 'n', xlab = "High pLI EM genes", ylab = "Difference in number\nof constrained domains")
abline(h=0, lty = "longdash",  col = rgb(0,0,0,0.5))
points((max(which(sorted_diffs < 0))+1):(max(which(sorted_diffs == 0))), 
       sorted_diffs[which(sorted_diffs == 0)], 
       col = 1, pch = 19, cex = 0.3)
points((max(which(sorted_diffs == 0))+1):(max(which(sorted_diffs > 0))), 
       sorted_diffs[which(sorted_diffs > 0)], 
       col = 2, pch = 19, cex = 0.3, xlim = c(0,145), ylim = c(-6,10))

legend <- legend("topleft", legend = c("more EM domains constrained", "equal number", 
                                       "more non-EM domains constrained"), pch = 19, 
                 col = c(2, 1, "royalblue1"), bty = 'n', cex = 0.62)
axis(1, at = c(1,70,140))
abline(h=0, lty = "longdash", col = rgb(0,0,0,0.3))
par(op)

dev.off()


quartz(file = "em_vs_other_domain_constraint_binary_2.pdf", width = 3, height = 3, type = "pdf")
op <- par(mar = c(3, 4.8, 1, 0) + 0.1)
barplot(c(length(which(em_domain_constraint_vector_high_pli == 1)), 
          length(which(em_domain_constraint_vector_high_pli == 0)),
          0), 
        names.arg = c("EM-specific\ndomains", "Other\nDomains", ""), 
        ylim = c(0, 100), yaxt = 'n', ylab = "Percentage of high\npLI genes",
        col = c(rep("white", 3)), las = 1, 
        cex.names = 0.75, xlim = c(0, 1.5), width = c(0.33,0.691,0.4), cex.lab = 0.8, 
        border = "white")
barplot(c(100*round(nonzero_em/144, 2), 
          100*round(zero_em/144, 2),
          0,
          100*round(nonzero_other/144, 2),
          100*round(zero_other/144, 2)), 
        names.arg = c("", "", "", 
                      "", ""), 
        ylim = c(0, 300), yaxt = 'n', ylab = "",
        col = c(rgb(1,0,0,0.7), "royalblue1", "white", rgb(1,0,0,0.7), "royalblue1"), las = 1, 
        cex.names = 1.25, xlim = c(0, 1.5), width = c(0.2, 0.2, 0.1, 0.2, 0.2), cex.lab = 2, 
        border = "white", add = TRUE)
legend <- legend("topright", legend = c("at least 1 domain\nconstrained", "", "no domains\nconstrained"), 
                 fill = c(rgb(1,0,0,0.7), "white", "royalblue1"), border = "white", bty = 'n', cex = 0.7)
axis(2, at = c(0, 50, 100), cex.axis = 0.7)
par(op)

dev.off()






########PRDM only
set_prdm <- em_domain_constraint_list_high_pli[1,grep("PRDM", colnames(em_domain_constraint_matrix_high_pli))]

c2h2_prdm <- other_domain_constraint_list_high_pli[grep("C2H2", 
                                                        rownames(other_domain_constraint_matrix_high_pli)), 
                                                   grep("PRDM", colnames(other_domain_constraint_matrix_high_pli))]

c2h2_prdm <- apply(c2h2_prdm, 2, function(xx) length(which(xx == 1)))

prdm <- data.frame(set = c(0,1,0,0), c2h2 = c(6,6,3,4))
rownames(prdm) <- c("PRDM16", "PRDM2", "PRDM1", "PRDM4")

quartz(file = "prdm_high_pli.pdf", height = 3.2, width = 3.7, type = "pdf")
par(mar = c(4.8,5.5,1,1)+0.1)
barplot(as.matrix(t(prdm)), beside = TRUE, col = c(rgb(1,0,0,0.7), "royalblue1"), 
        border = "white", las = 2, yaxt = 'n', ylab = "number of constrained\ndomains")

axis(2, at = c(0,3,6))
legend <- legend("topright", legend = c("SET domain", "C2H2 Zinc finger"), fill = c(rgb(1,0,0,0.7), 
                                                                                    "royalblue1"), border = "white", 
                 cex = 0.7, bty = "n")
par(op)

dev.off()


