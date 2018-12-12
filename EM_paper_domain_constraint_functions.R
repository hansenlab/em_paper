getDomainMultiplicityforGene <- function(gene_name, domains_vector){
  df <- domain_constraint[which(domain_constraint$gene_name %in% c(gene_name)), ]
  sapply(domains_vector, function(xx) {
    df <- df[which(df$domain_name %in% c(xx)), ]
    method_1 <- grep("domain_1", df$domain_number)
    method_2 <- grep("domain_2", df$domain_number)
    method_3 <- grep("domain_3", df$domain_number)
    method_4 <- grep("domain_4", df$domain_number)
    method_0 <- grep("domain_0", df$domain_number)
    
    if (!(length(method_1) %in% c(0))){
      output <- length(unique(df$domain_number[method_1])) + length(unique(df$domain_number[method_0]))
    } else {
      if (!(length(method_2) %in% c(0))){
        output <- length(unique(df$domain_number[method_2])) + length(unique(df$domain_number[method_0]))
      } else {
        if (!(length(method_3) %in% c(0))){
          output <- length(unique(df$domain_number[method_3])) + length(unique(df$domain_number[method_0]))
        } else {
          output <- length(unique(df$domain_number[method_4])) + length(unique(df$domain_number[method_0]))
        }
      }
    }
    output})
}



####the following is the main function that outputs the constraint of the domains of a given gene
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
          100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
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
            100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
            #weighted_constraint
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
              100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
              #weighted_constraint
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
                100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
                #weighted_constraint
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
                100*round(sum(df_above_90th_pct$bp_overlap)/sum(df$bp_overlap), 2)
                #weighted_constraint
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
