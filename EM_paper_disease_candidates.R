load(file = "data/epiGenesDF.rda")
load(file = "em_adjmatrix.rda")


no_disease <- epiGenesDF$Gene_name[which(epiGenesDF$is.neuro.associated == FALSE
                                         & epiGenesDF$is.ca.associated == FALSE
                                         & epiGenesDF$is.MDEM == FALSE)]

ajhg_neuro <- c("KDM2A", "KDM2B", "KDM3A", "KDM3B", 
                "KDM4A", "KDM4B", "JMJD1C", "KDM5A", "KDM6B", "KDM7A", "PHF2", #the hmt and hdm candidates for developmental disorders from the 2018 ajhg paper 
                "SUV39H1", "SUV39H2", "SETDB1", "KMT2B", "KMT2C", "KMT2E", "ASH1L", "KMT5A", "KMT5B", "PRDM2")


non_disease_associated <- epiGenesDF$Gene_name_exac[which(epiGenesDF$Gene_name %in% 
                                                            unique(c(no_disease[-which(no_disease %in% ajhg_neuro)])))]

non_disease_associated <- non_disease_associated[-which(non_disease_associated 
                                                        %in% epiGenesDF$Gene_name_exac[which(epiGenesDF$Chr %in% c("chrX", "chrY"))])]
#non_disease_associated_plis <- getpLIs(non_disease_associated, EMinput = TRUE, name_to_use = "exac")


#probably use this instead
non_disease_associated_plis <- sapply(non_disease_associated, 
                                      function(xx) epiGenesDF$exac.pLI[which(epiGenesDF$Gene_name_exac %in% c(xx))])


only_drivers <- epiGenesDF$Gene_name[which(epiGenesDF$is.ca.associated == TRUE & 
                                             epiGenesDF$is.MDEM == FALSE &
                                             epiGenesDF$is.neuro.associated == FALSE)]

only_driver_associated <- epiGenesDF$Gene_name_exac[which(epiGenesDF$Gene_name 
                                                          %in% only_drivers[-which(only_drivers %in% ajhg_neuro)])]

only_driver_associated <- only_driver_associated[-which(only_driver_associated 
                                                        %in% epiGenesDF$Gene_name_exac[which(epiGenesDF$Chr %in% c("chrX", "chrY"))])]


only_driver_associated_plis <- sapply(only_driver_associated, 
                                      function(xx) epiGenesDF$exac.pLI[which(epiGenesDF$Gene_name_exac %in% c(xx))])


discan1 <- non_disease_associated[which(non_disease_associated_plis > 0.9)]

discan2 <- only_driver_associated[which(only_driver_associated_plis > 0.9)]

discan <- epiGenesDF$Gene_name[which(epiGenesDF$Gene_name_exac %in% c(discan1, discan2))]


coexpression_status <- rep(NA, length(discan))

disease_candidates <- data.frame(gene_name = discan, coexpression_status = coexpression_status)


disease_candidates$coexpression_status[which(disease_candidates$gene_name 
                                             %in% epiGenesDF$Gene_name[which(epiGenesDF$Gene_name_gtex %in% rownames(adjmatrix)[1:74])])] <- "highly coexpressed"

disease_candidates$coexpression_status[which(disease_candidates$gene_name 
                                             %in% epiGenesDF$Gene_name[which(epiGenesDF$Gene_name_gtex %in% rownames(adjmatrix)[75:157])])] <- "coexpressed"

disease_candidates$coexpression_status[which(disease_candidates$gene_name 
                                             %in% epiGenesDF$Gene_name[which(epiGenesDF$Gene_name_gtex %in% rownames(adjmatrix)[158:270])])] <- "not coexpressed"



write_csv(disease_candidates, "EM_disease_candidates.csv")
