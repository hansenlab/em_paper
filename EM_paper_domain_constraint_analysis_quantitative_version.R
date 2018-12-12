epigenetic_domains <- c("SET domain", 
                        "GNAT domain", "Histone acetyltransferase domain, MYST-type", "Histone acetyltransferase Rtt109/CBP",
                        "C-5 cytosine methyltransferase", 
                        "JmjC domain", "FAD/NAD(P)-binding domain", 
                        "Histone deacetylase domain", "Sirtuin family, catalytic core domain", 
                        "2OGFeDO, oxygenase domain", "Chromo domain", 
                        "Zinc finger, PHD-type", "Tudor domain", 
                        "PWWP domain", "Protein ASX-like, PHD domain", 
                        "Bromo adjacent homology (BAH) domain", "ADD domain", 
                        "Mbt repeat", "Zinc finger, CW-type", "WD40 repeat", 
                        "Ankyrin repeat", "Bromodomain", "Zinc finger, CXXC-type",
                        "C-terminal domain", 
                        "Methyl-CpG DNA binding", "Helicase superfamily 1/2, ATP-binding domain", 
                        "Helicase, C-terminal", "Helicase/SANT-associated domain")

epigenetic_domains_for_labeling <- c("SET domain", "GNAT domain", 
                                     "Histone acetyltransferase domain, MYST-type", "Histone acetyltransferase Rtt109/CBP",
                                     "C-5 cytosine methyltransferase", "JmjC domain", "FAD/NAD(P)-binding domain", 
                                     "Histone deacetylase domain", "Sirtuin family, catalytic core domain", 
                                     "2OGFeDO, oxygenase domain", paste0("Chromo domain ", 1:2), 
                                     paste0("Zinc finger, PHD-type ", 1:7),  paste0("Tudor domain ", 1:8), 
                                     paste0("PWWP domain ", 1:2), "Protein ASX-like, PHD domain", 
                                     paste0("Bromo adjacent homology (BAH) domain ", 1:2), "ADD domain", 
                                     paste0("Mbt repeat ", 1:4), "Zinc finger, CW-type", paste0("WD40 repeat ", 1:8), 
                                     paste0("Ankyrin repeat ", 1:7), paste0("Bromodomain ", 1:6), paste0("Zinc finger, CXXC-type ", 1:3),
                                     "C-terminal domain", "Methyl-CpG DNA binding", "Helicase superfamily 1/2, ATP-binding domain", 
                                     "Helicase, C-terminal", "Helicase/SANT-associated domain")


domain_constraint <- as.data.frame(read.table("EMdomainsvCCR_2.bed", header = TRUE, 
                                              sep="\t", stringsAsFactors=FALSE, quote=""))

epigenetic_domains <- gsub("[ ]" , "_", epigenetic_domains)
epigenetic_domains <- gsub("[-]", "_", gsub("[,]", "", epigenetic_domains))

epigenes_included <- unique(domain_constraint$gene_name) #there are no epigenes_incuded on the X or Y chromosome
#epigenes_included_plis <- sapply(epigenes_included[-which(epigenes_included %in% epiGenesDF$Gene_name[which(epiGenesDF$Chr %in% c("chrX", "chrY"))])], function(xx) exactab$pLI[which(exactab$gene %in% c(xx))])
epigenes_included_plis <- sapply(epigenes_included, function(xx) epiGenesDF$exac.pLI[which(epiGenesDF$Gene_name %in% c(xx))])

epigenes_included_high_pli <- epigenes_included[which(epigenes_included_plis > 0.9)]
epigenes_included_low_pli <- epigenes_included[which(epigenes_included_plis < 0.1)]


###get a matrix where rows are genes and columns are the domains, and each entry is equal to the multiplicity of a given domain
###in a given gene
em_domain_multiplicity_matrix <- sapply(epigenes_included, function(xx) 
  getDomainMultiplicityforGene(xx, epigenetic_domains))
em_max_multiplicities <- rowMaxs(em_domain_multiplicity_matrix)
em_domains_appearing_multiple_times <- epigenetic_domains[which(em_max_multiplicities > 1)]



em_domain_constraint_list_high_pli <- sapply(epigenes_included_high_pli, function(xx)
  getDomainPercentilesforGene(xx, epigenetic_domains, 
                              em_domains_appearing_multiple_times, em_max_multiplicities))
em_domain_constraint_list_low_pli <- sapply(epigenes_included_low_pli, function(xx)
  getDomainPercentilesforGene(xx, epigenetic_domains, 
                              em_domains_appearing_multiple_times, em_max_multiplicities))

em_domain_constraint_matrix_high_pli <- as.matrix(as.data.frame(em_domain_constraint_list_high_pli))
em_domain_constraint_matrix_low_pli <- as.matrix(as.data.frame(em_domain_constraint_list_low_pli))


####now other domains
other_domains <- unique(domain_constraint$domain_name)
other_domains <- other_domains[-which(other_domains %in% epigenetic_domains)]

other_domain_multiplicity_matrix <- sapply(epigenes_included, function(xx) 
  getDomainMultiplicityforGene(xx, other_domains))
other_max_multiplicities <- rowMaxs(other_domain_multiplicity_matrix)
other_domains_appearing_multiple_times <- other_domains[which(other_max_multiplicities > 1)]


other_domain_constraint_list_high_pli <- sapply(epigenes_included_high_pli, function(xx)
  getDomainPercentilesforGene(xx, other_domains, 
                              other_domains_appearing_multiple_times, other_max_multiplicities))
other_domain_constraint_list_low_pli <- sapply(epigenes_included_low_pli, function(xx)
  getDomainPercentilesforGene(xx, other_domains, 
                              other_domains_appearing_multiple_times, other_max_multiplicities))

other_domain_constraint_matrix_high_pli <- as.matrix(as.data.frame(other_domain_constraint_list_high_pli))
other_domain_constraint_matrix_low_pli <- as.matrix(as.data.frame(other_domain_constraint_list_low_pli))



em_domain_constraint_vector_high_pli <- as.vector(em_domain_constraint_matrix_high_pli)
em_domain_constraint_vector_high_pli <- em_domain_constraint_vector_high_pli[-
                           which(em_domain_constraint_vector_high_pli %in% c(190))]


other_domain_constraint_vector_high_pli <- as.vector(other_domain_constraint_matrix_high_pli)
other_domain_constraint_vector_high_pli <- other_domain_constraint_vector_high_pli[-
                           which(other_domain_constraint_vector_high_pli %in% c(190))]

em_domain_constraint_vector_low_pli <- as.vector(em_domain_constraint_matrix_low_pli)
em_domain_constraint_vector_low_pli <- em_domain_constraint_vector_low_pli[-
                                 which(em_domain_constraint_vector_low_pli %in% c(190))]

other_domain_constraint_vector_low_pli <- as.vector(other_domain_constraint_matrix_low_pli)
other_domain_constraint_vector_low_pli <- other_domain_constraint_vector_low_pli[-
                                 which(other_domain_constraint_vector_low_pli %in% c(190))]



mean_constraint_of_em_domains_per_gene <- apply(em_domain_constraint_matrix_high_pli, 2, function(xx) mean(xx[-which(xx %in% c(190))]))
mean_constraint_of_other_domains_per_gene <- apply(other_domain_constraint_matrix_high_pli, 2, function(xx) mean(xx[-which(xx %in% c(190))]))

diffs_high_pli <- mean_constraint_of_em_domains_per_gene - mean_constraint_of_other_domains_per_gene


mean_constraint_of_em_domains_per_gene <- apply(em_domain_constraint_matrix_low_pli, 2, function(xx) mean(xx[-which(xx %in% c(190))]))
mean_constraint_of_other_domains_per_gene <- apply(other_domain_constraint_matrix_low_pli, 2, function(xx) mean(xx[-which(xx %in% c(190))]))

diffs_low_pli <- mean_constraint_of_em_domains_per_gene - mean_constraint_of_other_domains_per_gene



