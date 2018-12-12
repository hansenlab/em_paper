load(file = "mdem/epigenetic_list/objects/exactab.rda")
all_interpro <- read.delim("data/interpro_domain_names.dat", header=FALSE, stringsAsFactors = FALSE)
all_interpro$V1 <- as.character(all_interpro$V1) #interpro ID
all_interpro$V2 <- as.character(all_interpro$V2) #domain name
names(all_interpro) = c("ID", "domain")
load(file = "mdem/epigenetic_list/objects/epiGenesDF.rda")

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

###now exclude some domains because they overlap with the above epigenetic domains
domains_to_exclude <- c("SNF2-related, N-terminal domain", 
                        "P-loop containing nucleoside triphosphate hydrolase", "Zinc finger, PHD-finger", 
                        "Zinc finger, RING-type", "Zinc finger, FYVE/PHD-type", 
                        "Zinc finger, RING/FYVE/PHD-type", "Acyl-CoA N-acyltransferase",
                        "S-adenosyl-L-methionine-dependent methyltransferase", "DNA (cytosine-5)-methyltransferase 1, metazoa",
                        "Amine oxidase", "Histone deacetylase superfamily", "Histone deacetylase", 
                        "Sirtuin family", "Sirtuin, class I", "DHS-like NAD/FAD-binding domain", 
                        "Histone lysine-specific demethylase", "Chromo/chromo shadow domain",
                        "Chromo domain-like", "Chromo domain subgroup", "Agenet domain, plant type", 
                        "ASX-like protein 1", "ASX-like protein 2", "ASX-like protein 3", 
                        "Polycomb protein ASX/ASX-like", "DNA (cytosine-5)-methyltransferase 3-like", 
                        "DNA (cytosine-5)-methyltransferase 3B", "DNA (cytosine-5)-methyltransferase 3A", 
                        "WD40/YVTN repeat-like-containing domain", "WD40-repeat-containing domain", 
                        "Ankyrin repeat-containing domain", "DNA-binding domain", "Zinc finger, C2H2", 
                        "RNA binding activity-knot of a chromodomain", "Zinc finger, MYND-type")
interpro_ids_to_exclude <- all_interpro$ID[which(all_interpro$domain %in% domains_to_exclude)]



library(ensembldb)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
library(biomaRt)
library("GenomicFeatures")
library(Pbase)
library(Biostrings)
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                  path="/biomart/martservice", dataset="hsapiens_gene_ensembl")


###create a function that does this for each scanning method separately, so that I avoid overlaps between domain predictions
mapDomainsToCoords <- function(gene_name, method = c("smart", "pfam", "pfscan", "supefamily")){
  initial_gene_object <- Proteins(edb, filter = ~ genename == gene_name)
  txids <- acols(initial_gene_object)$tx_id
  txid <- txids[which(txids %in% gsub("[.].*$", "", exactab$transcript[which(exactab$gene %in% c(epiGenesDF$Gene_name_exac[which(epiGenesDF$Gene_name %in% c(gene_name))]))]))]
  gene_object <- Proteins(edb, filter = TxIdFilter(txid))
  gene_map <- mapToGenome(gene_object, edb, idType = "tx_id", id = "tx_id")
  domain_ranges <- ranges(gene_map[[1]])
  domain_scanner_method_id <- names(domain_ranges)
  
  
  ###now select only those ranges that correspond to the interpro domains included in the epiGenesDF column
  ##and only from the smart scanner
  domains_irange <- pranges(gene_object)$ProteinDomains[[1]]
  domains_source <- mcols(domains_irange)$protein_domain_source
  domains_interpro_id <- mcols(domains_irange)$interpro_accession
  
  
  selected_method_interpro_ids <- domains_interpro_id[which(domains_source %in% method)]
  selected_method_domains_method_id <- names(gene_map[[1]])[which(domains_source %in% method)]
  
  interpro_ids <- strsplit(epiGenesDF$Interpro.Domains[which(epiGenesDF$Gene_name %in% c(gene_name))], ";")[[1]]
  if (length(which(interpro_ids %in% interpro_ids_to_exclude)) > 0){
    interpro_ids <- interpro_ids[-which(interpro_ids %in% interpro_ids_to_exclude)]}
  
  if (!(length(which(selected_method_interpro_ids %in% interpro_ids)) %in% c(0))){
    ids_to_use <- selected_method_interpro_ids[which(selected_method_interpro_ids %in% interpro_ids)]
    indices_to_use1 <- which(domains_interpro_id %in% ids_to_use & domains_source %in% method)
    domain_names_to_use <- sapply(indices_to_use1, function(xx) all_interpro$domain[which(all_interpro$ID %in% 
                                                                                            domains_interpro_id[xx])])
    irange_to_use <- domains_irange[indices_to_use1]
    domain_coord_df <- data.frame(start(irange_to_use), end(irange_to_use), domain_names_to_use)
    
    colnames(domain_coord_df) <- c("amino_acid_start", "amino_acid_end", "domain_name")
    domain_coord_df <- domain_coord_df[order(domain_coord_df$amino_acid_start),]
    for (i in 1:(dim(domain_coord_df)[1])){
      if (method %in% "smart"){
        domain_coord_df$domain_number[i] <- paste0("domain_", 1, i)
      } else if (method %in% "pfam"){
        domain_coord_df$domain_number[i] <- paste0("domain_", 2, i)
      } else if (method %in% "pfscan"){
        domain_coord_df$domain_number[i] <- paste0("domain_", 3, i)
      } else {
        domain_coord_df$domain_number[i] <- paste0("domain_", 4, i)
      }
    }
    for (i in 1:(dim(domain_coord_df)[1])){
      domain_coord_df$gene_name[i] <- gene_name
    }
    domain_coord_df
    
    domain_genomic_coords <- domain_ranges[which(names(domain_ranges) %in% names(domains_irange)[c(indices_to_use1)])]
    domain_genomic_coords <- as.data.frame(domain_genomic_coords)
    if (epiGenesDF$Strand[which(epiGenesDF$Gene_name %in% c(gene_name))] %in% "-"){
      domain_genomic_coords <- domain_genomic_coords[order(domain_genomic_coords$start, decreasing = TRUE), ]
    }
    
    cumulative_number_of_aa_in_domains <- vector()
    cumulative_number_of_aa_in_domains[1] <- domain_coord_df$amino_acid_end[1] - domain_coord_df$amino_acid_start[1] + 1
    if (length(domain_coord_df$amino_acid_start) > 1){
      for (i in 2:length(domain_coord_df$amino_acid_start)){
        cumulative_number_of_aa_in_domains[i] <- domain_coord_df$amino_acid_end[i] - domain_coord_df$amino_acid_start[i] + 1 + cumulative_number_of_aa_in_domains[i-1]
      }
    }
    
    cumulative_number_of_aa_encoded <- vector()
    effective_width <- vector()
    cumulative_number_of_aa_encoded[1] <- (domain_genomic_coords$width[1]-(domain_genomic_coords$width[1]%%3))/3
    effective_width[1] <- domain_genomic_coords$width[1]
    if (dim(domain_genomic_coords)[1] > 1){
      for (i in 2:(dim(domain_genomic_coords)[1])){
        if (cumulative_number_of_aa_encoded[i-1] %in% cumulative_number_of_aa_in_domains){
          cumulative_number_of_aa_encoded[i] <- cumulative_number_of_aa_encoded[i-1] + (domain_genomic_coords$width[i]-(domain_genomic_coords$width[i]%%3))/3
          effective_width[i] <- domain_genomic_coords$width[i]
        } else {
          aminoacid_to_add <- ifelse(effective_width[i-1]%%3 %in% c(0), 0, 1)
          if (effective_width[i-1]%%3 %in% c(0)){
            cumulative_number_of_aa_encoded[i] <- cumulative_number_of_aa_encoded[i-1] + (domain_genomic_coords$width[i]-(domain_genomic_coords$width[i]%%3))/3 + aminoacid_to_add
            effective_width[i] <- domain_genomic_coords$width[i]
          } else if (effective_width[i-1]%%3 %in% c(1)){
            cumulative_number_of_aa_encoded[i] <- cumulative_number_of_aa_encoded[i-1] + ((domain_genomic_coords$width[i]-2)-((domain_genomic_coords$width[i]-2)%%3))/3 + aminoacid_to_add
            effective_width[i] <- domain_genomic_coords$width[i] - 2
          } else {
            cumulative_number_of_aa_encoded[i] <- cumulative_number_of_aa_encoded[i-1] + ((domain_genomic_coords$width[i]-1)-((domain_genomic_coords$width[i]-1)%%3))/3 + aminoacid_to_add
            effective_width[i] <- domain_genomic_coords$width[i] - 1
          }
        }
      }
    }
    
    indices_to_break <- which(cumulative_number_of_aa_encoded %in% cumulative_number_of_aa_in_domains)
    coord_list <- list()
    coord_list[[1]] <- domain_genomic_coords[1:indices_to_break[1],]
    if (length(domain_coord_df$amino_acid_start) > 1){
      for (i in 2:length(indices_to_break)){
        coord_list[[i]] <- domain_genomic_coords[(indices_to_break[i-1]+1):indices_to_break[i],]
      }
    }
    coord_list[[length(indices_to_break)+1]] <- domain_coord_df
    for (i in 1:(length(coord_list)-1)){
      coord_list[[i]]$domain_name <-  coord_list[[length(coord_list)]]$domain_name[i]
      coord_list[[i]]$domain_number <-  coord_list[[length(coord_list)]]$domain_number[i]
      coord_list[[i]]$gene_name <-  coord_list[[length(coord_list)]]$gene_name[i]
    }
    final_coord_list <- list(rbindlist(lapply(1:(length(coord_list)-1), function(xx) coord_list[[xx]])), 
                             coord_list[[length(coord_list)]])
    final_coord_list
  } else {
    coord_list <- list(list(), list())
    coord_list
  }
}  

