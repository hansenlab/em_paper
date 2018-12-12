#######Code that downloads good quality (defined as "optimal idr thresholded peaks") 
#######TF chip-seq data from ENCODE. Specifically we will download the data for K562



full_encode_metadata <- read_tsv("metadata.tsv")
full_encode_metadata$`Biosample treatments`[
  which(is.na(full_encode_metadata$`Biosample treatments`))] <- "none"


getChIPInfoDF <- function(cell_line){
  df <- data.frame(regulatory_factor = full_encode_metadata$`Experiment target`[
    which(full_encode_metadata$`Output type` == "optimal idr thresholded peaks" 
          & full_encode_metadata$`Biosample treatments` == "none" 
          & full_encode_metadata$`Biosample term name` == cell_line)], 
    accession = full_encode_metadata$`File accession`[
      which(full_encode_metadata$`Output type` == "optimal idr thresholded peaks" 
            & full_encode_metadata$`Biosample treatments` == "none" 
            & full_encode_metadata$`Biosample term name` == cell_line)], 
    url = full_encode_metadata$`File download URL`[
      which(full_encode_metadata$`Output type` == "optimal idr thresholded peaks" 
            & full_encode_metadata$`Biosample treatments` == "none" 
            & full_encode_metadata$`Biosample term name` == cell_line)],
    assembly = full_encode_metadata$Assembly[
      which(full_encode_metadata$`Output type` == "optimal idr thresholded peaks" 
            & full_encode_metadata$`Biosample treatments` == "none" 
            & full_encode_metadata$`Biosample term name` == cell_line)])
  
  df$regulatory_factor <- gsub("-human", "", df$regulatory_factor)
  df$regulatory_factor <- gsub("eGFP-", "", df$regulatory_factor)
  df$regulatory_factor <- gsub("FLAG-", "", df$regulatory_factor)
  df <- df[-which(duplicated(df$regulatory_factor)), ]
  df
}

k562_chip <- getChIPInfoDF("K562")
file_names_to_use <- sapply(k562_chip$regulatory_factor, function(xx) 
  paste0("encode_tf_chip/k562/", xx, ".bed"))
sapply(1:dim(k562_chip)[1], function(xx) download.file(k562_chip$url[xx], 
                                                       destfile = file_names_to_use[xx]))

k562_dfs_list <- lapply(1:dim(k562_chip)[1], function(xx) {
  df <- as.data.frame(
    read.table(file_names_to_use[xx], header = FALSE, 
               sep = "\t", stringsAsFactors = FALSE, quote = ""))
  df <- df[, c(1,2,3)]
  colnames(df) <- c("chr", "start", "stop")
  df$factor_name <- k562_chip$regulatory_factor[xx]
  df$assembly <- k562_chip$assembly[xx]
  df <- df[which(df$chr %in% unique(epiGenesDF$Chr)), ]
  df
})


k562_granges_list <- lapply(k562_dfs_list, function(xx) 
  GRanges(seqnames = xx$chr, IRanges(start = xx$start, 
                                     width = (xx$stop - xx$start+1)), 
          factor_name = xx$factor_name, assembly = xx$assembly))

grch38_chip_granges <- lapply(which(k562_chip$assembly == "GRCh38"), function(xx) k562_granges_list[[xx]])
hg19_chip_granges <- lapply(which(k562_chip$assembly == "hg19"), function(xx) k562_granges_list[[xx]])


converted_grch38_granges <- lapply(grch38_chip_granges, function(xx) {
  #path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  #ch = import.chain(path)
  converted_granges_list <- liftOver(xx, ch)
  converted_grange <- unlist(converted_granges_list)
  converted_grange$assembly <- "hg19"
  converted_grange
}) 

final_k562_granges_list <- list()
for (i in 1:length(hg19_chip_granges)){
  final_k562_granges_list[[i]] <- hg19_chip_granges[[i]]
}
for (i in 1:length(converted_grch38_granges)){
  final_k562_granges_list[[length(hg19_chip_granges)+i]] <- converted_grch38_granges[[i]]
}

tfs_k562 <- k562_chip$regulatory_factor
