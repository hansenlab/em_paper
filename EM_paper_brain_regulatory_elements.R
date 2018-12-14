####Regulatory regions around EM genes. First get 1Mb regions around the TSS
####of EM genes and then intersect with brain enhancers

names_of_converted_in_order <- epiGenesDF$Gene_name[-c(121, 188)] #converted refers to liftover from hg38 to hg19
strand_of_converted_in_order <- epiGenesDF$Strand[-c(121, 188)]

em_hg19_coords <- as.data.frame(read.table("em_hg19_coords.bed",header = FALSE, 
                                           sep="\t", stringsAsFactors=FALSE, quote=""))
colnames(em_hg19_coords) <- c("chr", "start", "end")

#add 1 to the start because ucsc works with 0 based coordinates
em_hg19_coords$start <- em_hg19_coords$start + 1

em_hg19_coords$gene_name <- names_of_converted_in_order
em_hg19_coords$strand <- strand_of_converted_in_order

em_plusminus1Mb_tss <- em_hg19_coords


em_plusminus1Mb_tss$start <- sapply(1:293, function(xx) {
  if (em_plusminus1Mb_tss$strand[xx] %in% c("+")){
    minus <- em_hg19_coords$start[xx] - 1000000
  } else if(em_plusminus1Mb_tss$strand[xx] %in% c("-")) {
    minus <- em_hg19_coords$end[xx] - 1000000
  }
  minus
})

em_plusminus1Mb_tss$end <- sapply(1:293, function(xx) {
  if (em_plusminus1Mb_tss$strand[xx] %in% c("+")){
    plus <- em_hg19_coords$start[xx] + 1000000
  } else if (em_plusminus1Mb_tss$strand[xx] %in% c("-")){
    plus <- em_hg19_coords$end[xx] + 1000000
  }
  plus
})


em_1Mb_granges <- makeGRangesFromDataFrame(em_plusminus1Mb_tss)



######shared enhancers across many brain regions
all_brain_enhancers <- read_csv('brain_shared_enhancers.csv')
all_brain_enhancers_granges <- makeGRangesFromDataFrame(all_brain_enhancers)

##all EM cortex
EM_all_brain_enhancer_overlaps <- all_brain_enhancers_granges[
  unique(queryHits(findOverlaps(all_brain_enhancers_granges, em_1Mb_granges)))]

EM_all_brain_regulatory_element_coords <- data.frame(chr = seqnames(EM_all_brain_enhancer_overlaps), 
                                                     start = start(EM_all_brain_enhancer_overlaps),
                                                     end = end(EM_all_brain_enhancer_overlaps))


write.table(EM_all_brain_regulatory_element_coords, file="EM_all_brain_regulatory_elements.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)

##highly coexp all brain
hc_all_brain_enhancer_overlaps <- all_brain_enhancers_granges[
  unique(queryHits(findOverlaps(all_brain_enhancers_granges, em_1Mb_granges[which(
    em_plusminus1Mb_tss$gene_name %in% epiGenesDF$Gene_name[which(epiGenesDF$Gene_name_gtex 
                                                                  %in% rownames(adjmatrix)[1:74])]
  )])))]

hc_all_brain_regulatory_element_coords <- data.frame(chr = seqnames(hc_all_brain_enhancer_overlaps), 
                                                     start = start(hc_all_brain_enhancer_overlaps),
                                                     end = end(hc_all_brain_enhancer_overlaps))


write.table(hc_all_brain_regulatory_element_coords, file="highly_coexp_all_brain_regulatory_elements.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)




