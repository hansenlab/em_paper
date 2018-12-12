library(readr)
gtex <- read_tsv("extdata/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz", skip=2)
gtex.df <- gtex[, 1:2]
colnames(gtex.df)[2]<-"gene_names"
colnames(gtex.df)[1]<-"ensembl"

gencode <- read_tsv("extdata/gencode.v19.genes.patched_contigs.gtf.gz", skip = 5, col_names = FALSE)
genes <- gencode[gencode[,3] == "transcript", ]
colnames(genes)[9]<-"X9"
genes$type <- gsub("\"", "", sub("gene_type ", "", sapply(strsplit(genes$X9, "; "), function(xx) xx[3])))
genes$id <- gsub("\"", "", sub("gene_id ", "", sapply(strsplit(genes$X9, "; "), function(xx) xx[1])))
table(genes$type)
head(genes)
ids <- genes$id[genes$type == "protein_coding"]





gtex <- gtex[,-(1:2)]
cSums <- colSums(gtex) / 10^7
wh.na <- which(is.na(cSums))
gtex.mat <- as.matrix(gtex[, -wh.na])
gtex.mat <- log2(sweep(gtex.mat + 1, 2, FUN = "/", cSums[-wh.na]))
gtex.cds.mat <- gtex.mat[gtex.df$ensembl %in% ids,]
gtex.cds.df <- gtex.df[gtex.df$ensembl %in% ids,]
save(gtex.cds.mat, gtex.cds.df, file = "extdata/gtex.cds.rda")

