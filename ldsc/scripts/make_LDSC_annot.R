# Make .annot file for use with ldsc containing Leandros's features.
# Peter Hickey
# 2018-11-28

# Setup ------------------------------------------------------------------------

library(rtracklayer)
library(readr)
library(here)

# Load Leandros's features -----------------------------------------------------

f <- list.files(
  path = here("data", "reldscoreregressionresponsetoreviewers"),
  pattern = "^.*\\.bed$",
  full.names = TRUE)
names(f) <- sub(".bed", "", basename(f))
x <- lapply(
  X = f,
  FUN = function(xx) {
    z <- import(xx)
    z <- sortSeqlevels(z)
    sort(z)
  })

# Make annotations -------------------------------------------------------------

categories <- x
saveRDS(categories, here("objects", "Leandros_features.rds"))

mclapply(1:22, function(sl) {
  message(sl)
  cds <- read_tsv(
    here(
      "extdata",
      "Phase1",
      "cell_type_groups",
      paste0("CNS.", sl, ".annot.gz")))
  cds_gr <- GRanges(paste0("chr", cds$CHR), IRanges(cds$BP, width = 1L))
  annot <- cds[, c("CHR", "BP", "SNP", "CM")]
  annot[names(categories)] <- mclapply(names(categories), function(cn) {
    stopifnot(isDisjoint(categories[[cn]]))
    as.integer(overlapsAny(cds_gr, categories[[cn]]))
  }, mc.cores = 4)
  # 'Marginal' annotation file
  mclapply(names(categories), function(cn) {
    fl <- here(
      "output",
      "ldsc",
      paste0(cn, ".Phase1.", sl, ".annot.gz"))
    write_tsv(annot[, c("CHR", "BP", "SNP", "CM", cn)], fl)
  }, mc.cores = 8)
}, mc.cores = 1)
