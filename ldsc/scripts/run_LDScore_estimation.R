# Run LD score estimation with Leandros's categories (made in make_LDSC_annot.R)
# Peter Hickey
# 2018-11-28

# Setup ------------------------------------------------------------------------

# NOTE: On JHPCE, requires `module load python/2.7.6`

args <- commandArgs(TRUE)
i <- as.integer(args[1])
message("i = ", i)

library(parallel)
library(here)

options("mc.cores" = 8)

# Load data --------------------------------------------------------------------

# NOTE: Should be able to re-use SLDSR results from BrainEpigenome paper.
categories <- readRDS(here("objects", "Leandros_features.rds"))

seqlevels <- 1:22

message("category = ", names(categories)[i])

# Run LD Score estimation ------------------------------------------------------

lapply(names(categories)[i], function(cn) {
  message(cn)
  mclapply(seqlevels, function(sl) {
    cmd <- paste0("python ",
                  "/users/phickey/software/ldsc/ldsc.py ",
                  "--l2 ",
                  "--bfile ../extdata/Phase1/1000G_plinkfiles/1000G.mac5eur.",
                  sl, " ",
                  "--ld-wind-cm 1 ",
                  "--annot ../output/ldsc/", cn, ".Phase1.", sl, ".annot.gz ",
                  "--out ../output/ldsc/", cn, ".Phase1.", sl, " ",
                  "--print-snps ../extdata/Phase1/hapmap3_snps/hm.", sl, ".snp")
    print(cmd)
    system(cmd)
  }, mc.cores = 10)
})
