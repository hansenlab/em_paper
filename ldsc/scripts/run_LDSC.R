# Run LDSC with Leandros's categories (made in make_LDSC_annot.R)
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

categories <- readRDS(here("objects", "Leandros_features.rds"))
message("category = ", names(categories)[i])

gwasss <- list.files(
  path = here("extdata", "munge_sumstats", "Phase1"),
  full.names = TRUE,
  pattern = glob2rx("*.sumstats.gz"))


# Run LDSC (adjusting for baseline) --------------------------------------------

lapply(names(categories)[i], function(cn) {
  mclapply(gwasss, function(x) {
    bn <- sub(".sumstats.gz", "", basename(x))
      cmd <- paste0(
        "python ",
        "/users/phickey/software/ldsc/ldsc.py ",
        "--h2 ", x, " ",
        "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
        "--ref-ld-chr ../output/ldsc/", cn, ".Phase1.,",
        "../extdata/Phase1/baseline/baseline. ",
        "--overlap-annot ",
        "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
        "--out ../output/ldsc/", cn, ".", bn, ".Phase1 ",
        "--print-coefficients")
      print(cmd)
      system(cmd)
  }, mc.cores = getOption("mc.cores"))
})
