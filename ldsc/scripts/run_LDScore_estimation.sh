# Run LDScore estimation.
# Peter Hickey
# 2018-10-31

# SGE variables  ---------------------------------------------------------------

#$ -l mem_free=4G
#$ -l h_vmem=4.1G
#$ -m n
#$ -l h_fsize=500G
#$ -l cegs
#$ -pe local 8

# Load modules -----------------------------------------------------------------

module load conda_R/3.5.x
module load python/2.7.6

# Run R script -----------------------------------------------------------------

Rscript run_LDScore_estimation.R ${SGE_TASK_ID}
