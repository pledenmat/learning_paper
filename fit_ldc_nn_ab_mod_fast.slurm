#!/bin/bash -l
#SBATCH --time=02:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=5g
#SBATCH --error=slurm-%j_%a.e
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pierre.ledenmat@kuleuven.be
#SBATCH --account=lp_desenderlab
#SBATCH --array=1-74

cd $VSC_DATA
source activate rmodel
Rscript --vanilla fit_ldc_nn_ab_mod_fast.R ${SLURM_ARRAY_TASK_ID}
