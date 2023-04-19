#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=1g
#SBATCH --error=%x-%j-%a.e
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pierre.ledenmat@kuleuven.be
#SBATCH --account=lp_desenderlab
#SBATCH --array=1-74

cd $PBS_O_WORKDIR
source activate rmodel
Rscript --vanilla fit_ldc_nn_ab_mod_alpha.R ${SLURM_ARRAY_TASK_ID}
