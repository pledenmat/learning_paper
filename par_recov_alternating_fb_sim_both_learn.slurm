#!/bin/bash -l
#SBATCH --time=72:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=5g
#SBATCH --cluster=genius
#SBATCH --error=slurm-%j_%a.e
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pierre.ledenmat@kuleuven.be
#SBATCH --account=lp_desenderlab
#SBATCH --array=1-2

cd $VSC_DATA
export PATH=${VSC_DATA}/miniconda3/bin:${PATH}
source activate rmodel
Rscript --vanilla par_recov_alternating_fb_sim_both_learn_mean_ev.R ${SLURM_ARRAY_TASK_ID}
