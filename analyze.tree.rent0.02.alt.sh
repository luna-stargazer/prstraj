#!/bin/bash
#SBATCH --job-name=neutral_analyze_tree_0_rent
#SBATCH --output=/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/neutral_sel_0.02/analyze.tree.alt.rent.out
#SBATCH --error=/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/neutral_sel_0.02/analyze.tree.rent.alt.err
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --account=lc_mde
#SBATCH --partition=cmb
#SBATCH --mem=6000
#SBATCH --mem-per-cpu=4000

/home/cmb-00/mde/lindadin/conda/envs/r_env/bin/Rscript /home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/analyze.trees.alt.R 0.02 0 rent
