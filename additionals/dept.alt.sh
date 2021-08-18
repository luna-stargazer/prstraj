#!/bin/bash
#SBATCH --job-name=getDeptv2.alt
#SBATCH --output=/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/getDeptv2.alt.out
#SBATCH --error=/home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/getDeptv2.alt.err
#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --account=lc_mde
#SBATCH --partition=cmb
#SBATCH --mem=6000
#SBATCH --mem-per-cpu=2000

/home/cmb-00/mde/lindadin/conda/envs/r_env/bin/Rscript /home/cmb-00/mde/lindadin/abbeyroad/rhps_coalescent/maintext_sims_rent_061318/getDeptv2.alt.R
