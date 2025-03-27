#!/bin/bash

# SLURM options:

#SBATCH --job-name=gmao25    # Job name
#SBATCH --output=merra2_%j.log   # Standard output and error log
#SBATCH --partition=lsst,htc               # Partition choice (htc by default)
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --mem=2000                    # Memory in MB per default
#SBATCH --time=0-05:00:00             # Max time limit = 7 days
#SBATCH --mail-user=sylvie.dagoret-campagne@ijclab.in2p3.fr          # Where to send the e-mail
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)


# Commands to be submitted:
export YEAR="2025"
export OBSERVATORY="lsst" 
source /pbs/throng/lsst/users/dagoret/desc/${YEAR}/setup_anaconda3-py311.sh
python ScanMERRA2_inst1_2d_asm_Nx_M2I1NXASMOneYear_obs_year.py -y ${YEAR} -o ${OBSERVATORY}
