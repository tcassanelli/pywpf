#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=15:00:00
#SBATCH --job-name PCAF
#SBATCH --output=/scratch/v/vanderli/cassane/output/pcaf_ex_%j.txt
#SBATCH --mail-type=ALL

module load texlive
export MPLCONFIGDIR="$MPLCONFIGDIR:$SCRATCH" 
cd $SCRATCH

# INPUT 
# file="pypcaf_run.py"
file="pypcaf_run_sim.py"

directory="/home/v/vanderli/cassane/PCA-Folding/scripts/"
path2file=$directory$file

# python environment
APY3="/home/v/vanderli/cassane/anaconda3/bin/python"

$APY3 $path2file