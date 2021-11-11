#!/bin/bash
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1:00:00
#SBATCH --job-name PCAF-debug
#SBATCH --output=/scratch/v/vanderli/cassane/output/pcaf_debug_%j.txt
#SBATCH --mail-type=FAIL

module load gcc/8.3.0
module load openmpi/3.1.3
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

mpirun -np 4 $APY3 $path2file