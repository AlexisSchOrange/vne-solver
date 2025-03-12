#!/bin/bash
#SBATCH --job-name=colge-different-k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=alexis.schneider@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=colge-different-k.txt
#SBATCH --error=colge-different-kerror.txt


# Run the job
srun julia whichkscript.jl