#!/bin/bash
#SBATCH --job-name=msprime
#SBATCH --output=logs/msprime.out
#SBATCH --error=logs/msprime.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=40000
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load python/3.10.5
#source ../myVilma/bin/activate 

echo "SLURM_JOBID="$SLURM_JOBID
rm -r .snakemake/
cat snakefile_main
snakemake --unlock -s snakefile_main
snakemake --profile cluster-setup/ -s snakefile_main --rerun-incomplete
