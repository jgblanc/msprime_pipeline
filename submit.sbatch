#!/bin/bash
#SBATCH --job-name=main
#SBATCH --output=logs/main.out
#SBATCH --error=logs/main.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=40000
#SBATCH --partition=tier1q

module load gcc/12.1.0
module load python/3.10.5
#source ../myVilma/bin/activate 

echo "SLURM_JOBID="$SLURM_JOBID
rm -r .snakemake/
cat snakefile_main
snakemake --unlock -s snakefile_main
snakemake --profile cluster-setup/ -s snakefile_main --rerun-incomplete
