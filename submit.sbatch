#!/bin/bash
#SBATCH --job-name=main
#SBATCH --output=logs/main.out
#SBATCH --error=logs/main.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=40000
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load python/3.10.5
module load plink/2.0
module load gcta/1.94.1
module load R
module load regenie
#source ../myVilma/bin/activate 

echo "SLURM_JOBID="$SLURM_JOBID
rm -r .snakemake/
cat snakefile_InUKBB
snakemake --unlock -s snakefile_InUKBB
snakemake --profile cluster-setup/ -s snakefile_InUKBB --rerun-incomplete
