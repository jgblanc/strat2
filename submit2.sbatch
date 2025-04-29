#!/bin/bash
#SBATCH --job-name=main2
#SBATCH --output=logs/main2.out
#SBATCH --error=logs/main2.err
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=40000
#SBATCH --partition=tier1q

module load gcc/12.1.0
module load python/3.10.5
module load plink/2.0
module load gcta/1.94.1
module load R
#source ../myVilma/bin/activate 

echo "SLURM_JOBID="$SLURM_JOBID
rm -r .snakemake/
cat snakefile_main2
snakemake --unlock -s snakefile_main2
snakemake --profile cluster-setup/ -s snakefile_main2 --rerun-incomplete
