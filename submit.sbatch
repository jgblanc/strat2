#!/bin/bash
#SBATCH --job-name=HGDP
#SBATCH --output=logs/HGDP.out
#SBATCH --error=logs/HGDP.err
#SBATCH --time=10:00:00
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
#rm -r .snakemake/
cat snakefile_HGDP1KG
snakemake --unlock -s snakefile_HGDP1KG
snakemake --profile cluster-setup/ -s snakefile_HGDP1KG --rerun-incomplete
