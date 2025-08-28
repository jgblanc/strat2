#!/bin/bash
pfile_path=$1
pheno_path=$2
outfile=$3
overlap_snps=$4

plink2 \
  --pfile $pfile_path \
  --extract $overlap_snps \
  --glm omit-ref allow-no-covars \
  --pheno $pheno_path \
  --pheno-name tvec \
  --out $outfile
