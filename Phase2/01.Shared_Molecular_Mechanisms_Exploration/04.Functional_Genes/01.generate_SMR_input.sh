#!/bin/bash
#PBS -N SMR
#PBS -l nodes=cu04:ppn=8
#PBS -e SMR.err
#PBS -o SMR.log
#PBS -q batch
#PBS -V

GWAS_sumstats_dir='/public/jiangjw/GWAS_sumstats/EUR' #modify into your path
result_dir='/public/jiangjw/02.anxiety_ADs_EUR/09.SMR/Github_code_test' #modify into your path

phenos=("Anxiety2021" "DiagnosedAR")

#### generate SMR input ####
for j in ${phenos[@]};
do
	cat ${GWAS_sumstats_dir}/${j}.clean.txt | cut -f 3-10 | awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$7,$4,$5,$6,$8}' | sed '1d' | sed "1i SNP\tA1\tA2\tfreq\tb\tse\tp\tn" > ${result_dir}/${j}_SMR.clean.txt
done

