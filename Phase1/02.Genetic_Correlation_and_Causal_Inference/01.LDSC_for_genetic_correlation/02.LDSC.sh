#!/bin/bash
#PBS -S /bin/bash
#PBS -N ldsc
#PBS -l nodes=cu03:ppn=1
#PBS -e ldsc.err
#PBS -o ldsc.log
#PBS -q batch
#PBS -V

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

source ~/miniconda3/bin/activate ldsc

ref_dir='/public/jiangjw/reference/eur_w_ld_chr/' #modify into your path
result_dir='/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC/test' #modify into your path
munge_results_dir='/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC/test' #modify into your path
run_dir='/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC/ldsc-master' #modify into your path

pheno1="Anxiety2021"
pheno2="DiagnosedAR"
 
for i in ${pheno1[@]}
do
	for j in ${pheno2[@]}
	do  
		python ${run_dir}/ldsc.py \
   		--rg ${munge_results_dir}/${i}.sumstats.gz,${munge_results_dir}/${j}.sumstats.gz \
   		--ref-ld-chr ${ref_dir} \
   		--w-ld-chr ${ref_dir} \
   		--out ${result_dir}/ldsc_${i}_${j}
	done
done

