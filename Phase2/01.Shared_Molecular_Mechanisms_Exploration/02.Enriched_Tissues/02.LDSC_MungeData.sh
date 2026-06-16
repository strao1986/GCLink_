#/bin/bash
#PBS -S /bin/bash
#PBS -N ldsc
#PBS -l nodes=cu04:ppn=5
#PBS -e ldsc.err
#PBS -o ldsc.log
#PBS -q batch
#PBS -V

export OMP_NUM_THREADS=5
export OPENBLAS_NUM_THREADS=5
export MKL_NUM_THREADS=5
export VECLIB_MAXIMUM_THREADS=5
export NUMEXPR_NUM_THREADS=5

source ~/miniconda3/bin/activate ldsc

ref_dir='/public/jiangjw/reference' #modify into your path
GWAS_sumstats_path='/public/jiangjw/GWAS_sumstats/EUR' #modify into your path
result_dir='/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC/test' #modify into your path
run_dir='/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC/ldsc-master' #modify into your path

Disease_phenos=("Anxiety2021" "DiagnosedAR")

#### Munge Data ####
for i in ${Disease_phenos[@]}
do 
	python ${run_dir}/munge_sumstats.py \
	--sumstats ${GWAS_sumstats_path}/${i}.clean.txt \
	--N-col N \
	--out ${result_dir}/${i} \
	--merge-alleles ${ref_dir}/w_hm3.snplist \
	--chunksize 500000 # if this option not set, the program may not run.
done


