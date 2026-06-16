#!/bin/bash

#### Reference file generate ####
run_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA' #modify into your path
reference_dir1='/mnt/e/Projects_documents/reference/g1000_eur' #modify into your path
reference_dir2='/mnt/e/Projects_documents/reference' #modify into your path
output_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/Github_GCLink/03.Shared_Molecular_Mechanisms_Exploration/02.Shared_Genetic_Architecture/04.Specific_Cell_Types' #modify into your path

${run_dir}/magma \
	--annotate window=10 \
	--snp-loc ${reference_dir1}/g1000_eur.bim \
	--gene-loc ${reference_dir2}/NCBI37.3.gene.loc \
	--out ${output_dir}/g1000.eur.magma.10k.snp2gene


#### run MAGMA program ####
GWAS_sumstats_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/Github_GCLink/GWAS_Summary_Statistics_Example_Preprocessed/Diseases' #modify into your path
phenos=("Anxiety2021" "DiagnosedAR")
   
for i in ${phenos[@]} 
do
	${run_dir}/magma \
	--bfile ${reference_dir1}/g1000_eur \
   	--pval ${GWAS_sumstats_dir}/${i}.clean.txt \
	use=rsid,p \
	ncol=N \
	--gene-annot ${output_dir}/g1000.eur.magma.10k.snp2gene.genes.annot \
	--out ${output_dir}/${i}.35UP.10DOWN.genes
done

