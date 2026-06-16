library(data.table)
library(dplyr)

TissueIndex_dir <- "/public/jiangjw/reference/LDSC_SEG_annotation_geneset/Multi_tissue_gene_expr_1000Gv3_ldscores/" #modify into your path
Output_dir <- "/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC/LDSC_SEG/Github_code_test/" #modify into your path

tissue_index <- fread("/public/jiangjw/reference/LDSC_SEG_annotation_geneset/tissue_index.csv",header=F) %>% as.data.frame()
colnames(tissue_index) <- c("tissue_name","index")
tissue_index$tissue_name <- gsub(" ","_",tissue_index$tissue_name)

ldcts_file <- data.frame()
for(i in 1:nrow(tissue_index))
{
	ldcts_file[i,1] <- tissue_index$tissue_name[i]
	ldcts_file[i,2] <- paste0(TissueIndex_dir,"GTEx.",i,".",",",TissueIndex_dir,"GTEx.control.")
}

fwrite(ldcts_file,paste0(Output_dir,"Multi_tissue_gene_expr.ldcts"),sep="\t",col.names=F)
