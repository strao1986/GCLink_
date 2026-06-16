library(data.table)
library(dplyr)

LDSC_SEG_Results_dir <- "/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC/LDSC_SEG/Github_code_test/" #modify into your path
phenos <- c("Anxiety2021","DiagnosedAR")

result_summary <- function(phenos)
{
	tmp <- fread(paste0(LDSC_SEG_Results_dir,phenos,"_multitissue_gene_expr.cell_type_results.txt"),header=T) %>% as.data.frame()
	tmp <- arrange(tmp,Coefficient_P_value)
	tmp$fdr_adj_Pvalue <- p.adjust(tmp$Coefficient_P_value,method="fdr")
	tmp$trait <- phenos
	tmp <- select(tmp,trait,Name,Coefficient,Coefficient_std_error,Coefficient_P_value,fdr_adj_Pvalue)
	colnames(tmp)[2] <- "tissue"
	return(tmp)
}

info <- data.frame()
for(i in 1:length(phenos))
{
	if(i == 1)
	{
		info <- result_summary(phenos[i])	
	}
	else
	{
		info <- rbind(info,result_summary(phenos[i]))
	}
}
info_p005 <- subset(info,Coefficient_P_value < 0.05)

fwrite(info,"results_summary_full.txt",sep="\t")
fwrite(info_p005,"results_summary_p005.txt",sep="\t")
