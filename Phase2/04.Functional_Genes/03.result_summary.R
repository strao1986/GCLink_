library(data.table)
library(dplyr)

SMR_result_dir <- "/public/jiangjw/02.anxiety_ADs_EUR/09.SMR/Github_code_test/" #modify into your path
ResultSummary_dir <- "/public/jiangjw/02.anxiety_ADs_EUR/09.SMR/Github_code_test/" #modify into your path

pheno1 <- "Anxiety2021"
pheno2 <- "DiagnosedAR"
tissue <- "Brain_Frontal_Cortex_BA9"

result_summary3 <- function(anxiety_dat,allergy_dat,tissue)
{
        anxiety_tmp <- fread(paste0(SMR_result_dir,anxiety_dat,"_",tissue,".smr")) %>% as.data.frame() %>% select(probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,p_HEIDI,nsnp_HEIDI)
        allergy_tmp <- fread(paste0(SMR_result_dir,allergy_dat,"_",tissue,".smr")) %>% as.data.frame() %>% select(probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,p_HEIDI,nsnp_HEIDI)
        
	anxiety_tmp_sig <- anxiety_tmp
        allergy_tmp_sig <- allergy_tmp
	for(i in 5:9)
	{
		colnames(anxiety_tmp_sig)[i] <- paste0(colnames(anxiety_tmp_sig)[i],"_trait1")
		colnames(allergy_tmp_sig)[i] <- paste0(colnames(allergy_tmp_sig)[i],"_trait2")
	}
        if(nrow(anxiety_tmp_sig)==0 | nrow(allergy_tmp_sig)==0)
        {
                return("can not be intersected")
        }
        else
        {
                merge_gene <- merge(anxiety_tmp_sig,allergy_tmp_sig,by=c("probeID","ProbeChr","Gene","Probe_bp"))
        }

        if(nrow(merge_gene) == 0)
        {
                return("can not be intersected")
        }
        else
        {
                merge_gene$tissue <- tissue
                merge_gene$trait1 <- anxiety_dat
                merge_gene$trait2 <- allergy_dat

                merge_gene$FDR_trait1 <- p.adjust(merge_gene$p_SMR_trait1, method = "fdr")
		merge_gene$FDR_trait2 <- p.adjust(merge_gene$p_SMR_trait2, method = "fdr")

                merge_gene <- subset(merge_gene,p_SMR_trait1<0.05 & p_SMR_trait2<0.05 & p_HEIDI_trait1>0.05 & p_HEIDI_trait2>0.05)

                if(nrow(merge_gene) == 0)
                {
                        return("can not be intersected")
                }
                else
                {
                merge_gene <- select(merge_gene,trait1,trait2,probeID,ProbeChr,Gene,Probe_bp,b_SMR_trait1,se_SMR_trait1,p_SMR_trait1,FDR_trait1,p_HEIDI_trait1,nsnp_HEIDI_trait1,b_SMR_trait2,se_SMR_trait2,p_SMR_trait2,FDR_trait2,p_HEIDI_trait2,nsnp_HEIDI_trait2,tissue)
                return(merge_gene)
                }
        }
}



info <- data.frame()
for(k in 1:length(pheno1))
{
	for(i in 1:length(pheno2))
	{
		for(j in 1:length(tissue))
		{
			temp <- result_summary3(pheno1[k],pheno2[i],tissue[j])
			if(class(temp) == "character")
			{
				next
			}
			else
			{
				if(nrow(info) == 0)
				{
					info <- temp
				}
				else
				{
					info <- rbind(info,temp)
				}
			}
		}
	}
}
fwrite(info,paste0(ResultSummary_dir,"SMR_ResultSummary_p005.txt"),sep="\t")
