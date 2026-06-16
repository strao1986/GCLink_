library(data.table)
library(dplyr)
library(TwoSampleMR)

#### MR results summary ####
clump_result_dir <- "/public/jiangjw/GWAS_sumstats/EUR/clump_result/02.p5e6/" #modify into your path
MRsummaryFunc_dir <- "/public/jiangjw/02.anxiety_ADs_EUR/01.MR/Github_code_test/"
source(paste0(MRsummaryFunc_dir,"MRsummary.func.R")) #the most suitable MR method selection function

### forward MR results summary ###
result_dir <- "/public/jiangjw/02.anxiety_ADs_EUR/01.MR/Github_code_test/forward_MRresult/" #modify into your path

exp_pheno <- c("DiagnosedAR","Anxiety2021")
out_pheno <- c("BASO","CD19","CD3","CD4CD8ratio","CD4","CD56","CD8","EOS","HB","HT","LYMPH","MCHC","MCH","MCV","MONO","NEUT","PLT","RBC","WBC")

exptran <- exp_pheno
outtran <- out_pheno
## summary of significant SNPs and employed SNPs in MR analysis (MRA) ##
info_clumped <- data.frame()
for(i in 1:length(exptran))
{
  clumped_snp <- fread(paste0(result_dir,exptran[i],"-",outtran[1],"/clump^exp_clean.csv"))
  info_clumped[i,1] <- exptran[i]
  info_clumped[i,2] <- nrow(clumped_snp)
}
colnames(info_clumped) <- c("exposure","clumped_snp")

## summary of sensitive analysis results (heterogeneity, pleiotropy and directionality test) ##
info_HPD <- data.frame()
for(i in 1:length(exptran))
{
  for(j in 1:length(outtran))
  {
  	if(length(list.files(paste0(result_dir,exptran[i],"-",outtran[j]),pattern = "heterogeneity")) == 0)
  	{
  		next
  	}
  	mr_result <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/MRresult.csv"))
  	heterogeneity <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/heterogeneity.csv"))
  	pleiotropy <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/pleiotropy.csv"))
  	directionality <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/directionality_test.csv"))
  	temp <- MR_summary(mr_result,heterogeneity,pleiotropy,directionality)
  	temp1 <- data.frame(temp$summary)

  	if(nrow(info_HPD) == 0)
  	{
    		info_HPD <- temp1
  	}
  	else
  	{
    		info_HPD <- rbind(info_HPD,temp1)
  	}


  }
}
## summary results when number of iv is 1 ##
for(i in 1:length(exptran))
{
  for(j in 1:length(outtran))
  {
  	if(length(list.files(paste0(result_dir,exptran[i],"-",outtran[j]),pattern = "heterogeneity")) == 0 && length(list.files(paste0(datapath,exptran[i],"-",outtran[j]),pattern = "MRresult.csv")) == 1)
  	{
    		mr_result <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/MRresult.csv"))
    		Wald <- subset(mr_result,method == "Wald ratio")
    		temp2 <- data.frame("exposure"=Wald$exposure,"outcome"=Wald$outcome,"No..of.SNPs.in.MRA"=Wald$nsnp,"Beta"=Wald$b,"Se"=Wald$se,"MR.pvalue"=Wald$pval,"Method"=Wald$method,"heterogeneity"="-","pleiotropy"="-","directionality"="-","OR"=Wald$or,"LCI95"=Wald$or_lci95,"UCI95"=Wald$or_uci95)
    		info_HPD <- rbind(info_HPD,temp2)
  	}	
  }
}

info_full <- merge(info_clumped,info_HPD,by='exposure')
info_full$beta_se <- paste0(sprintf("%0.3f", info_full$Beta),"(",sprintf("%0.3f", info_full$Se),")")
info_full$or_95CI <- paste0(sprintf("%0.3f", info_full$OR),"(",sprintf("%0.3f", info_full$LCI95),",",sprintf("%0.3f", info_full$UCI95),")")
info_full <- select(info_full,exposure,outcome,clumped_snp,No..of.SNPs.in.MRA,MR.pvalue,beta_se,Method,heterogeneity,pleiotropy,directionality,everything())
info_full <- arrange(info_full,exposure,outcome)

fwrite(info_full,paste0(result_dir,"MR_ResultsSummary.csv"))
fwrite(subset(info_full,MR.pvalue<0.05),paste0(result_dir,"MR_ResultsSummary_p005.csv"))

## MR results and sensitive analysis results visualization ##
for(i in 1:length(exptran))
{
	for(j in 1:length(outtran))
	{
		result_path_plot <- paste0(result_dir,as.character(exptran[i]),"-",outtran[j])
  		setwd(result_path_plot)
  		dat_plot <- read.csv("dat_plot.csv")
  		mr_result_plot <- read.csv("MRresult_plot.csv")
  		res_single_plot <- read.csv("res_single_plot.csv")
  		leave_one_out_plot <- read.csv("leave_one_out.csv")

		# scatter plot #
		tmp_method_choose <- subset(info_full,exposure==exptran[i] & outcome==outtran[j])
		if(nrow(tmp_method_choose) == 2)
                {
                        pdf(file=paste0("dot_",tmp_method_choose$Method[1],".pdf"))
			print(mr_scatter_plot(subset(mr_result_plot,method == tmp_method_choose$Method[1]),dat_plot))
			dev.off()

			pdf(file=paste0("dot_",tmp_method_choose$Method[2],".pdf"))
                        print(mr_scatter_plot(subset(mr_result_plot,method == tmp_method_choose$Method[2]),dat_plot))
			dev.off()
                }
                else
                {
                        pdf(file='dot.pdf')
			print(mr_scatter_plot(subset(mr_result_plot,method == tmp_method_choose$Method),dat_plot))
			dev.off()
                }

		# forest plot #
  		pdf(file='forest.pdf', height=15,width=8)
  		print(mr_forest_plot(res_single_plot))
  		dev.off()
		
		# leave one out testing visualization #
  		pdf(file='leave_one_out.pdf', height=15,width=8)
  		print(mr_leaveoneout_plot(leave_one_out_plot))
  		dev.off()

		# funnel plot #
  		pdf(file='funnel.pdf', height=6,width=8)
  		print(mr_funnel_plot(res_single_plot))
  		dev.off()
	}
}


### reverse MR results summary ###
result_dir <- "/public/jiangjw/02.anxiety_ADs_EUR/01.MR/Github_code_test/reverse_MRresult/" #modify into your path

exp_pheno <- c("BASO","CD19","CD3","CD4CD8ratio","CD4","CD56","CD8","EOS","HB","HT","LYMPH","MCHC","MCH","MCV","MONO","NEUT","PLT","RBC","WBC")
out_pheno <- c("DiagnosedAR","Anxiety2021")

exptran <- exp_pheno
outtran <- out_pheno
## summary of significant SNPs and employed SNPs in MR analysis (MRA) ##
info_clumped <- data.frame()
for(i in 1:length(exptran))
{
  clumped_snp <- fread(paste0(result_dir,exptran[i],"-",outtran[1],"/clump^exp_clean.csv"))
  info_clumped[i,1] <- exptran[i]
  info_clumped[i,2] <- nrow(clumped_snp)
}
colnames(info_clumped) <- c("exposure","clumped_snp")

## summary of sensitive analysis results (heterogeneity, pleiotropy and directionality test) ##
info_HPD <- data.frame()
for(i in 1:length(exptran))
{
  for(j in 1:length(outtran))
  {
        if(length(list.files(paste0(result_dir,exptran[i],"-",outtran[j]),pattern = "heterogeneity")) == 0)
        {
                next
        }
        mr_result <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/MRresult.csv"))
        heterogeneity <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/heterogeneity.csv"))
        pleiotropy <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/pleiotropy.csv"))
        directionality <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/directionality_test.csv"))
        temp <- MR_summary(mr_result,heterogeneity,pleiotropy,directionality)
        temp1 <- data.frame(temp$summary)

        if(nrow(info_HPD) == 0)
        {
                info_HPD <- temp1
        }
        else
        {
                info_HPD <- rbind(info_HPD,temp1)
        }


  }
}

## summary results when number of iv is 1 ##
for(i in 1:length(exptran))
{
  for(j in 1:length(outtran))
  {
        if(length(list.files(paste0(result_dir,exptran[i],"-",outtran[j]),pattern = "heterogeneity")) == 0 && length(list.files(paste0(datapath,exptran[i],"-",outtran[j]),pattern = "MRresult.csv")) == 1)
        {
                mr_result <- fread(paste0(result_dir,exptran[i],"-",outtran[j],"/MRresult.csv"))
                Wald <- subset(mr_result,method == "Wald ratio")
                temp2 <- data.frame("exposure"=Wald$exposure,"outcome"=Wald$outcome,"No..of.SNPs.in.MRA"=Wald$nsnp,"Beta"=Wald$b,"Se"=Wald$se,"MR.pvalue"=Wald$pval,"Method"=Wald$method,"heterogeneity"="-","pleiotropy"="-","directionality"="-","OR"=Wald$or,"LCI95"=Wald$or_lci95,"UCI95"=Wald$or_uci95)
                info_HPD <- rbind(info_HPD,temp2)
        }
  }
}

info_full <- merge(info_clumped,info_HPD,by='exposure')
info_full$beta_se <- paste0(sprintf("%0.3f", info_full$Beta),"(",sprintf("%0.3f", info_full$Se),")")
info_full$or_95CI <- paste0(sprintf("%0.3f", info_full$OR),"(",sprintf("%0.3f", info_full$LCI95),",",sprintf("%0.3f", info_full$UCI95),")")
info_full <- select(info_full,exposure,outcome,clumped_snp,No..of.SNPs.in.MRA,MR.pvalue,beta_se,Method,heterogeneity,pleiotropy,directionality,everything())
info_full <- arrange(info_full,exposure,outcome)

fwrite(info_full,paste0(result_dir,"MR_ResultsSummary.csv"))
fwrite(subset(info_full,MR.pvalue<0.05),paste0(result_dir,"MR_ResultsSummary_p005.csv"))

## MR results and sensitive analysis results visualization ##
for(i in 1:length(exptran))
{
        for(j in 1:length(outtran))
        {
                result_path_plot <- paste0(result_dir,as.character(exptran[i]),"-",outtran[j])
                setwd(result_path_plot)
                dat_plot <- read.csv("dat_plot.csv")
                mr_result_plot <- read.csv("MRresult_plot.csv")
                res_single_plot <- read.csv("res_single_plot.csv")
                leave_one_out_plot <- read.csv("leave_one_out.csv")

                # scatter plot #
                tmp_method_choose <- subset(info_full,exposure==exptran[i] & outcome==outtran[j])
                if(nrow(tmp_method_choose) == 2)
                {
                        pdf(file=paste0("dot_",tmp_method_choose$Method[1],".pdf"))
                        print(mr_scatter_plot(subset(mr_result_plot,method == tmp_method_choose$Method[1]),dat_plot))
                        dev.off()

                        pdf(file=paste0("dot_",tmp_method_choose$Method[2],".pdf"))
                        print(mr_scatter_plot(subset(mr_result_plot,method == tmp_method_choose$Method[2]),dat_plot))
                        dev.off()
                }
                else
                {
                        pdf(file='dot.pdf')
                        print(mr_scatter_plot(subset(mr_result_plot,method == tmp_method_choose$Method),dat_plot))
                        dev.off()
                }

                # forest plot #
                pdf(file='forest.pdf', height=15,width=8)
                print(mr_forest_plot(res_single_plot))
                dev.off()

                # leave one out testing visualization #
                pdf(file='leave_one_out.pdf', height=15,width=8)
                print(mr_leaveoneout_plot(leave_one_out_plot))
                dev.off()

                # funnel plot #
                pdf(file='funnel.pdf', height=6,width=8)
                print(mr_funnel_plot(res_single_plot))
                dev.off()
        }
}

