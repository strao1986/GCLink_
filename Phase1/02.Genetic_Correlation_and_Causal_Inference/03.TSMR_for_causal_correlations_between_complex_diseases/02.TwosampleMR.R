library(data.table)
library(TwoSampleMR)
library(mr.raps)

#### TSMR function ####
mr_pip <- function(exposure, outcome, iv_plink, exp_trait, out_trait, pvalue, respath,samplesize_exp,samplesize_out){
  
  tmp_exp <- fread(exposure) %>% as.data.frame()  
  tmp_exp$F <- (tmp_exp$beta^2)/(tmp_exp$se^2)
  tmp_exp <- subset(tmp_exp,F>10)
  tmp_out <- fread(outcome) %>% as.data.frame() 
  tmp_iv <- fread(iv_plink) %>% as.data.frame() 
  
  tmp_exp1 <- tmp_exp[which(tmp_exp$p < pvalue), ]   
  ### set result path ###
  tmp_respath <- paste0(respath, '/', exp_trait, '-', out_trait) 
  dir.create(tmp_respath, recursive = T)
  ### input exposure data ###
  exp_dat <- format_data(    
    type = 'exposure',       
    tmp_exp1,
    snp_col = "rsid",
    chr_col = "chr",
    pos_col = "pos",
    effect_allele_col = "A1",  
    other_allele_col = "A2",  
    #eaf_col = "eaf",
    beta_col = "beta",
    se_col = "se",
    pval_col = "p",
    samplesize_col = "N")
  
  iv <- as.data.frame(tmp_iv$SNP)  
  colnames(iv) <- 'SNP'   
  
  exp_iv <- merge(iv, exp_dat, by = "SNP")   
  exp_iv$exposure <- exp_trait  
  
  ### write down IV information ###
  write.csv(exp_iv, paste0(tmp_respath, '/', "clump^exp_clean.csv"), row.names = F)  
  
  ### get IVs in  outcome data ###
  merge <- merge(exp_iv, tmp_out, by.x = "SNP", by.y = "rsid")
  
  if (nrow(merge) != 0) {
    out_dat <- format_data(   
      type = "outcome",        
      merge,
      snp_col = "SNP",     
      chr_col = "chr",
      pos_col = "pos",
      effect_allele_col = "A1",
      other_allele_col = "A2",  
      #eaf_col = "eaf",
      beta_col = "beta",
      se_col = "se",
      pval_col = "p",
      samplesize_col = "N")
    
    out_dat$outcome <- out_trait 
    dat_raw <- harmonise_data(exp_iv, out_dat, action = 2)   
    dat <- as.data.frame(dat_raw[which(dat_raw$mr_keep == "TRUE"),]) 
    
    directionality_test <- directionality_test(dat) 
    
    if (length(rownames(dat)) != 0){
      ### perform MR analysis ###
      mr_results <- mr(dat, method_list = c('mr_wald_ratio', 'mr_egger_regression',
                                            'mr_weighted_median', 'mr_ivw', 'mr_simple_mode',
                                            'mr_weighted_mode', 'mr_ivw_fe'))  
      write.csv(mr_results, paste0(tmp_respath, "/", "MRresult_plot.csv"), row.names = F)
      write.csv(dat, paste0(tmp_respath, "/", "dat_plot.csv"), row.names = F)
      
      mr_raps <- as.data.frame(mr.raps( dat$beta.exposure, dat$beta.outcome, 
                                        dat$se.exposure, dat$se.outcome))		
      
      ### add cols so that we can combine mr & raps result ###
      mr_results$naive.se <- "NA"
      mr_results$chi.sq.test <- "NA"
      
      mr_raps$outcome <- out_trait
      mr_raps$exposure <- exp_trait
      mr_raps$id.exposure <- dat[1,which(colnames(dat) == "id.exposure")]
      mr_raps$id.outcome <- dat[1,which(colnames(dat) == "id.outcome")]
      mr_raps$method <- "Robust Adjusted Profile Score"
      mr_raps$nsnp <- length(which(dat$mr_keep) == "TRUE")
      
      ### correct col locations by mr_results ###
      mr_raps <- mr_raps[,c(8, 9, 6, 7, 10, 11, 1, 2, 3, 4, 5)]  
      colnames(mr_raps) <- c("id.exposure", "id.outcome",
                             "outcome", "exposure", "method",
                             "nsnp", "b", "se", "pval",
                             "naive.se", "chi.sq.test")
      mr_results <- rbind(mr_results, mr_raps) 
      mr_results <- generate_odds_ratios(mr_results)
      
      ### single snp analysis ###
      res_single <- mr_singlesnp(dat) 
      write.csv(res_single, paste0(tmp_respath, "/", "res_single_plot.csv"), row.names = F)
      ### output result ###
      write.csv(out_dat, paste0(tmp_respath, "/", "outcome.csv"), row.names = F)
      write.csv(dat_raw, paste0(tmp_respath, "/", "harmonize.csv"), row.names = F)
      write.csv(mr_results, paste0(tmp_respath, "/", "MRresult.csv"), row.names = F)
      
      ### sensitive analysis ###
      if (length(which(dat$mr_keep) == "TRUE") > 1){
        heter_dat <- mr_heterogeneity(dat)  
        pleio_dat <- mr_pleiotropy_test(dat) 
        leave_one_dat <- mr_leaveoneout(dat) 
        
        write.csv(heter_dat, paste0(tmp_respath, "/", "heterogeneity.csv"), row.names = F)
        write.csv(pleio_dat, paste0(tmp_respath, "/", "pleiotropy.csv"), row.names = F)
        write.csv(leave_one_dat, paste0(tmp_respath, "/", "leave_one_out.csv"), row.names = F)
        write.csv(res_single, paste0(tmp_respath, "/", "single_snp.result.csv"), row.names = F)
        write.csv(directionality_test, paste0(tmp_respath, "/", "directionality_test.csv"), row.names = F)
      }
    }
  }						
}

#### call funcyion ####
clump_result_dir <- "/public/jiangjw/GWAS_sumstats/EUR/clump_result/01.F10/" #modify into your path
GWAS_sumstats_dir <- "/public/jiangjw/GWAS_sumstats/EUR/" #modify into your path
forward_MRresult_dir <- "/public/jiangjw/02.anxiety_ADs_EUR/01.MR/Github_code_test/forward_MRresult/" #modify into your path
reverse_MRresult_dir <- "/public/jiangjw/02.anxiety_ADs_EUR/01.MR/Github_code_test/reverse_MRresult/" #modify into your path

### forward MR analysis (regarding AR as an exposure and AD as an outcome) ###
exp_pheno <- "DiagnosedAR"
out_pheno <- "Anxiety2021"

for(i in 1:length(exp_pheno))
{
	for(j in 1:length(out_pheno))
	{
		        tmp_exp <- paste0(GWAS_sumstats_dir,exp_pheno[i],".clean.txt")
        		tmp_out <- paste0(GWAS_sumstats_dir,out_pheno[j],".clean.txt")
			tmp_iv <- paste0(clump_result_dir,exp_pheno[i],".LD_result.clumped")
			mr_pip(exposure = tmp_exp,    
         			outcome = tmp_out,  
         			iv_plink = tmp_iv,
         			exp_trait = exp_pheno[i],
         			out_trait = out_pheno[j],
         			pvalue = 1,
         			respath = forward_MRresult_dir)
	}
}

### reverse MR analysis (regarding AD as an exposure and AR as an outcome) ###
exp_pheno <- "Anxiety2021"
out_pheno <- "DiagnosedAR"

for(i in 1:length(exp_pheno))
{
        for(j in 1:length(out_pheno))
        {               
                        tmp_exp <- paste0(GWAS_sumstats_dir,exp_pheno[i],".clean.txt")
                        tmp_out <- paste0(GWAS_sumstats_dir,out_pheno[j],".clean.txt")
                        tmp_iv <- paste0(clump_result_dir,exp_pheno[i],".LD_result.clumped")
                        mr_pip(exposure = tmp_exp,
                                outcome = tmp_out,     
                                iv_plink = tmp_iv,     
                                exp_trait = exp_pheno[i],
                                out_trait = out_pheno[j],
                                pvalue = 1,
                                respath = reverse_MRresult_dir)
        }
}


