library(data.table)
library(dplyr)

rm(list = ls())
gc()

EpiData_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/Github_GCLink/01.Phenotypic_Association_Analysis/" #modify into your path

#### logistic regression ####

Epi_AD <- "AD_EUR"
Epi_Allergy <- "AR_EUR"

logistic_regression <- function(anxiety_phenotype,allergy_phenotype,response_status)
{
  logistic_regression_resultsummary <- data.frame()
  temp_anxiety <- fread(paste0(EpiData_dir,anxiety_phenotype,".txt"))
  temp_allergy <- fread(paste0(EpiData_dir,allergy_phenotype,".txt"))
  logistic_data <- merge(x=temp_anxiety,y=temp_allergy,by=c("eid","sex","age"))
  if(response_status == 1)  #response_status == 1 means regarding anxiety as response variable
  {
    logistic_model <- glm(anxiety_status ~ sex + age + allergy_status, data = logistic_data, family = binomial)
  }
  
  if(response_status == 2)  #response_status == 2 means regarding allergic diseases as response variable
  {
    logistic_model <- glm(allergy_status ~ sex + age + anxiety_status, data = logistic_data, family = binomial)
  }
  OR_CI <- data.frame(exp(cbind(OR=coef(logistic_model),confint(logistic_model,level = 0.95))))
  info <- data.frame(coef(summary(logistic_model)))
  colnames(OR_CI) <- c("OR","95%_Lower","95%_Upper")
  logistic_regression_resultsummary[1,1] <- anxiety_phenotype
  logistic_regression_resultsummary[1,2] <- allergy_phenotype
  logistic_regression_resultsummary[1,3] <- OR_CI$OR[4]
  logistic_regression_resultsummary[1,4] <- OR_CI$`95%_Lower`[4]
  logistic_regression_resultsummary[1,5] <- OR_CI$`95%_Upper`[4]
  logistic_regression_resultsummary[1,6] <- info$z.value[4]
  logistic_regression_resultsummary[1,7] <- info$Pr...z..[4]
  return(logistic_regression_resultsummary)
}

####response status = 1 
logistic_responAnxiety_resultsummary_sum <- data.frame()
for(i in Epi_AD)
{
  for(j in Epi_Allergy)
  {
    if(nrow(logistic_responAnxiety_resultsummary_sum) == 0)
    {
      logistic_responAnxiety_resultsummary_sum <- logistic_regression(anxiety_phenotype = i,allergy_phenotype = j,response_status = 1)
    }
    else
    {
      logistic_responAnxiety_resultsummary_sum <- rbind(logistic_responAnxiety_resultsummary_sum,logistic_regression(anxiety_phenotype = i,allergy_phenotype = j,response_status = 1))
    }
  }
}
colnames(logistic_responAnxiety_resultsummary_sum) <- c("anxiety_phenotype","allergy_phenotype","OR(allergy)","95% Lower","95% Upper","z value","Pr(>|z|)")
#fwrite(logistic_responAnxiety_resultsummary_sum,"logistic_responAnxiety_summary_EUR20240820.txt",sep = "\t")

####response status = 2 
logistic_responAllergy_resultsummary_sum <- data.frame()
for(i in Epi_AD)
{
  for(j in Epi_Allergy)
  {
    if(nrow(logistic_responAllergy_resultsummary_sum) == 0)
    {
      logistic_responAllergy_resultsummary_sum <- logistic_regression(anxiety_phenotype = i,allergy_phenotype = j,response_status = 2)
    }
    else
    {
      logistic_responAllergy_resultsummary_sum <- rbind(logistic_responAllergy_resultsummary_sum,logistic_regression(anxiety_phenotype = i,allergy_phenotype = j,response_status = 2))
    }
  }
}
colnames(logistic_responAllergy_resultsummary_sum) <- c("anxiety_phenotype","allergy_phenotype","OR(anxiety)","95% Lower","95% Upper","z value","Pr(>|z|)")
#fwrite(logistic_responAllergy_resultsummary_sum,"logistic_responAllergy_summary.txt",sep = "\t")