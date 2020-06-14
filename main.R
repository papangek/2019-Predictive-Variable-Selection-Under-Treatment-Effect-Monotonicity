#################################################
# Example for running info theoretic approaches 
#################################################

setwd("your/working/directory")
output_path = "your/output/folder"

covariates = 50
observations = c(100,250,500,750,1000)
experiments = 200
theta = 1
models = c('M1')

source("generate_data.R")
source("info_theoretic_criteria.R")

library(infotheo)
library(entropy)
library(MASS)

for (m in models){
  print(m)
for (i in 1:length(observations)){
  print(observations[i])
  dataset = Get.Data(observations[i],covariates,theta,2,m)
  K = length(dataset$predictive)
  ranking_info_plus = matrix(data=NA, nrow = experiments, ncol = K)
  ranking_pmi = matrix(data=NA, nrow = experiments, ncol = K)
  ranking_jpmi = matrix(data=NA, nrow = experiments, ncol = K)

  for (exp in 1:experiments){
    print("Experiment...")
    print(exp)
    
    bins = sample(c(2,3,4,5),1)
    
    dataset = Get.Data(observations[i],covariates,theta,bins,m)
    
    res <- INFOplus(dataset$data,dataset$labels,dataset$treatment,K)
    ranking_info_plus[exp,] = res$ranking_scores
    
    res <- PRED_PMI(dataset$data,dataset$labels,dataset$treatment,K)
    ranking_pmi[exp,] = res$ranking_scores
    
    res <- PRED_JPMI(dataset$data,dataset$labels,dataset$treatment,K)
    ranking_jpmi[exp,] = res$ranking_scores
    
  }
  
  path = paste(output_path,toString(m),"/INFOplus_obs_",toString(observations[i]),".txt",sep="")
  write.table(ranking_info_plus,file=path,row.names=FALSE,col.names=FALSE,sep=',')
  
  path = paste(output_path,toString(m),"/PRED_PMI_obs_",toString(observations[i]),".txt",sep="")
  write.table(ranking_pmi,file=path,row.names=FALSE,col.names=FALSE,sep=',')
  
  path = paste(output_path,toString(m),"/PRED_JPMI_obs_",toString(observations[i]),".txt",sep="")
  write.table(ranking_jpmi,file=path,row.names=FALSE,col.names=FALSE,sep=',')
}
}
