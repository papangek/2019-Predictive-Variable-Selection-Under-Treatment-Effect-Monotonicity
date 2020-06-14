##################################################################
########### Example for running various methods on Jobs
##################################################################

rm(list=ls()) # clear the workspace
setwd("your/working/directory")
input_data_link = "jobs/data/directory/nsw.csv"
output_path = "your/output/folder"

##################################
##### Load libraries/sources #####
##################################
source("vt_pred.R")
source("info_theoretic_criteria.R")
library(infotheo)
library(entropy)

# Libraries for VT
library(aVirtualTwins)
library(randomForest, verbose = F)

##############################################
##### Pre-process data   #####
##############################################
data_original = read.csv(input_data_link, header = TRUE)

t = data_original$treat
age = data_original$age
educ = data_original$education
black = data_original$black
hispanic = data_original$hispanic
married = data_original$married
nodegree = data_original$nodegree
re75 = data_original$re75
re78 = data_original$re78
u75 = as.numeric(re75==0)
yf = as.numeric(re78 > re75)

### Equal width discretisation
bins = 5
educ_disc = t(infotheo::discretize( educ, disc="equalwidth", nbins= bins))
age_disc = t(infotheo::discretize( age, disc="equalwidth", nbins= bins))
re75_disc = t(infotheo::discretize( age, disc="equalwidth", nbins= bins))

num_examples = nrow(data_original)
data = matrix(data=NA, nrow = num_examples, ncol = 8)
data[,1] = age_disc
data[,2] = educ_disc
data[,3] = black
data[,4] = hispanic
data[,5] = married
data[,6] = nodegree
data[,7] = u75
data[,8] = re75_disc
data_frame = as.data.frame(data)

num_features = 8
number_datasets = 500

#########################################################################
##### Initial arrays where the rankings for each run will be stored #####
#########################################################################
ranking_infoplus <-array(0,dim=c(number_datasets,num_features))
ranking_pmi <-array(0,dim=c(number_datasets,num_features))
VT_s <-array(0,dim=c(number_datasets,num_features))
ranking_jpmi <-array(0,dim=c(number_datasets,num_features))

#################################################
##### For loop to generate sample datasets  #####
#################################################
for (index_dataset in 1:number_datasets){
  print(sprintf(" %d bootrstrap out of %d", index_dataset,number_datasets))
  
  boots <- sample(num_examples, num_examples, replace=TRUE)
  
  #########################################################
  #### Variable importance using Information-Theoretic ####
  #########################################################
  
  print("Information Theoretic")
  res <- INFOplus(data_frame[boots,],yf[boots],t[boots],num_features)
  ranking_info_plus[exp,] = res$ranking_scores
  res <- PRED_PMI(data_frame[boots,],yf[boots],t[boots],num_features)
  ranking_pmi[exp,] = res$ranking_scores
  res <- PRED_JPMI(data_frame[boots,],yf[boots],t[boots],num_features)
  ranking_jpmi[exp,] = res$ranking_scores
  
  #########################################
  #### Variable importance using VT### ####
  #########################################
  print("Virtual Twins")
  dataset_frame_vt <- data.frame(dataset.labels=yf[boots],dataset.treatment=t[boots],dataset.data = lapply(data_frame[boots,], as.factor))
  VT_s[index_dataset,] <- VT.Predictive.one(dataset_frame_vt)$ranking_scores # Returns the variable importance using VT
}

path = paste(output_path,toString(m),"/INFOplus_obs_",toString(observations[i]),".txt",sep="")
write.table(ranking_info_plus,file=path,row.names=FALSE,col.names=FALSE,sep=',')
path = paste(output_path,toString(m),"/PRED_PMI_obs_",toString(observations[i]),".txt",sep="")
write.table(ranking_pmi,file=path,row.names=FALSE,col.names=FALSE,sep=',')
path = paste(output_path,toString(m),"/PRED_JPMI_obs_",toString(observations[i]),".txt",sep="")
write.table(ranking_jpmi,file=path,row.names=FALSE,col.names=FALSE,sep=',')
path = paste(output_path,toString(m),"/VT_s_obs_",toString(observations[i]),".txt",sep="")
write.table(VT_s,file=path,row.names=FALSE,col.names=FALSE,sep=',')

# Check population CR 
B = 1000
ATE_boot = rep(0,B)
for (b in 1:B){
  boots = sample(num_examples, num_examples, replace=TRUE)
  x_boot = data[boots,]
  yf_boot = yf[boots]
  t_boot = t[boots]
  idx_treated = which(t_boot == 1)
  idx_control = which(t_boot == 0)
  ATE_boot[b] = mean(yf_boot[idx_treated])/mean(yf_boot[idx_control])
}
hist(ATE_boot)
mean(ATE_boot)
quantile(ATE_boot,probs=c(0.025,0.975))

# Check degree subgroups
B = 1000
feat_uq = unique(data[,6])
ATE_degree = matrix(data=NA,nrow = length(feat_uq),ncol= B)
num_examples = nrow(data)
for(b in 1:B){
  boots = sample(num_examples, num_examples, replace=TRUE)
  x_boot = data[boots,]
  yf_boot = yf[boots]
  t_boot = t[boots]
  idx_treated = which(t_boot == 1)
  idx_control = which(t_boot == 0)
for (uq_val in 1:length(feat_uq)){
  idx = which(x_boot[,6]==feat_uq[uq_val])
  ATE_degree[uq_val,b] = mean(yf_boot[intersect(idx_treated,idx)])/mean(yf_boot[intersect(idx_control,idx)])
}
}
mean(ATE_degree[1,])
quantile(ATE_degree[1,],probs=c(0.025,0.975))
mean(ATE_degree[2,])
quantile(ATE_degree[2,],probs=c(0.025,0.975))

# Check degree/education subgroup
B = 1000
ATE_sg = matrix(data=NA,nrow = 2,ncol= B)
for (b in 1:B){
  boots = sample(num_examples, num_examples, replace=TRUE)
  x_boot = data[boots,]
  yf_boot = yf[boots]
  t_boot = t[boots]
  feat1_1 = which(x_boot[,2]==1)
  feat1_4 = which(x_boot[,2]==4)
  feat1_5 = which(x_boot[,2]==5)
  feat2_0 = which(x_boot[,6]==0)
  feat2_1 = which(x_boot[,6]==1)
  idx1 = intersect(feat1_1,feat2_1)
  idx3 = intersect(feat1_4,feat2_1)
  idx = union(idx1,idx3)
  sg1 = idx
  sg2 = setdiff(c(1:num_examples),sg1) 
  idx_treated = which(t_boot == 1)
  idx_control = which(t_boot == 0)
  ATE_sg[1,b] = mean(yf_boot[intersect(idx_treated,sg1)])/mean(yf_boot[intersect(idx_control,sg1)])
  ATE_sg[2,b] = mean(yf_boot[intersect(idx_treated,sg2)])/mean(yf_boot[intersect(idx_control,sg2)])
}
mean(ATE_sg[1,])
quantile(ATE_sg[1,],probs=c(0.025,0.975))
mean(ATE_sg[2,])
quantile(ATE_sg[2,],probs=c(0.025,0.975))

