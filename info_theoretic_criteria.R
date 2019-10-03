EPSILON = 0.001
CONST = -100000

######################################
############## INFO  #################
#####################################
INFO <- function(data,labels,treatment){
  
  num_features <- ncol(data)
  mi_scores <- 	rep(0, num_features)
  for (index_feature in 1:num_features){ 
    mi_scores[index_feature] =  condinformation_normalised(treatment,labels,data[,index_feature])# condinformation(treatment,labels,data[,index_feature],method="shrink")  
  }
  
  sorted_scores <- sort(mi_scores, decreasing=T,method='shell',index.return=TRUE) 
  ranking_scores <- sort(sorted_scores$ix, decreasing=F,method='shell',index.return=TRUE)
  results <- list("scores" = mi_scores, "ranking" = sorted_scores$ix, "ranking_scores" = ranking_scores$ix)
  
  return(results)
  
}
#######################################
############# INFO+ ###################
#######################################
INFOplus <- function(data,labels,treatment, top_k){
  num_features <- ncol(data)
  mi_scores <- 	rep(0, num_features)
  ranking_scores <- 	rep(0, num_features)
  ranking <- 	rep(0, num_features)
  selected_features <- 0   
  VT.First <-  INFO(data,labels,treatment)
  selected_features[1]<-VT.First$ranking[1]
  ranking_scores[selected_features[1]] <- 1
  mi_scores[selected_features[1]] <- VT.First$scores[VT.First$ranking[1]]
  
  if (top_k == 1){
    ranking_scores[ranking_scores==0] <- (top_k+1):num_features
    results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores)
    
    return(results)
  }
  
  not_selected_features <- setdiff(1:num_features,selected_features)
  score_per_feature <- array(0,dim=c(1,num_features))
  
  score_per_feature[selected_features[1]]<-NA
  count_cmi <- num_features
  for (count in 2:top_k){
    for (index_feature_ns in 1:length(not_selected_features)){
      conditioning_features <- do.call(interaction,data[,c(not_selected_features[index_feature_ns], selected_features[count-1])])
      score_per_feature[not_selected_features[index_feature_ns]] <-  score_per_feature[not_selected_features[index_feature_ns]] + condinformation_normalised(treatment,labels,conditioning_features)# condinformation(treatment,labels,conditioning_features,method="shrink") 
      count_cmi <- count_cmi+1
      
    }
    
    selected_features[count] <- which.max(score_per_feature) 
    ranking_scores[selected_features[count]] <- count
    
    mi_scores[selected_features[count]] <-  score_per_feature[selected_features[count]]
    score_per_feature[selected_features[count]]<-NA
    not_selected_features <- setdiff(1:num_features,selected_features)
    
  }
  
  ranking_scores[ranking_scores==0] <- (top_k+1):num_features
  results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores,"count_cmi" = count_cmi)
  
  return(results)
}

condinformation_normalised <- function(treatment, labels, features)
{
  
  cmi_normalised <- infotheo::condinformation(treatment,labels,features,method="shrink")   / sqrt(infotheo::condentropy(treatment, features, method="shrink")*infotheo::condentropy(labels, features, method="shrink"))
  return(cmi_normalised)
}

mutualinformation_normalised <- function(labels, features)
{
  mi_normalised <- infotheo::mutinformation(labels,features,method = "shrink") / sqrt(infotheo::entropy(features,method="shrink")*infotheo::entropy(labels,method="shrink"))
  return(mi_normalised)
}

############################################################
##### Calculates the average PMI for each x ################
############################################################
APMI <- function(features,labels,treatment){
  
  N = length(labels)
  
  ### p(T)
  freq_T = table(treatment)
  p_T = freqs.shrink(freq_T,verbose=FALSE)
  
  ### p(X)
  freq_X = table(features)
  p_X = freqs.shrink(freq_X,verbose=FALSE)
  
  ### p(Y,T)
  YT = interaction(labels,treatment)
  freq_YT = table(YT)
  p_YT = freqs.shrink(freq_YT,verbose = FALSE)
  
  ### p(X,T)
  XT = interaction(features,treatment)
  freq_XT = table(XT)
  p_XT = freqs.shrink(freq_XT,verbose = FALSE)
  
  ### p(Y,X,T)
  YXT = interaction(labels,features,treatment)
  freq_YXT = table(YXT)
  p_YXT = freqs.shrink(freq_YXT,verbose = FALSE)
  
  unique_X = unique(features)
  unique_Y = unique(labels)
  unique_T = unique(treatment)
  
  p_Y_given_X_Tis1 = matrix(, nrow = length(unique_Y), ncol = length(unique_X))
  rownames(p_Y_given_X_Tis1) = unique_Y
  colnames(p_Y_given_X_Tis1) = unique_X
  p_Y_given_X_Tis0 = matrix(, nrow = length(unique_Y), ncol = length(unique_X))
  rownames(p_Y_given_X_Tis0) = unique_Y
  colnames(p_Y_given_X_Tis0) = unique_X
  p_YX_Tis1 = matrix(, nrow = length(unique_Y), ncol = length(unique_X))
  rownames(p_YX_Tis1) = unique_Y
  colnames(p_YX_Tis1) = unique_X
  p_YX_Tis0 = matrix(, nrow = length(unique_Y), ncol = length(unique_X))
  rownames(p_YX_Tis0) = unique_Y
  colnames(p_YX_Tis0) = unique_X
  
  ratio_outcome_1 = c()
  ratio_outcome_0 = c()
  test = c()
  for (j in unique_X){
    p_Y_given_X_Tis1["1",toString(j)] = p_YXT[paste("1",toString(j),'1',sep=".")]/p_XT[paste(toString(j),'1',sep=".")] 
    p_Y_given_X_Tis0["1",toString(j)] = p_YXT[paste("1",toString(j),'0',sep=".")]/p_XT[paste(toString(j),'0',sep=".")]
    p_YX_Tis1["1",toString(j)] = p_Y_given_X_Tis1["1",toString(j)]*(p_XT[paste(toString(j),'1',sep=".")]/p_T['1'])
    p_YX_Tis0["1",toString(j)] = p_Y_given_X_Tis0["1",toString(j)]*(p_XT[paste(toString(j),'0',sep=".")]/p_T['0'])
    ratio_outcome_1[toString(j)] = (p_YX_Tis1['1',toString(j)]/(p_YT[paste('1','1',sep=".")]/p_T["1"]))*log(p_Y_given_X_Tis1['1',toString(j)]/(p_YT[paste('1','1',sep=".")]/p_T["1"])) 
    ratio_outcome_0[toString(j)] = (p_YX_Tis0['1',toString(j)]/(p_YT[paste('1','0',sep=".")]/p_T["0"]))*log(p_Y_given_X_Tis0['1',toString(j)]/(p_YT[paste('1','0',sep=".")]/p_T["0"])) 
    test[toString(j)] = abs(p_Y_given_X_Tis1['1',toString(j)] - p_Y_given_X_Tis0['1',toString(j)] - (p_YT[paste('1','1',sep=".")]/p_T["1"] - p_YT[paste('1','0',sep=".")]/p_T["0"]))
  }
  
  use = T
  if (sum(test) < EPSILON)
    use = F
  
  results <- list("use" = use,"ratio_outcome_1" = ratio_outcome_1,"ratio_outcome_0" = ratio_outcome_0)
  
  return(results) 
}

PRED_PMIM <- function(data,labels,treatment){
  
  num_features <- ncol(data)
  mi_scores <- 	rep(0, num_features)
  mi_scores_corrected <- rep(0, num_features)
  use_in_sum <- rep(T, num_features)
  for (index_feature in 1:num_features){ 
    results_temp = APMI(data[,index_feature],labels,treatment)
    score_temp = sum(results_temp$ratio_outcome_1) - sum(results_temp$ratio_outcome_0)
    mi_scores[index_feature] =  score_temp
    mi_scores_corrected[index_feature] = score_temp
    if (results_temp$use == F){
      mi_scores_corrected[index_feature] = CONST 
      use_in_sum[index_feature] = F
    }
  }

  sorted_scores <- sort(mi_scores_corrected, decreasing=T,method='shell',index.return=TRUE) 
  ranking_scores <- sort(sorted_scores$ix, decreasing=F,method='shell',index.return=TRUE)
  results <- list("use" = use_in_sum, "scores" = mi_scores, "ranking" = sorted_scores$ix, "ranking_scores" = ranking_scores$ix)
  
  return(results)
}

PRED_PMI <- function(data,labels,treatment,top_k){
  
  num_features <- ncol(data)
  mi_scores <- 	rep(0, num_features)
  ranking_scores <- 	rep(NA, num_features)
  ranking <- 	rep(0, num_features)
  selected_features <- 0   
  INFOdiff.First <-  PRED_PMIM(data,labels,treatment)
  selected_features[1]<-INFOdiff.First$ranking[1]
  ranking_scores[selected_features[1]] <- 1
  mi_scores[selected_features[1]] <- INFOdiff.First$scores[INFOdiff.First$ranking[1]]
  
  if (top_k == 1){
    ranking_scores[is.na(ranking_scores)] <- (top_k+1):num_features
    results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores)
    return(results)
  }
  
  not_selected_features <- setdiff(1:num_features,selected_features)
  score_per_feature <- array(NA,dim=c(1,num_features))
  score_per_feature[selected_features[1]]<-NA
  
  for (count in 2:top_k){
    S = length(selected_features)
    for (index_feature_ns in 1:length(not_selected_features)){
      if (INFOdiff.First$use==T){
        score_per_feature[not_selected_features[index_feature_ns]] = (1-S)*INFOdiff.First$scores[not_selected_features[index_feature_ns]]
      }
      else{
        score_per_feature[not_selected_features[index_feature_ns]] = 0
      }
      counter = 0
      for (index_feature_s in 1:length(selected_features)){
        conditioning_features <- do.call(interaction,data[,c(not_selected_features[index_feature_ns], selected_features[index_feature_s])])
        results_temp = APMI(conditioning_features,labels,treatment)
        score_temp = sum(results_temp$ratio_outcome_1) - sum(results_temp$ratio_outcome_0) 
        if (results_temp$use == F){ #there is no treatment difference - do not account this in the score
          score_temp = 0
        }
        score_per_feature[not_selected_features[index_feature_ns]] <-  score_per_feature[not_selected_features[index_feature_ns]] + score_temp 
      }
      
    }
    
    selected_features[count] <- which.max(score_per_feature) 
    ranking_scores[selected_features[count]] <- count
    mi_scores[selected_features[count]] <-  score_per_feature[selected_features[count]]
    score_per_feature[selected_features[count]]<-NA
    not_selected_features <- setdiff(1:num_features,selected_features)
  }
  
  ranking_scores[is.na(ranking_scores)] <- (top_k+1):num_features
  results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores)
  
  return(results)
}

PRED_JPMI <- function(data,labels,treatment, top_k){
  
  num_features <- ncol(data)
  mi_scores <- 	rep(0, num_features)
  ranking_scores <- 	rep(NA, num_features)
  ranking <- 	rep(0, num_features)
  selected_features <- 0   
  INFOdiff.First <-  PRED_PMIM(data,labels,treatment)
  selected_features[1]<-INFOdiff.First$ranking[1]
  ranking_scores[selected_features[1]] <- 1
  mi_scores[selected_features[1]] <- INFOdiff.First$scores[INFOdiff.First$ranking[1]]
  
  if (top_k == 1){
    ranking_scores[is.na(ranking_scores)] <- (top_k+1):num_features
    results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores)
    
    return(results)
  }
  
  not_selected_features <- setdiff(1:num_features,selected_features)
  score_per_feature <- array(NA,dim=c(1,num_features))
  score_per_feature[selected_features[1]]<-NA
  
  for (count in 2:top_k){
    S = length(selected_features)
    for (index_feature_ns in 1:length(not_selected_features)){
      score_per_feature[not_selected_features[index_feature_ns]] = 0
      counter = 0
      for (index_feature_s in 1:length(selected_features)){
        conditioning_features <- do.call(interaction,data[,c(not_selected_features[index_feature_ns], selected_features[index_feature_s])])
        results_temp = APMI(conditioning_features,labels,treatment)
        score_temp = sum(results_temp$ratio_outcome_1) - sum(results_temp$ratio_outcome_0) 
        if (results_temp$use == F){ #there is no treatment difference - do not account this in the score
          score_temp = 0
        }
        score_per_feature[not_selected_features[index_feature_ns]] <-  score_per_feature[not_selected_features[index_feature_ns]] + score_temp 
      }
      
    }
    
    selected_features[count] <- which.max(score_per_feature) 
    ranking_scores[selected_features[count]] <- count
    mi_scores[selected_features[count]] <-  score_per_feature[selected_features[count]]
    score_per_feature[selected_features[count]]<-NA
    not_selected_features <- setdiff(1:num_features,selected_features)
  }
  
  ranking_scores[is.na(ranking_scores)] <- (top_k+1):num_features
  results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores)
  
  return(results)
}

PRED_JPMI.efficient <- function(data,labels,treatment, top_k){
  
  num_features <- ncol(data)
  mi_scores <- 	rep(0, num_features)
  ranking_scores <- 	rep(0, num_features)
  ranking <- 	rep(0, num_features)
  selected_features <- 0   
  INFOdiff.First <-  PRED_PMIM(data,labels,treatment)
  selected_features[1]<-INFOdiff.First$ranking[1]
  ranking_scores[selected_features[1]] <- 1
  mi_scores[selected_features[1]] <- INFOdiff.First$scores[INFOdiff.First$ranking[1]]
  
  if (top_k == 1){
    ranking_scores[ranking_scores==0] <- (top_k+1):num_features
    results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores)
    
    return(results)
  }
  
  not_selected_features <- setdiff(1:num_features,selected_features)
  score_per_feature <- array(0,dim=c(1,num_features))
  score_per_feature[selected_features[1]]<-NA
  
  for (count in 2:top_k){
    S = length(selected_features)
    for (index_feature_ns in 1:length(not_selected_features)){
      conditioning_features <- do.call(interaction,data[,c(not_selected_features[index_feature_ns], selected_features[count-1])])
      results_temp = APMI(conditioning_features,labels,treatment)
      score_temp = sum(results_temp$ratio_outcome_1) - sum(results_temp$ratio_outcome_0)
      if (results_temp$use == F){ #there is no treatment difference - do not account this in the score
        score_temp = 0
      }
      score_per_feature[not_selected_features[index_feature_ns]] <-  score_per_feature[not_selected_features[index_feature_ns]] + score_temp
    }
    
    selected_features[count] <- which.max(score_per_feature) 
    ranking_scores[selected_features[count]] <- count
    mi_scores[selected_features[count]] <-  score_per_feature[selected_features[count]]
    score_per_feature[selected_features[count]]<-NA
    not_selected_features <- setdiff(1:num_features,selected_features)
  }
  
  ranking_scores[ranking_scores==0] <- (top_k+1):num_features
  results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores)
  
  return(results)
}

PRED_PMI.efficient <- function(data,labels,treatment, top_k){
  
  num_features <- ncol(data)
  mi_scores <- 	rep(0, num_features)
  ranking_scores <- 	rep(0, num_features)
  ranking <- 	rep(0, num_features)
  selected_features <- 0   
  
  INFOdiff.First <-  PRED_PMIM(data,labels,treatment)
  
  selected_features[1]<-INFOdiff.First$ranking[1]
  ranking_scores[selected_features[1]] <- 1
  mi_scores[selected_features[1]] <- INFOdiff.First$scores[INFOdiff.First$ranking[1]]
  
  if (top_k == 1){
    ranking_scores[ranking_scores==0] <- (top_k+1):num_features
    results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores)
    
    return(results)
  }
  
  not_selected_features <- setdiff(1:num_features,selected_features)
  score_per_feature <- array(0,dim=c(1,num_features))
  score_per_feature[selected_features[1]]<-NA
  
  for (count in 2:top_k){
    S = length(selected_features)
    for (index_feature_ns in 1:length(not_selected_features)){
      if (INFOdiff.First$use==T & count > 2){
        score_per_feature[not_selected_features[index_feature_ns]] =  score_per_feature[not_selected_features[index_feature_ns]] - (3-count)*INFOdiff.First$scores[not_selected_features[index_feature_ns]] + (2-count)*INFOdiff.First$scores[not_selected_features[index_feature_ns]]
      }
      conditioning_features <- do.call(interaction,data[,c(not_selected_features[index_feature_ns], selected_features[count-1])])
      results_temp = APMI(conditioning_features,labels,treatment)
      score_temp = sum(results_temp$ratio_outcome_1) - sum(results_temp$ratio_outcome_0)
      if (results_temp$use == F){ 
        score_temp = 0
      }
      score_per_feature[not_selected_features[index_feature_ns]] <-  score_per_feature[not_selected_features[index_feature_ns]] + score_temp
    }
    
    selected_features[count] <- which.max(score_per_feature) 
    ranking_scores[selected_features[count]] <- count
    mi_scores[selected_features[count]] <-  score_per_feature[selected_features[count]]
    score_per_feature[selected_features[count]]<-NA
    not_selected_features <- setdiff(1:num_features,selected_features)
  }
  
  ranking_scores[ranking_scores==0] <- (top_k+1):num_features
  results <- list("scores" = mi_scores, "ranking" = selected_features, "ranking_scores" = ranking_scores)
  
  return(results)
}