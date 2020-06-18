###########################################################
################ GENERATE DATA ############################
###########################################################
library(MASS)

Get.Data <- function(sample_size,num_features,bins,model){
  
  predictive = c()
  prognostic = c()
  
  if(model=='M1')
  { 
    correl = 0.7
    sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
    sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
    sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
    diag(sigma) = 1
    x = mvrnorm(sample_size, rep(0, num_features), sigma)
    logit_y0 = -1 + 0.5*x[,1] + 0.5*x[,2] - 0.5*x[,3] + 0.5*x[,2]*x[,3]
    logit_y1 = logit_y0 + 0.1 + 0.9*(x[,1]>0)*(x[,2]<0)

    treatment <- rbinom(sample_size,1,0.5)
    predictive = c(1,2)
    prognostic = c(3)
    
    for (index_feature in 1:num_features){ 
      x[,index_feature] = t(infotheo::discretize( x[,index_feature], disc="equalwidth", nbins= bins))
    }
  }
  
  if(model=='M2')
  { 
    correl = 0.7
    sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
    sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
    sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
    diag(sigma) = 1
    x = mvrnorm(sample_size, rep(0, num_features), sigma)
    logit_y0 = -1 + 0.5*x[,1] + 0.5*x[,2] - 0.5*x[,3] + 0.5*x[,2]*x[,3]
    logit_y1 = logit_y0 + 0.1 + 0.9*(x[,1]>=-0.545)*(x[,2]<=0.545)
    treatment <- rbinom(sample_size,1,0.5)
    predictive = c(1,2)
    prognostic = c(3)
    
    for (index_feature in 1:num_features){ 
      x[,index_feature] = t(infotheo::discretize( x[,index_feature], disc="equalwidth", nbins= bins))
    }
  }
  
  if (model=='M3')
  {
    correl = 0.3
    sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
    sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
    sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
    diag(sigma) = 1
    x = mvrnorm(sample_size, rep(0, num_features), sigma)
    for (index_feature in 1:num_features){
      x[,index_feature] = 1*x[,index_feature]>0
    }
    
    p0_temp = 0.3 +0.05*x[,1]*x[,2] + 0.05*x[,3]*x[,4] + 0.05*x[,5]*x[,6]
    logit_y0 = log(p0_temp/(1-p0_temp))
    p1_temp = 0.3 + 0.35*x[,1]*x[,2] + 0.25*x[,3]*x[,4] + 0.05*x[,5]*x[,6]
    logit_y1 = log(p1_temp/(1-p1_temp))
    
    treatment <- rbinom(sample_size,1,0.5)
    predictive = c(1,2,3,4)
    prognostic = c(3)
  }
  
  
  if (model=='M4')
  {
    correl = 0.3
    sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
    sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
    sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
    diag(sigma) = 1
    x = mvrnorm(sample_size, rep(0, num_features), sigma)
    for (index_feature in 1:num_features){ 
      x[,index_feature] = 1*x[,index_feature]>0
    }
    
    p0_temp = 0.3 + 0.05*x[,1]*x[,2]*x[,3] + 0.05*x[,4]*x[,5]*x[,6]
    logit_y0 = log(p0_temp/(1-p0_temp))
    p1_temp = 0.3 + 0.35*x[,1]*x[,2]*x[,3] + 0.25*x[,4]*x[,5]*x[,6] 
    logit_y1 = log(p1_temp/(1-p1_temp))
    treatment <- rbinom(sample_size,1,0.5)
    predictive = c(1,2,3,4,5,6)
    prognostic = c(3)
  }
  
  
  if(model=='M5')
  { 
    correl = 0.7
    sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
    sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
    sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
    diag(sigma) = 1
    x = mvrnorm(sample_size, rep(0, num_features), sigma)
    logit_y0 = x[,1] + x[,2] + x[,3] + x[,4] + x[,5] + x[,6] 
    logit_y1 = logit_y0 + 5*(x[,7]>0)*(x[,8]<0)
    treatment <- rbinom(sample_size,1,0.5)
    predictive = c(7,8)
    prognostic = c(1,2,3,4,5,6)
    
    for (index_feature in 1:num_features){ 
      x[,index_feature] = t(infotheo::discretize( x[,index_feature], disc="equalwidth", nbins= bins))
    }
  }
  
  if(model=='M6')
  { 
    correl = 0.7
    sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
    sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
    sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
    diag(sigma) = 1
    x = mvrnorm(sample_size, rep(0, num_features), sigma)
    logit_y0 = x[,1] + x[,2] + x[,3] + x[,4] + x[,5] + x[,6]
    logit_y1 = logit_y0 + 5*(x[,7]>-0.545)*(x[,8]<0.545)*(x[,9]>0) 
    treatment <- rbinom(sample_size,1,0.5)
    predictive = c(7,8,9)
    prognostic = c(1,2,3,4,5,6)
    
    for (index_feature in 1:num_features){ 
      x[,index_feature] = t(infotheo::discretize( x[,index_feature], disc="equalwidth", nbins= bins))
    }
  }
  
  if(model=='M7') # To test biases of INFO+ wrt the treatment assignment
  { 
    correl = 0
    sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
    sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
    sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
    diag(sigma) = 1
    x = mvrnorm(sample_size, rep(0, num_features), sigma)
    logit_y0 = x[,1] + x[,2] + x[,3] + x[,4] + x[,5] + x[,6]
    logit_y1 = logit_y0 + 1*(x[,1]>0)
    logit_pt = 0.5*x[,1]
    pt = 1/(1+exp(-logit_pt))
    treatment <- rbinom(sample_size,1,pt)
    predictive = c(1)
    prognostic = c(1,2,3,4,5,6)
    
    for (index_feature in 1:num_features){ 
      x[,index_feature] = t(infotheo::discretize( x[,index_feature], disc="equalwidth", nbins= bins))
    }
  }
  
  if(model=='M8') # To test biases of INFO+ wrt the treatment assignment
  { 
    correl = 0
    sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
    sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
    sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
    diag(sigma) = 1
    x = mvrnorm(sample_size, rep(0, num_features), sigma)
    logit_y0 = x[,1] + x[,2] + x[,3] + x[,4] + x[,5] + x[,6]
    logit_y1 = logit_y0 + 1*(x[,1]>0)
    logit_pt = 0.5*x[,2]
    pt = 1/(1+exp(-logit_pt))
    treatment <- rbinom(sample_size,1,pt)
    predictive = c(1)
    prognostic = c(1,2,3,4,5,6)
    
    for (index_feature in 1:num_features){ 
      x[,index_feature] = t(infotheo::discretize( x[,index_feature], disc="equalwidth", nbins= bins))
    }
  }
  
  if(model=='M9') # To test whether it distinguishes between positive and negative effects
  { 
  correl = 0.7
  sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
  sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
  sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
  diag(sigma) = 1
  x = mvrnorm(sample_size, rep(0, num_features), sigma)
  logit_y0 = x[,1] + x[,2] + x[,3] + x[,4] + x[,5] + x[,6] + x[,7] + x[,8]
  logit_y1 = logit_y0 -5*(x[,1]>0)*(x[,2]>0) + 3*(x[,3]>0)*(x[,4]>0)
  treatment <- rbinom(sample_size,1,0.5)
  predictive = c(1,2,3,4)
  prognostic = c(5,6,7,8)
  
  for (index_feature in 1:num_features){ 
    x[,index_feature] = t(infotheo::discretize( x[,index_feature], disc="equalwidth", nbins= bins))
  }
  }
  
  if(model=='M10') # To test whether it distinguishes between positive and negative effects
  { 
    correl = 0.7
    sigma = matrix(rep(0, num_features*num_features), nrow = num_features, ncol = num_features)
    sigma[seq(1,num_features,2),seq(1,num_features,2)] <-correl 
    sigma[seq(2,num_features,2),seq(2,num_features,2)] <-correl 
    diag(sigma) = 1
    x = mvrnorm(sample_size, rep(0, num_features), sigma)
    alpha = mvrnorm(sample_size, 0, 0.25)
    beta = mvrnorm(sample_size, 0, 0.25)
    logit_y0 = x[,1] + x[,2] + x[,3] + x[,4] + x[,5] + x[,6] + x[,7] + x[,8]
    logit_y1 = logit_y0 -1*(x[,1]>0)*(x[,2]>0) + 3*(x[,3]>0)*(x[,4]>0)
    treatment <- rbinom(sample_size,1,0.5)
    predictive = c(1,2,3,4)
    prognostic = c(5,6,7,8)
    
    for (index_feature in 1:num_features){ 
      x[,index_feature] = t(infotheo::discretize( x[,index_feature], disc="equalwidth", nbins= bins))
    }
  }
  
  pY0 = 1/(1+exp(-logit_y0))
  pY1 = 1/(1+exp(-logit_y1))  
  y0 =  rbinom(sample_size,1,pY0);
  y1 =  rbinom(sample_size,1,pY1);
  idx_treated = treatment==1
  idx_control = treatment==0
  yf = y1
  yf[idx_control] = y0[idx_control]
  data = as.data.frame.matrix(x)
  
  synthetic_dataset <- c()
  synthetic_dataset$data <- data
  synthetic_dataset$treatment <- treatment
  synthetic_dataset$labels <- yf
  synthetic_dataset$prognostic <- prognostic
  synthetic_dataset$predictive <- predictive  
  synthetic_dataset$y1 <- y1
  synthetic_dataset$y0 <- y0
  
  return(synthetic_dataset)
  
}
