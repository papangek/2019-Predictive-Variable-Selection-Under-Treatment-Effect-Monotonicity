#############################################
################ Virtual Twins ##############
#############################################
library(aVirtualTwins)
library(randomForest, verbose = F)

VT.Predictive <- function(dataset_frame){
  # Uses T-learner
  vt.o <- vt.data(dataset_frame, "dataset.labels", "dataset.treatment", interactions = T)
  model.rf.trt1 <- randomForest(x = vt.o$getX(trt = 1),
                                y = vt.o$getY(trt = 1),ntree=1000)
  model.rf.trt0 <- randomForest(x = vt.o$getX(trt = 0),
                                y = vt.o$getY(trt = 0),ntree=1000)
  vt.doublef.rf <- vt.forest("double",
                             vt.data = vt.o, 
                             model_trt1 = model.rf.trt1, 
                             model_trt0 = model.rf.trt0)
  
  z <- vt.doublef.rf$difft
  model.rf.z <- randomForest(dataset_frame[,3:length(dataset_frame) ], y = z,ntree=1000)
  sorted_scores <- sort(model.rf.z$importance, decreasing=T,method='shell',index.return=TRUE) 
  ranking_scores <- sort(sorted_scores$ix, decreasing=F,method='shell',index.return=TRUE)
  results <- list("scores" = as.vector(t(model.rf.z$importance)), "ranking" = sorted_scores$ix, "ranking_scores" = ranking_scores$ix)
  return(results)
  
}

VT.Predictive.one <- function(dataset_frame){
  # Uses S-learner
  vt.o <- vt.data(dataset_frame, "dataset.labels", "dataset.treatment", interactions = T)
  model.rf <- randomForest(x = vt.o$getX(interactions = T),
                           y = vt.o$getY(),
                           ntree = 1000)
  vt.f.rf <- vt.forest("one", vt.data = vt.o, model = model.rf, interactions = T)
  z <- vt.f.rf$difft
  model.rf.z <- randomForest(dataset_frame[,3:length(dataset_frame) ], y = z,ntree=1000)
  sorted_scores <- sort(model.rf.z$importance, decreasing=T,method='shell',index.return=TRUE) 
  ranking_scores <- sort(sorted_scores$ix, decreasing=F,method='shell',index.return=TRUE)
  results <- list("scores" = as.vector(t(model.rf.z$importance)), "ranking" = sorted_scores$ix, "ranking_scores" = ranking_scores$ix)
  return(results)
  
}


