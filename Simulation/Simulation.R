set.seed(314)
library(MASS)
library(pbmcapply)
library(ggplot2)
source('../Recursive_Partitioning/InteractionTree.R')
GeneratePhenotypes <- function(gamma_int,parts = 10000,MAF = 0.3,vars=100,gamma_0 = 0.3,abs_only = F,perc_effect=0.01,g,z,epsx){
  #Effect size of moderating variants.
  if(gamma_int == 0){
    gamma = rep(0,vars)
  }else{
    if(abs_only){
      gamma = c(rep(0,(1-perc_effect)*vars),rnorm(perc_effect*vars, gamma_int, 0.01))
    }else{
      gamma = c(rep(0,(1-perc_effect)*vars),rnorm(floor(perc_effect/2*vars), gamma_int, 0.01),rnorm(ceiling(perc_effect/2*vars), -gamma_int, 0.01))
    }
  }

  #Moderating effect of each individual
  effect = g%*%gamma
  
  x = (gamma_0 + effect)*z + epsx
  
  #Set Names
  snp_names <- paste0('rs',1:vars)
  part_IDs <- paste0('Indv',1:parts)
  x <- as.vector(x)
  names(x) <- part_IDs
  z <- as.vector(z)
  names(z) <- part_IDs
  colnames(g) <- snp_names
  rownames(g) <- part_IDs
  covariates=matrix(ncol = 0,nrow = length(x))
  
  names(gamma) <- snp_names
  effect <- as.vector(effect)
  names(effect) <- part_IDs
  
  
  return(list(phenotypes=x,dosageTarget=z,dosageMatrix=g,covariates=covariates,metadata = list(gamma=gamma,effect=effect,gamma_int=gamma_int,gamma_0=gamma_0)))
}
RunLmPrediction <- function(sim_obj,train_set,test_set,metadata){
  int_train <- lapply(1:ncol(sim_obj$dosageMatrix[train_set,]),function(i) summary(lm('y~x*g',data.frame(y=sim_obj$phenotypes[train_set],
                                                                                                         x=sim_obj$dosageTarget[train_set],
                                                                                                         g=sim_obj$dosageMatrix[train_set,i]))))
  p_values <- sapply(int_train,function(x) x$coefficients[4,"Pr(>|t|)"])
  top_int <- which.min(p_values)
  
  mdl_int <- lm('y~ 0 + x + x:g',data.frame(y=sim_obj$phenotypes[train_set],
                                   x=sim_obj$dosageTarget[train_set],
                                   g=sim_obj$dosageMatrix[train_set,top_int]))
  pred_delta_int <- predict.lm(mdl_int,newdata = data.frame(y=sim_obj$phenotypes[test_set],
                                                           x=sim_obj$dosageTarget[test_set],
                                                           g=sim_obj$dosageMatrix[test_set,top_int])) / sim_obj$dosageTarget[test_set]
  
  mdl_no_int <- lm('y~ 0 + x',data.frame(y=sim_obj$phenotypes[train_set],
                                        x=sim_obj$dosageTarget[train_set],
                                        g=sim_obj$dosageMatrix[train_set,top_int]))
  pred_delta_no_int <- predict.lm(mdl_no_int,newdata = data.frame(y=sim_obj$phenotypes[test_set],
                                                                 x=sim_obj$dosageTarget[test_set],
                                                                 g=sim_obj$dosageMatrix[test_set,top_int])) / sim_obj$dosageTarget[test_set]
  
  true_delta <- as.vector(metadata$metadata$gamma_0 + sim_obj$dosageMatrix[test_set,]%*%metadata$metadata$gamma)
  test_set_data <- sapply(sim_obj,function(x) SubsetData(x,test_set))
  rmse_int <- sqrt(mean((true_delta - pred_delta_int)^2)) 
  rmse_no_int <- sqrt(mean((true_delta - pred_delta_no_int)^2)) 
  return(list(rmse_int=rmse_int,rmse_no_int=rmse_no_int))
}
ConstructSimRF <- function(train_set_dat,test_set_dat,n_bootstrap,minSize = 500,n_features = 75){
  pred_delta <- pbmclapply(1:n_bootstrap,function(i){
    #Train tree using bootstrap sample
    boostrap_samples <- sample(1:length(train_set_dat$phenotypes),size = length(train_set_dat$phenotypes) * 2/3)
    boostrap_set_dat <- lapply(train_set_dat,function(x) SubsetData(x,boostrap_samples))
    #Construct boostrap tree
    boostrap_tree <- ConstructTree(data = boostrap_set_dat,minSize = minSize,n_features = n_features)
    boostrap_tree_party <- party(boostrap_tree,as.data.frame(boostrap_set_dat$dosageMatrix))
    #Predict using testing set
    beta_pred_string <- nodeapply(boostrap_tree_party,ids = as.numeric(predict(boostrap_tree_party,as.data.frame(test_set_dat$dosageMatrix))),FUN = function(n) n$info)
    beta_pred <- sapply(beta_pred_string,function(x) as.numeric(unlist(strsplit(x = unlist(x),split = "  "))[1]))
    return(beta_pred)
  },mc.cores = 25)
  return(colMeans(do.call(rbind,pred_delta)))
}

CompareModels <- function(gamma_int,perc_effect,min_node_size,parts,train_set,test_set,abs_only=T,g,z,epsx){
  beta_0p2 = GeneratePhenotypes(gamma_int=gamma_int,parts = length(train_set) + length(test_set),MAF = 0.3,vars=100,gamma_0 = 0.3,abs_only=abs_only,perc_effect=perc_effect,g,z,epsx)
  beta_0p2_dat <- beta_0p2[!names(beta_0p2) == "metadata"]
  metadata <- beta_0p2[names(beta_0p2) == "metadata"]
  pred <- RunLmPrediction(beta_0p2_dat,train_set,test_set,metadata)
  
  train_set_dat <- lapply(beta_0p2_dat,function(x) SubsetData(x,train_set))
  test_set_dat <- lapply(beta_0p2_dat,function(x) SubsetData(x,test_set))
  
  test_set_true_delta <- as.vector(metadata$metadata$gamma_0 + test_set_dat$dosageMatrix%*%metadata$metadata$gamma)
  test_set_pred <- ConstructSimRF(train_set_dat,test_set_dat,25,minSize=min_node_size)
  rmse_RF <- sqrt(mean((test_set_true_delta - test_set_pred)^2))
  return(c(pred,list(rmse_RF=rmse_RF)))
}
parts <- 500000
vars <- 100
MAF = 0.3
train_set <- sample(1:parts,size = 2/3 * parts)
test_set <- setdiff(1:parts,train_set)

g = array(rbinom(parts*vars, 2, MAF), dim=c(parts, vars))    # Moderating variants
z = rnorm(parts)                                             # Pharmacomimetic score
epsx = rnorm(parts) #Noise term for risk factor

test_cond <- expand.grid(gamma_int=c(0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2),perc_effect=c(0.06))
mdl_results <- lapply(1:nrow(test_cond),function(i) CompareModels(test_cond$gamma_int[i],test_cond$perc_effect[i],5000,parts,train_set,test_set,abs_only=T,g,z,epsx))
rmse_RF <- sapply(mdl_results,function(x) x$rmse_RF)
rmse_int <- sapply(mdl_results,function(x) x$rmse_int)
rmse_no_int <- sapply(mdl_results,function(x) x$rmse_no_int)
test_df <- cbind(rbind(test_cond,test_cond,test_cond))
test_df$RMSE <- c(sapply(mdl_results,function(x) x$rmse_no_int),sapply(mdl_results,function(x) x$rmse_int),sapply(mdl_results,function(x) x$rmse_RF))
test_df$Model <- c(rep("Main Effect",nrow(test_cond)),rep("Main Effect and Top Moderator",nrow(test_cond)),rep("RFIT",nrow(test_cond)))

# p_5000 <- ggplot2::ggplot(test_df,aes(x=gamma_int,y=RMSE,color = Model)) + geom_point() + geom_line() +
#   xlab(expression(Moderating~Effect~(gamma[int]))) + ylab('RMSE') + scale_color_discrete(name = 'Method') + ggtitle('a)')

save.image('/mnt/data3/mxu/InteractionTree/SimResults_5000.rda')

test_cond <- expand.grid(gamma_int=c(0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2),perc_effect=c(0.06))
mdl_results <- lapply(1:nrow(test_cond),function(i) CompareModels(test_cond$gamma_int[i],test_cond$perc_effect[i],5000,parts,train_set,test_set,abs_only=F,g,z,epsx))
rmse_RF <- sapply(mdl_results,function(x) x$rmse_RF)
rmse_int <- sapply(mdl_results,function(x) x$rmse_int)
rmse_no_int <- sapply(mdl_results,function(x) x$rmse_no_int)
test_df <- cbind(rbind(test_cond,test_cond,test_cond))
test_df$RMSE <- c(sapply(mdl_results,function(x) x$rmse_no_int),sapply(mdl_results,function(x) x$rmse_int),sapply(mdl_results,function(x) x$rmse_RF))
test_df$Model <- c(rep("Main Effect",nrow(test_cond)),rep("Main Effect and Top Moderator",nrow(test_cond)),rep("RFIT",nrow(test_cond)))
# 
# p_5000Neg <- ggplot2::ggplot(test_df,aes(x=gamma_int,y=RMSE,color = Model)) + geom_point() + geom_line() +
#   xlab(expression(Moderating~Effect~(gamma[int]))) + ylab('RMSE') + scale_color_discrete(name = 'Method') + ggtitle('b)')

save.image('/mnt/data3/mxu/InteractionTree/SimResults_5000Neg.rda')

load('/mnt/data3/mxu/InteractionTree/SimResults_5000.rda')
test_df$Model <- c(rep("Single Value",nrow(test_cond)),rep("Single Moderating Variant",nrow(test_cond)),rep("RFIT",nrow(test_cond)))
p_5000 <- ggplot2::ggplot(test_df,aes(x=gamma_int,y=RMSE,color = Model)) + geom_point() + geom_line() +
  xlab(expression(Moderating~Effect~(gamma[int]))) + ylab('RMSE') + scale_color_discrete(name = 'Method') + ggtitle('a)')
load('/mnt/data3/mxu/InteractionTree/SimResults_5000Neg.rda')
test_df$Model <- c(rep("Single Value",nrow(test_cond)),rep("Single Moderating Variant",nrow(test_cond)),rep("RFIT",nrow(test_cond)))
p_5000Neg <- ggplot2::ggplot(test_df,aes(x=gamma_int,y=RMSE,color = Model)) + geom_point() + geom_line() +
  xlab(expression(Magnitude~of~Moderating~Effect~("|"~gamma[int]~"|"))) + ylab('RMSE') + scale_color_discrete(name = 'Method') + ggtitle('b)')
ggarrange(p_5000,p_5000Neg,common.legend = T)
