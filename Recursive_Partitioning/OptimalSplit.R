library(RcppEigen)
library(partykit)
#Calculate interaction induced by split, by comparing the two resulting splits. 
CalculateSplitInteraction <- function(dosageTargetVector,dosageVector,phenotypes,covariates,split){
  #indicator of whether sample meets split criteria
  splitIndicator <- ifelse(dosageVector %in% split,1,0)
  #fit linear model to test significance of split-treatment interactions.
  mdlMat <- cbind('Intercept' = rep(1,length(dosageVector)),'treatment' = as.vector(dosageTargetVector),'split' = splitIndicator,'int' = as.vector(dosageTargetVector) * splitIndicator,as.matrix(covariates))
  #Fit linear model and calculate stats
  fit <- RcppEigen::fastLmPure(y = phenotypes,X = mdlMat)
  coeff <- fit$coefficients
  se <- fit$se
  t <- coeff/se
  G_s <- t['int'] ^ 2
  return(list(G_s = as.numeric(G_s),splitIndicator=splitIndicator))
}
#Find optimal splits from candidate splits
FindOptimalSplit <- function(dosageTarget,dosageMatrix,phenotypes,covariates){
  #Find best split out of all candidate splits
  splits <- list(0,2)
  bestSplitInfo <- list(rs=NA,split=NA,indicator=c(),G_s=-Inf)
  for(i in 1:ncol(dosageMatrix)){
    uniqueValues <- unique(dosageMatrix[,i])
    allowableSplits <- splits[sapply(splits,function(x) x %in% uniqueValues)]
    if(length(uniqueValues) == 1 | length(allowableSplits) == 0){
      next
    }
    for(j in 1:length(allowableSplits)){
      splitInfo <- CalculateSplitInteraction(dosageTarget,dosageMatrix[,i],phenotypes,covariates,allowableSplits[[j]])
      comparison <- splitInfo$G_s > bestSplitInfo$G_s
      if(comparison){
        bestSplitInfo$G_s <- splitInfo$G_s
        bestSplitInfo$rs <- colnames(dosageMatrix)[i]
        bestSplitInfo$split <- allowableSplits[[j]]
        bestSplitInfo$indicator <- splitInfo$splitIndicator
      }
    }
  }
  return(bestSplitInfo)
}
#Construct leaf node by fitting main effect model of target.
ConstructLeafNode <- function(dosageTarget,phenotypes,covariates,nodeID){
  mdlMat <- cbind('Intercept' = rep(1,length(dosageTarget)),'treatment' = as.vector(dosageTarget),as.matrix(covariates))
  fit <- RcppEigen::fastLmPure(y = phenotypes,X = mdlMat)
  node <- partykit::partynode(id = nodeID,info = paste0('beta=',signif(fit$coeff['treatment'],5),'|n=',length(dosageTarget)))
  return(node)
}

SubsetData <- function(data,subset){
  if(is.null(nrow(data))){
    data[subset]
  }else if(nrow(data) == 1){
    data[1,subset]
  }else{
    data[subset,]
  }
}
#Construct Tree
Split <- function(node,data,minSize){
  #Seperate data into subgroups for right node and left node
  if(node$split == 0){
    leftData <- data; rightData <- data
    leftIndicator <- as.logical(node$indicator)
    rightIndicator <- !leftIndicator
  }else if(node$split == 2){
    leftData <- data; rightData <- data
    rightIndicator <- as.logical(node$indicator)
    leftIndicator <- !rightIndicator
  }
  leftData <- lapply(data,function(x) SubsetData(x,leftIndicator))
  rightData <- lapply(data,function(x) SubsetData(x,rightIndicator))
  
  #Process Left node
  nodeIndex <<- nodeIndex + 1
  if(length(leftData$dosageTarget) < minSize){
    leftNode <- ConstructLeafNode(leftData$dosageTarget,leftData$phenotypes,leftData$covariates,nodeIndex)
  }else{
    leftBestSplit <- FindOptimalSplit(leftData$dosageTarget,leftData$dosageMatrix,leftData$phenotypes,leftData$covariates)
    leftSplitObj <- partysplit(which(colnames(leftData$dosageMatrix) == leftBestSplit$rs),
                               breaks = ifelse(leftBestSplit$split==0,0,2),
                               right = ifelse(leftBestSplit$split==0,T,F))
    leftNode <- partynode(nodeIndex,split=leftSplitObj,kids = Split(leftBestSplit,leftData,minSize),info=leftBestSplit$G_s)
  }

  
  #Process Right node
  nodeIndex <<- nodeIndex + 1
  if(length(rightData$dosageTarget) < minSize){
    rightNode <- ConstructLeafNode(rightData$dosageTarget,rightData$phenotypes,rightData$covariates,nodeIndex)
  }else{
    rightBestSplit <- FindOptimalSplit(rightData$dosageTarget,rightData$dosageMatrix,rightData$phenotypes,rightData$covariates)
    rightSplitObj <- partysplit(which(colnames(rightData$dosageMatrix) == rightBestSplit$rs),
                               breaks = ifelse(rightBestSplit$split==0,0,2),
                               right = ifelse(rightBestSplit$split==0,T,F))
    rightNode <- partynode(nodeIndex,split=rightSplitObj,kids = Split(rightBestSplit,rightData,minSize),info=rightBestSplit$G_s)
  }
  return(list(leftNode,rightNode))
}
ConstructTree <- function(data,minSize){
  rootBestSplit <- FindOptimalSplit(data$dosageTarget,data$dosageMatrix,data$phenotypes,data$covariates)
  splitObj <- partysplit(which(colnames(data$dosageMatrix) == rootBestSplit$rs),
                         breaks = ifelse(rootBestSplit$split==0,0,2),
                         right = ifelse(rootBestSplit$split==0,T,F))
  #Keep track of index used in nodes
  nodeIndex <<- 1L
  nodeObj <- partynode(nodeIndex,split=splitObj,kids = Split(rootBestSplit,data,minSize),info=splitObj$G_s)
  return(nodeObj)
}


pn <- ConstructTree(test,10000)
py <- party(pn,as.data.frame(test$dosageMatrix))
# split_test <- FindOptimalSplit(test$dosageTarget,test$dosageMatrix,test$phenotypesAndCov)
# pn <- partynode(1L,split=sp,kids = ConstructLeafNode(test$dosageTarget,test$phenotypesAndCov,split_test$indicator,c(2L,3L)))