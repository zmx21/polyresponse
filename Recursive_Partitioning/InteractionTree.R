library(RcppEigen)
library(partykit)
#Calculate interaction according to specified indeicator vector
CalculateInteraction <- function(dosageTargetVector,splitIndicator,phenotypes,covariates){
  #fit linear model to test significance of split-treatment interactions.
  mdlMat <- cbind('Intercept' = rep(1,length(dosageTargetVector)),'treatment' = as.vector(dosageTargetVector),'split' = splitIndicator,'int' = as.vector(dosageTargetVector) * splitIndicator,as.matrix(covariates))
  #Fit linear model and calculate stats
  fit <- RcppEigen::fastLmPure(y = phenotypes,X = mdlMat)
  coeff <- fit$coefficients
  se <- fit$se
  t <- coeff/se
  G_s <- t['int'] ^ 2
  return(as.numeric(G_s))
}
#Calculate interaction induced by split, by comparing the two resulting splits. 
CalculateSplitInteraction <- function(dosageTargetVector,dosageVector,phenotypes,covariates,split){
  #indicator of whether sample meets split criteria
  splitIndicator <- ifelse(dosageVector %in% split,1,0)
  G_s <- CalculateInteraction(dosageTargetVector,splitIndicator,phenotypes,covariates)
  return(list(G_s = G_s,splitIndicator=splitIndicator))
}
#Find optimal splits from candidate splits
FindOptimalSplit <- function(dosageTarget,dosageMatrix,phenotypes,covariates,n_features,min_size){
  #Find best split out of all candidate splits
  splits <- list(0,2)
  bestSplitInfo <- list(rs=NA,split=NA,indicator=c(),G_s=-Inf)
  #Sample random set of features (random forest framework)
  if(n_features < ncol(dosageMatrix)){
    randomFeatures <- sample(1:ncol(dosageMatrix),size=n_features,replace = F)
  }else if(n_features == ncol(dosageMatrix)){
    randomFeatures <- 1:ncol(dosageMatrix)
  }else{
    stop('number of features exceed availiable features')
  }
  dosageMatrix <- dosageMatrix[,randomFeatures]
  for(i in 1:ncol(dosageMatrix)){
    uniqueValues <- unique(dosageMatrix[,i])
    allowableSplits <- splits[sapply(splits,function(x) x %in% uniqueValues)]
    if(length(uniqueValues) == 1 | length(allowableSplits) == 0){
      next
    }
    for(j in 1:length(allowableSplits)){
      splitInfo <- CalculateSplitInteraction(dosageTarget,dosageMatrix[,i],phenotypes,covariates,allowableSplits[[j]])
      if(is.na(splitInfo$G_s) | any(table(splitInfo$splitIndicator) < min_size)){
        next
      }
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
  node <- partykit::partynode(id = nodeID,info = paste0(signif(fit$coeff['treatment'],4),"  ",length(dosageTarget)))
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
Split <- function(node,data,minSize,n_features){
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
    leftBestSplit <- FindOptimalSplit(leftData$dosageTarget,leftData$dosageMatrix,leftData$phenotypes,leftData$covariates,n_features,minSize)
    #If no allowable splits remaining, create terminal node
    if(is.na(leftBestSplit$rs)){
      leftNode <- ConstructLeafNode(leftData$dosageTarget,leftData$phenotypes,leftData$covariates,nodeIndex)
    }else{
      leftSplitRSIndex <- which(colnames(leftData$dosageMatrix) == leftBestSplit$rs)
      leftSplitObj <- partysplit(leftSplitRSIndex,
                                 breaks = ifelse(leftBestSplit$split==0,0,2),
                                 right = ifelse(leftBestSplit$split==0,T,F))
      leftNode <- partynode(nodeIndex,split=leftSplitObj,kids = Split(leftBestSplit,leftData,minSize,n_features),info=leftBestSplit$rs)
    }
  }

  
  #Process Right node
  nodeIndex <<- nodeIndex + 1
  if(length(rightData$dosageTarget) < minSize){
    rightNode <- ConstructLeafNode(rightData$dosageTarget,rightData$phenotypes,rightData$covariates,nodeIndex)
  }else{
    rightBestSplit <- FindOptimalSplit(rightData$dosageTarget,rightData$dosageMatrix,rightData$phenotypes,rightData$covariates,n_features,minSize)
    #If no allowable splits remaining, create terminal node
    if(is.na(rightBestSplit$rs)){
      rightNode <- ConstructLeafNode(rightData$dosageTarget,rightData$phenotypes,rightData$covariates,nodeIndex)
    }else{
      rightSplitRSIndex <- which(colnames(rightData$dosageMatrix) == rightBestSplit$rs)
      rightSplitObj <- partysplit(rightSplitRSIndex,
                                  breaks = ifelse(rightBestSplit$split==0,0,2),
                                  right = ifelse(rightBestSplit$split==0,T,F))
      rightNode <- partynode(nodeIndex,split=rightSplitObj,kids = Split(rightBestSplit,rightData,minSize,n_features),info=rightBestSplit$rs)
    }
  }
  return(list(leftNode,rightNode))
}
ConstructTree <- function(data,minSize,n_features){
  rootBestSplit <- FindOptimalSplit(data$dosageTarget,data$dosageMatrix,data$phenotypes,data$covariates,n_features,minSize)
  splitObj <- partysplit(which(colnames(data$dosageMatrix) == rootBestSplit$rs),
                         breaks = ifelse(rootBestSplit$split==0,0,2),
                         right = ifelse(rootBestSplit$split==0,T,F))
  #Keep track of index used in nodes
  nodeIndex <<- 1L
  nodeObj <- partynode(nodeIndex,split=splitObj,kids = Split(rootBestSplit,data,minSize,n_features),info=rootBestSplit$rs)
  return(nodeObj)
}
#Predict beta coefficient of target SNP/Gene, given a genotype vector of a sample, based on the constructed tree
PredictBeta <- function(genotype,tree){
  leafNodeId <- predict(tree,genotype)
  info <- nodeapply(tree,ids = as.numeric(leafNodeId),FUN = function(n) n$info)
  beta <- as.numeric(unlist(strsplit(x = unlist(info),split = "  "))[1])
  return(beta)
}
#Calculate interaction of each split, by sending a sample genotype data down the tree.
CalculateTotalInteraction <- function(genotypeData,tree,dosageTargetVector,phenotypes,covariates){
  terminalNodes <- nodeids(tree,terminal = T)
  internalNodes <- setdiff(nodeids(tree),terminalNodes)
  leafNodeAssignment <- predict(tree,genotypeData)
  G_s <- rep(NA,length(internalNodes))
  for(i in 1:length(internalNodes)){
    #Get information of daughter nodes
    daughterInfo <- nodeapply(tree,ids = internalNodes[i],FUN = function(n) n$kids)
    #Obtain subgroups belonging to left subtree and right subtree
    leftDaughter <- as.numeric(unlist(daughterInfo[[1]][1])['id'])
    leftSubTreeLeaves <- intersect(as.numeric(names(tree[[leftDaughter]])),terminalNodes)
    leftSubGroup <- names(leafNodeAssignment)[leafNodeAssignment %in% leftSubTreeLeaves]
    
    rightDaughter <- as.numeric(unlist(daughterInfo[[1]][2])['id'])
    rightSubTreeLeaves <- intersect(as.numeric(names(tree[[rightDaughter]])),terminalNodes)
    rightSubGroup <- names(leafNodeAssignment)[leafNodeAssignment %in% rightSubTreeLeaves]
    
    #Keep data included in both subtree
    overallSubGroup <- c(leftSubGroup,rightSubGroup)
    dosageSubset <- dosageTargetVector[overallSubGroup]
    #Assign indicator to seperate subtree
    splitIndicator <- rep(0,length(dosageSubset))
    splitIndicator[names(dosageSubset) %in% leftSubGroup] <- 1
    G_s[i] <- CalculateInteraction(dosageSubset,splitIndicator,phenotypes[overallSubGroup],covariates[overallSubGroup,])
  }
  names(G_s) <- as.character(unlist(nodeapply(tree,ids=internalNodes,FUN = function(n) n$info)))
  return(G_s)
}

# pn <- ConstructTree(test,50000,5)
# py <- party(pn,as.data.frame(test$dosageMatrix))
# plot(py)
# int <- CalculateTotalInteraction(as.data.frame(test$dosageMatrix),py,test$dosageTarget,test$phenotypes,test$covariates)
