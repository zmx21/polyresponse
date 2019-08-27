RunEnrichR <- function(p_thresh,file_name){
  library(dplyr)
  resultPath <- '~/bsu_scratch/LDL_Project_Data_Aug2019/Interaction_Data/annovar/'

  #Read in annovar results
  annovarResult <- data.table::fread(paste0(resultPath,file_name))
  annovarResultFilt <- dplyr::filter(annovarResult,V8 < p_thresh) %>% dplyr::select(type=V1,genes=V2,chr=V3,pos=V4,p=V8)
  annovarResultFilt$score <- (-log10(annovarResultFilt$p) / max(-log10(annovarResultFilt$p)))
  
  #Seperate according to comma
  allGenes <- annovarResultFilt$genes
  allGenes <- unlist(lapply(allGenes,function(x) strsplit(x,',',fixed = T)))
  #Remove none entries and NM entries
  removeEntries <- sapply(allGenes,function(x) grepl('NONE',x)) | sapply(allGenes,function(x) grepl('NM_',x))
  allGenes <- allGenes[!removeEntries]
  #Remove dist brackets
  allGenes <- sapply(allGenes,function(x) unlist(strsplit(x,'(',fixed = T))[1],USE.NAMES = F)
  #Remove repeating elements
  allGenes <- unique(allGenes)
  
  weights <- rep(NA,length(allGenes))
  for(i in 1:length(allGenes)){
    curGene <- allGenes[i]
    matchingEntries <- sapply(annovarResultFilt$genes,function(x) grepl(curGene,x))
    lowestP <- which.min(annovarResultFilt$p[matchingEntries])
    weights[i] <- annovarResultFilt$score[i]
  }
  
  
  gene_df <- data.frame(gene=allGenes,weights=weights)
  write(as.character(gene_df$gene),'~/bsu_scratch/LDL_Project_Data_Aug2019/Interaction_Data/annovar/parsed_genes.txt',sep = '\n')
  
  library(enrichR)
  dbs <- c('GO_Molecular_Function_2015','Reactome_2016')
  enrichr_results <- enrichr(gene_df,dbs)
  
  mol_func <- enrichr_results$GO_Molecular_Function_2015[1:5,]
  reactome <- enrichr_results$Reactome_2016[1:5,]
  return(list(mol_func=mol_func,reactome=reactome,gene_df=gene_df,annovarResultFilt=annovarResultFilt))
}


DrawDendroGraphGOMolFunc <- function(data,gene_df,annovarResultFilt){
  data_genes <- lapply(data$Genes,function(x) unlist(strsplit(x,';')))
  data$Term <- sapply(data$Term,function(x) gsub(x = x,pattern = ' (',replacement = '\n(',fixed = T))
  data$Term <- sapply(1:length(data$Term),function(i) paste0(data$Term[i],'\n','p=',signif(data$P.value[i],3)))
  
  #Set up order of genes
  includedGenes <- c()
  for(i in 1:nrow(data)){
    curGenes <- data_genes[[i]]
    includedGenes <- c(includedGenes,setdiff(curGenes,includedGenes))
  }
  includedGenes <- sapply(includedGenes,function(x) paste0(x,' (p=',as.character(signif(10^(-1*gene_df$weights[which(gene_df$gene==x)]*max(-log10(annovarResultFilt$p))),3)),')'))
  df <- data.frame(genes=factor(includedGenes,levels = rev(includedGenes)))
  for(i in 1:nrow(data)){
    curGo <- data$Term[i]
    curGene <- data_genes[[i]]
    geneInclVector <- ifelse(includedGenes %in% sapply(curGene,function(x) paste0(x,' (p=',as.character(signif(10^(-1*gene_df$weights[which(gene_df$gene==x)]*max(-log10(annovarResultFilt$p))),3)),')')),'T','F')
    curDf <- as.data.frame(geneInclVector)
    colnames(curDf) <- curGo
    df <- cbind(df,curDf)
  }
  
  library(ggplot2); library(reshape2)
  melted_df <- melt(df, id.var = 'genes')
  p <- ggplot(melted_df, aes(variable, genes)) + geom_tile(aes(fill = value),
                                                            colour = "black",show.legend = F,size=1.1) + scale_fill_manual(values=c("grey", "red")) + 
    theme(axis.text.x = element_text(size = 12,angle = 60,hjust = 1,vjust=1,family='Myriad Pro'),axis.title.x = element_text(size=14,family='Myriad Pro'),axis.title.y = element_text(size=14,family='Myriad Pro'),axis.text.y = element_text(size=12,family='Myriad Pro')) + xlab('') + ylab('Gene\n')
  return(p)
}
DrawDendroGraphReactome<- function(data,gene_df,annovarResultFilt){
  data_genes <- lapply(data$Genes,function(x) unlist(strsplit(x,';')))
  data$Term <- sapply(data$Term,function(x) gsub(x = x,pattern = ' (',replacement = '\n(',fixed = T))
  data$Term <- sapply(1:length(data$Term),function(i) paste0(data$Term[i],'\n','p=',signif(data$P.value[i],3)))
  
  #Set up order of genes
  includedGenes <- c()
  for(i in 1:nrow(data)){
    curGenes <- data_genes[[i]]
    includedGenes <- c(includedGenes,setdiff(curGenes,includedGenes))
  }
  includedGenes <- sapply(includedGenes,function(x) paste0(x,' (p=',as.character(signif(10^(-1*gene_df$weights[which(gene_df$gene==x)]*max(-log10(annovarResultFilt$p))),3)),')'))
  df <- data.frame(genes=factor(includedGenes,levels = rev(includedGenes)))
  for(i in 1:nrow(data)){
    curGo <- data$Term[i]
    curGo <- unlist(strsplit(curGo,split = '_',fixed = T))[1]
    curGo <- gsub(pattern = 'of',replacement = 'of\n',x = curGo,fixed = T)
    curGo <- gsub(pattern = 'and',replacement = 'and\n',x = curGo,fixed=T)
    curGo <- gsub(pattern = ' in ',replacement = ' in\n',x = curGo,fixed = T)
    curGo <- gsub(pattern = ' by ',replacement = ' by\n',x = curGo,fixed = T)
    curGo <- gsub(pattern = ' causes ',replacement = ' causes\n',x = curGo,fixed = T)
    
    curGene <- data_genes[[i]]
    geneInclVector <- ifelse(includedGenes %in% sapply(curGene,function(x) paste0(x,' (p=',as.character(signif(10^(-1*gene_df$weights[which(gene_df$gene==x)]*max(-log10(annovarResultFilt$p))),3)),')')),'T','F')
    curDf <- as.data.frame(geneInclVector)
    colnames(curDf) <- curGo
    df <- cbind(df,curDf)
  }
  
  library(ggplot2); library(reshape2)
  melted_df <- melt(df, id.var = 'genes')
  p <- ggplot(melted_df, aes(variable, genes)) + geom_tile(aes(fill = value),
                                                           colour = "white",show.legend = F) + scale_fill_manual(values=c("grey", "red")) + 
    theme(axis.text.x = element_text(size = 12,angle = 60,hjust = 1,vjust=1,family='Myriad Pro'),axis.title.x = element_text(size=14,family='Myriad Pro'),axis.title.y = element_text(size=14,family='Myriad Pro'),axis.text.y = element_text(size=12,family='Myriad Pro')) + xlab('Reactome Pathway') + ylab('Gene')
  return(p)
}


HMGCR_res <- RunEnrichR(1e-5,'HMGCR_int_indep_parsed.variant_function')
HMGCR_mol_func <- DrawDendroGraphGOMolFunc(HMGCR_res$mol_func,HMGCR_res$gene_df,HMGCR_res$annovarResultFilt)
HMGCR_reactome <- DrawDendroGraphReactome(HMGCR_res$reactome,HMGCR_res$gene_df,HMGCR_res$annovarResultFilt)

# PCSK9_res <- RunEnrichR(5e-6,'PCSK9_int_parsed.variant_function')
# PCSK9_mol_func <- DrawDendroGraphGOMolFunc(PCSK9_res$mol_func,PCSK9_res$gene_df,PCSK9_res$annovarResultFilt)
# PCSK9_reactome <- DrawDendroGraphReactome(PCSK9_res$reactome,PCSK9_res$gene_df,PCSK9_res$annovarResultFilt)
# save.image(file='~/bsu_scratch/LDL_Project_Data/int_variant_anno_5e-6.RData')

# ff <- factor(as.matrix(mol_func_df[,2:ncol(mol_func_df)]),
#              levels = c('T','F'),
#              labels = c('T','F'))
# fx <- matrix(as.numeric(ff),ncol = ncol(mol_func_df) - 1)
# 
# col <- c('T'='red',F='white')
# imgflip<-function(x) {t(x[nrow(x):1,])}
# image(imgflip(fx),
#       breaks=(1:(nlevels(ff)+1))-.5,
#       col=col[levels(ff)],
#       xaxt="n", yaxt="n"
# )
# axis(2, at=seq(0,1,length.out=nrow(mol_func_df)), labels=rev(mol_func_df$genes), las=2)
# axis(3, at=seq(0,1,length.out=length(names(mol_func_df))-1), labels=names(mol_func_df)[-1],las=2)


