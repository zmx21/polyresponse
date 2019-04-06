library(dplyr)
resultPath <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/annovar/'
p_thresh <- 1e-4

#Read in annovar results
annovarResult <- data.table::fread(paste0(resultPath,'HMGCR_int_parsed.variant_function'))
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
data.table::fwrite(gene_df,'~/bsu_scratch/LDL_Project_Data/Interaction_Data/annovar/parsed_genes.txt',sep = ',',col.names = F)

library(enrichR)
dbs <- c('GO_Molecular_Function_2015','Reactome_2016')
enrichr_results <- enrichr(gene_df,dbs)

mol_func <- enrichr_results$GO_Molecular_Function_2015[1:5,]
mol_func_genes <- lapply(mol_func$Genes,function(x) unlist(strsplit(x,';')))
mol_func$Term <- sapply(mol_func$Term,function(x) gsub(x = x,pattern = ' (',replacement = '\n(',fixed = T))
mol_func$Term <- sapply(1:length(mol_func$Term),function(i) paste0(mol_func$Term[i],'\n','p=',signif(mol_func$P.value[i],3)))

#Set up order of genes
includedGenes <- c()
for(i in 1:nrow(mol_func)){
  curGenes <- mol_func_genes[[i]]
  includedGenes <- c(includedGenes,setdiff(curGenes,includedGenes))
}
includedGenes <- sapply(includedGenes,function(x) paste0(x,' (p=',as.character(signif(10^(-1*gene_df$weights[which(gene_df$gene==x)]*max(-log10(annovarResultFilt$p))),3)),')'))
mol_func_df <- data.frame(genes=factor(includedGenes,levels = rev(includedGenes)))
for(i in 1:nrow(mol_func)){
  curGo <- mol_func$Term[i]
  curGene <- mol_func_genes[[i]]
  geneInclVector <- ifelse(includedGenes %in% sapply(curGene,function(x) paste0(x,' (p=',as.character(signif(10^(-1*gene_df$weights[which(gene_df$gene==x)]*max(-log10(annovarResultFilt$p))),3)),')')),'T','F')
  curDf <- as.data.frame(geneInclVector)
  colnames(curDf) <- curGo
  mol_func_df <- cbind(mol_func_df,curDf)
}

library(ggplot2); library(reshape2)
melted_mol_func <- melt(mol_func_df, id.var = 'genes')
ggplot(melted_mol_func, aes(variable, genes)) + geom_tile(aes(fill = value),
                                                colour = "white",show.legend = F) + scale_fill_manual(values=c("grey", "red")) + 
  theme(axis.text.x = element_text(size = 12,angle = 55,hjust = 1),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),axis.text.y = element_text(size=12)) + xlab('GO Molecular Process') + ylab('Gene')
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


