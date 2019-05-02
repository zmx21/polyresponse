library(dplyr)
library(dbplyr)
library(data.table)
library(RSQLite)
interaction_path <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/HMGCR_LDL_known.txt'
#interaction_path <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/PCSK9_LDL_known.txt'

interaction_results <- data.table::fread(interaction_path)

WriteAnnovarInput <- function(){
  #Filter variants based on threshold
  interaction_results <- dplyr::filter(interaction_results,p_int < 1e-2) %>% dplyr::arrange(p_int)
  interaction_results <- dplyr::filter(interaction_results,rsid %in% colnames(readRDS('~/bsu_scratch/LDL_Project_Data/Random_Forest/rs12916_rs17238484_rs5909_rs2303152_rs10066707_rs2006760_LDLdirect/data_p_1e-04.rds')$dosageMatrix))
  
  #Connect to rsid annotation database. 
  anno_sql_name<- "all_snp_stats.sqlite"
  path <- '~/bsu_scratch/SQL/'
  curwd <- getwd()
  setwd(path)
  anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
  anno_db <- dplyr::tbl(anno_con,'all_snp_stats')
  #Get all rsids meeting criteria
  anno <- anno_db %>% dplyr::filter(rsid %in% interaction_results$rsid) %>% collect()
  dbDisconnect(anno_con)
  interaction_results <- interaction_results %>% dplyr::left_join(anno,by=c('rsid'='rsid')) 
  chr_pos <- interaction_results %>% dplyr::select(chromosome,position)
  chr_pos <- cbind(chr_pos,data.frame(position2=chr_pos$position))
  output <- cbind(chr_pos,dplyr::select(interaction_results,alleleA,alleleB,p_int))
  output$p_int <- as.character(output$p_int)
  
  data.table::fwrite(output,'~/bsu_scratch/LDL_Project_Data/Interaction_Data/annovar/HMGCR_int_indep_parsed',col.names = F,sep='\t')
  #data.table::fwrite(output,'~/bsu_scratch/LDL_Project_Data/Interaction_Data/annovar/PCSK9_int_parsed',col.names = F,sep='\t')
}

GetAllInteractionMetadata <- function(){
  #Filter variants based on threshold
  interaction_results <- dplyr::filter(interaction_results,p_int < 1e-5) %>% dplyr::arrange(p_int)
  
  #Connect to rsid annotation database. 
  anno_sql_name<- "all_snp_stats.sqlite"
  path <- '~/bsu_scratch/SQL/'
  curwd <- getwd()
  setwd(path)
  anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
  anno_db <- dplyr::tbl(anno_con,'all_snp_stats')
  
  #Get all rsids meeting criteria
  anno <- anno_db %>% dplyr::filter(rsid %in% interaction_results$rsid) %>% collect()
  dbDisconnect(anno_con)
  interaction_results <- interaction_results %>% dplyr::left_join(anno,by=c('rsid'='rsid')) 
  interaction_results <- interaction_results %>% dplyr::select(rsid,'p-value'=p_int,'Chr'=chromosome,'Position'=position,'Minor Allele'=minor_allele,'Major Allele'=major_allele,'MAF' = minor_allele_frequency)
  interaction_results$`p-value` <- signif(interaction_results$`p-value`,4)
  interaction_results$MAF <- signif(as.numeric(interaction_results$MAF),4)
  

  #Read in annovar results
  resultPath <- '~/bsu_scratch/LDL_Project_Data/Interaction_Data/annovar/'
  file_name <- 'HMGCR_int_parsed.variant_function'
  annovarResult <- data.table::fread(paste0(resultPath,file_name))
  annovarResultFilt <- annovarResult %>% dplyr::filter(V8 < 1e-5) %>% dplyr::arrange(V8) 
  annovarResultFilt <- annovarResultFilt %>% dplyr::select(type=V1,genes=V2,chr=V3,pos=V4,p=V8)
  annovarResultFilt <- cbind(annovarResultFilt,interaction_results)
  data.table::fwrite(annovarResultFilt,file='~/bsu_scratch/LDL_Project_Data/HMGCR_Interaction_Stats.csv',sep = ',')
  
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
  
  geneRegionAnnovar <- annovarResultFilt[0,]
  for(i in 1:length(allGenes)){
    curGene <- allGenes[i]
    matchingEntries <- sapply(annovarResultFilt$genes,function(x) grepl(curGene,x))
    lowestP <- which.min(annovarResultFilt$p[matchingEntries])
    geneRegionAnnovar <- rbind(geneRegionAnnovar,annovarResultFilt[which(matchingEntries)[lowestP],])
  }
  geneRegionAnnovar$indepGene <- allGenes
  geneRegionAnnovar <- geneRegionAnnovar[,!(colnames(geneRegionAnnovar) %in% c('genes','type','chr','pos','p'))]
  
  geneRegionAnnovarUnique <- geneRegionAnnovar[0,]
  uniqueRS <- unique(geneRegionAnnovar$rsid)
  genes <- rep('',length(uniqueRS))
  
  for(j in 1:length(uniqueRS)){
    curRS <- uniqueRS[j]
    matchingEntries <- dplyr::filter(geneRegionAnnovar,rsid == curRS)
    genes[j] <- paste(matchingEntries$indepGene,collapse = ';')
    geneRegionAnnovarUnique <- rbind(geneRegionAnnovarUnique,matchingEntries[1,])
  }
  geneRegionAnnovarUnique$indepGene <- genes
  data.table::fwrite(geneRegionAnnovarUnique[order(geneRegionAnnovarUnique$`p-value`),],file='~/bsu_scratch/LDL_Project_Data/HMGCR_Interaction_Stats_Indep.csv',sep = ',')
  
}














