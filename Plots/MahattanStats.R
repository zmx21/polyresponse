source('~/MRC_BSU_Internship_LDL/Gene_Phenotype_Association/GetSNPsOfGene.R')
source('~/MRC_BSU_Internship_LDL/Gene_Phenotype_Association/GenePhenotypeAssociation.R')
includedSNPS <- c('rs11206510','rs2479409','rs2149041','rs2479394','rs10888897','rs7552841','rs562556')
geneName <- 'PCSK9'
# includedSNPS <- c('rs12916,rs17238484,rs5909,rs2303152,rs10066707,rs2006760')
# geneName <- 'HMGCR'
phenotype <- 'LDLdirect'

gene_snps <- CollectBetaCoeff(geneName,100000,100000,phenotype,0.01,0.5,12) %>% 
  dplyr::select('MarkerName'=rsid,'P-value'=p) %>% dplyr::filter(!is.na('P.value'))

gene_snps$Included <- gene_snps$MarkerName %in% includedSNPS
data.table::fwrite(gene_snps %>% dplyr::filter(!is.na('P.value')),
                   paste0('~/bsu_scratch/LDL_Project_Data/Locus_Zoom/',geneName,'_SNPs.txt'),row.names = F,sep = '\t')

path <-  '~/bsu_scratch/LDL_Project_Data/Genotype_Data/'
bgen_file_prefix <- 'ukb_imp_chr#_HRConly'

refSNP <- includedSNPS[1] #gene_snps$MarkerName[which.min(gene_snps$`P.value`)]
refdosageMatrix <- as.vector(LoadBgen(path,bgen_file_prefix,refSNP))

otherSNPS <- setdiff(gene_snps$MarkerName,refSNP)
r2 <- pbmcapply::pbmclapply(otherSNPS,function(x){
  tryCatch({
    curDosage <- as.vector(LoadBgen(path,bgen_file_prefix,x));
    cor(refdosageMatrix,curDosage) ^ 2
  },error=function(e){
    NA
  }
  )
},mc.cores = 12)

r2_df <- data.frame(snp1=otherSNPS,snp2=rep(refSNP,length(otherSNPS)),dprime=rep(0,length(otherSNPS)),rsquare=unlist(r2))
data.table::fwrite(r2_df,file=paste0('~/bsu_scratch/LDL_Project_Data/Locus_Zoom/',geneName,'_r2.txt'),row.names = F,sep = ' ',quote=F)
