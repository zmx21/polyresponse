source('~/MRC_BSU_Internship/Gene_Phenotype_Association/GetSNPsOfGene.R')
source('~/MRC_BSU_Internship/Gene_Phenotype_Association/GenePhenotypeAssociation.R')
CACNA1D_SNPs <- CollectBetaCoeff('CACNA1D',50000,50000,'sbp',0.01,0.5,16) %>% dplyr::select('MarkerName'=rsid,'P.value'=p) %>% dplyr::filter(!is.na('P.value'))

includedSNPS <- c('rs3821843','rs7340705','rs113210396','rs312487','rs11719824','rs3774530','rs3821856')
CACNA1D_SNPs$Included <- CACNA1D_SNPs$MarkerName %in% includedSNPS
data.table::fwrite(CACNA1D_SNPs %>% dplyr::filter(!is.na('P.value')),'~/bsu_scratch/Locus_Zoom/CACNA1D_SNPs.txt',row.names = F,sep = '\t')

path <-  '~/bsu_scratch/UKB_Data/'
bgen_file_prefix <- 'ukb_imp_chr#_HRConly'

refSNP <- CACNA1D_SNPs$MarkerName[which.min(CACNA1D_SNPs$`P.value`)]
refdosageMatrix <- as.vector(LoadBgen(path,bgen_file_prefix,refSNP))

otherSNPS <- setdiff(CACNA1D_SNPs$MarkerName,refSNP)
r2 <- pbmcapply::pbmclapply(otherSNPS,function(x){
  tryCatch({
    curDosage <- as.vector(LoadBgen(path,bgen_file_prefix,x));
    cor(refdosageMatrix,curDosage) ^ 2
  },error=function(e){
    NA
  }
  )
},mc.cores = 16)

r2_df <- data.frame(snp1=otherSNPS,snp2=rep(refSNP,length(otherSNPS)),dprime=rep(0,length(otherSNPS)),rsquare=unlist(r2))
data.table::fwrite(r2_df,file='~/bsu_scratch/Locus_Zoom/CACNA1D_r2.txt',row.names = F,sep = ' ',quote=F)
