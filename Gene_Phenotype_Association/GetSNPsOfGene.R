####################################################################################
#Collects all SNPs within a gene region, their rsids and positions
#Input: Gene name and upstream/downstream distance
#Output: Info of all SNPs in the region, as a data frame
####################################################################################
library(RSQLite)
library(dplyr)
library(dbplyr)
library(biomaRt)
#Get all SNPs in a gene region, upstream_dist and downstream_dist from the start site and stop site of the gene respectively. 
GetFlankingSNPs <- function(chr,gene_position_start,gene_position_end,upstream_dist,downstream_dist){
  #Connect to annotation database
  anno_sql_name<- "all_snp_stats.sqlite"
  path <- '~/bsu_scratch/SQL/'
  setwd(path)
  anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
  anno_db <- tbl(anno_con,'all_snp_stats')
  
  #Start of flanking region
  start_of_region <- gene_position_start - upstream_dist
  #End of flanking region
  end_of_region <- gene_position_end + downstream_dist
  
  snps_within_region <- dplyr::filter(anno_db,as.numeric(chromosome)==chr & as.numeric(position) >= start_of_region & as.numeric(position) <= end_of_region) %>% dplyr::select(chromosome,position,rsid,minor_allele_frequency,info,alleleA,alleleB,alleleA_frequency,alleleB_frequency) %>% collect()
  RSQLite::dbDisconnect(anno_con)
  return(snps_within_region)   
}
#Get genomic position, given a gene name. 
GetGenePosition <- function(gene_name){
  #UKB data is based on the GRCh37 reference genome. All SNP positions are based on this. 
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  #Get gene chr and position information
  gene_pos <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters = 'hgnc_symbol', values =gene_name, mart = ensembl)
  return(gene_pos)
}
#Main function, calls GetGenePosition() to get genomic position, 
#then calls GetFlankingSNPs() to get all SNPs in the gene region
AllSNPsOfGene <- function(gene_name,upstream_dist,downstream_dist){
  print('Getting Gene Pos')
  gene_pos <- GetGenePosition(gene_name)
  print('Getting Flanking SNPs')
  flanking_snps <- GetFlankingSNPs(gene_pos$chromosome_name,gene_pos$start_position,gene_pos$end_position,upstream_dist,downstream_dist)
  return(flanking_snps)
}