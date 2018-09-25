RunDosageCalculation <- function(path,bgen_file_prefix,chr,out_path){
  chunkSize=50
  
  library(data.table)
  library(pbmcapply)
  source('/home/zmx21/MRC_BSU_Internship/Target_Gene_Interactions/LoadBgen.R')
  
  #Generate prefix
  bgen_file_prefix <- gsub(pattern = '#',replacement = chr,x = bgen_file_prefix)
  
  #Find all rsids, and generate chunks to read.
  print('Loading rsID')
  allRSIds <- FindAllRSIds(path,bgen_file_prefix)
  allRSIds <- unique(allRSIds$rsid)
  rsIDChunks <- split(allRSIds,seq(length(allRSIds)-1)%/%chunkSize)
  
  #Run chunks in parallel
  print(paste0('Writing Chunks: ',length(rsIDChunks),' chunks'))
  # pb <- progress_bar$new(
  #   format = " writing [:bar] :percent eta: :eta",
  #   total = length(rsIDChunks), clear = FALSE,force = T,show_after = 0)
  junk <- pbmclapply(1:length(rsIDChunks),function(i) {
    # pb$tick()
    currentRSIdChunk <- rsIDChunks[[i]]
    genotype_data <- LoadBgen(path,bgen_file_prefix,currentRSIdChunk)
    
    #Calculate allele dosage based on genotype probability
    alleleProbMatrix <- genotype_data$data
    dosageMatrix <- alleleProbMatrix[,,'g=1'] + 2*alleleProbMatrix[,,'g=2']
    dosageMatrix <- round(dosageMatrix,2)
    data.table::fwrite(as.data.frame(dosageMatrix),
                       file = paste0(out_path,'chunk',i,'.csv'),
                       append = F,quote = F,row.names = T,col.names = F)
  },mc.cores=16,ignore.interactive = T)
}
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  #supply default values
  path <-  '/mrc-bsu/scratch/zmx21/UKB_Data/'
  sample_file_prefix <- 'ukbb_eur_all_sbp'
  bgen_file_prefix <- 'ukb_imp_chr#_HRConly'
  chr <- '10'
  out_path <- '~/bsu_scratch/UKB_Data/dosage_matrix/'
  
}else if (length(args)!=5){
  stop("You need to supply:\n",
       "# 1: Input Path\n",
       '# 2: Sample File Prefix\n',
       "# 3: Bgen File Rrefix\n",
       "# 4: Chromosome to Analyze\n",
       "# 5: Out Path\n",
       "Exiting...", call.=FALSE)
  
}else{
  print('All Arguments Supplied')
  str <- c("# 1: Input Path:",
           '# 2: Sample File Prefix:',
           "# 3: Bgen File Prefix:",
           "# 4: Chromosome to Analyze:",
           "# 5: Out Path:")
  variables <- c('path','sample_file_prefix','bgen_file_prefix','chr','out_path')
  for(i in 1:length(args)){
    eval(parse(text=paste0(variables[i],'=',"'",args[[i]],"'")))
    print(paste0(str[i],"   ",args[i]))
  }
  system(paste0('mkdir -p ',paste0(out_path,'/chr',chr,'/')))
  system(paste0('chmod a+rwx ',paste0(out_path,'/chr',chr,'/')))
  
}
RunDosageCalculation(path,bgen_file_prefix,chr,paste0(out_path,'/chr',chr,'/'))
  