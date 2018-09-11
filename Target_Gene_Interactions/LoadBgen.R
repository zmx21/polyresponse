sample_file_prefix <- 'ukbb_eur_all_sbp'
bgen_file_prefix <- 'sbp_chr10'
library(rbgen)
#Load sample file
samplesTbl <- read.table(file = paste0(sample_file_prefix,'.sample'),header = F,stringsAsFactors = F)
colnames(samplesTbl) <- samplesTbl[1,]#Add header as first row
samplesTbl <- samplesTbl[-c(1,2),]#Remove first two rows

#Load bgen file
genotype_data <- rbgen::bgen.load(filename = paste0(bgen_file_prefix,'.bgen'),
                                  ranges = data.frame(chromosome=10,start=1,end=1.35e8))
