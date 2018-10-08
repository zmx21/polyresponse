library(httr)
library(jsonlite)
library(xml2)

# GetPairwiseLD <- function(SNP1,SNP2){
#   #Use GRCh37 reference genome
#   server <- "http://grch37.rest.ensembl.org/"
#   #Generate query
#   ext <- paste0("/ld/human/pairwise/",SNP1,"/",SNP2,"?population_name=1000GENOMES:phase_3:ALL")
#   #Sent query as GET request
#   r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
#   stop_for_status(r)
#   
#   result <- head(fromJSON(toJSON(content(r))))
#   return(as.numeric(unlist(result$r2,use.names = F)))
# }
GetPairwiseLD <- function(SNP1,SNP2){
  result <- system(command = paste0("curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldpair?var1=",SNP1,"&var2=",SNP2,"&pop=EUR'"),intern = T,ignore.stdout = T,ignore.stderr = T)
  r2 <- result[grepl('R2',result)]
  r2 <- as.numeric(unlist(strsplit(r2,split = ':'))[2])
  return(r2)
}
