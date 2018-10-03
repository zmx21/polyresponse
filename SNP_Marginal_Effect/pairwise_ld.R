library(httr)
library(jsonlite)
library(xml2)

GetPairwiseLD <- function(SNP1,SNP2){
  #Use GRCh37 reference genome
  server <- "http://grch37.rest.ensembl.org/"
  #Generate query
  ext <- paste0("/ld/human/pairwise/",SNP1,"/",SNP2,"?population_name=1000GENOMES:phase_3:EUR")
  #Sent query as GET request
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  
  result <- head(fromJSON(toJSON(content(r))))
  return(as.numeric(unlist(result$r2,use.names = F)))
}
