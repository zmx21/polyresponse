library(qqman)
library(RSQLite)
library(dplyr)
library(pbmcapply)
#Connect to rsid annotation database. 
anno_sql_name<- "all_snp_stats.sqlite"
path <- '~/bsu_scratch/SQL/'
setwd(path)
anno_con <- RSQLite::dbConnect(SQLite(), dbname = anno_sql_name)
anno_db <- tbl(anno_con,'all_snp_stats')


#Read in results
results_tbl <- data.table::fread(paste0('~/bsu_scratch/UKB_Data/rs1262894_sbp_eur_age_sex_bmi_PC_med5/all_chr.txt'),header = T)
results_tbl$fdr <- p.adjust(results_tbl$p_int,method = 'fdr')

sig_results <- results_tbl %>% dplyr::filter(fdr<0.9)
#Join with annotation information
sig_results_annotation <- anno_db %>% dplyr::select(rsid,chromosome,position) %>% dplyr::filter(rsid %in% sig_results$rsid) %>% collect()
sig_results <- dplyr::left_join(sig_results,sig_results_annotation,by=c('rsid'='rsid'))

#Mahattan plot
manhattan(data.frame(BP=as.numeric(sig_results$position),CHR=as.numeric(sig_results$chromosome),P=sig_results$p_int,SNP=sig_results$rsid),suggestiveline = F,genomewideline = -log10(5e-8),annotateTop = T,ylab='-log(p-value)',annotatePval = 1e-8)

chr12 <- sig_results %>% dplyr::filter(as.numeric(chromosome)==12)
manhattan(data.frame(BP=as.numeric(chr12$position),CHR=as.numeric(chr12$chromosome),P=chr12$p_int,SNP=chr12$rsid),suggestiveline = F,genomewideline = -log10(5e-8),annotateTop = T,ylab='-log(p-value)',annotatePval = 1e-8,main='chr12')

#QQPlot
#Convert p values to chi-sq.
# PToChiSq <- function(p,df){
#   return(qt(p/2,df) ^2)
# }
# results_tbl_anno <- anno_db %>% dplyr::select(rsid,chromosome) %>% dplyr::filter(rsid %in% results_tbl$rsid) %>% collect()


qq.chisq <- function (x, df = 1, x.max, main = "QQ plot", sub = paste("Expected distribution: chi-squared (", 
                                                          df, " df)", sep = ""), xlab = "Expected", ylab = "Observed", 
          conc = c(0.025, 0.975), overdisp = FALSE, trim = 0.5, slope.one = FALSE, 
          slope.lambda = FALSE, pvals = FALSE, thin = c(0.25, 50), 
          oor.pch = 24, col.shade = "gray",ylim=c(),xlim=c(), ...) 
{
  shade <- function(x1, y1, x2, y2, color = col.shade) {
    n <- length(x2)
    polygon(c(x1, x2[n:1]), c(y1, y2[n:1]), border = NA, 
            col = color)
  }
  obsvd <- sort(x, na.last = NA)
  N <- length(obsvd)
  top <- obsvd[N]
  if (missing(x.max)) {
    Np <- N
  }
  else {
    Np <- sum(obsvd <= abs(x.max))
    if (Np < N || x.max < 0) 
      top <- abs(x.max)
  }
  if (Np == 0) 
    stop("Nothing to plot")
  if (df == 2) {
    expctd <- 2 * cumsum(1/(N:1))
  }
  else {
    expctd <- qchisq(p = (1:N)/(N + 1), df = df)
  }
  if (!is.null(conc)) {
    if (conc[1] > 0) {
      e.low <- qchisq(p = qbeta(conc[1], 1:N, N:1), df = df)
    }
    else {
      e.low <- rep(0, N)
    }
    if (conc[2] < 1) {
      e.high <- qchisq(p = qbeta(conc[2], 1:N, N:1), df = df)
    }
    else {
      e.high <- 1.1 * rep(max(x), N)
    }
  }
  right <- expctd[N]
  par(mar = c(5, 4, 4, 2) + 0.1, las = 1)
  if (pvals) {
    mlp <- floor(-log10(pchisq(top, df = df, lower.tail = FALSE)))
    if (mlp > 0) {
      gap <- ceiling(mlp/5)
      lp.vals <- seq(mlp, 1, -gap)
      chi2.vals <- qchisq(10^(-lp.vals), df = df, lower.tail = FALSE)
      par(mar = c(5, 4, 4, 4) + 0.1, las = 1)
    }
    else pvals <- FALSE
  }
  if(length(xlim)==0 & length(ylim)==0){
    plot(c(0, right), c(0, top), type = "n", xlab = xlab, ylab = ylab, 
         main = main, sub = sub)
    
  }else if(length(xlim)==0){
    plot(c(0, right), c(0, top), type = "n", xlab = xlab, ylab = ylab, 
         main = main, sub = sub,ylim=ylim)
    
  }else if(length(ylim)==0){
    plot(c(0, right), c(0, top), type = "n", xlab = xlab, ylab = ylab, 
         main = main, sub = sub,xlim=xlim)
    
  }else{
    plot(c(0, right), c(0, top), type = "n", xlab = xlab, ylab = ylab, 
         main = main, sub = sub,ylim=ylim,xlim=xlim)
  }
  if (pvals) {
    nvals <- length(lp.vals)
    for (i in 1:nvals) axis(side = 4, at = chi2.vals[i], 
                            labels = substitute(10^{
                              a
                            }, list(a = -lp.vals[i])), xaxt = "n")
    mtext("P-value", side = 4, line = 3, las = 0, padj = 0)
  }
  if (is.na(thin[1])) {
    show <- 1:Np
  }
  else if (length(thin) != 2 || thin[1] < 0 || thin[1] > 1 || 
           thin[2] < 1) {
    warning("invalid thin parameter; no thinning carried out")
    show <- 1:Np
  }
  else {
    space <- right * thin[1]/floor(thin[2])
    iat <- round((N + 1) * pchisq(q = (1:floor(thin[2])) * 
                                    space, df = df))
    if (max(iat) > thin[2]) 
      show <- unique(c(iat, (1 + max(iat)):Np))
    else show <- 1:Np
  }
  Nu <- floor(trim * N)
  if (Nu > 0) 
    lambda <- mean(obsvd[1:Nu])/mean(expctd[1:Nu])
  if (!is.null(conc)) {
    if (Np < N) 
      vert <- c(show, (Np + 1):N)
    else vert <- show
    if (overdisp) 
      shade(expctd[vert], lambda * e.low[vert], expctd[vert], 
            lambda * e.high[vert])
    else shade(expctd[vert], e.low[vert], expctd[vert], e.high[vert])
  }
  points(expctd[show], obsvd[show], ...)
  if (Np < N) {
    over <- (Np + 1):N
    points(expctd[over], rep(top, N - Np), pch = oor.pch)
  }
  line.types <- c("solid", "dashed", "dotted")
  key <- NULL
  txt <- NULL
  if (slope.one) {
    key <- c(key, line.types[1])
    txt <- c(txt, "y = x")
    abline(a = 0, b = 1, lty = line.types[1])
  }
  if (slope.lambda && Nu > 0) {
    key <- c(key, line.types[2])
    txt <- c(txt, paste("y = ", format(lambda, digits = 4), 
                        "x", sep = ""))
    if (!is.null(conc)) {
      if (Np < N) 
        vert <- c(show, (Np + 1):N)
      else vert <- show
    }
    abline(a = 0, b = lambda, lty = line.types[2])
  }
  if (!is.null(key)) 
    legend(0, top, legend = txt, lty = key)
  c(N = N, omitted = N - Np, lambda = lambda)
}

# rand_sample <- sample(1:nrow(results_tbl),size = 1000000,replace = F)
qq.chisq(results_tbl$t_int^2,main = 'rs1262894 all chr (Age,Sex,BMI,5PC)')

# results_tbl$chi_sq <- unlist(pbmclapply(results_tbl$p,function(x) PToChiSq(x,460209),mc.cores = 30),use.names = F)
# qq.chisq(results_tbl$chi_sq,main = 'rs1262894 chr10 EUR (Age,Sex,BMI)')
# 
# maf_info_filtered <-dplyr::filter(anno_db,rsid %in% results_tbl$rsid & minor_allele_frequency > 0.05 & info > 0.5) %>% dplyr::select(rsid) %>% collect() %>% dplyr::left_join(results_tbl,by=c('rsid'='rsid'))
# qq.chisq(maf_info_filtered$chi_sq,main = 'rs1262894 chr10 EUR MAF and Info Filtered')
# 

# qq.chisq(results_tbl$t_int^2,main = 'rs1262894 chr10 EUR \n (Age,Sex,BMI,5 PC)')
