library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(readxl)
library(OUTRIDER)
library(splines)
library(ggplot2)
library(ggsci)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install('OUTRIDER')

## Our NFR datset
ods <- readRDS("/storage/liu/home/u247123/drop_analysis/drop_batch_2023Mar/drop-udn-20230314-ultima/output/processed_results/aberrant_expression/v39/outrider/RNA/ods.Rds")

###test data provided by OUTRIDER

baseURL <- paste0("https://static-content.springer.com/esm/","art%3A10.1038%2Fncomms15824/MediaObjects/")
count_URL <- paste0(baseURL, "41467_2017_BFncomms15824_MOESM390_ESM.txt")
anno_URL <- paste0(baseURL, "41467_2017_BFncomms15824_MOESM397_ESM.txt")
ctsTable <- read.table(count_URL, sep="\t") %>% filter(rowSums(.) > 1000)
annoTable <- read.table(anno_URL, sep="\t", header=TRUE)
annoTable$sampleID <- annoTable$RNA_ID
ods <- OutriderDataSet(countData=ctsTable, colData=annoTable)

#ods <- filterExpression(ods,minCounts=TRUE, filterGenes=TRUE)
## preprocessing only for test data from matrix
ods <- plotCountCorHeatmap(ods, colGroups=c("SEX"),
                           normalized=FALSE, nRowCluster=4)
## our data
ods <- plotCountCorHeatmap(ods, colGroups=c("bin"),
                           normalized=FALSE, nRowCluster=4)

### autoencoder correction
ods <- estimateSizeFactors(ods)
ods_autoencoder <- controlForConfounders(ods, q=21, iterations=5)
hist(theta(ods_autoencoder))
# test data
plotCountCorHeatmap(ods_autoencoder, colGroups=c("SEX"),
                    normalized=TRUE, nRowCluster=4)
# our data
plotCountCorHeatmap(ods_autoencoder, colGroups=c("bin"),
                    normalized=TRUE, nRowCluster=4)
## power analysis
plotExpressionRank(ods_autoencoder, "AGRN", basePlot=TRUE)
plotPowerAnalysis(ods_autoencoder)


##customized power analysis
parametricDispersionFit <- function (means, disps){
  coefs <- c(0.1, 1)
  iter <- 0
  while (TRUE) {
    residuals <- disps/(coefs[1] + coefs[2]/means)
    good <- which((residuals > 1e-04) & (residuals < 15))
    suppressWarnings({
      fit <- glm(disps[good] ~ I(1/means[good]), 
                 family=Gamma(link="identity"), start=coefs)
    })
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if (!all(coefs > 0)){
      warning("Parametric dispersion fit failed.", 
              " Using last working coefficients:", 
              paste0(round(oldcoefs, 3), sep=", "))
      coefs <- oldcoefs
      break
    }
    if ((sum(log(coefs/oldcoefs)^2) < 1e-06) & fit$converged) 
      break
    iter <- iter + 1
    if (iter > 100) {
      warning("Dispersion fit did not converge after 100 ",
              "iterations. We stopped here.")
      break
    }
  }
  names(coefs) <- c("asymptDisp", "extraPois")
  ans <- function(q) coefs[1] + coefs[2]/q
  attr(ans, "coefficients") <- coefs
  ans
}

getDispEstsData <- function(ods, mu=NULL){
  if(is.null(theta(ods))){
    stop('Please fit the ods first. ods <- fit(ods)')
  }
  odsMu <- rowMeans(counts(ods, normalized=TRUE))
  if(is.null(mu)){
    mu <- odsMu
  }
  theta <- theta(ods)
  xidx <- 10^(seq.int(max(-5,log10(min(mu))-1), log10(max(mu))+0.1, 
                      length.out = 500))
  
  # fit DESeq2 parametric Disp Fit
  fit <- parametricDispersionFit(mu, 1/theta)
  pred <- fit(xidx)
  return(list(
    mu=mu,
    disp=theta,
    xpred=xidx,
    ypred=pred,
    fit=fit
  ))
}

plotPowerAnalysis_1 <- function(ods){
  dispfit <-getDispEstsData(ods)
  m <- 10^seq.int(0,4.5,length.out = 1E4)
  d <- 1/dispfit$fit(m)
  dt<-rbindlist(lapply(c(0,0.2,0.5,2,5), function(frac)
    data.table(mean=m, disp=d, frac=frac,
               pVal=pmin(0.5, pnbinom(round(frac * m), mu = m, size=d),
                         1 - pnbinom(round(frac * m), mu = m, size=d) +
                           dnbinom(round(frac * m), mu = m, size=d))
    )))
  
  dt[,negLog10pVal:=-log10(pVal)]
  dt[,Fraction:=as.factor(frac)]
  dt[,ExprType:= ifelse(frac<1, 'Downregulation', 'Overexpression')]
  ggplot(dt, aes(mean, negLog10pVal, col=Fraction, linetype=ExprType)) +
    geom_smooth(method=lm, formula = y ~ bs(x, 10), se = FALSE) +
    scale_x_log10(breaks=c(1,10,100,1000,5000,10000)) +
    labs(x="Gene read counts", y='-log10(P-value)',color='Expression level',linetype='Expression change') + 
    ylim(0,15) + 
    geom_hline(yintercept = 5.585,linetype = "dashed", color = "black", linewidth = 0.6)+
    theme_classic()
}

plotPowerAnalysis_ashg <- function(ods){
  dispfit <-getDispEstsData(ods)
  m <- 10^seq.int(0,4.5,length.out = 1E4)
  d <- 1/dispfit$fit(m)
  dt<-rbindlist(lapply(c(0,0.2,0.5,2,5), function(frac)
    data.table(mean=m, disp=d, frac=frac,
               pVal=pmin(0.5, pnbinom(round(frac * m), mu = m, size=d),
                         1 - pnbinom(round(frac * m), mu = m, size=d) +
                           dnbinom(round(frac * m), mu = m, size=d))
    )))
  
  dt[,negLog10pVal:=-log10(pVal)]
  dt[,Fraction:=as.factor(frac)]
  dt[,ExprType:= ifelse(frac<1, 'Downregulation', 'Overexpression')]
  ggplot(dt, aes(mean, negLog10pVal, col=Fraction, linetype=ExprType)) +
    geom_smooth(method=lm, formula = y ~ bs(x, 10), se = FALSE) +
    scale_x_log10(breaks=c(1,10,100,1000,5000,10000)) +
    labs(x="Gene read counts", y='-log10(P-value)',color='Expression level',linetype='Expression change') + 
    ylim(0,15) + 
    geom_hline(yintercept = 5.585,linetype = "dashed", color = "black", linewidth = 0.6)+
    xlab(NULL)+
    theme_classic()+
    theme(text = element_text(size = 23),
          axis.text.x = element_text(angle = 45, hjust = 1))
}
plotPowerAnalysis_ashg(ods_autoencoder)

## get specific value
dispfit <-getDispEstsData(ods_autoencoder)
m <- 10^seq.int(0,4.5,length.out = 1E4)
d <- 1/dispfit$fit(m)
dt<-rbindlist(lapply(c(0,0.1,0.2,0.3,0.5, 2,5,10), function(frac)
  data.table(mean=m, disp=d, frac=frac,
             pVal=pmin(0.5, pnbinom(round(frac * m), mu = m, size=d),
                       1 - pnbinom(round(frac * m), mu = m, size=d) +
                         dnbinom(round(frac * m), mu = m, size=d))
  )))

