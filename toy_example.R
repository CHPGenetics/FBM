args = commandArgs(trailingOnly=TRUE)

require(MCMCpack)
library(doParallel)
library(foreach)
library(cluster)
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(scales)

wd <- args[1]
setwd(wd)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste0(args[2],"/FBM.cpp")) #put in the path of source code FBM.cpp here

# Note ##
# Order of elements in the data list must be followed for Rcpp to process: e.g. missing indices be before observed indices.

data <- list()
data$Y <- read.table("example_pheno.txt",header=T)
rownames(data$Y) <- data$Y[,1]
data$Y <- as.numeric(data$Y[,-1])

data$geneObs <- read.table("example_gene_expression.txt",header=T)
rownames(data$geneObs) <- data$geneObs[,1]
data$geneObs <- as.matrix(data$geneObs[,-1])

data$methyObs <- read.table("example_methylation_m_value.txt",header=T)
rownames(data$methyObs) <- data$methyObs[,1]
data$methyObs <- as.matrix(data$methyObs[,-1])


data$mapMethy<- read.table("example_map_methylation_to_gene.txt", header=T)
rownames(data$mapMethy) <- data$mapMethy[,1]
data$mapMethy <- as.numeric(data$mapMethy[,-1])


data$C <- read.table("example_clinical.txt", header=T)
rownames(data$C) <- data$C[,1]
data$C <- as.matrix(data$C[,-1])

data$missGeneIdx<- read.table("example_observed_gene_sample_indices.txt", header=T)
rownames(data$missGeneIdx) <- data$missGeneIdx[,1]
data$missGeneIdx <- as.numeric(data$missGeneIdx[,-1])

data$obsGeneIdx <- data$missGeneIdx
data$missGeneIdx <- setdiff( (1:length(data$Y)), data$obsGeneIdx)


data$missMethyIdx<- read.table("example_observed_methylation_sample_indices.txt", header=T)
rownames(data$missMethyIdx) <- data$missMethyIdx[,1]
data$missMethyIdx <- as.numeric(data$missMethyIdx[,-1])

data$obsMethyIdx <- data$missMethyIdx
data$missMethyIdx <- setdiff( (1:length(data$Y)), data$obsMethyIdx)



results <- FBMcpp(data, 1000, 15217, "Full")

plot(1:1000, results$gammaM[,1], type="l",col="firebrick", ylab="Gene Effect Estimates",xlab="Iterations")
lines(1:1000, results$gammaM[,2], col="dodgerblue")

