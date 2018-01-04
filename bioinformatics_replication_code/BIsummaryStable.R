args = commandArgs(trailingOnly=TRUE)

require(MCMCpack)
library(doParallel)
library(foreach)
library(cluster)
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(scales)

#Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
#sourceCpp("~/BI/BI.cpp")

localTest <- 0
singleTest <- 0
multipleReplicates <- 0
fullRun <- 0

if(args[7] == "local"){
    localTest <- 1
} else if (args[7] == "single"){
    singleTest <- 1
} else if (args[7] == "multiple"){
    multipleReplicates <- 1
} else if (args[7] == "full"){
    fullRun <- 1
}

Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

simulateData <- function(nSample, nGene, nMethy, nC, nCausalGene, exchangeable=FALSE, MCAR, epsilon, epsilonGeneMbar, thetaGene, thetaMethy, seed){
    set.seed(seed)
    
    gene <- matrix(NA, nSample,nGene) # NxJ
    methy <- matrix(NA, nSample, nMethy) # NxK
    geneMTrue <- matrix(NA, nSample, nGene)

    if(nC == 0){
        print("Warning: no clinical variables, dummy variable will be created")
    }
    C <- matrix(NA, nSample,nC) #clinical # NxL
    
    if(nMethy < nGene){
        stop("Not all methylation sites mapped to gene")
    }
    if(thetaGene+thetaMethy > 0.99){
        stop("Missing proportion greater than 99%.")
    }
    mapMethy <- sort(c(1:nGene, sample(1:nGene,nMethy-nGene,replace=TRUE))) #mapping methlyation to gene (sorted)

    missGene <- rep(NA, nSample) # 1 means missing, note it means this sample missed entire gene exprs
    
    nuGene <- runif(nGene, -5, 5) #mean, exchangeable
    etaGene <- runif(nGene, 0.5,2) #sd, exchangeable
    
                                        #epsilonGeneMbar <- 1.2^epsilon #we generate geneMbar ~ N(0, epsilonGeneMbar) is sigmak
                                        #  epsilonGeneMbar <- epsilon
                                        #  nuMethy <- runif(nMethy, 0.7, 0.9)
                                        #  etaMethy <- runif(nMethy, 0.15,0.25)
                                        #nuMethy <- rep(0.8, nMethy) #We shift methy value to mean zero first
                                        #    etaMethy <- rep(0.2, nMethy)
    nuMethy <- rep(0, nMethy) #we assume all mean and var of methylation are same
    etaMethy <- rep(1, nMethy)

    nuC <- rep(0,nC)
    etaC <- rep(2, nC)

    if(nC!=0){
        for(l in 1:nC){
            C[,l] <- rnorm(nSample,nuC[l],etaC[l])
        }
    }else{
        nC <-1
        C <- matrix(0, nSample, nC)
    }
    
    for(j in 1:nMethy){
                                        #methy[,j] <- rtruncnorm(nSample,a=0, b=1, mean=nuMethy[j], etaMethy[j])
        methy[,j] <- rnorm(nSample,mean=nuMethy[j], etaMethy[j]) #fix
                                        #methy[,j] <- rtruncnorm(nSample,a=-0.8, b=0.8, mean=nuMethy[j], etaMethy[j])
    }
    
                                        # Model true parameters and true features
    if(nCausalGene <0){
        causalGeneMIdx <- sort(sample(1:nGene, round(nGene/10)))
        causalGeneMbarIdx <- sort(sample(1:nGene, round(nGene/10)))
    }else{
        causalGeneMIdx <- sort(sample(1:nGene, nCausalGene))
        causalGeneMbarIdx <- sort(sample(1:nGene, nCausalGene))
    }
    
    omegaTrue <- rep(5, nGene) #for any methylation towards gene
    
    gammaTrue <- c(0, 10, 10, 10) # beta0, betaGeneM, betaGeneMbar, betaC(any clinical)
                                        #  gammaTrue <- c(0, 2, 2/(epsilon), 2)

    effectiveMethy <- rep(0, nMethy) #Adding sparsity, 0 means sparse, 1 means methy is effetive on gene.
    for (k in 1:nGene){
        effectiveMethy[which(mapMethy==k)][1]=1 #make sure at least one methy is effective per gene.
    }
                                        #    effectiveMethy[which(effectiveMethy==0)] <- rbinom(nMethy-nGene, 1, 1) #0% sparsity on the rest of the methy.
    effectiveMethy <- rep(1, nMethy) # We use all methy    
    
    if(exchangeable==TRUE){
        for (k in 1:nGene){
            gene[,k] <- rnorm(nSample, nuGene[k], etaGene[k])
        }
    } else {
        for (k in 1:nGene){
            geneMTrue[,k] <- omegaTrue[k] * apply(as.matrix(methy[,which(mapMethy==k & effectiveMethy==1)]), 1, sum) 
            tmpGeneMbarTrue <- rnorm(nSample, 0, epsilonGeneMbar) # when geneMbar is purely residual
                                        #            tmpGeneMbarTrue <- rnorm(nSample, 0.5*apply(as.matrix(C),1,sum), epsilonGeneMbar) #fix, when geneMbar is an dependent on C term
                                        #            omega0tmp <- rnorm(1,0,1)
                                        #            tmpGeneMbarTrue <- rnorm(nSample, omega0tmp*apply(as.matrix(C),1,sum), epsilonGeneMbar)
            gene[,k] <- geneMTrue[,k] + tmpGeneMbarTrue
        }
    }

    tmpMean  <- gammaTrue[1] + gammaTrue[2]*apply(as.matrix(geneMTrue[,causalGeneMIdx]),1,sum) + gammaTrue[3]*apply( as.matrix(gene[,causalGeneMbarIdx]-geneMTrue[,causalGeneMbarIdx]),1,sum) + gammaTrue[4]*apply(as.matrix(C),1,sum)
    
    Y <- rnorm(nSample, tmpMean, epsilon)

    if(nSample < 2){
        stop("More than 1 sample is needed")
    }
    if(MCAR==TRUE){
        obsGeneIdx <- 0
        while(length(obsGeneIdx) <2){
            missGene <- rbinom(nSample, 1, thetaGene)
            missGeneIdx <- which(missGene==1)
            obsGeneIdx <- setdiff(1:nSample, missGeneIdx)
        }        
        missMethy <- rep(0,nSample)
        missMethy[obsGeneIdx] <- rbinom(length(obsGeneIdx), 1, thetaMethy/(1-thetaGene))
    } else {# phenotype is a reasonable variable associated with missingness
        psiTrue <- c(NA, NA)
        psiTrue[2] <- 1/mean(Y)
        psiTrue[1] <- qnorm(thetaGene) - 1
        missGene <- rbinom(nSample, 1, pnorm(psiTrue[1] + psiTrue[2]*Y))
        missGeneIdx <- which(missGene==1)
        obsGeneIdx <- setdiff(1:nSample, missGeneIdx)

        psiTrue[1] <- qnorm(thetaMethy) - 1
        missMethy <- rep(0,nSample)
        missMethy[obsGeneIdx] <- rbinom(length(obsGeneIdx), 1,  min(1,pnorm(psiTrue[1] + psiTrue[2]*Y)/(1-thetaGene)))
    }
    missMethyIdx <- which(missMethy==1)
    obsMethyIdx <- setdiff(1:nSample, missMethyIdx)

    if( length(obsMethyIdx) <2){
        stop("More than 1 sample is needed")
    }
    geneObs <- as.matrix(gene[obsGeneIdx,])
    methyObs <- as.matrix(methy[obsMethyIdx,])
    colnames(methy) <- paste("methy", seq(1:nMethy), ".",seq(1:nMethy), ".gene", mapMethy, ".",seq(1:nMethy) , sep="" )
    rownames(methy) <- paste("sample",seq(1:nSample),sep="")
    colnames(geneObs) <- paste("gene",1:nGene, sep="")
    rownames(geneObs) <- paste("sample",obsGeneIdx,sep="")

    colnames(C) <- paste("clinical",1:nC, sep="")
    rownames(C) <- paste("sample",seq(1:nSample),sep="")
    
                                        #write.csv(methy, paste(wd, "/methyData.csv", sep=""), row.names=F,quote=F )
                                        #write.csv(geneObs, paste(wd, "/geneData.csv",sep=""),row.names=F, quote=F )
                                        #write.csv(C, paste(wd, "/clinicalData.csv", sep=""), row.names=F, quote=F )
    write.csv(methy, paste(wd, "/methyData.csv", sep=""), quote=F ) # in future write out methyObs
    write.csv(geneObs, paste(wd, "/geneData.csv",sep=""), quote=F )
    write.csv(C, paste(wd, "/clinicalData.csv", sep=""), quote=F )
    print("Done simulating data")
    
    return(list(Y=Y, geneObs=geneObs, methyObs=methyObs, mapMethy=mapMethy, C=C, missGeneIdx=missGeneIdx, obsGeneIdx=obsGeneIdx, missMethyIdx=missMethyIdx, obsMethyIdx=obsMethyIdx, causalGeneMIdx=causalGeneMIdx, causalGeneMbarIdx=causalGeneMbarIdx, geneTrue=gene, geneMTrue=geneMTrue, methyTrue=methy, gammaTrue=gammaTrue, omegaTrue=omegaTrue, effectiveMethy=effectiveMethy))

}

BFDR <- function(Ipm, n, nGene){ #Ipm is the feature selection posterior mean.
    BFDRtmp <- matrix(NA,nGene,2) #t, BFDR(t)
    PMtmp <- rep(NA, nGene)
    qtmp <- rep(NA, nGene)
                                        #   PMtmp <- 2*( 1- pt( abs((est[[r]]$gammaMEst) / (est[[r]]$gammaMSd / sqrt(all_n[n]))), all_n[n]-1) )
                                        #PMtmp <- 2*( 1- pt( abs(mean) / (sd / sqrt(n)), n-1) )
    PMtmp <- 1-Ipm 
    sortedPM <- sort(PMtmp)
    for(k in 1:nGene){
        BFDRtmp[k,1] <- sortedPM[k]
        BFDRtmp[k,2] <- sum(sortedPM[1:k])/k
    }
    for(k in 1:nGene){
        qtmp[k] <- min( BFDRtmp[ which(BFDRtmp[,1] >= PMtmp[k]), 2] )
    }

    return(list(q=qtmp))
}
#######################################################
### Bayesian Imputation (BI) with gene/methy missing ##
#######################################################
                                        #now using Rcpp





nGene <- as.numeric(args[8])
nCausalGene <-  as.numeric(args[9])
nMethy <-  as.numeric(args[10])
nC <- 2
exchangeable <- FALSE
MCAR <- TRUE

if(localTest){
    wd <- "~/Research/Bayesian_Missing/BINew/localTest"
    setwd(wd)
    nSample <- 200
    epsilon <- 0.1
    epsilonGeneMbar <- 2
    nItr <- 5000
    seed <- 2334
    thetaGene <- 0.1
    thetaMethy <- 0.1
    
    data <- simulateData(nSample, nGene, nMethy,nC, nCausalGene, exchangeable, MCAR, epsilon, epsilonGeneMbar, thetaGene, thetaMethy, seed)
    
    cl <- makeCluster(1)
    registerDoParallel(cl)
    getDoParWorkers()
    
    result <- BI(data, nItr, seed, "BI")
    
    stopCluster(cl)

    est <- estimating(data,result,wd)
}

if(singleTest){
    wd <- args[1]
    setwd(wd)
    
    nItr <- as.numeric(args[4])
    seed <- as.numeric(args[3])
    thetaGene <- as.numeric(args[5])

    thetaMethy <- as.numeric(args[6])

    cl <- makeCluster(as.numeric(args[2]))
    registerDoParallel(cl)
    getDoParWorkers()

    nSample <- 200
    epsilon <- 0.5
    epsilonGeneMbar <- 2
    
    data <- simulateData(nSample, nGene, nMethy,nC, nCausalGene, exchangeable, MCAR, epsilon, epsilonGeneMbar, thetaGene, thetaMethy, seed)
    print(paste("Start single Test. With ","thetaGene: ", thetaGene, ", nSample: ", nSample, ", number of complete case: ", length(data$obsGeneIdx), ", number of missing case: ", length(data$missGeneIdx), sep=" "))
    
    result <- BI(data, nItr, seed, "BI")
    
    stopCluster(cl)

    est <- estimating(data,result,wd)
}

if( fullRun | multipleReplicates ){
    wd <- args[1]
    setwd(wd)

    noCFlag <- 0 #flag for no clinical variable input
    nCtmp <- nC
    if(nC==0){
        noCFlag=1
        nCtmp <- 1
    }
    nItr <- as.numeric(args[4])
    seed <- as.numeric(args[3])
    thetaGene <- as.numeric(args[5])
    thetaMethy <- as.numeric(args[6])

                                        #    cl <- makeCluster(as.numeric(args[2]))
                                        #    registerDoParallel(cl)

                                        #    all_n <- c(20,30,50,100,200,500)
#    all_n <- c(50,100,200,500)
    all_n <- c(50,100,200,500)
    if(args[12] == "strong"){
        epsilon <- 0.01
    } else if( args[12]=="mid"){
        epsilon <- 0.5
    } else if( args[12]=="weak"){
        epsilon <- 5
    }else if( args[12]=="midstrong"){
        epsilon <- 0.1
    }else if( args[12]=="superweak"){
        epsilon <- 8
    }else if( args[12]=="ultraweak"){
        epsilon <- 10
    }

    epsilonGeneMbar <- 2

    if( multipleReplicates ){
        all_n <- 100
    }

    parNames <- c( expression(paste(gamma^M," causal")),  expression(paste(gamma^bar(M), " causal")),  expression(paste(gamma^M, " non-causal")),  expression(paste(gamma^bar(M), " non-causal")),  expression(gamma^C),  expression(omega), expression(hat(Y)))
    nPar <- length(parNames) #  gammaM_causal, gammaMbar_causal,gammaM_nc, gammaMbar_nc, gammaC, omegaMean, yHat

    sensSpecNames <- c( expression(paste(I^M,"sens")),  expression(paste(I^M,"spec")),  expression(paste(I^bar(M),"sens")),expression(paste(I^bar(M),"spec")), expression("Overall sens"),  expression("Overall spec"))

    methods <- c("CC", "BI", "Full")
    
    nReplicates <-  as.numeric(args[11]) #50 replicates (MCMC chains)

                                        #    colors <- c(colors()[461],  colors()[282], colors()[555], colors()[610])
    colors<- c(colors()[26],  colors()[261], colors()[35], colors()[614],  colors()[621]) #blue, grey, red, green, yellow.
    
    tmp <- rep(NA, nReplicates*nCausalGene*length(all_n))
    gammaMAll <- array(tmp, c(length(all_n),nReplicates, nCausalGene )) # #n/e,replicates, ngenes  so: 5,20,nCausalGene
    gammaMbarAll <- array(tmp, c(length(all_n),nReplicates, nCausalGene))
    gammaMAll_sd<- array(tmp, c(length(all_n),nReplicates, nCausalGene )) # #n/e,replicates, ngenes  so: 5,20,nCausalGene
    gammaMbarAll_sd<- array(tmp, c(length(all_n),nReplicates, nCausalGene))

    gammaMAll_CC <- array(tmp, c(length(all_n),nReplicates, nCausalGene )) # #n/e,replicates, ngenes  so: 5,20,nCausalGene
    gammaMbarAll_CC <- array(tmp, c(length(all_n),nReplicates, nCausalGene))
    gammaMAll_sd_CC <- array(tmp, c(length(all_n),nReplicates, nCausalGene )) # #n/e,replicates, ngenes  so: 5,20,nCausalGene
    gammaMbarAll_sd_CC <- array(tmp, c(length(all_n),nReplicates, nCausalGene))

    gammaMAll_Full <- array(tmp, c(length(all_n),nReplicates, nCausalGene )) # #n/e,replicates, ngenes  so: 5,20,nCausalGene
    gammaMbarAll_Full <- array(tmp, c(length(all_n),nReplicates, nCausalGene))
    gammaMAll_sd_Full <- array(tmp, c(length(all_n),nReplicates, nCausalGene )) # #n/e,replicates, ngenes  so: 5,20,nCausalGene
    gammaMbarAll_sd_Full <- array(tmp, c(length(all_n),nReplicates, nCausalGene))

    tmp <- rep(NA, nReplicates*(nGene-nCausalGene)*length(all_n)) 
    gammaMAll_nc<- array(tmp, c(length(all_n),nReplicates, nGene-nCausalGene )) # #none causal
    gammaMbarAll_nc<- array(tmp, c(length(all_n),nReplicates,  nGene-nCausalGene))
    gammaMAll_sd_nc<- array(tmp, c(length(all_n),nReplicates, nGene-nCausalGene )) # #none causal
    gammaMbarAll_sd_nc<- array(tmp, c(length(all_n),nReplicates,  nGene-nCausalGene))

    gammaMAll_nc_CC<- array(tmp, c(length(all_n),nReplicates, nGene-nCausalGene )) # #none causal
    gammaMbarAll_nc_CC<- array(tmp, c(length(all_n),nReplicates,  nGene-nCausalGene))
    gammaMAll_sd_nc_CC<- array(tmp, c(length(all_n),nReplicates, nGene-nCausalGene )) # #none causal
    gammaMbarAll_sd_nc_CC<- array(tmp, c(length(all_n),nReplicates,  nGene-nCausalGene))

    gammaMAll_nc_Full<- array(tmp, c(length(all_n),nReplicates, nGene-nCausalGene )) # #none causal
    gammaMbarAll_nc_Full<- array(tmp, c(length(all_n),nReplicates,  nGene-nCausalGene))
    gammaMAll_sd_nc_Full<- array(tmp, c(length(all_n),nReplicates, nGene-nCausalGene )) # #none causal
    gammaMbarAll_sd_nc_Full<- array(tmp, c(length(all_n),nReplicates,  nGene-nCausalGene))


    tmp <- rep(NA, nReplicates*nCtmp*length(all_n)) 
    gammaCAll <- array(tmp, c(length(all_n),nReplicates, nCtmp))
    gammaCAll_sd <- array(tmp, c(length(all_n),nReplicates, nCtmp))
    gammaCAll_CC <- array(tmp, c(length(all_n),nReplicates, nCtmp))
    gammaCAll_sd_CC <- array(tmp, c(length(all_n),nReplicates, nCtmp))
    gammaCAll_Full <- array(tmp, c(length(all_n),nReplicates, nCtmp))
    gammaCAll_sd_Full <- array(tmp, c(length(all_n),nReplicates, nCtmp))

    tmp <- rep(NA, nReplicates*nGene*length(all_n)) 
    omegaAll <- array(tmp, c(length(all_n),nReplicates, 1))  #we only summarize on the mean of omega
    omegaAll_sd <- array(tmp, c(length(all_n),nReplicates, 1))
    omegaAll_CC <- array(tmp, c(length(all_n),nReplicates, 1))
    omegaAll_sd_CC <- array(tmp, c(length(all_n),nReplicates, 1))
    omegaAll_Full <- array(tmp, c(length(all_n),nReplicates, 1))
    omegaAll_sd_Full <- array(tmp, c(length(all_n),nReplicates, 1))
    
    tmp <- rep(NA, nReplicates*1*length(all_n)) 
    yHatAll <- array(tmp, c(length(all_n),nReplicates, 1))  #different from all the other, stores MSE(yHat-yTrue) for each replicate
    yHatAll_CC <- array(tmp, c(length(all_n),nReplicates, 1)) 
    yHatAll_Full <- array(tmp, c(length(all_n),nReplicates, 1)) 

    tableMSE <- matrix(NA, 8, 3) #rows are: gammaM_causal, gammaMbar_causal,gammaM_nc, gammaMbar_nc, gammaC, omegaMean, yHat
                                        # MSE mean(se) CC BI BIp
    tableCI <- matrix(NA, 6, 3)        # coverage(credible interval width)
    tableFS <- matrix(NA, 4, 3)        # Feature selection: sens spec

                                        # coverage(interval width)(BI),  coverage(credible interval width) (CC)
    tableMSE[1,] <- c("\ ", "MSE(SE) CC", "BI ")
    tableMSE[,1] <- c("\ ", "$\\gamma^M_c$", "$\\gamma^{\\overline M}_c$", "$\\gamma^M_{nc}$", "$\\gamma^{\\overline M}_{nc}$", "$\\gamma^C$" ,"$\\omega$", "$\\hat{Y}$")

    tableCI[1,] <- c("\ ", "Coverage(Width) CC", "BI")
    tableCI[,1] <- c("\ ", "$\\gamma^M_c$", "$\\gamma^{\\overline M}_c$", "$\\gamma^M_{nc}$", "$\\gamma^{\\overline M}_{nc}$", "$\\gamma^C$" )

    tableFS[1,] <- c("\ ", "Sens(Spec) CC", "BI ")
    tableFS[,1] <- c("\ ", "$\\gamma^M$", "$\\gamma^{\\overline M}$", "Overall")   

                                        #   as for IM, IMbar, the sensitivity and specificity are plotted out directly

                                        # the 4 in confusion matrix correspond to
                                        #    T F
                                        #Pos 1 3
                                        #Neg 2 4
    IMConfusion <- matrix(0, length(all_n),4)
    IMbarConfusion <- matrix(0, length(all_n),4)
    IMConfusion_CC<- matrix(0, length(all_n),4) 
    IMbarConfusion_CC<- matrix(0, length(all_n),4)
    IMConfusion_Full<- matrix(0, length(all_n),4) 
    IMbarConfusion_Full<- matrix(0, length(all_n),4)

    sensSpec <- matrix(NA, length(all_n), length(sensSpecNames)) #IM sens, IM spec, IMbar sens, IMbar spec, overall
    sensSpec_CC<- matrix(NA, length(all_n), length(sensSpecNames)) #IM sens, IM spec, IMbar sens, IMbar spec
    sensSpec_Full<- matrix(NA, length(all_n), length(sensSpecNames)) #IM sens, IM spec, IMbar sens, IMbar spec

                                        # for a different FDR cutoff
    IMConfusion2 <- matrix(0, length(all_n),4)
    IMbarConfusion2 <- matrix(0, length(all_n),4)
    IMConfusion2_CC<- matrix(0, length(all_n),4) 
    IMbarConfusion2_CC<- matrix(0, length(all_n),4)
    IMConfusion2_Full<- matrix(0, length(all_n),4) 
    IMbarConfusion2_Full<- matrix(0, length(all_n),4)

    IMConfusionAll <- matrix(0, length(all_n),4)
    IMbarConfusionAll <- matrix(0, length(all_n),4)
    IMConfusionAll_CC<- matrix(0, length(all_n),4) 
    IMbarConfusionAll_CC<- matrix(0, length(all_n),4)
    IMConfusionAll_Full<- matrix(0, length(all_n),4) 
    IMbarConfusionAll_Full<- matrix(0, length(all_n),4)



    sensSpec2 <- matrix(NA, length(all_n), length(sensSpecNames)) #IM sens, IM spec, IMbar sens, IMbar spec, overall
    sensSpec2_CC<- matrix(NA, length(all_n), length(sensSpecNames)) #IM sens, IM spec, IMbar sens, IMbar spec
    sensSpec2_Full<- matrix(NA, length(all_n), length(sensSpecNames)) #IM sens, IM spec, IMbar sens, IMbar spec


    
    truthAll <- rep(NA, nPar)

    meanAll <- matrix(NA, nPar, length(all_n))
    tmp <- rep(NA, nReplicates*nPar*length(all_n)) 
    biasAll <- array(tmp, c(nPar, length(all_n), 2)) # Used to be   biasAll <- matrix(NA, nPar, length(all_n)),  now store mean(bias) and sd(bias)
    mse <- array(tmp, c(nPar, length(all_n), 2))#used to be mse <- matrix(NA,nPar, length(all_n)), now stores mean(MSE) and sd(MSE)
    coverage <- matrix(0, nPar-2, length(all_n))
    width <- matrix(0, nPar-2, length(all_n))
    
    meanAllCC <- matrix(NA, nPar, length(all_n))
    biasAllCC <- array(tmp, c(nPar, length(all_n), 2))
    mseCC <- array(tmp, c(nPar, length(all_n), 2))
    coverageCC <- matrix(0, nPar-2, length(all_n))
    widthCC <- matrix(0, nPar-2, length(all_n))

    meanAllFull <- matrix(NA, nPar, length(all_n))
    biasAllFull <- array(tmp, c(nPar, length(all_n), 2))
    mseFull <- array(tmp, c(nPar, length(all_n), 2))
    coverageFull <- matrix(0, nPar-2, length(all_n))
    widthFull <- matrix(0, nPar-2, length(all_n))


    bfdr <- matrix(NA, 2, length(all_n)) #M, Mbar
    bfdrCC <- matrix(NA, 2, length(all_n))
    bfdrFull <- matrix(NA, 2, length(all_n))
    
    bi <- list()
    biCC<- list()
    biFull <- list()
    est <- list()
    estCC <- list()
    estFull <- list()
    data <- list()
    dataCC <- list()
    dataFull <- list()
    testData <- list() # same test data for BI, CC and FUll
    # ROC, Youden, FDR, and TP
    pdf(paste(wd,"/ROC.pdf",sep=""))
    par(mfrow=c(1,1),mar=c(5.1,4.1,1.9,2.1),pin=c(4,4))
    
    seg <- 200 #200 points (bins) for 1-spec
    for (n in 1:length(all_n)){
      fdr = matrix(NA,length(all_n),2) #0.05 and 0.1
      fdrCC = matrix(NA, length(all_n),2)
      fdrFull = matrix(NA, length(all_n),2)
      
      sens <- list() #changed
      sensCC <- list()
      sensFull <- list()
      spec <- list() #store 1-spec
      specCC <- list()
      specFull <- list()
      
      tp <- list()
      tpCC <- list() #changed
      tpFull <- list()
      
      roc = matrix(NA,seg+2,4) #sens lower, mean, upper, 1-spec mean (S.E instead of sd) #changed
      rocCC = matrix(NA, seg+2,4)
      rocFull = matrix(NA, seg+2,4)
      
      roc[seg+1,] = c(0,0,0,0)
      roc[seg+2,] = c(1,1,1,1)
      rocCC[seg+1,] = c(0,0,0,0) #could also be like 0.1,0,0,0,0,0
      rocCC[seg+2,] = c(1,1,1,1)
      rocFull[seg+1,] = c(0,0,0,0) #could also be like 0.1,0,0,0,0,0
      rocFull[seg+2,] = c(1,1,1,1)
      
      youden = rep(NA,seg)
      youdenCC = rep(NA, seg)
      youdenFull = rep(NA, seg)
      
      youden[c(1,seg)] = c(0,0)
      youdenCC[c(1,seg)] = c(0,0)     
      youdenFull[c(1,seg)] = c(0,0)
      
        for(r in 1:nReplicates){
            data[[r]] <- simulateData(all_n[n], nGene, nMethy,nC, nCausalGene, exchangeable, MCAR, epsilon,epsilonGeneMbar, thetaGene, thetaMethy, seed+r)
            dataCC[[r]] <- data[[r]]
            dataCC[[r]]$missGeneIdx <- integer(0)
            dataCC[[r]]$missMethyIdx <- integer(0)
            missJointIdx <- sort(c( data[[r]]$missGeneIdx, data[[r]]$missMethyIdx))
            obsJointIdx <- setdiff(1:all_n[n], missJointIdx)
            dataCC[[r]]$obsMethyIdx <- 1:length(obsJointIdx)            
            dataCC[[r]]$obsGeneIdx <- 1:length(obsJointIdx)
            dataCC[[r]]$Y <- dataCC[[r]]$Y[obsJointIdx]
            dataCC[[r]]$methyObs <- as.matrix(dataCC[[r]]$methyTrue[obsJointIdx,])
            dataCC[[r]]$geneObs <- as.matrix(dataCC[[r]]$geneTrue[obsJointIdx,])
            dataCC[[r]]$C <- as.matrix(dataCC[[r]]$C[obsJointIdx,])
            dataCC[[r]]$geneMTrue <- as.matrix(dataCC[[r]]$geneMTrue[obsJointIdx,])
            dataCC[[r]]$geneTrue <- as.matrix(dataCC[[r]]$geneTrue[obsJointIdx,])
            
            dataFull[[r]] <- data[[r]]
            dataFull[[r]]$missGeneIdx <- integer(0)
            dataFull[[r]]$missMethyIdx <- integer(0)
            dataFull[[r]]$obsMethyIdx <- 1:length(data[[r]]$Y)            
            dataFull[[r]]$obsGeneIdx <- 1:length(data[[r]]$Y)
            dataFull[[r]]$methyObs <- as.matrix(data[[r]]$methyTrue)
            dataFull[[r]]$geneObs <- as.matrix(data[[r]]$geneTrue)
        }

#        for(r in 1:nReplicates){
#            set.seed(seed)
#            bi[[r]] <- BIcpp(data[[r]], nItr, seed, "BI")
#        }
        for (r in 1:nReplicates){
            wdSub <- paste(wd,"/n",n,"_r",r,sep="")
#            dir.create(wdSub)
#            est[[r]] <- estimating(data[[r]],bi[[r]],wdSub,testData[[r]])
#            plotting(data[[r]], bi[[r]], wdSub)
            load(paste(wdSub,"/run.RData",sep=""))
            gammaMAll[n, r, ] <- est[[r]]$gammaMEst[data[[r]]$causalGeneMIdx]
            gammaMbarAll[n, r, ] <- est[[r]]$gammaMbarEst[data[[r]]$causalGeneMbarIdx]
            gammaMAll_nc[n, r, ] <- est[[r]]$gammaMEst[-data[[r]]$causalGeneMIdx]
            gammaMbarAll_nc[n, r, ] <- est[[r]]$gammaMbarEst[-data[[r]]$causalGeneMbarIdx]
            gammaCAll[n, r, ] <- est[[r]]$gammaCEst
            omegaAll[n, r, ] <- est[[r]]$omegaMeanEst

            gammaMAll_sd[n, r, ] <- est[[r]]$gammaMSd[data[[r]]$causalGeneMIdx]
            gammaMbarAll_sd[n, r, ] <- est[[r]]$gammaMbarSd[data[[r]]$causalGeneMbarIdx]
            gammaMAll_sd_nc[n, r, ] <- est[[r]]$gammaMSd[-data[[r]]$causalGeneMIdx]
            gammaMbarAll_sd_nc[n, r, ] <- est[[r]]$gammaMbarSd[-data[[r]]$causalGeneMbarIdx]
            gammaCAll_sd[n, r, ] <- est[[r]]$gammaCSd
            omegaAll_sd[n, r, ] <- est[[r]]$omegaSd


            yHatAll[n,r,] <- est[[r]]$MSEyHatEst #MSE(yHat-yTrue)

            cut <- 0.05  #BFDR 0.05
                                        #          qtmp <- BFDR(est[[r]]$IMEst, all_n[n], nGene)$q
            qtmp <- est[[r]]$qM
            qtmpM <- qtmp
            IMConfusion[n,1] <- IMConfusion[n,1] + length(which(qtmp[data[[r]]$causalGeneMIdx] < cut))
            IMConfusion[n,2] <- IMConfusion[n,2] + length(which(qtmp[data[[r]]$causalGeneMIdx] >= cut))
            IMConfusion[n,3] <- IMConfusion[n,3] + length(which(qtmp[-data[[r]]$causalGeneMIdx] < cut))
            IMConfusion[n,4] <- IMConfusion[n,4] + length(which(qtmp[-data[[r]]$causalGeneMIdx] >= cut))
            
            cut <- 0.1   #BFDR 0.1
            IMConfusion2[n,1] <- IMConfusion2[n,1] + length(which(qtmp[data[[r]]$causalGeneMIdx] < cut))
            IMConfusion2[n,2] <- IMConfusion2[n,2] + length(which(qtmp[data[[r]]$causalGeneMIdx] >= cut))
            IMConfusion2[n,3] <- IMConfusion2[n,3] + length(which(qtmp[-data[[r]]$causalGeneMIdx] < cut))
            IMConfusion2[n,4] <- IMConfusion2[n,4] + length(which(qtmp[-data[[r]]$causalGeneMIdx] >= cut))

            cut <- 0.05
                                        #            qtmp <- BFDR(est[[r]]$IMbarEst, all_n[n], nGene)$q
                        qtmp <- est[[r]]$qMbar
            qtmpMbar <- qtmp
            IMbarConfusion[n,1] <- IMbarConfusion[n,1] + length(which(qtmp[data[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion[n,2] <- IMbarConfusion[n,2] + length(which(qtmp[data[[r]]$causalGeneMbarIdx] >= cut))
            IMbarConfusion[n,3] <- IMbarConfusion[n,3] + length(which(qtmp[-data[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion[n,4] <- IMbarConfusion[n,4] + length(which(qtmp[-data[[r]]$causalGeneMbarIdx] >= cut))
            cut <- 0.1
            IMbarConfusion2[n,1] <- IMbarConfusion2[n,1] + length(which(qtmp[data[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion2[n,2] <- IMbarConfusion2[n,2] + length(which(qtmp[data[[r]]$causalGeneMbarIdx] >= cut))
            IMbarConfusion2[n,3] <- IMbarConfusion2[n,3] + length(which(qtmp[-data[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion2[n,4] <- IMbarConfusion2[n,4] + length(which(qtmp[-data[[r]]$causalGeneMbarIdx] >= cut))

            ######## calculating sens, spec, TP/T
            qtmpMsort <- sort(qtmpM)
            qtmpMidx <- order(qtmpM)
            qtmpMbarsort <- sort(qtmpMbar)
            qtmpMbaridx <- order(qtmpMbar)
            
            tp[[r]] <- rep(0,2*nGene)
            sens[[r]] <- tp[[r]]
            spec[[r]] <- tp[[r]]
            tpTmp <- matrix(0,2,nGene) #M and Mbar
            specTmp <- tpTmp
            
            tpTmp[1,1] <- 1*(qtmpMidx[1] %in% data[[r]]$causalGeneMIdx)
            tpTmp[2,1] <- 1*(qtmpMbaridx[1] %in% data[[r]]$causalGeneMbarIdx)
            tp[[r]][1:2] <- c(tpTmp[1,1], tpTmp[1,1]+tpTmp[2,1])
            sens[[r]][1:2] <- tp[[r]][1:2]
            
            specTmp[1,1] <- abs(1*(qtmpMidx[1] %in% data[[r]]$causalGeneMIdx)-1)
            specTmp[2,1] <- abs(1*(qtmpMbaridx[1] %in% data[[r]]$causalGeneMbarIdx)-1)
            spec[[r]][1:2] <- c(specTmp[1,1], specTmp[1,1]+specTmp[2,1])
            
            for (i in 2:nGene){
              tpTmp[1,i] <- tpTmp[1,i-1] + 1*(qtmpMidx[i] %in% data[[r]]$causalGeneMIdx)
              tpTmp[2,i] <- tpTmp[2,i-1] + 1*(qtmpMbaridx[i] %in% data[[r]]$causalGeneMbarIdx)
              specTmp[1,i] <- specTmp[1,i-1] + abs(1*(qtmpMidx[i] %in% data[[r]]$causalGeneMIdx)-1)
              specTmp[2,i] <- specTmp[2,i-1] + abs(1*(qtmpMbaridx[i] %in% data[[r]]$causalGeneMbarIdx)-1)
            }
            for (i in 2:nGene){
              if(tpTmp[1,i] == tpTmp[1,i-1]){
                tp[[r]][2*i-1] <- tp[[r]][2*i-2]
              }else{
                tp[[r]][2*i-1] <- tp[[r]][2*i-2]+1
              }
              if(tpTmp[2,i] == tpTmp[2,i-1]){
                tp[[r]][2*i] <- tp[[r]][2*i-1]
              }else{
                tp[[r]][2*i] <- tp[[r]][2*i-1]+1
              }
              if(specTmp[1,i] == specTmp[1,i-1]){
                spec[[r]][2*i-1] <- spec[[r]][2*i-2]
              }else{
                spec[[r]][2*i-1] <- spec[[r]][2*i-2]+1
              }
              if(specTmp[2,i] == specTmp[2,i-1]){
                spec[[r]][2*i] <- spec[[r]][2*i-1]
              }else{
                spec[[r]][2*i] <- spec[[r]][2*i-1]+1
              }
            }
            tp[[r]] <- tp[[r]]/(2*length(data[[r]]$causalGeneMbarIdx))
            sens[[r]] <- tp[[r]]
            spec[[r]] <- spec[[r]] /(2*nGene-2*length(data[[r]]$causalGeneMbarIdx))
            #finished finding TP/T and sens, spec
        }

#        for(r in 1:nReplicates){
#            set.seed(seed)
#            biCC[[r]] <- BIcpp(dataCC[[r]], nItr, seed, "CC")
#        }
        for (r in 1:nReplicates){
            wdSub <- paste(wd,"/n",n,"_r",r,"_CC",sep="")
#            dir.create(wdSub)
 #           estCC[[r]] <- estimating(dataCC[[r]],biCC[[r]],wdSub,testData[[r]])
  #          plotting(dataCC[[r]], biCC[[r]], wdSub)
            load(paste(wdSub,"/run.RData",sep=""))

            gammaMAll_CC[n, r, ] <- estCC[[r]]$gammaMEst[dataCC[[r]]$causalGeneMIdx]
            gammaMbarAll_CC[n, r, ] <- estCC[[r]]$gammaMbarEst[dataCC[[r]]$causalGeneMbarIdx]
            gammaMAll_nc_CC[n, r, ] <- estCC[[r]]$gammaMEst[-dataCC[[r]]$causalGeneMIdx]
            gammaMbarAll_nc_CC[n, r, ] <- estCC[[r]]$gammaMbarEst[-dataCC[[r]]$causalGeneMbarIdx]
            gammaCAll_CC[n, r, ] <- estCC[[r]]$gammaCEst
            omegaAll_CC[n, r, ] <- estCC[[r]]$omegaMeanEst

            gammaMAll_sd_CC[n, r, ] <- estCC[[r]]$gammaMSd[dataCC[[r]]$causalGeneMIdx]
            gammaMbarAll_sd_CC[n, r, ] <- estCC[[r]]$gammaMbarSd[dataCC[[r]]$causalGeneMbarIdx]
            gammaMAll_sd_nc_CC[n, r, ] <- estCC[[r]]$gammaMSd[-dataCC[[r]]$causalGeneMIdx]
            gammaMbarAll_sd_nc_CC[n, r, ] <- estCC[[r]]$gammaMbarSd[-dataCC[[r]]$causalGeneMbarIdx]
            gammaCAll_sd_CC[n, r, ] <- estCC[[r]]$gammaCSd
            omegaAll_sd_CC[n, r, ] <- estCC[[r]]$omegaSd

            yHatAll_CC[n,r,] <- estCC[[r]]$MSEyHatEst 

            cut <- 0.05  #BFDR 0.05
                                        #            qtmp <- BFDR(estCC[[r]]$IMEst, all_n[n], nGene)$q
            qtmp <- estCC[[r]]$qM
            qtmpM <- qtmp
            IMConfusion_CC[n,1] <- IMConfusion_CC[n,1] + length(which(qtmp[dataCC[[r]]$causalGeneMIdx] < cut))
            IMConfusion_CC[n,2] <- IMConfusion_CC[n,2] + length(which(qtmp[dataCC[[r]]$causalGeneMIdx] >= cut))
            IMConfusion_CC[n,3] <- IMConfusion_CC[n,3] + length(which(qtmp[-dataCC[[r]]$causalGeneMIdx] < cut))
            IMConfusion_CC[n,4] <- IMConfusion_CC[n,4] + length(which(qtmp[-dataCC[[r]]$causalGeneMIdx] >= cut))
            
            cut <- 0.1   #BFDR 0.1
            IMConfusion2_CC[n,1] <- IMConfusion2_CC[n,1] + length(which(qtmp[dataCC[[r]]$causalGeneMIdx] < cut))
            IMConfusion2_CC[n,2] <- IMConfusion2_CC[n,2] + length(which(qtmp[dataCC[[r]]$causalGeneMIdx] >= cut))
            IMConfusion2_CC[n,3] <- IMConfusion2_CC[n,3] + length(which(qtmp[-dataCC[[r]]$causalGeneMIdx] < cut))
            IMConfusion2_CC[n,4] <- IMConfusion2_CC[n,4] + length(which(qtmp[-dataCC[[r]]$causalGeneMIdx] >= cut))

            cut <- 0.05
                                        #            qtmp <- BFDR(estCC[[r]]$IMbarEst, all_n[n], nGene)$q
            qtmp <- estCC[[r]]$qMbar
            qtmpMbar <- qtmp
            IMbarConfusion_CC[n,1] <- IMbarConfusion_CC[n,1] + length(which(qtmp[dataCC[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion_CC[n,2] <- IMbarConfusion_CC[n,2] + length(which(qtmp[dataCC[[r]]$causalGeneMbarIdx] >= cut))
            IMbarConfusion_CC[n,3] <- IMbarConfusion_CC[n,3] + length(which(qtmp[-dataCC[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion_CC[n,4] <- IMbarConfusion_CC[n,4] + length(which(qtmp[-dataCC[[r]]$causalGeneMbarIdx] >= cut))
            cut <- 0.1
            IMbarConfusion2_CC[n,1] <- IMbarConfusion2_CC[n,1] + length(which(qtmp[dataCC[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion2_CC[n,2] <- IMbarConfusion2_CC[n,2] + length(which(qtmp[dataCC[[r]]$causalGeneMbarIdx] >= cut))
            IMbarConfusion2_CC[n,3] <- IMbarConfusion2_CC[n,3] + length(which(qtmp[-dataCC[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion2_CC[n,4] <- IMbarConfusion2_CC[n,4] + length(which(qtmp[-dataCC[[r]]$causalGeneMbarIdx] >= cut))
            
            ######## calculating sens, spec, TP/T
            qtmpMsort <- sort(qtmpM)
            qtmpMidx <- order(qtmpM)
            qtmpMbarsort <- sort(qtmpMbar)
            qtmpMbaridx <- order(qtmpMbar)
            
            tpCC[[r]] <- rep(0,2*nGene)
            sensCC[[r]] <- tpCC[[r]]
            specCC[[r]] <- tpCC[[r]]
            tpTmp <- matrix(0,2,nGene) #M and Mbar
            specTmp <- tpTmp
            
            tpTmp[1,1] <- 1*(qtmpMidx[1] %in% dataCC[[r]]$causalGeneMIdx)
            tpTmp[2,1] <- 1*(qtmpMbaridx[1] %in% dataCC[[r]]$causalGeneMbarIdx)
            tpCC[[r]][1:2] <- c(tpTmp[1,1], tpTmp[1,1]+tpTmp[2,1])
            sensCC[[r]][1:2] <- tpCC[[r]][1:2]
            
            specTmp[1,1] <- abs(1*(qtmpMidx[1] %in% dataCC[[r]]$causalGeneMIdx)-1)
            specTmp[2,1] <- abs(1*(qtmpMbaridx[1] %in% dataCC[[r]]$causalGeneMbarIdx)-1)
            specCC[[r]][1:2] <- c(specTmp[1,1], specTmp[1,1]+specTmp[2,1])
            
            for (i in 2:nGene){
              tpTmp[1,i] <- tpTmp[1,i-1] + 1*(qtmpMidx[i] %in% dataCC[[r]]$causalGeneMIdx)
              tpTmp[2,i] <- tpTmp[2,i-1] + 1*(qtmpMbaridx[i] %in% dataCC[[r]]$causalGeneMbarIdx)
              specTmp[1,i] <- specTmp[1,i-1] + abs(1*(qtmpMidx[i] %in% dataCC[[r]]$causalGeneMIdx)-1)
              specTmp[2,i] <- specTmp[2,i-1] + abs(1*(qtmpMbaridx[i] %in% dataCC[[r]]$causalGeneMbarIdx)-1)
            }
            for (i in 2:nGene){
              if(tpTmp[1,i] == tpTmp[1,i-1]){
                tpCC[[r]][2*i-1] <- tpCC[[r]][2*i-2]
              }else{
                tpCC[[r]][2*i-1] <- tpCC[[r]][2*i-2]+1
              }
              if(tpTmp[2,i] == tpTmp[2,i-1]){
                tpCC[[r]][2*i] <- tpCC[[r]][2*i-1]
              }else{
                tpCC[[r]][2*i] <- tpCC[[r]][2*i-1]+1
              }
              if(specTmp[1,i] == specTmp[1,i-1]){
                specCC[[r]][2*i-1] <- specCC[[r]][2*i-2]
              }else{
                specCC[[r]][2*i-1] <- specCC[[r]][2*i-2]+1
              }
              if(specTmp[2,i] == specTmp[2,i-1]){
                specCC[[r]][2*i] <- specCC[[r]][2*i-1]
              }else{
                specCC[[r]][2*i] <- specCC[[r]][2*i-1]+1
              }
            }
            tpCC[[r]] <- tpCC[[r]]/(2*length(dataCC[[r]]$causalGeneMbarIdx))
            sensCC[[r]] <- tpCC[[r]]
            specCC[[r]] <- specCC[[r]] /(2*nGene-2*length(dataCC[[r]]$causalGeneMbarIdx))
            #finished finding TP/T and sens, spec
        }

#      for(r in 1:nReplicates){
#            set.seed(seed)
#            biFull[[r]] <- BIcpp(dataFull[[r]], nItr, seed, "BI")
#        }
        for (r in 1:nReplicates){
            wdSub <- paste(wd,"/n",n,"_r",r,"_Full",sep="")
  #          dir.create(wdSub)
 #           estFull[[r]] <- estimating(dataFull[[r]],biFull[[r]],wdSub, testData[[r]])
  #          plotting(dataFull[[r]], biFull[[r]], wdSub)
            load(paste(wdSub,"/run.RData",sep=""))
            gammaMAll_Full[n, r, ] <- estFull[[r]]$gammaMEst[dataFull[[r]]$causalGeneMIdx]
            gammaMbarAll_Full[n, r, ] <- estFull[[r]]$gammaMbarEst[dataFull[[r]]$causalGeneMbarIdx]
            gammaMAll_nc_Full[n, r, ] <- estFull[[r]]$gammaMEst[-dataFull[[r]]$causalGeneMIdx]
            gammaMbarAll_nc_Full[n, r, ] <- estFull[[r]]$gammaMbarEst[-dataFull[[r]]$causalGeneMbarIdx]
            gammaCAll_Full[n, r, ] <- estFull[[r]]$gammaCEst
            omegaAll_Full[n, r, ] <- estFull[[r]]$omegaMeanEst

            gammaMAll_sd_Full[n, r, ] <- estFull[[r]]$gammaMSd[dataFull[[r]]$causalGeneMIdx]
            gammaMbarAll_sd_Full[n, r, ] <- estFull[[r]]$gammaMbarSd[dataFull[[r]]$causalGeneMbarIdx]
            gammaMAll_sd_nc_Full[n, r, ] <- estFull[[r]]$gammaMSd[-dataFull[[r]]$causalGeneMIdx]
            gammaMbarAll_sd_nc_Full[n, r, ] <- estFull[[r]]$gammaMbarSd[-dataFull[[r]]$causalGeneMbarIdx]
            gammaCAll_sd_Full[n, r, ] <- estFull[[r]]$gammaCSd
            omegaAll_sd_Full[n, r, ] <- estFull[[r]]$omegaSd

            yHatAll_Full[n,r,] <- estFull[[r]]$MSEyHatEst 

            cut <- 0.05  #BFDR 0.05
  #          qtmp <- BFDR(estFull[[r]]$IMEst, all_n[n], nGene)$q
            qtmp <- estFull[[r]]$qM
            qtmpM <- qtmp
            IMConfusion_Full[n,1] <- IMConfusion_Full[n,1] + length(which(qtmp[dataFull[[r]]$causalGeneMIdx] < cut))
            IMConfusion_Full[n,2] <- IMConfusion_Full[n,2] + length(which(qtmp[dataFull[[r]]$causalGeneMIdx] >= cut))
            IMConfusion_Full[n,3] <- IMConfusion_Full[n,3] + length(which(qtmp[-dataFull[[r]]$causalGeneMIdx] < cut))
            IMConfusion_Full[n,4] <- IMConfusion_Full[n,4] + length(which(qtmp[-dataFull[[r]]$causalGeneMIdx] >= cut))
            
            cut <- 0.1   #BFDR 0.1
            IMConfusion2_Full[n,1] <- IMConfusion2_Full[n,1] + length(which(qtmp[dataFull[[r]]$causalGeneMIdx] < cut))
            IMConfusion2_Full[n,2] <- IMConfusion2_Full[n,2] + length(which(qtmp[dataFull[[r]]$causalGeneMIdx] >= cut))
            IMConfusion2_Full[n,3] <- IMConfusion2_Full[n,3] + length(which(qtmp[-dataFull[[r]]$causalGeneMIdx] < cut))
            IMConfusion2_Full[n,4] <- IMConfusion2_Full[n,4] + length(which(qtmp[-dataFull[[r]]$causalGeneMIdx] >= cut))

            cut <- 0.05
                                        #            qtmp <- BFDR(estFull[[r]]$IMbarEst, all_n[n], nGene)$q
                        qtmp <- estFull[[r]]$qMbar
                        qtmpMbar <- qtmp

            IMbarConfusion_Full[n,1] <- IMbarConfusion_Full[n,1] + length(which(qtmp[dataFull[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion_Full[n,2] <- IMbarConfusion_Full[n,2] + length(which(qtmp[dataFull[[r]]$causalGeneMbarIdx] >= cut))
            IMbarConfusion_Full[n,3] <- IMbarConfusion_Full[n,3] + length(which(qtmp[-dataFull[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion_Full[n,4] <- IMbarConfusion_Full[n,4] + length(which(qtmp[-dataFull[[r]]$causalGeneMbarIdx] >= cut))
            cut <- 0.1
            IMbarConfusion2_Full[n,1] <- IMbarConfusion2_Full[n,1] + length(which(qtmp[dataFull[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion2_Full[n,2] <- IMbarConfusion2_Full[n,2] + length(which(qtmp[dataFull[[r]]$causalGeneMbarIdx] >= cut))
            IMbarConfusion2_Full[n,3] <- IMbarConfusion2_Full[n,3] + length(which(qtmp[-dataFull[[r]]$causalGeneMbarIdx] < cut))
            IMbarConfusion2_Full[n,4] <- IMbarConfusion2_Full[n,4] + length(which(qtmp[-dataFull[[r]]$causalGeneMbarIdx] >= cut))
            
            ######## calculating sens, spec, TP/T
            qtmpMsort <- sort(qtmpM)
            qtmpMidx <- order(qtmpM)
            qtmpMbarsort <- sort(qtmpMbar)
            qtmpMbaridx <- order(qtmpMbar)
            
            tpFull[[r]] <- rep(0,2*nGene)
            sensFull[[r]] <- tpFull[[r]]
            specFull[[r]] <- tpFull[[r]]
            tpTmp <- matrix(0,2,nGene) #M and Mbar
            specTmp <- tpTmp
            
            tpTmp[1,1] <- 1*(qtmpMidx[1] %in% dataFull[[r]]$causalGeneMIdx)
            tpTmp[2,1] <- 1*(qtmpMbaridx[1] %in% dataFull[[r]]$causalGeneMbarIdx)
            tpFull[[r]][1:2] <- c(tpTmp[1,1], tpTmp[1,1]+tpTmp[2,1])
            sensFull[[r]][1:2] <- tpFull[[r]][1:2]
            
            specTmp[1,1] <- abs(1*(qtmpMidx[1] %in% dataFull[[r]]$causalGeneMIdx)-1)
            specTmp[2,1] <- abs(1*(qtmpMbaridx[1] %in% dataFull[[r]]$causalGeneMbarIdx)-1)
            specFull[[r]][1:2] <- c(specTmp[1,1], specTmp[1,1]+specTmp[2,1])
            
            for (i in 2:nGene){
              tpTmp[1,i] <- tpTmp[1,i-1] + 1*(qtmpMidx[i] %in% dataFull[[r]]$causalGeneMIdx)
              tpTmp[2,i] <- tpTmp[2,i-1] + 1*(qtmpMbaridx[i] %in% dataFull[[r]]$causalGeneMbarIdx)
              specTmp[1,i] <- specTmp[1,i-1] + abs(1*(qtmpMidx[i] %in% dataFull[[r]]$causalGeneMIdx)-1)
              specTmp[2,i] <- specTmp[2,i-1] + abs(1*(qtmpMbaridx[i] %in% dataFull[[r]]$causalGeneMbarIdx)-1)
            }
            for (i in 2:nGene){
              if(tpTmp[1,i] == tpTmp[1,i-1]){
                tpFull[[r]][2*i-1] <- tpFull[[r]][2*i-2]
              }else{
                tpFull[[r]][2*i-1] <- tpFull[[r]][2*i-2]+1
              }
              if(tpTmp[2,i] == tpTmp[2,i-1]){
                tpFull[[r]][2*i] <- tpFull[[r]][2*i-1]
              }else{
                tpFull[[r]][2*i] <- tpFull[[r]][2*i-1]+1
              }
              if(specTmp[1,i] == specTmp[1,i-1]){
                specFull[[r]][2*i-1] <- specFull[[r]][2*i-2]
              }else{
                specFull[[r]][2*i-1] <- specFull[[r]][2*i-2]+1
              }
              if(specTmp[2,i] == specTmp[2,i-1]){
                specFull[[r]][2*i] <- specFull[[r]][2*i-1]
              }else{
                specFull[[r]][2*i] <- specFull[[r]][2*i-1]+1
              }
            }
            tpFull[[r]] <- tpFull[[r]]/(2*length(dataFull[[r]]$causalGeneMbarIdx))
            sensFull[[r]] <- tpFull[[r]]
            specFull[[r]] <- specFull[[r]] /(2*nGene-2*length(dataFull[[r]]$causalGeneMbarIdx))
            #finished finding TP/T and sens, spec
        }

        sensSpec[n,1] <- IMConfusion[n,1]/(IMConfusion[n,1]+IMConfusion[n,2])
        sensSpec[n,2] <- IMConfusion[n,4]/(IMConfusion[n,3]+IMConfusion[n,4])
        sensSpec[n,3] <- IMbarConfusion[n,1]/(IMbarConfusion[n,1]+IMbarConfusion[n,2])
        sensSpec[n,4] <- IMbarConfusion[n,4]/(IMbarConfusion[n,3]+IMbarConfusion[n,4])
        sensSpec[n,5] <- (IMConfusion[n,1] + IMbarConfusion[n,1]) / (IMConfusion[n,1]+IMConfusion[n,2] + IMbarConfusion[n,1]+IMbarConfusion[n,2])
        sensSpec[n,6] <- (IMConfusion[n,4] + IMbarConfusion[n,4]) / (IMConfusion[n,3]+IMConfusion[n,4] + IMbarConfusion[n,3]+IMbarConfusion[n,4])

        sensSpec[which(is.na(sensSpec))] <- 0.5
        
        sensSpec_CC[n,1] <- IMConfusion_CC[n,1]/(IMConfusion_CC[n,1]+IMConfusion_CC[n,2])
        sensSpec_CC[n,2] <- IMConfusion_CC[n,4]/(IMConfusion_CC[n,3]+IMConfusion_CC[n,4])
        sensSpec_CC[n,3] <- IMbarConfusion_CC[n,1]/(IMbarConfusion_CC[n,1]+IMbarConfusion_CC[n,2])
        sensSpec_CC[n,4] <- IMbarConfusion_CC[n,4]/(IMbarConfusion_CC[n,3]+IMbarConfusion_CC[n,4])
        sensSpec_CC[n,5] <- (IMConfusion_CC[n,1] + IMbarConfusion_CC[n,1]) / (IMConfusion_CC[n,1]+IMConfusion_CC[n,2] + IMbarConfusion_CC[n,1]+IMbarConfusion_CC[n,2])
        sensSpec_CC[n,6] <- (IMConfusion_CC[n,4] + IMbarConfusion_CC[n,4]) / (IMConfusion_CC[n,3]+IMConfusion_CC[n,4] + IMbarConfusion_CC[n,3]+IMbarConfusion_CC[n,4])

        sensSpec_CC[which(is.na(sensSpec_CC))] <- 0.5

              sensSpec_Full[n,1] <- IMConfusion_Full[n,1]/(IMConfusion_Full[n,1]+IMConfusion_Full[n,2])
        sensSpec_Full[n,2] <- IMConfusion_Full[n,4]/(IMConfusion_Full[n,3]+IMConfusion_Full[n,4])
        sensSpec_Full[n,3] <- IMbarConfusion_Full[n,1]/(IMbarConfusion_Full[n,1]+IMbarConfusion_Full[n,2])
        sensSpec_Full[n,4] <- IMbarConfusion_Full[n,4]/(IMbarConfusion_Full[n,3]+IMbarConfusion_Full[n,4])
        sensSpec_Full[n,5] <- (IMConfusion_Full[n,1] + IMbarConfusion_Full[n,1]) / (IMConfusion_Full[n,1]+IMConfusion_Full[n,2] + IMbarConfusion_Full[n,1]+IMbarConfusion_Full[n,2])
        sensSpec_Full[n,6] <- (IMConfusion_Full[n,4] + IMbarConfusion_Full[n,4]) / (IMConfusion_Full[n,3]+IMConfusion_Full[n,4] + IMbarConfusion_Full[n,3]+IMbarConfusion_Full[n,4])

        sensSpec_Full[which(is.na(sensSpec_Full))] <- 0.5

        
                                        # different cutoff
        sensSpec2[n,1] <- IMConfusion2[n,1]/(IMConfusion2[n,1]+IMConfusion2[n,2])
        sensSpec2[n,2] <- IMConfusion2[n,4]/(IMConfusion2[n,3]+IMConfusion2[n,4])
        sensSpec2[n,3] <- IMbarConfusion2[n,1]/(IMbarConfusion2[n,1]+IMbarConfusion2[n,2])
        sensSpec2[n,4] <- IMbarConfusion2[n,4]/(IMbarConfusion2[n,3]+IMbarConfusion2[n,4])
        sensSpec2[n,5] <- (IMConfusion2[n,1] + IMbarConfusion2[n,1]) / (IMConfusion2[n,1]+IMConfusion2[n,2] + IMbarConfusion2[n,1]+IMbarConfusion2[n,2])
        sensSpec2[n,6] <- (IMConfusion2[n,4] + IMbarConfusion2[n,4]) / (IMConfusion2[n,3]+IMConfusion2[n,4] + IMbarConfusion2[n,3]+IMbarConfusion2[n,4])
        
        sensSpec2_CC[n,1] <- IMConfusion2_CC[n,1]/(IMConfusion2_CC[n,1]+IMConfusion2_CC[n,2])
        sensSpec2_CC[n,2] <- IMConfusion2_CC[n,4]/(IMConfusion2_CC[n,3]+IMConfusion2_CC[n,4])
        sensSpec2_CC[n,3] <- IMbarConfusion2_CC[n,1]/(IMbarConfusion2_CC[n,1]+IMbarConfusion2_CC[n,2])
        sensSpec2_CC[n,4] <- IMbarConfusion2_CC[n,4]/(IMbarConfusion2_CC[n,3]+IMbarConfusion2_CC[n,4])
        sensSpec2_CC[n,5] <- (IMConfusion2_CC[n,1] + IMbarConfusion2_CC[n,1]) / (IMConfusion2_CC[n,1]+IMConfusion2_CC[n,2] + IMbarConfusion2_CC[n,1]+IMbarConfusion2_CC[n,2])
        sensSpec2_CC[n,6] <- (IMConfusion2_CC[n,4] + IMbarConfusion2_CC[n,4]) / (IMConfusion2_CC[n,3]+IMConfusion2_CC[n,4] + IMbarConfusion2_CC[n,3]+IMbarConfusion2_CC[n,4])

      sensSpec2_Full[n,1] <- IMConfusion2_Full[n,1]/(IMConfusion2_Full[n,1]+IMConfusion2_Full[n,2])
        sensSpec2_Full[n,2] <- IMConfusion2_Full[n,4]/(IMConfusion2_Full[n,3]+IMConfusion2_Full[n,4])
        sensSpec2_Full[n,3] <- IMbarConfusion2_Full[n,1]/(IMbarConfusion2_Full[n,1]+IMbarConfusion2_Full[n,2])
        sensSpec2_Full[n,4] <- IMbarConfusion2_Full[n,4]/(IMbarConfusion2_Full[n,3]+IMbarConfusion2_Full[n,4])
        sensSpec2_Full[n,5] <- (IMConfusion2_Full[n,1] + IMbarConfusion2_Full[n,1]) / (IMConfusion2_Full[n,1]+IMConfusion2_Full[n,2] + IMbarConfusion2_Full[n,1]+IMbarConfusion2_Full[n,2])
        sensSpec2_Full[n,6] <- (IMConfusion2_Full[n,4] + IMbarConfusion2_Full[n,4]) / (IMConfusion2_Full[n,3]+IMConfusion2_Full[n,4] + IMbarConfusion2_Full[n,3]+IMbarConfusion2_Full[n,4])

        sensSpec2[which(is.na(sensSpec2))] <- 0.5
        sensSpec2_CC[which(is.na(sensSpec2_CC))] <- 0.5
  sensSpec2_Full[which(is.na(sensSpec2_Full))] <- 0.5

  ########PLOT
  auc=0
  aucCC=0
  aucFull=0
  
  roc[,4] <- c(seq(0,1,by=1/(seg+1))[2:(seg+1)],0,1)
  rocCC[,4] <- roc[,4]
  rocFull[,4] <- roc[,4]
  
  bins <- matrix(0,2,seg) #lower upper of each bin of spec
  incremental <- 1/(seg+1)
  bins[1,1] <- incremental/2
  bins[2,1] <- bins[1,1]+incremental
  
  sensSpecAll <- list()
  sensSpecAllCC <- list()
  sensSpecAllFull <- list()
  tpAll <- rep(0,2*nGene)
  tpAllCC <- tpAll
  tpAllFull <- tpAll
  
  for(i in 2:seg){
    bins[,i] <- bins[,(i-1)]+incremental
  }
  sensSpecAll[[seg]] <- integer(0)
  sensSpecAllCC[[seg]] <- integer(0)
  sensSpecAllFull[[seg]] <- integer(0)
  
  for(r in 1:nReplicates){
    tpAll <- tpAll + tp[[r]]
    tpAllCC <- tpAllCC + tpCC[[r]]
    tpAllFull <- tpAllFull + tpFull[[r]]
    
    binStart=1
    for(g in which(spec[[r]]<=bins[2,seg])){
      binEnd=binStart
      while( spec[[r]][g] > bins[2,binEnd] && binEnd<seg ){
        binEnd = binEnd+1
      }
      sensSpecAll[[binEnd]] <- c(sensSpecAll[[binEnd]], sens[[r]][g])
      binStart = binEnd # speed up computing by changing starting point
    }
    
    binStart=1
    for(g in which(specCC[[r]]<=bins[2,seg])){
      binEnd=binStart
      while( specCC[[r]][g] > bins[2,binEnd] && binEnd<seg ){
        binEnd = binEnd+1
      }
      sensSpecAllCC[[binEnd]] <- c(sensSpecAllCC[[binEnd]], sensCC[[r]][g])
      binStart = binEnd # speed up computing by changing starting point
    }
    
    binStart=1
    for(g in which(specFull[[r]]<=bins[2,seg])){
      binEnd=binStart
      while( specFull[[r]][g] > bins[2,binEnd] && binEnd<seg ){
        binEnd = binEnd+1
      }
      sensSpecAllFull[[binEnd]] <- c(sensSpecAllFull[[binEnd]], sensFull[[r]][g])
      binStart = binEnd # speed up computing by changing starting point
    }
  }
  tpAll <- tpAll/nReplicates
  tpAllCC <- tpAllCC/nReplicates
  tpAllFull <- tpAllFull/nReplicates
  
  if(is.null(sensSpecAll[[1]])){
    sensSpecAll[[1]] <- 0
    roc[1,1:3] <- 0
  }else{
    tmpMean <- mean(sensSpecAll[[1]])
    tmpSe <- sd(sensSpecAll[[1]])/sqrt(length(sensSpecAll[[1]]))
    if(length(sensSpecAll[[1]])==1){
      tmpSe <- 0
    }
    roc[1,1] <- max(tmpMean - tmpSe, 0)
    roc[1,2] <- tmpMean
    roc[1,3] <- min(tmpMean + tmpSe, 1)
  }
  for(g in 2:seg){
    if(length(sensSpecAll[[g]])==0){
      roc[g,1:3] <- roc[g-1,1:3]
    }else{
      tmpMean <- mean(sensSpecAll[[g]])
      tmpSe <- sd(sensSpecAll[[g]])/sqrt(length(sensSpecAll[[g]]))
      if(length(sensSpecAll[[g]])==1){
        tmpSe <- 0
      }
      roc[g,1] <- max(tmpMean - tmpSe, 0)
      roc[g,2] <- tmpMean
      roc[g,3] <- min(tmpMean + tmpSe, 1)
    }
  }
  ##CC
  if(is.null(sensSpecAllCC[[1]])){
    sensSpecAllCC[[1]] <- 0
    rocCC[1,1:3] <- 0
  }else{
    tmpMean <- mean(sensSpecAllCC[[1]])
    tmpSe <- sd(sensSpecAllCC[[1]])/sqrt(length(sensSpecAllCC[[1]]))
    if(length(sensSpecAllCC[[1]])==1){
      tmpSe <- 0
    }
    rocCC[1,1] <- max(tmpMean - tmpSe, 0)
    rocCC[1,2] <- tmpMean
    rocCC[1,3] <- min(tmpMean + tmpSe, 1)
  }
  for(g in 2:seg){
    if(length(sensSpecAllCC[[g]])==0){
      rocCC[g,1:3] <- rocCC[g-1,1:3]
    }else{
      tmpMean <- mean(sensSpecAllCC[[g]])
      tmpSe <- sd(sensSpecAllCC[[g]])/sqrt(length(sensSpecAllCC[[g]]))
      if(length(sensSpecAllCC[[g]])==1){
        tmpSe <- 0
      }
      rocCC[g,1] <- max(tmpMean - tmpSe, 0)
      rocCC[g,2] <- tmpMean
      rocCC[g,3] <- min(tmpMean + tmpSe, 1)
    }
  }
  
  ##Full
  if(is.null(sensSpecAllFull[[1]])){
    sensSpecAllFull[[1]] <- 0
    rocFull[1,1:3] <- 0
  }else{
    tmpMean <- mean(sensSpecAllFull[[1]])
    tmpSe <- sd(sensSpecAllFull[[1]])/sqrt(length(sensSpecAllFull[[1]]))
    if(length(sensSpecAllFull[[1]])==1){
      tmpSe <- 0
    }
    rocFull[1,1] <- max(tmpMean - tmpSe, 0)
    rocFull[1,2] <- tmpMean
    rocFull[1,3] <- min(tmpMean + tmpSe, 1)
  }
  for(g in 2:seg){
    if(length(sensSpecAllFull[[g]])==0){
      rocFull[g,1:3] <- rocFull[g-1,1:3]
    }else{
      tmpMean <- mean(sensSpecAllFull[[g]])
      tmpSe <- sd(sensSpecAllFull[[g]])/sqrt(length(sensSpecAllFull[[g]]))
      if(length(sensSpecAllFull[[g]])==1){
        tmpSe <- 0
      }
      rocFull[g,1] <- max(tmpMean - tmpSe, 0)
      rocFull[g,2] <- tmpMean
      rocFull[g,3] <- min(tmpMean + tmpSe, 1)
    }
  }
  
  maxIdx <- min( round(max(c(which(tpAll==1)[1],which(tpAllCC==1)[1],which(tpAllFull==1)[1]))*1.2), 2*nGene)
  
#  plot(0:maxIdx, c(0,tpAll[1:maxIdx]), main="TP/T vs P", xlab="Positive",ylab="Sensitivity",sub=paste("nSample",all_n[n]),col=colors[1],type="l")
      plot(0:maxIdx, c(0,tpAll[1:maxIdx]), main="TP/T vs P", xlab="Positive",ylab="Sensitivity",col=colors[1],type="l")

  abline(h=1, col=colors[2],lty=2)
  lines(0:maxIdx, c(0,tpAllCC[1:maxIdx]), col=colors[3],type="l")
  lines(0:maxIdx, c(0,tpAllFull[1:maxIdx]), col=colors[4],type="l")
  
  roc <- roc[order(roc[,4]),]
  rocCC <- rocCC[order(rocCC[,4]),]
      rocFull <- rocFull[order(rocFull[,4]),]
 
  youden <- roc[,2]-roc[,4]
  youdenCC <- rocCC[,2]-rocCC[,4]
  youdenFull <- rocFull[,2]-rocFull[,4]
  
  for( i in 1:(seg-1)){
    auc = auc + (roc[i+1,4] - roc[i,4]) * (roc[i+1,2]+roc[i,2])/2
    aucCC = aucCC + (rocCC[i+1,4] - rocCC[i,4]) * (rocCC[i+1,2]+rocCC[i,2])/2
    aucFull = aucFull + (rocFull[i+1,4] - rocFull[i,4]) * (rocFull[i+1,2]+rocFull[i,2])/2
  }

      
  auc <- signif(auc,4)
  aucCC <- signif(aucCC,4)
  aucFull <- signif(aucFull,4)
  
  #plotting bands of ROC
  #polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
  #x <- 1:10
  #y <- c(2,4,6,8,7,12,14,16,18,20)
  #lo <- loess(y~x)
  #plot(x,y)
  #lines(predict(lo), col='red', lwd=2)
  plot(0:1,0:1,main="ROC",
#       xlab="1-Specificity",ylab="Sensitivity",xaxt="n",sub=paste("nSample",all_n[n]),col=colors[2],type="l",lty=2)
       xlab="1-Specificity",ylab="Sensitivity",xaxt="n",col=colors[2],type="l",lty=2)

  axis(1,at=c(0,0.25,0.5,0.75,1),label=c(0,0.25,0.5,0.75,1),cex.axis=1)
  
  polygon(c(rev(roc[,4]),roc[,4]), c(rev(roc[,1]), roc[,3]), col = alpha(colors[1],0.2), border=NA)
  polygon(c(rev(rocCC[,4]),rocCC[,4]), c(rev(rocCC[,1]), rocCC[,3]), col = alpha(colors[3],0.2), border=NA)
  polygon(c(rev(rocFull[,4]),rocFull[,4]), c(rev(rocFull[,1]), rocFull[,3]), col = alpha(colors[4],0.2), border=NA)
  
  lines(roc[,4], roc[,2], col=colors[1],type="l")
  lines(rocCC[,4], rocCC[,2], col=colors[3],type="l")
  lines(rocFull[,4], rocFull[,2], col=colors[4],type="l")
  
  points(roc[which(youden==max(youden))[1],4], roc[which(youden==max(youden))[1],2], pch=8, cex=1.5) #changed
  points(rocCC[which(youdenCC==max(youdenCC))[1],4], rocCC[which(youdenCC==max(youdenCC))[1],2], pch=8, cex=1.5)
  points(rocFull[which(youdenFull==max(youdenFull))[1],4], rocFull[which(youdenFull==max(youdenFull))[1],2], pch=8, cex=1.5)
  legend("bottomright", bty = "n", c(paste("CC(AUC=",aucCC,")",sep=""), paste("BI(AUC=",auc,")",sep=""), paste("Full(AUC=",aucFull,")",sep="")), cex=1,  col=colors[c(3,1,4)], lty=c(1,1,1))
  
  plot(c(0,1),c(0,1)) #changed 
  text(min(max(roc[which(youden==max(youden))[1],4],0.1),0.9), min(max(roc[which(youden==max(youden))[1],2]*1.1,0.1),0.9), labels = signif(max(youden),3), col=colors[1])
  text(min(max(rocCC[which(youdenCC==max(youdenCC))[1],4],0.1),0.9),min(max(rocCC[which(youdenCC==max(youdenCC))[1],2]*1.1,0.1),0.9), labels = signif(max(youdenCC),3), col=colors[3])
  text(min(max(rocFull[which(youdenFull==max(youdenFull))[1],4],0.1),0.9),min(max(rocFull[which(youdenFull==max(youdenFull))[1],2]*1.1,0.1),0.9), labels = signif(max(youdenFull),3), col=colors[4])
  
  legend("bottomright", bty = "n", c(paste("CC(AUC=",aucCC,")",sep=""), paste("BI(AUC=",auc,")",sep=""), paste("Full(AUC=",aucFull,")",sep="")), cex=1,  col=colors[c(3,1,4)], lty=c(1,1,1))
  
  ## tp <- tp[order(tp[,1]),]
  ## tpCC <- tpCC[order(tpCC[,1]),]
  ##       tpFull <- tpFull[order(tpFull[,1]),]
  
  ## plot(tp[,1],tp[,2],main="TP/P",
  ##      xlab="Positive",ylab="TruePositive",sub=paste("nSample",all_n[n]),col=colors[1],type="l")
  ## lines(tpCC[,1], tpCC[,2], col=colors[3], type="l")
  ## lines(tpFull[,1], tpFull[,2], col=colors[4], type="l")
  
  # don't save some of the data that are not useful for plotting
  estCC$omegaEst <- integer(0)
  est$omegaEst <- integer(0)
  data$geneObs <- integer(0)
  data$methyObs <- integer(0)
  data$geneTrue <- integer(0)
  data$methyTrue <- integer(0)
  dataCC$geneObs <- integer(0)
  dataCC$methyObs <- integer(0)
  dataCC$geneTrue <- integer(0)
  dataCC$methyTrue <- integer(0)
  
  estFull$omegaEst <- integer(0)
  dataFull$geneObs <- integer(0)
  dataFull$methyObs <- integer(0)
  dataFull$geneTrue <- integer(0)
  dataFull$methyTrue <- integer(0)
      save(dataCC,dataFull,data, estCC,estFull,est,all_n, auc, aucFull,aucCC, roc, rocCC, rocFull,youden, youdenFull,youdenCC,file=paste(wd, "/data_n",n,".RData", sep=""))
    }
    dev.off()
                                        #fix, do we use a random gamma, or we average them?
                                        #for now we average them
    truthAll[1] <- data[[1]]$gammaTrue[2]
    truthAll[2] <- data[[1]]$gammaTrue[3]
    truthAll[3] <- 0
    truthAll[4] <- 0
    truthAll[5] <- data[[1]]$gammaTrue[4]
    truthAll[6] <- data[[1]]$omegaTrue[1]
    truthAll[7] <- 0 #the est is already yHat-yTrue

    for (n in 1:length(all_n)){
        tmpAll <- list()
        tmpAll[[1]] <- as.matrix(gammaMAll[n,,])
        tmpAll[[2]] <- as.matrix(gammaMbarAll[n,,])
        tmpAll[[3]] <- as.matrix(gammaMAll_nc[n,,])
        tmpAll[[4]] <- as.matrix(gammaMbarAll_nc[n,,])
        tmpAll[[5]] <- as.matrix(gammaCAll[n,,])
        tmpAll[[6]] <- as.matrix(omegaAll[n,,])

        tmpAllCC <- list()
        tmpAllCC[[1]] <- as.matrix(gammaMAll_CC[n,,])
        tmpAllCC[[2]] <- as.matrix(gammaMbarAll_CC[n,,])
        tmpAllCC[[3]] <- as.matrix(gammaMAll_nc_CC[n,,])
        tmpAllCC[[4]] <- as.matrix(gammaMbarAll_nc_CC[n,,])
        tmpAllCC[[5]] <- as.matrix(gammaCAll_CC[n,,])
        tmpAllCC[[6]] <- as.matrix(omegaAll_CC[n,,])

                tmpAllFull <- list()
        tmpAllFull[[1]] <- as.matrix(gammaMAll_Full[n,,])
        tmpAllFull[[2]] <- as.matrix(gammaMbarAll_Full[n,,])
        tmpAllFull[[3]] <- as.matrix(gammaMAll_nc_Full[n,,])
        tmpAllFull[[4]] <- as.matrix(gammaMbarAll_nc_Full[n,,])
        tmpAllFull[[5]] <- as.matrix(gammaCAll_Full[n,,])
        tmpAllFull[[6]] <- as.matrix(omegaAll_Full[n,,])

        tmpAll_sd <- list()
        tmpAll_sd[[1]] <- as.matrix(gammaMAll_sd[n,,])
        tmpAll_sd[[2]] <- as.matrix(gammaMbarAll_sd[n,,])
        tmpAll_sd[[3]] <- as.matrix(gammaMAll_sd_nc[n,,])
        tmpAll_sd[[4]] <- as.matrix(gammaMbarAll_sd_nc[n,,])
        tmpAll_sd[[5]] <- as.matrix(gammaCAll_sd[n,,])
        tmpAll_sd[[6]] <- as.matrix(omegaAll_sd[n,,])

        tmpAll_sdCC <- list()
        tmpAll_sdCC[[1]] <- as.matrix(gammaMAll_sd_CC[n,,])
        tmpAll_sdCC[[2]] <- as.matrix(gammaMbarAll_sd_CC[n,,])
        tmpAll_sdCC[[3]] <- as.matrix(gammaMAll_sd_nc_CC[n,,])
        tmpAll_sdCC[[4]] <- as.matrix(gammaMbarAll_sd_nc_CC[n,,])
        tmpAll_sdCC[[5]] <- as.matrix(gammaCAll_sd_CC[n,,])
        tmpAll_sdCC[[6]] <- as.matrix(omegaAll_sd_CC[n,,])

        tmpAll_sdFull <- list()
        tmpAll_sdFull[[1]] <- as.matrix(gammaMAll_sd_Full[n,,])
        tmpAll_sdFull[[2]] <- as.matrix(gammaMbarAll_sd_Full[n,,])
        tmpAll_sdFull[[3]] <- as.matrix(gammaMAll_sd_nc_Full[n,,])
        tmpAll_sdFull[[4]] <- as.matrix(gammaMbarAll_sd_nc_Full[n,,])
        tmpAll_sdFull[[5]] <- as.matrix(gammaCAll_sd_Full[n,,])
        tmpAll_sdFull[[6]] <- as.matrix(omegaAll_sd_Full[n,,])

        upper <- rep(NA, nPar-2) #omega and Yhat have no CI
        lower <- rep(NA, nPar-2)
        for(r in 1:nReplicates){
            for(i in 1:(nPar-2)){
                upper[i] <- tmpAll[[i]][r,1] + 1.96*tmpAll_sd[[i]][r,1] #only take the 1st one, might waste some but should still be unbiased. Furthermore, we aren't interested in coverage for now.
                lower[i] <- tmpAll[[i]][r,1] - 1.96*tmpAll_sd[[i]][r,1] 

                coverage[i,n] <- coverage[i,n] + sum( upper[i] >= truthAll[i] & lower[i] <= truthAll[i])
                width[i,n] <- width[i,n] + upper[i] - lower[i]
            }
        }

        coverage[1,n] <- coverage[1,n] / (length(data[[1]]$causalGeneMIdx)*nReplicates)
        coverage[2,n] <- coverage[2,n] / (length(data[[1]]$causalGeneMbarIdx)*nReplicates)
        coverage[3,n] <- coverage[3,n] / ( (nGene - length(data[[1]]$causalGeneMIdx))*nReplicates)
        coverage[4,n] <- coverage[4,n] / ( (nGene - length(data[[1]]$causalGeneMbarIdx))*nReplicates)
        coverage[5,n] <- coverage[5,n] / (nCtmp*nReplicates)

        width[,n] <- width[,n] / nReplicates

        upper <- rep(NA, nPar-2) #omega and Yhat have no CI
        lower <- rep(NA, nPar-2)
        for(r in 1:nReplicates){
            for(i in 1:(nPar-2)){
                upper[i] <- tmpAllCC[[i]][r,1] + 1.96*tmpAll_sdCC[[i]][r,1] #only take the 1st one, might waste some but should still be unbiased. Furthermore, we aren't interested in coverage for now.
                lower[i] <- tmpAllCC[[i]][r,1] - 1.96*tmpAll_sdCC[[i]][r,1] 

                coverageCC[i,n] <- coverageCC[i,n] + sum( upper[i] >= truthAll[i] & lower[i] <= truthAll[i])
                widthCC[i,n] <- widthCC[i,n] + upper[i] - lower[i]

                                upper[i] <- tmpAllFull[[i]][r,1] + 1.96*tmpAll_sdFull[[i]][r,1] #only take the 1st one, might waste some but should still be unbiased. Furthermore, we aren't interested in coverage for now.
                lower[i] <- tmpAllFull[[i]][r,1] - 1.96*tmpAll_sdFull[[i]][r,1] 

                coverageFull[i,n] <- coverageFull[i,n] + sum( upper[i] >= truthAll[i] & lower[i] <= truthAll[i])
                widthFull[i,n] <- widthFull[i,n] + upper[i] - lower[i]

            }
        }

        coverageCC[1,n] <- coverageCC[1,n] / (length(dataCC[[1]]$causalGeneMIdx)*nReplicates)
        coverageCC[2,n] <- coverageCC[2,n] / (length(dataCC[[1]]$causalGeneMbarIdx)*nReplicates)
        coverageCC[3,n] <- coverageCC[3,n] / ( (nGene - length(dataCC[[1]]$causalGeneMIdx))*nReplicates)
        coverageCC[4,n] <- coverageCC[4,n] / ( (nGene - length(dataCC[[1]]$causalGeneMbarIdx))*nReplicates)
        coverageCC[5,n] <- coverageCC[5,n] / (nCtmp*nReplicates)
        
        widthCC[,n] <- widthCC[,n] / nReplicates

                coverageFull[1,n] <- coverageFull[1,n] / (length(dataFull[[1]]$causalGeneMIdx)*nReplicates)
        coverageFull[2,n] <- coverageFull[2,n] / (length(dataFull[[1]]$causalGeneMbarIdx)*nReplicates)
        coverageFull[3,n] <- coverageFull[3,n] / ( (nGene - length(dataFull[[1]]$causalGeneMIdx))*nReplicates)
        coverageFull[4,n] <- coverageFull[4,n] / ( (nGene - length(dataFull[[1]]$causalGeneMbarIdx))*nReplicates)
        coverageFull[5,n] <- coverageFull[5,n] / (nCtmp*nReplicates)
        
        widthFull[,n] <- widthFull[,n] / nReplicates

        coverage <- signif(coverage, 3)
        coverageCC <- signif(coverageCC, 3)
        coverageFull <- signif(coverageFull, 3)

        width <- signif(width, 3)
        widthCC <- signif(widthCC, 3)
        widthFull <- signif(widthFull, 3)


################################

                                        #Now we are no longer using multiple imputation to calculate all sd
                                        #sdAll[1,n] <- sqrt( (1+1/nReplicates)*(sd(as.numeric(gammaMAll[n,,]))^2) + mean(as.numeric(gammaMAll_sd[n,,])^2))
        
        tmpNPar <- c(nCausalGene, nCausalGene, nGene-nCausalGene, nGene-nCausalGene, nCtmp)

        
        
        for (i in 1:nPar){
            if(i==nPar){
                mse[i,n,1] <- signif( mean(yHatAll[n,,1]), 3)
                mseCC[i,n,1] <- signif( mean(yHatAll_CC[n,,1]), 3)
                mseFull[i,n,1] <- signif( mean(yHatAll_Full[n,,1]), 3)
                mse[i,n,2] <- signif( sd(yHatAll[n,,1])/sqrt(nReplicates), 3)
                mseCC[i,n,2] <- signif( sd(yHatAll_CC[n,,1])/sqrt(nReplicates), 3)
                mseFull[i,n,2] <- signif( sd(yHatAll_Full[n,,1])/sqrt(nReplicates), 3)
                biasAll[i,n,1] <- 0
                biasAllCC[i,n,1] <- 0
                biasAllFull[i,n,1] <- 0
                biasAll[i,n,2] <- 0
                biasAllCC[i,n,2] <- 0
                biasAllFull[i,n,2] <- 0
            }else{
                meanAll[i,n] <- mean(as.numeric(tmpAll[[i]]))
                biasAll[i,n,2] <- signif(sd(as.numeric(tmpAll[[i]])),3)            

                meanAllCC[i,n] <- mean(as.numeric(tmpAllCC[[i]]))
                biasAllCC[i,n,2] <- signif(sd(as.numeric(tmpAllCC[[i]])),3)            

                meanAllFull[i,n] <- mean(as.numeric(tmpAllFull[[i]]))
                biasAllFull[i,n,2] <- signif(sd(as.numeric(tmpAllFull[[i]])),3)            


                mse[i,n,1] <- signif(  mean((as.numeric(tmpAll[[i]])-truthAll[i])^2 + sd(as.numeric(tmpAll[[i]]))^2) , 3)
                mseCC[i,n,1] <- signif(  mean((as.numeric(tmpAllCC[[i]])-truthAll[i])^2 + sd(as.numeric(tmpAllCC[[i]]))^2) , 3)
                mseFull[i,n,1] <- signif(  mean((as.numeric(tmpAllFull[[i]])-truthAll[i])^2 + sd(as.numeric(tmpAllFull[[i]]))^2) , 3)

                mse[i,n,2] <- signif(  sd((as.numeric(tmpAll[[i]])-truthAll[i])^2 + sd(as.numeric(tmpAll[[i]]))^2)/sqrt(nCtmp*nReplicates) , 3)
                mseCC[i,n,2] <- signif(  sd((as.numeric(tmpAllCC[[i]])-truthAll[i])^2 + sd(as.numeric(tmpAllCC[[i]]))^2)/sqrt(nCtmp*nReplicates) , 3)
                mseFull[i,n,2] <- signif(  sd((as.numeric(tmpAllFull[[i]])-truthAll[i])^2 + sd(as.numeric(tmpAllFull[[i]]))^2)/sqrt(nCtmp*nReplicates) , 3)

                biasAll[i,n,1] <- signif(meanAll[i,n] - truthAll[i] , 3)
                biasAllCC[i,n,1] <- signif(meanAllCC[i,n] - truthAll[i]  , 3)
                biasAllFull[i,n,1] <- signif(meanAllFull[i,n] - truthAll[i]  , 3)
            }
            
            tableMSE[i+1,2] <- paste(mseCC[i,n,1],"(", mseCC[i,n,2] ,")",sep="")
                                        #            if(abs(biasAll[i,n,1]) > 0.5*sdAll[i,n]){
                                        #               table[i+1,2] <- paste( "\\textbf{", biasAll[i,n,1],"}(", signif(sqrt(mse[i,n,1]),3) ,")",sep="")
                                        #          }
            tableMSE[i+1,3] <- paste(mse[i,n,1],"(", mse[i,n,2]  ,")",sep="")
                                        #         if(abs(biasAllCC[i,n,1]) > 0.5*sdAllCC[i,n]){
                                        #            table[i+1,3] <- paste( "\\textbf{", biasAllCC[i,n,1],"}(", signif(sqrt(mseCC[i,n,1]),3) ,")",sep="")
                                        #       }
            if(i<=5){
                tableCI[i+1,2] <- paste(coverageCC[i,n],"(", widthCC[i,n] ,")",sep="")
                tableCI[i+1,3] <- paste(coverage[i,n],"(", width[i,n] ,")",sep="")
            }

        }
        
        tableFS[2,2] <- paste(signif(sensSpec_CC[n,1],3),"(", signif(sensSpec_CC[n,2],3) ,")",sep="") #gammaM Sens spec
        tableFS[3,2] <- paste(signif(sensSpec_CC[n,3],3),"(", signif(sensSpec_CC[n,4],3) ,")",sep="") #gammaMbar
        tableFS[4,2] <- paste(signif(sensSpec_CC[n,5],3),"(", signif(sensSpec_CC[n,6],3) ,")",sep="") #overall

        tableFS[2,3] <- paste(signif(sensSpec[n,1],3),"(", signif(sensSpec[n,2],3) ,")",sep="") #gammaM Sens spec
        tableFS[3,3] <- paste(signif(sensSpec[n,3],3),"(", signif(sensSpec[n,4],3) ,")",sep="") #gammaMbar
        tableFS[4,3] <- paste(signif(sensSpec[n,5],3),"(", signif(sensSpec[n,6],3) ,")",sep="") #overall

        write.table(tableMSE, paste(wd, "/MSE_n", n, sep=""), sep=" & ", quote=F, col.names=F, row.names=F, eol="\\\\\n")
        write.table(tableCI, paste(wd, "/CI_n", n, sep=""), sep=" & ", quote=F, col.names=F, row.names=F, eol="\\\\\n")
        write.table(tableFS, paste(wd, "/FS_n", n, sep=""), sep=" & ", quote=F, col.names=F, row.names=F, eol="\\\\\n")


    }

                                        #  save.image(paste(wd,"/all_n.RData",sep=""))
                                        #    save(dataCC,data, estCC,est, biasAll, biasAllCC,  gammaMAll, gammaMbarAll,gammaMAll_nc, gammaMbarAll_nc, gammaCAll, omegaAll, yHatAll, IMConfusion, IMbarConfusion, sensSpec, sensSpec_CC,  sensSpec2, sensSpec2_CC, mse,mseCC, file=paste(wd, "/allData.RData", sep=""))

#############################
####### plotting ############
#############################    

    pdf(paste(wd,"/CI.pdf",sep=""))
    par(mfrow=c(2,2),mar=c(5.1,4.1,1.9,2.1))
    barEpsilon=0.04
    x<-seq(1:length(all_n))

    for(i in 1:(nPar-1)){ #yhat has no bias-sd plot
        ylimTmp <- max(biasAllFull[i,,1], biasAllCC[i,,1],  biasAll[i,,1] ) + max(biasAllCC[i,,2],biasAllFull[i,,2],biasAll[i,,2])
        ylimTmp <- ylimTmp+abs(ylimTmp)*0.2
        
        ylimTmpNeg <- min(biasAllCC[i,,1],biasAllFull[i,,1],biasAll[i,,1])-max(biasAllCC[i,,2],biasAllFull[i,,2],biasAll[i,,2])
        ylimTmpNeg<-ylimTmpNeg - abs(ylimTmp)*0.1
        
        plot(x, biasAll[i,,1], ylim=c(ylimTmpNeg,ylimTmp),xlab="N",main=parNames[i],
             ylab="Bias",xaxt="n",col=colors[1],type="l")
        axis(1,at=(1:length(all_n)),label=all_n,las=2,cex.axis=0.8)
        segments(x,biasAll[i,,1]-biasAll[i,,2], x, biasAll[i,,1]+biasAll[i,,2],col=colors[1])
        segments(x-barEpsilon, biasAll[i,,1]-biasAll[i,,2], x+barEpsilon, biasAll[i,,1]-biasAll[i,,2],col=colors[1] )
        segments(x-barEpsilon, biasAll[i,,1]+biasAll[i,,2], x+barEpsilon, biasAll[i,,1]+biasAll[i,,2],col=colors[1] )

        lines(x+barEpsilon, biasAllCC[i,,1],col=colors[3],type="l")
        segments(x+barEpsilon,biasAllCC[i,,1]-biasAllCC[i,,2], x+barEpsilon, biasAllCC[i,,1]+biasAllCC[i,,2],col=colors[3])
        segments(x, biasAllCC[i,,1]-biasAllCC[i,,2], x+2*barEpsilon, biasAllCC[i,,1]-biasAllCC[i,,2],col=colors[3] )
        segments(x, biasAllCC[i,,1]+biasAllCC[i,,2], x+2*barEpsilon, biasAllCC[i,,1]+biasAllCC[i,,2],col=colors[3] )

                lines(x+barEpsilon, biasAllFull[i,,1],col=colors[4],type="l")
        segments(x+barEpsilon,biasAllFull[i,,1]-biasAllFull[i,,2], x+barEpsilon, biasAllFull[i,,1]+biasAllFull[i,,2],col=colors[4])
        segments(x, biasAllFull[i,,1]-biasAllFull[i,,2], x+2*barEpsilon, biasAllFull[i,,1]-biasAllFull[i,,2],col=colors[4] )
        segments(x, biasAllFull[i,,1]+biasAllFull[i,,2], x+2*barEpsilon, biasAllFull[i,,1]+biasAllFull[i,,2],col=colors[4] )

        abline(0,0,col=colors[2],lty=2)
                                        #      legend(1, ylimTmp ,bty = "n", sixMethods, cex=0.8,
        if(i==1){
            legend("topright", bty = "n", methods, cex=0.8, 
                   col=colors[c(3,1,4)], lty=c(1,1,1))
        }
    }
    save(mse, mseCC,mseFull, all_n,file=paste(wd, "/mse.RData", sep=""))
    
    for(i in 1:(nPar)){ 
        ylimTmp <- max(mseFull[i,,1],mseCC[i,,1],mse[i,,1]) + max(mseCC[i,,2],mseFull[i,,2],mse[i,,2])
        ylimTmp <- ylimTmp*1.2
        
        ylimTmpNeg <- min(mseCC[i,,1],mseFull[i,,1],mse[i,,1])-max(mseCC[i,,2],mseFull[i,,2],mse[i,,2])
        ylimTmpNeg<-ylimTmpNeg*1.1
        
        plot(x, mse[i,,1], ylim=c(ylimTmpNeg,ylimTmp),xlab="N",main=parNames[i],
             ylab="MSE",xaxt="n",col=colors[1],type="l")
        axis(1,at=(1:length(all_n)),label=all_n,las=2,cex.axis=0.8)
        segments(x,mse[i,,1]-mse[i,,2], x, mse[i,,1]+mse[i,,2],col=colors[1])
        segments(x-barEpsilon, mse[i,,1]-mse[i,,2], x+barEpsilon, mse[i,,1]-mse[i,,2],col=colors[1] )
        segments(x-barEpsilon, mse[i,,1]+mse[i,,2], x+barEpsilon, mse[i,,1]+mse[i,,2],col=colors[1] )

        lines(x+barEpsilon, mseCC[i,,1],col=colors[3],type="l")
        segments(x+barEpsilon,mseCC[i,,1]-mseCC[i,,2], x+barEpsilon, mseCC[i,,1]+mseCC[i,,2],col=colors[3])
        segments(x, mseCC[i,,1]-mseCC[i,,2], x+2*barEpsilon, mseCC[i,,1]-mseCC[i,,2],col=colors[3] )
        segments(x, mseCC[i,,1]+mseCC[i,,2], x+2*barEpsilon, mseCC[i,,1]+mseCC[i,,2],col=colors[3] )

        lines(x+barEpsilon, mseFull[i,,1],col=colors[4],type="l")
        segments(x+barEpsilon,mseFull[i,,1]-mseFull[i,,2], x+barEpsilon, mseFull[i,,1]+mseFull[i,,2],col=colors[4])
        segments(x, mseFull[i,,1]-mseFull[i,,2], x+2*barEpsilon, mseFull[i,,1]-mseFull[i,,2],col=colors[4] )
        segments(x, mseFull[i,,1]+mseFull[i,,2], x+2*barEpsilon, mseFull[i,,1]+mseFull[i,,2],col=colors[4] )

        if(i==nPar){
            abline(epsilon^2,0,col=colors[2],lty=2) #mse of yhat should shrink to sigma^2
        }
        if(i==1){
            legend("topright", bty = "n", methods, cex=0.8, 
                   col=colors[c(3,1,4)], lty=c(1,1,1))
        }
    }
    
    for(i in 1:(length(sensSpecNames))){
        ylimTmp <- max(sensSpec[,i], sensSpec_Full[,i], sensSpec_CC[,i])
        ylimTmp <- ylimTmp*1.2
        ylimTmp <- min(ylimTmp, 1)
        ylimTmpNeg <- min(sensSpec[,i], sensSpec_Full[,i], sensSpec_CC[,i])
        ylimTmpNeg <- ylimTmpNeg*0.8
        ylimTmpNeg <- min(ylimTmpNeg,0.5)
        
        expression(paste(gamma^M," causal"))
        
        plot(x, sensSpec[,i], xlab="N", main=sensSpecNames[i],ylim=c(ylimTmpNeg,ylimTmp), ylab=NA,xaxt="n",col=colors[1], type="l",sub="BFDR=0.05")
        axis(1,at=(1:length(all_n)),label=all_n,las=2,cex.axis=0.8)
                                        #  legend("bottomright", "FDR=0.05", bty="n")
        
        lines(x, sensSpec_CC[,i],col=colors[3],type="l")
        lines(x, sensSpec_Full[,i],col=colors[4],type="l")
        abline(1,0,col=colors[2],lty=2)
    }

    for(i in 1:(length(sensSpecNames))){
        ylimTmp <- max(sensSpec2[,i], sensSpec2_Full[,i], sensSpec2_CC[,i])
        ylimTmp <- ylimTmp*1.2
        ylimTmp <- min(ylimTmp, 1)
        ylimTmpNeg <- min(sensSpec2[,i], sensSpec2_Full[,i], sensSpec2_CC[,i])
        ylimTmpNeg <- ylimTmpNeg*0.8
        ylimTmpNeg <- min(ylimTmpNeg,0.5)
        
        expression(paste(gamma^M," causal"))
        
        plot(x, sensSpec2[,i], xlab="N", main=sensSpecNames[i], ylim=c(ylimTmpNeg,ylimTmp), ylab=NA,xaxt="n",col=colors[1], type="l",sub="BFDR=0.1")
        axis(1,at=(1:length(all_n)),label=all_n,las=2,cex.axis=0.8)
                                        #        legend("topright", "FDR=0.1", bty="n")
        lines(x, sensSpec2_CC[,i],col=colors[3],type="l")
        lines(x, sensSpec2_Full[,i],col=colors[4],type="l")
        abline(1,0,col=colors[2],lty=2)
    }
    
    dev.off()
    
    #plot MSE(yhat) and feature selections separatedly
    pdf(paste(wd,"/MSEyHat.pdf",sep=""))
    par(mfrow=c(1,1),mar=c(5.1,4.1,1.9,2.1),pin=c(4,4))
    i=nPar
    xZoom <- list()
    xZoom[[1]] <- seq(1:length(all_n))
    if(length(all_n)==3){
            xZoom[[2]] <- 1:2
            xZoom[[3]] <- 2:3
    }else if(length(all_n)==4){
        xZoom[[2]] <- 2:4
        xZoom[[3]] <- 1:2
        xZoom[[4]] <- 2:3
        xZoom[[5]] <- 3:4
    }
    for(x in xZoom){
    ylimTmp <- max(mseFull[i,x,1],mseCC[i,x,1],mse[i,x,1]) + max(mseCC[i,x,2],mseFull[i,x,2],mse[i,x,2])
    ylimTmp <- ylimTmp*1.2
    
    ylimTmpNeg <- min(mseCC[i,x,1],mseFull[i,x,1],mse[i,x,1])-max(mseCC[i,x,2],mseFull[i,x,2],mse[i,x,2])
    ylimTmpNeg<-ylimTmpNeg*1.1
    
    plot(x, mse[i,x,1], ylim=c(ylimTmpNeg,ylimTmp),xlab="N",main=parNames[i],
         ylab="MSE",xaxt="n",col=colors[1],type="l")
    axis(1,at=x,label=all_n[x],las=2,cex.axis=0.8)
    segments(x,mse[i,x,1]-mse[i,x,2], x, mse[i,x,1]+mse[i,x,2],col=colors[1])
    segments(x-barEpsilon, mse[i,x,1]-mse[i,x,2], x+barEpsilon, mse[i,x,1]-mse[i,x,2],col=colors[1] )
    segments(x-barEpsilon, mse[i,x,1]+mse[i,x,2], x+barEpsilon, mse[i,x,1]+mse[i,x,2],col=colors[1] )
    
    lines(x+barEpsilon, mseCC[i,x,1],col=colors[3],type="l")
    segments(x+barEpsilon,mseCC[i,x,1]-mseCC[i,x,2], x+barEpsilon, mseCC[i,x,1]+mseCC[i,x,2],col=colors[3])
    segments(x, mseCC[i,x,1]-mseCC[i,x,2], x+2*barEpsilon, mseCC[i,x,1]-mseCC[i,x,2],col=colors[3] )
    segments(x, mseCC[i,x,1]+mseCC[i,x,2], x+2*barEpsilon, mseCC[i,x,1]+mseCC[i,x,2],col=colors[3] )

        lines(x+barEpsilon, mseFull[i,x,1],col=colors[4],type="l")
    segments(x+barEpsilon,mseFull[i,x,1]-mseFull[i,x,2], x+barEpsilon, mseFull[i,x,1]+mseFull[i,x,2],col=colors[4])
    segments(x, mseFull[i,x,1]-mseFull[i,x,2], x+2*barEpsilon, mseFull[i,x,1]-mseFull[i,x,2],col=colors[4] )
    segments(x, mseFull[i,x,1]+mseFull[i,x,2], x+2*barEpsilon, mseFull[i,x,1]+mseFull[i,x,2],col=colors[4] )

    abline(epsilon^2,0,col=colors[2],lty=2) #mse of yhat should shrink to sigma^2
    legend("topright", bty = "n", methods, cex=0.8, 
           col=colors[c(3,1,4)], lty=c(1,1,1))
    }    
    
    dev.off()
    
    
    pdf(paste(wd,"/FS.pdf",sep=""))
    par(mfrow=c(2,2),mar=c(5.1,4.1,1.9,2.1))
    tmpLength = length(sensSpecNames)

#    x<-seq(1:length(all_n))
    for(x in xZoom){
    
    for(i in (tmpLength-1):tmpLength){
      ylimTmp <- max(sensSpec[x,i], sensSpec_Full[x,i], sensSpec_CC[x,i] )
      ylimTmp <- ylimTmp*1.2
      ylimTmp <- min(ylimTmp, 1)
      ylimTmpNeg <- min(sensSpec[x,i], sensSpec_Full[x,i], sensSpec_CC[x,i])
      ylimTmpNeg <- ylimTmpNeg*0.8
      ylimTmpNeg <- min(ylimTmpNeg,0.5)
      
      plot(x, sensSpec[x,i], xlab="N", main=sensSpecNames[i],ylim=c(ylimTmpNeg,ylimTmp), ylab=NA,xaxt="n",col=colors[1], type="l",sub="BFDR=0.05")
      axis(1,at=x,label=all_n[x],las=2,cex.axis=0.8)
      
      lines(x, sensSpec_CC[x,i],col=colors[3],type="l")
      lines(x, sensSpec_Full[x,i],col=colors[4],type="l")
      if(i==(tmpLength-1)){
        legend("bottomright", bty = "n", methods, cex=0.8, 
               col=colors[c(3,1,4)], lty=c(1,1,1))
      }
      
      abline(1,0,col=colors[2],lty=2)
    }
    
    for(i in (tmpLength-1):tmpLength){
      ylimTmp <- max(sensSpec2[x,i], sensSpec2_Full[x,i], sensSpec2_CC[x,i])
      ylimTmp <- ylimTmp*1.2
      ylimTmp <- min(ylimTmp, 1)
      ylimTmpNeg <- min(sensSpec2[x,i], sensSpec2_Full[x,i], sensSpec2_CC[x,i])
      ylimTmpNeg <- ylimTmpNeg*0.8
      ylimTmpNeg <- min(ylimTmpNeg,0.5)
      
      plot(x, sensSpec2[x,i], xlab="N", main=sensSpecNames[i], ylim=c(ylimTmpNeg,ylimTmp), ylab=NA,xaxt="n",col=colors[1], type="l",sub="BFDR=0.1")
      axis(1,at=x,label=all_n[x],las=2,cex.axis=0.8)
      #        legend("topright", "FDR=0.1", bty="n")
      lines(x, sensSpec2_Full[x,i],col=colors[4],type="l")
            lines(x, sensSpec2_CC[x,i],col=colors[3],type="l")

      abline(1,0,col=colors[2],lty=2)
    }
}   
    dev.off()
        
    print("Done for varying nSample")
    
    
}



