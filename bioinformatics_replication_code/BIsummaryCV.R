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


    methods <- c("CC", "BI", "BIgene")
    
    nReplicates <-  as.numeric(args[11]) #50 replicates (MCMC chains)

                                        #    colors <- c(colors()[461],  colors()[282], colors()[555], colors()[610])
    colors<- c(colors()[26],  colors()[261], colors()[35], colors()[614],  colors()[621]) #blue, grey, red, green, yellow.
    
    
    tmp <- rep(NA, nReplicates*1*length(all_n)) 
    yHatAll <- array(tmp, c(length(all_n),nReplicates, 1))  #different from all the other, stores MSE(yHat-yTrue) for each replicate
    yHatAll_CC <- array(tmp, c(length(all_n),nReplicates, 1)) 
    yHatAll_Gene <- array(tmp, c(length(all_n),nReplicates, 1)) 

    mse <- array(tmp, c(1, length(all_n), 2))#used to be mse <- matrix(NA,nPar, length(all_n)), now stores mean(MSE) and sd(MSE)
    mseCC <- array(tmp, c(1, length(all_n), 2))
    mseGene <- array(tmp, c(1, length(all_n), 2))

    
    bi <- list()
    biCC<- list()
    biGene <- list()
    est <- list()
    estCC <- list()
    estGene <- list()
    data <- list()
    dataCC <- list()
    dataGene <- list()
    testData <- list() # same test data for BI, CC and FUll

    decisionCV <- matrix(NA,4, nReplicates) 

    decisionTrue <- matrix(NA, 4, nReplicates)

    decision3by3 <- list()
    decision3by3adjust <- list()
    
    
    for (n in 1:length(all_n)){
        for (r in 1:nReplicates){
            wdSub <- paste(wd,"/n",n,"_r",r,sep="")
            load(paste(wdSub,"/run.RData",sep=""))
            yHatAll[n,r,] <- est[[r]]$MSEyHatEst #MSE(yHat-yTrue)
        }

        for (r in 1:nReplicates){
            wdSub <- paste(wd,"/n",n,"_r",r,"_CC",sep="")
            load(paste(wdSub,"/run.RData",sep=""))
            yHatAll_CC[n,r,] <- estCC[[r]]$MSEyHatEst 
        }

        if(thetaGene !=0 & thetaMethy !=0){
            for (r in 1:nReplicates){
                wdSub <- paste(wd,"/n",n,"_r",r,"_Gene",sep="")
                load(paste(wdSub,"/run.RData",sep=""))
                yHatAll_Gene[n,r,] <- estGene[[r]]$MSEyHatEst
            }
        } else {
            for (r in 1:nReplicates){
                yHatAll_Gene[n,r,] <- 0
            }
        }

        for (r in 1:nReplicates){
            wdSub <- paste(wd,"/n",n,"_r",r,"_CV",sep="")
            load(paste(wdSub,"/run.RData",sep=""))
            decisionCV[n,r] <- which(decisions == max(decisions))[1]
        }

        save(estCC,estGene,est,file=paste(wd, "/data_n",n,".RData", sep=""))
    }
                                        #fix, do we use a random gamma, or we average them?

    badRound <- integer(0) #unconverged rounds by Geweke diagnostics
    nReplicatesActual <- nReplicates-length(badRound)
#    nReplicatesActual <- nReplicates
    
    decisionErr <- rep(0, 4) #for 50 100 200 500
    decisionErrAdjust <- rep(0,4)
    tolerance <- 0.01
    toleranceBase <- 200
    newStat <- matrix(0,4,nReplicatesActual)


    for (n in 1:length(all_n)){                
        i=1
        ## mse[i,n,1] <- signif( mean(yHatAll[n,-badRound,1]), 3)
        ## mseCC[i,n,1] <- signif( mean(yHatAll_CC[n,-badRound,1]), 3)
        ## mseGene[i,n,1] <- signif( mean(yHatAll_Gene[n,-badRound,1]), 3)
        ## mse[i,n,2] <- signif( sd(yHatAll[n,-badRound,1])/sqrt(nReplicates), 3)
        ## mseCC[i,n,2] <- signif( sd(yHatAll_CC[n,-badRound,1])/sqrt(nReplicates), 3)
        ## mseGene[i,n,2] <- signif( sd(yHatAll_Gene[n,-badRound,1])/sqrt(nReplicates), 3)

        mse[i,n,1] <- signif( mean(yHatAll[n,,1]), 3)
        mseCC[i,n,1] <- signif( mean(yHatAll_CC[n,,1]), 3)
        mseGene[i,n,1] <- signif( mean(yHatAll_Gene[n,,1]), 3)
        mse[i,n,2] <- signif( sd(yHatAll[n,,1])/sqrt(nReplicates), 3)
        mseCC[i,n,2] <- signif( sd(yHatAll_CC[n,,1])/sqrt(nReplicates), 3)
        mseGene[i,n,2] <- signif( sd(yHatAll_Gene[n,,1])/sqrt(nReplicates), 3)
        # we don't actually need to use this average mse

        if(thetaGene==0 | thetaMethy==0){
            mseGene[i,n,1] = 2* max(mseCC[i,n,1], mse[i,n,1])
        }

        decision3by3[[n]] <- matrix(0, 3, 3)
        decision3by3adjust[[n]] <- matrix(0, 3, 3)
        for(r in 1:nReplicates){
            if(thetaGene==0 | thetaMethy==0)
                threeMSE <- c(yHatAll_CC[n,r,1], yHatAll[n,r,1])
            else
                threeMSE <- c(yHatAll_CC[n,r,1], yHatAll[n,r,1], yHatAll_Gene[n,r,1])

            decisionTrue[n,r] <- which(threeMSE==min(threeMSE))
        }


        for(r in 1:nReplicates){
            decision3by3[[n]][decisionCV[n,r], decisionTrue[n,r]] <-  decision3by3[[n]][decisionCV[n,r], decisionTrue[n,r]] + 1
            decision3by3adjust[[n]][decisionCV[n,r], decisionTrue[n,r]] <-  decision3by3adjust[[n]][decisionCV[n,r], decisionTrue[n,r]] + 1
            
            if( decisionCV[n,r] != decisionTrue[n,r]){
                newStat[n,r] <- ( sqrt(threeMSE[decisionCV[n,r]]) - sqrt(threeMSE[decisionTrue[n,r]]) ) / ( sqrt(threeMSE[decisionTrue[n,r]]) + toleranceBase )
                if(newStat[n,r] < tolerance)
                    decision3by3adjust[[n]][decisionCV[n,r], decisionTrue[n,r]] <-  decision3by3adjust[[n]][decisionCV[n,r], decisionTrue[n,r]] - 1                                 
            }
        }
    }

    for (n in 1:length(all_n)){
        for(i in 1:3){
            decisionErr[n] <- decisionErr[n] + decision3by3[[n]][i,i]
            decisionErrAdjust[n] <- decisionErrAdjust[n] + decision3by3adjust[[n]][i,i]
        }
        decisionErr[n]  = 1 - decisionErr[n] / sum(decision3by3[[n]])
        decisionErrAdjust[n]  = 1 - decisionErrAdjust[n] / sum(decision3by3adjust[[n]])
    }

    #print(decisionTrue)
    print(decision3by3)
    print(decision3by3adjust)
    print(decisionErr)
    print(decisionErrAdjust)

    save(mse, mseCC,mseGene,  decisionTrue, decisionCV,decisionErr, decisionErrAdjust, newStat, yHatAll, yHatAll_CC, yHatAll_Gene, decision3by3, decision3by3adjust, newStat, file=paste(wd, "/mse.RData", sep=""))     
    print("Done for varying nSample")

    for(n in 1:4){
        write.table(decision3by3[[n]], paste(wd,"/decision",n, sep=""), quote=F, sep=" & ")
        write.table(decision3by3adjust[[n]], paste(wd,"/decisionAdjust",n, sep=""), quote=F, sep=" & ")
    }
    write.table(decisionErrAdjust, paste(wd,"/adjusterr", sep=""), quote=F, sep=" & ")
    write.table(decisionErr, paste(wd,"/err", sep=""), quote=F, sep=" & ")

}



