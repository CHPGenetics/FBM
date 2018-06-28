args = commandArgs(trailingOnly=TRUE)

require(MCMCpack)
library(doParallel)
library(foreach)
library(cluster)
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(scales)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("~/fbm/FBM.cpp") #put in the path of source code FBM.cpp here

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
    
    gammaTrue <- c(0, 5, 5, 5) # beta0, betaGeneM, betaGeneMbar, betaC(any clinical)
                                        #  gammaTrue <- c(0, 2, 2/(epsilon), 2)

    effectiveMethy <- rep(0, nMethy) #Adding sparsity, 0 means sparse, 1 means methy is effetive on gene.
    for (k in 1:nGene){
        effectiveMethy[which(mapMethy==k)][1]=1 #make sure at least one methy is effective per gene.
    }
                                        #    effectiveMethy[which(effectiveMethy==0)] <- rbinom(nMethy-nGene, 1, 1) #0% sparsity on the rest of the methy.
                                        #effectiveMethy <- rep(0, nMethy)
    
    effectiveMethy <- rep(1, nMethy)     #fix: if gene not causal, omega=0;
                                        #   for(k in 1:nGene){
                                        #       if(k %in%causalGeneMIdx){
                                        #           effectiveMethy[which(mapMethy==k)] <- 1
                                        #       }
                                        #      if(k %in%causalGeneMbarIdx){
                                        #         effectiveMethy[which(mapMethy==k)] <- 1
                                        #    }
                                        #}

    if(exchangeable==TRUE){
        for (k in 1:nGene){
            gene[,k] <- rnorm(nSample, nuGene[k], etaGene[k])
        }
    } else {
                                        #now generate gene with covariance strucuture
        for(k in 1:nGene){
            geneMTrue[,k] <- omegaTrue[k] * apply(as.matrix(methy[,which(mapMethy==k & effectiveMethy==1)]), 1, sum)
        }
        r <- 0.2 #0.2 or 0.8
        covMat <- matrix(epsilonGeneMbar^2, nGene, nGene)
        tmpFill <- epsilonGeneMbar^2
        for(k1 in 1:nGene){
            for(k2 in 1:nGene){
                covMat[k1, k2] <- tmpFill*r^(abs(k1-k2))
            }
        }

        for(n in 1:nSample){
            gene[n,] <- mvrnorm(1, geneMTrue[n,], Sigma=covMat)
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
    ## colnames(methy) <- paste("methy", seq(1:nMethy), ".",seq(1:nMethy), ".gene", mapMethy, ".",seq(1:nMethy) , sep="" )
    ## rownames(methy) <- paste("sample",seq(1:nSample),sep="")
    ## colnames(geneObs) <- paste("gene",1:nGene, sep="")
    ## rownames(geneObs) <- paste("sample",obsGeneIdx,sep="")

    ## colnames(C) <- paste("clinical",1:nC, sep="")
    ## rownames(C) <- paste("sample",seq(1:nSample),sep="")
    
    if(is.null(colnames(methy))){
        colnames(methy) <- paste("methy", seq(1:nMethy), ".",seq(1:nMethy), ".gene", mapMethy, ".",seq(1:nMethy) , sep="" )
    }
    if(is.null(rownames(methy))){
        rownames(methy) <- paste("sample",seq(1:nSample),sep="")
    }
    if(is.null(colnames(geneObs))){
        colnames(geneObs) <- paste("gene",1:nGene, sep="")
    }
    if(is.null(rownames(geneObs))){
        rownames(geneObs) <- paste("sample",obsGeneIdx,sep="")
    }
    if(is.null(colnames(C))){
        colnames(C) <- paste("clinical",1:nC, sep="")
    }
    if(is.null(rownames(C))){
        rownames(C) <- paste("sample",seq(1:nSample),sep="")
    }
    
                                        #write.csv(methy, paste(wd, "/methyData.csv", sep=""), row.names=F,quote=F )
                                        #write.csv(geneObs, paste(wd, "/geneData.csv",sep=""),row.names=F, quote=F )
                                        #write.csv(C, paste(wd, "/clinicalData.csv", sep=""), row.names=F, quote=F )
    write.csv(methy, paste(wd, "/methyData.csv", sep=""), quote=F ) # in future write out methyObs
    write.csv(geneObs, paste(wd, "/geneData.csv",sep=""), quote=F )
    write.csv(C, paste(wd, "/clinicalData.csv", sep=""), quote=F )
    print("Done simulating data")
    
    return(list(Y=Y, geneObs=geneObs, methyObs=methyObs, mapMethy=mapMethy, C=C, missGeneIdx=missGeneIdx, obsGeneIdx=obsGeneIdx, missMethyIdx=missMethyIdx, obsMethyIdx=obsMethyIdx, causalGeneMIdx=causalGeneMIdx, causalGeneMbarIdx=causalGeneMbarIdx, geneTrue=gene, geneMTrue=geneMTrue, methyTrue=methy, gammaTrue=gammaTrue, omegaTrue=omegaTrue, effectiveMethy=effectiveMethy))

}

#######################################################
### Bayesian Imputation (BI) with gene/methy missing ##
#######################################################
                                        #now using Rcpp
plotting <- function(data, result, wd){
    opar <- par()
    
    pdf(paste(wd, "/convergence.pdf",sep=""))
    par(mfrow=c(3,2))
    
    for (k in 1:(min(nGene,50))){
        for(j in which(data$mapMethy==k)){
            plot(result$omega[1:nItr,j],type="l", ylab=paste("omega",j,k),xlab="Iterations")
            abline(h=data$omegaTrue[1]*data$effectiveMethy[j], col=3)
        }
    }
    
    for (k in c(1:(min(nGene,30)),data$causalGeneMIdx) ){
        plot(result$gammaM[1:nItr,k],type="l", ylab=paste("gammaM",k),xlab="Iterations")
        if(k %in% data$causalGeneMIdx ){
            abline(h=data$gammaTrue[2], col=3)  
        }else{
            abline(h=0, col=3)  
        }
    }
    
    for (k in  c(1:(min(nGene,30)),data$causalGeneMbarIdx) ){
        plot(result$gammaMbar[1:nItr,k],type="l", ylab=paste("gammaMbar",k),xlab="Iterations")
        if(k %in% data$causalGeneMbarIdx ){
            abline(h=data$gammaTrue[3], col=3)  
        }else{
            abline(h=0, col=3)  
        }
    }
    if(nC!=0){
        for (k in 1:nC){
            plot(result$gammaC[1:nItr,k],type="l", ylab=paste("gammaC",k),xlab="Iterations")
            abline(h=data$gammaTrue[4], col=3)  
        }
    }
    
    plot(result$sigma[1:nItr],type="l", ylab="sigma",xlab="Iterations")
    abline(h=epsilon, col=3)
    
    for (k in 1:(min(nGene,50))){
        plot(result$sigmak[1:nItr,k],type="l", ylab=paste("sigmak",k),xlab="Iterations")
    }
    
    plot(result$tauM[1:nItr],type="l", ylab=paste("tauM",k),xlab="Iterations")
    abline(h=1, col=3)
    
    plot(result$tauMbar[1:nItr],type="l", ylab=paste("tauMbar",k),xlab="Iterations")
    abline(h=1, col=3)

    for (j in 1:(min(nMethy,100))){
        plot(result$IOmega[1:nItr,j],type="l", ylab=paste("IOmega ",j,xlab="Iterations"))
        abline(h=1, col=3)
    }

    for (k in 1:(min(nGene,50))){
        plot(result$IM[1:nItr,k],type="l", ylab=paste("IM ",k,xlab="Iterations"))
        if(k %in% data$causalGeneMIdx){
            abline(h=1, col=3)
        }else{
            abline(h=0, col=3)
        }
    }
    for (k in 1:(min(nGene,50))){
        plot(result$IMbar[1:nItr,k],type="l", ylab=paste("IMbar ",k,xlab="Iterations"))
        if(k %in% data$causalGeneMbarIdx){
            abline(h=1, col=3)
        }else{
            abline(h=0, col=3)
        }
    }
    dev.off()
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
        qtmp[k] <- min( BFDRtmp[ which(BFDRtmp[,1] > PMtmp[k]), 2],1 ) #fix: >= or > ?
    }

    return(list(q=qtmp))
}

estimating <- function(data, result,wd, testData){
    nItrActual <- result$nItrActual
    nBurninActual <- result$nBurninActual
                                        #  nBurninActual <- round(nItr/3)
    nGene <- dim(data$geneObs)[2]    
    nC <- dim(data$C)[2]
    gammaMTrue <- rep(0, nGene)
    gammaMTrue[data$causalGeneMIdx] <- data$gammaTrue[2]
    gammaMbarTrue <- rep(0, nGene)
    gammaMbarTrue[data$causalGeneMbarIdx] <- data$gammaTrue[3]
    gammaCTrue <- rep(data$gammaTrue[4], nC)
    
    IMEst <- rep(0,nGene)
    IMbarEst <- rep(0,nGene)
    gammaMEst <- rep(0,nGene)
    gammaMbarEst <- rep(0,nGene)
    gammaMSd <- rep(0,nGene)
    gammaMbarSd <- rep(0,nGene)
    gammaCEst <- rep(0,nC)
    gammaCSd <- rep(0,nC) 
    
    for(k in 1:nGene){
        tmpData <- as.numeric(result$gammaM[nBurninActual:nItrActual,k])
        gammaMEst[k] <- signif(median(tmpData),4)
        gammaMSd[k] <- signif(sd(tmpData),4)
        IMEst[k] <- mean(result$IM[nBurninActual:nItrActual,k])

        tmpData <- as.numeric(result$gammaMbar[nBurninActual:nItrActual,k])
        gammaMbarEst[k] <- signif(median(tmpData),4)
        gammaMbarSd[k] <- signif(sd(tmpData),4)
        IMbarEst[k] <- mean(result$IMbar[nBurninActual:nItrActual,k])

    }
    
    for(k in 1:nC){
        tmpData <- as.numeric(result$gammaC[nBurninActual:nItrActual,k])
        gammaCEst[k] <- signif(median(tmpData),4)
        gammaCSd[k] <- signif(sd(tmpData),4)
    }
    
    
    table <- matrix(NA, 1+2*nGene+nC, 6)
    table[1,] <- c("\ ","geneName","Estimate", "Truth", "Selected","q")
    table[2:(1+2*nGene+nC),1] <- rep("\ ", 2*nGene+nC )
    table[2:(1+2*nGene+nC),2] <- c(colnames(data$geneObs),colnames(data$geneObs), colnames(data$C))
    table[2:(1+2*nGene+nC),3] <- c(gammaMEst, gammaMbarEst, gammaCEst)
    table[2:(1+2*nGene+nC),4] <- c(gammaMTrue, gammaMbarTrue, gammaCTrue)
    table[2:(1+2*nGene+nC),5] <- c(IMEst, IMbarEst, rep("\ ", nC))

    qM <- BFDR(IMEst, dim(data$geneObs)[1], dim(data$geneObs)[2] )$q
    qMbar <- BFDR(IMbarEst, dim(data$geneObs)[1], dim(data$geneObs)[2] )$q

    table[2:(1+2*nGene+nC),6] <- c(qM, qMbar, rep("\ ", nC))

                                        #we have only nMethy meaningful omega
                                        # without loss of generalizability, we only use the first omega
    sigmaEst <- signif(median(result$sigma[nBurninActual:nItrActual]),4)
    tauMEst <- signif(median(result$tauM[nBurninActual:nItrActual]),4)
    tauMbarEst <- signif(median(result$tauMbar[nBurninActual:nItrActual]),4)
    sigmakEst <- signif(apply(result$sigmak[nBurninActual:nItrActual,],2,median),4)
    
    yHatEstTmp <- as.matrix(testData$geneMTrue[,testData$causalGeneMIdx])%*%gammaMEst[data$causalGeneMIdx] + as.matrix(testData$geneTrue[,testData$causalGeneMbarIdx]-testData$geneMTrue[,testData$causalGeneMbarIdx])%*%gammaMbarEst[data$causalGeneMbarIdx] + as.matrix(testData$C)%*%gammaCEst
    + as.matrix(testData$geneMTrue[,-testData$causalGeneMIdx])%*%gammaMEst[-data$causalGeneMIdx] + as.matrix(testData$geneTrue[,-testData$causalGeneMbarIdx]-testData$geneMTrue[,-testData$causalGeneMbarIdx])%*%gammaMbarEst[-data$causalGeneMbarIdx]

    MSEyHatEst <- yHatEstTmp-testData$Y
    dftmp <- ( length(testData$Y) - 2*nCausalGene - nC)
    if(dftmp <= 1){
        dftmp=1
    }
    MSEyHatEst <- sum(MSEyHatEst^2)/ dftmp
    
    omegaEst <- signif(apply(result$omega[nBurninActual:nItrActual,],2,median),4)
    
    omegaMeanEst <- median(omegaEst)
    omegaSd <- sd(omegaEst)
    
    write.table(table, paste(wd, "/estimatesMain",sep=""),
                sep="&", quote=F, col.names=F, row.names=F, eol="\\\\\n")
    
    cat("mean(omega)\n",file=paste(wd, "/estimatesOther", sep=""),append=T)
    cat(omegaMeanEst, file=paste(wd, "/estimatesOther", sep=""),append=T, sep=",")
    cat("\nsigmaEst,tauMEst,tauMbarEst\n",file=paste(wd, "/estimatesOther", sep=""),append=T)
    cat(c(sigmaEst,tauMEst,tauMbarEst),file=paste(wd, "/estimatesOther", sep=""),append=T, sep=",")
    cat("\nsigmakEst\n",file=paste(wd, "/estimatesOther", sep=""),append=T)
    cat(sigmakEst, file=paste(wd, "/estimatesOther", sep=""),append=T, sep=",")
    
    return(list(gammaMEst=gammaMEst, gammaMbarEst=gammaMbarEst, gammaCEst=gammaCEst, omegaEst=omegaEst,omegaMeanEst=omegaMeanEst, IMEst=IMEst,IMbarEst=IMbarEst, sigmaEst=sigmaEst, tauMEst=tauMEst, tauMbarEst=tauMbarEst, sigmakEst=sigmakEst, MSEyHatEst = MSEyHatEst, gammaMSd=gammaMSd, gammaMbarSd=gammaMbarSd, gammaCSd=gammaCSd, omegaSd=omegaSd, qM=qM, qMbar=qMbar))
}


nGene <- as.numeric(args[8])
nCausalGene <-  as.numeric(args[9])
nMethy <-  as.numeric(args[10])
nC <- 2
exchangeable <- FALSE
MCAR <- TRUE

if( fullRun | multipleReplicates ){
    wd <- args[1]
dir.create(wd)
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

    parNames <- c( expression(paste(gamma^M," causal")),  expression(paste(gamma^bar(M), " causal")),  expression(paste(gamma^M, " non-causal")),  expression(paste(gamma^bar(M), " non-causal")),  expression(gamma^C),  expression(omega), expression(hat(Y)))
    nPar <- length(parNames) #  gammaM_causal, gammaMbar_causal,gammaM_nc, gammaMbar_nc, gammaC, omegaMean, yHat

    sensSpecNames <- c( expression(paste(I^M,"sens")),  expression(paste(I^M,"spec")),  expression(paste(I^bar(M),"sens")),expression(paste(I^bar(M),"spec")), expression("Overall sens"),  expression("Overall spec"))

    methods <- c("CC", "BI", "Full")
    
    nReplicates <-  as.numeric(args[11])
    r <- as.numeric(args[13])

    colors<- c(colors()[26],  colors()[261], colors()[35], colors()[614],  colors()[621]) #blue, grey, red, green, yellow.
    
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
        data[[r]] <- simulateData(all_n[n], nGene, nMethy,nC, nCausalGene, exchangeable, MCAR, epsilon,epsilonGeneMbar, thetaGene, thetaMethy, seed+r)
        testData[[r]] <- simulateData(1000, nGene, nMethy,nC, nCausalGene, exchangeable, MCAR, epsilon,epsilonGeneMbar, 0, 0, seed-r)
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


                                        #FULL
        if(thetaGene==0.5 & thetaMethy==0){
            set.seed(seed)
          print(paste("Full",Sys.time())) #test
          biFull[[r]] <- FBMcpp(dataFull[[r]], nItr, seed, "Full")
          print(paste("Full",Sys.time())) #test
          wdSub <- paste(wd,"/n",n,"_r",r,"_Full",sep="")
            dir.create(wdSub)
            estFull[[r]] <- estimating(dataFull[[r]],biFull[[r]],wdSub, testData[[r]])
            plotting(dataFull[[r]], biFull[[r]], wdSub)
            save(estFull, file=paste(wdSub,"/run.RData",sep=""))
        }  #fix: no need to run BI full data for more than one settings.
        
                                        #BI
        set.seed(seed)
        print(paste("FBM",Sys.time())) #test
        bi[[r]] <- FBMcpp(data[[r]], nItr, seed, "BI")
        print(paste("FBM",Sys.time())) #test
        wdSub <- paste(wd,"/n",n,"_r",r,sep="")
        dir.create(wdSub)
        est[[r]] <- estimating(data[[r]],bi[[r]],wdSub,testData[[r]])
        plotting(data[[r]], bi[[r]], wdSub)
        save(est, file=paste(wdSub,"/run.RData",sep=""))

                                        #CC

        set.seed(seed)
        biCC[[r]] <- FBMcpp(dataCC[[r]], nItr, seed, "CC")
        wdSub <- paste(wd,"/n",n,"_r",r,"_CC",sep="")
        dir.create(wdSub)
        estCC[[r]] <- estimating(dataCC[[r]],biCC[[r]],wdSub,testData[[r]])
        plotting(dataCC[[r]], biCC[[r]], wdSub)
        save(estCC, file=paste(wdSub,"/run.RData",sep=""))


    }
}

