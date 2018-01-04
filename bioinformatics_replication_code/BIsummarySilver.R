args = commandArgs(trailingOnly=TRUE)

require(MCMCpack)
library(doParallel)
library(foreach)
library(cluster)
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(scales)
library(coda)
#Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
#sourceCpp("~/BI/BI.cpp")
qcutAll <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)

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


#######################################################
### Bayesian Imputation (BI) with gene/methy missing ##
#######################################################
#now using Rcpp


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




nGene <- as.numeric(args[8])
nCausalM <-  as.numeric(args[9])
nCausalMbar <-  as.numeric(args[10])
nMethy <-  as.numeric(args[11])
nC <- 2
exchangeable <- FALSE
MCAR <- TRUE

nGene<-1000
nMethy <- 2619
#nCausalM <- 10 #fix
#nCausalMbar <- 10 

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
  all_n <- c(460)
  
  if(args[13] == "strong"){
    epsilon <- 0.01
  } else if( args[13]=="mid"){
    epsilon <- 0.5
  } else if( args[13]=="weak"){
    epsilon <- 2
  }else if( args[13]=="midstrong"){
    epsilon <- 0.1
  }else if( args[13]=="superweak"){
    epsilon <- 3
  }else if( args[13]=="ultraweak"){
    epsilon <- 5
  }
  
  epsilonGeneMbar <- 2
  
  if( multipleReplicates ){
    all_n <- 100
  }
  load("~/BI/nasal/nasalProcessed.rdata")
  geneNames <- colnames(gene)
  
  load(paste(wd, "/dataSilver.rdata",sep=""))
  nGene= dim(dataSilver$geneMTrue)[2]
  nMethy = dim(dataSilver$methyTrue)[2]
  nCausalM <- length(dataSilver$causalGeneMIdx)
  nCausalMbar <-  length(dataSilver$causalGeneMbarIdx)
  
  parNames <- c( expression(paste(gamma^M," causal")),  expression(paste(gamma^bar(M), " causal")),  expression(paste(gamma^M, " non-causal")),  expression(paste(gamma^bar(M), " non-causal")),  expression(gamma^C),  expression(omega), expression(hat(Y)))
  nPar <- length(parNames) #  gammaM_causal, gammaMbar_causal,gammaM_nc, gammaMbar_nc, gammaC, omegaMean, yHat
  
  sensSpecNames <- c( expression(paste(I^M,"sens")),  expression(paste(I^M,"spec")),  expression(paste(I^bar(M),"sens")),expression(paste(I^bar(M),"spec")), expression("Overall sens"),  expression("Overall spec"))
  
  methods <- c("CC", "BI","Full","BIgene")
  
  nReplicates <-  as.numeric(args[12]) #50 replicates (MCMC chains)
  
  qcut <- as.numeric(args[14])
  #    colors <- c(colors()[461],  colors()[282], colors()[555], colors()[610])
  colors<- c(colors()[26],  colors()[261], colors()[35], colors()[614],  colors()[621], colors()[121]) #blue, grey, red, green, yellow. light blue[6]
  
  print( nReplicates)
  print(paste("CausalM",nCausalM))
  print(paste("CausalMbar",nCausalMbar))
  tmp <- rep(NA, nReplicates*length(qcutAll)*length(all_n)) 
  yHatAll <- array(tmp, c(length(all_n),nReplicates, length(qcutAll)))  #different from all the other, stores MSE(yHat-yTrue) for each replicate. Store different for each qcut
  yHatAll_CC <- array(tmp, c(length(all_n),nReplicates, length(qcutAll))) 
  yHatAll_Full <- array(tmp, c(length(all_n),nReplicates, length(qcutAll))) 
  yHatAll_Gene <- array(tmp, c(length(all_n),nReplicates, length(qcutAll))) 
  
  
  tmp <- rep(NA, length(all_n)*2*length(qcutAll))
  mse <- array(tmp, c(length(all_n), 2, length(qcutAll)))
  mseCC <- array(tmp, c(length(all_n), 2, length(qcutAll)))
  mseFull <- array(tmp, c(length(all_n), 2, length(qcutAll)))
  mseGene <- array(tmp, c(length(all_n), 2, length(qcutAll)))
  
  bi <- list()
  biCC<- list()
  biFull <- list()
  biGene <- list()
  est <- list()
  estCC <- list()
  estFull <- list()
  estGene <- list()
  
  pdf(paste(wd,"/TP.pdf",sep=""))
  par(mfrow=c(1,1),mar=c(5.1,4.1,1.9,2.1),pin=c(4,4))
  
  repLoop <- 1:nReplicates 
  #    if(thetaMethy!=0 & thetaGene!=0) ### manually remove unconverged replications by Geweke
  #    repLoop <- (1:nReplicates)[-2]
  
  set.seed(seed)
  for (n in 1:length(all_n)){
    
    #        qMSilverIdx <- 1:10
    #        qMbarSilverIdx <- 1:10 #fix      
    #qMSilverIdx <- which(estSilver$qM < qcut)
    #qMbarSilverIdx <- which(estSilver$qMbar < qcut)
    
    #qMSilverIdx <- order(estSilver$qM)
    #qMbarSilverIdx <- order(estSilver$qMbar)
    
    #Instead of comparing with silver standard, we now compare with the results from each of the full training set:
    
    qSilverIdx <- list()
    geneNamesSilver <- list()
    tp <- rep(0, nGene) #combining M and Mbar together , so top 1 p is 1M+1Mbar
    tpCC <- rep(0, nGene)
    tpGene <- rep(0, nGene)
    tpCV <- rep(0,nGene)
    
    pp <- matrix(0,10, nGene) #combining M and Mbar together , so top 1 p is 1M+1Mbar
    ppCC <- matrix(0,10, nGene)
    ppGene <- matrix(0,10, nGene)
    ppCV <- matrix(0,10,nGene)
    decisionFinal <- rep(1,length(repLoop))
    for (r in repLoop){
      #            wdSub <- paste(wd,"/n",n,"_r",r,"_Full",sep="")
      #            load(paste(wdSub,"/run.RData",sep=""))
      #            qSilverIdx[[r]] <- order(estFull[[r]]$q)
      #            qSilverIdx[[r]] <- order(estSilver$q) #fix: use silver instead of full for each fold
      #            geneNamesSilver[[r]] <- geneNames[order(estFull[[r]]$q)]
      geneNamesSilver[[r]] <- geneNames[order(estSilver$q)]
      #            yHatAll_Full[n,r,] <- as.numeric(estFull[[r]]$MSEyHatEst) #MSE(yHat-yTrue)
    }
    for (r in repLoop){
      wdSub <- paste(wd,"/n",n,"_r",r,"_CV",sep="")
      load(paste(wdSub,"/run.RData",sep=""))
      decisionFinal[r] <- which(decisions == max(decisions))[1]   
    }
    decisions <- decisionFinal
    repTmp <- 0
    repTmpCV <- 0 #repTmpCV
    for (r in repLoop){
      wdSub <- paste(wd,"/n",n,"_r",r,sep="")
      load(paste(wdSub,"/run.RData",sep=""))
      
      yHatAll[n,r,] <- as.numeric(est[[r]]$MSEyHatEst) #MSE(yHat-yTrue)
      
      #            qMIdx <- order(est[[r]]$qM)
      #            qMbarIdx <- order(est[[r]]$qMbar)
      shuffleIdx <- sample(1:nGene)
      #          shuffleIdx <- 1:nGene
      #            qIdx <- order(qShuffled)
      geneNamesShuffle <- (geneNames[shuffleIdx])[ order(est[[r]]$q[shuffleIdx]) ]
      for(i in 1:nGene){
        #                tp[i] = tp[i] + length( which(qMIdx[1:i] %in% qMSilverIdx[1:i])) +  length( which(qMbarIdx[1:i] %in% qMbarSilverIdx[1:i])) #intead of comparing
        tp[i] = tp[i] + length( which(geneNamesShuffle[1:i] %in% geneNamesSilver[[r]][1:i]) )
      }
      repTmp <- repTmp +1
      
      
      for(i in 1:nGene){
        #                tp[i] = tp[i] + length( which(qMIdx[1:i] %in% qMSilverIdx[1:i])) +  length( which(qMbarIdx[1:i] %in% qMbarSilverIdx[1:i])) #intead of comparing
        tpCV[i] = tpCV[i] + length( which(geneNamesShuffle[1:i] %in% geneNamesSilver[[r]][1:i]) )
      }
      repTmpCV <- repTmpCV + 1
      
      
      for(i in 1:nGene){
        pp[r,i] = length( which(geneNamesShuffle[1:i] %in% geneNamesSilver[[r]][1:i]) )
      }
    }
    tp <- tp/repTmp
    tpCV <- tpCV /repTmpCV
    pp[r,] <- pp[r,]
    
    
    repTmp <- 0
    for (r in repLoop){
      wdSub <- paste(wd,"/n",n,"_r",r,"_CC",sep="")
      load(paste(wdSub,"/run.RData",sep=""))
      
      yHatAll_CC[n,r,] <- as.numeric(estCC[[r]]$MSEyHatEst) 
      
      #            qMIdx <- order(estCC[[r]]$qM)
      #            qMbarIdx <- order(estCC[[r]]$qMbar)
      shuffleIdx <- sample(1:nGene)
      geneNamesShuffle <- (geneNames[shuffleIdx])[ order(estCC[[r]]$q[shuffleIdx]) ]
      
      for(i in 1:nGene){
        tpCC[i] = tpCC[i] + length( which(geneNamesShuffle[1:i] %in% geneNamesSilver[[r]][1:i]) )
      }
      repTmp <- repTmp +1
      for(i in 1:nGene){
        ppCC[r,i] = length( which(geneNamesShuffle[1:i] %in% geneNamesSilver[[r]][1:i]) )
      }
      
    }
    tpCC <- tpCC/repTmp
    
    
    ppCC[r,] <- ppCC[r,]
    
    if(thetaGene !=0 & thetaMethy!=0){
      repTmp <- 0
      for (r in repLoop){
        wdSub <- paste(wd,"/n",n,"_r",r,"_Gene",sep="")
        load(paste(wdSub,"/run.RData",sep=""))
        yHatAll_Gene[n,r,] <- as.numeric(estGene[[r]]$MSEyHatEst) 
        #            qMIdx <- order(estCC[[r]]$qM)
        #            qMbarIdx <- order(estCC[[r]]$qMbar)
        #                qIdx <- order(estGene[[r]]$q)
        shuffleIdx <- sample(1:nGene)  
        geneNamesShuffle <- (geneNames[shuffleIdx])[ order(estGene[[r]]$q[shuffleIdx]) ]
        
        for(i in 1:nGene){
          tpGene[i] = tpGene[i] + length( which(geneNamesShuffle[1:i] %in% geneNamesSilver[[r]][1:i]) )
        }
        repTmp <- repTmp+1
        
        for(i in 1:nGene){
          ppGene[r,i] = length( which(geneNamesShuffle[1:i] %in% geneNamesSilver[[r]][1:i]) )
        }
      }
      #            print(tpGene)
      #print(repTmp)
      tpGene <- tpGene/repTmp
      ppGene[r,] <- ppGene[r,]
      
      
    }
    
    
    #        maxIdx <- min( max( min(which( tp == 1)), min(which( tpCC == 1)) )*1.4, nGene)
    
    maxIdx <- c(nGene,10, 20, 50, 100, 200, 300, 500, length(which(estSilver$q !=1 )))
    
    for(m in 1:length(maxIdx)){
      ylimTmp <- max(tp[maxIdx[m]], tpCC[maxIdx[m]])*1.2
      
      plot(0:maxIdx[m], c(0,tp[1:maxIdx[m]]), ylim=c(0, ylimTmp), main=NA, xlab="Positive features",ylab="Overlapped features",sub=paste("nSample",all_n[n]),col=colors[1],type="l")
      #  abline(h=1, col=colors[2],lty=2)
      curve(0.001*x^2, 0,maxIdx[m], add = TRUE, col = colors[2], lwd=0.7)
      
      lines(0:maxIdx[m], c(0,tpCC[1:maxIdx[m]]), col=colors[3],type="l")
      
      
      if(thetaGene !=0 & thetaMethy!=0){
        lines(0:maxIdx[m], c(0,tpGene[1:maxIdx[m]]), col=colors[6],type="l")
        legend("bottomright", bty = "n", c("CC","BI","BIgene"), cex=1,  col=colors[c(3,1,6)], lty=c(1,1))
      } else{
        legend("bottomright", bty = "n", c("CC","BI"), cex=1,  col=colors[c(3,1)], lty=c(1,1))
      }
      
      r = repLoop[1]
      plot(0:maxIdx[m], c(0,pp[r,1:maxIdx[m]]), ylim=c(0, ylimTmp), main=NA, xlab="Positive features",ylab="Overlapped features",sub=paste("nSample",all_n[n]),col=colors()[180+r],type="l")
      lines(0:maxIdx[m], c(0,ppCC[r,1:maxIdx[m]]), col=colors()[180+r],lty=2)
      if(thetaGene !=0 & thetaMethy!=0)
        lines(0:maxIdx[m], c(0,ppGene[r,1:maxIdx[m]]), col=colors()[180+r],lty=3)
      #                legend("bottomright", bty = "n", c("CC","BI"), cex=1,  col=colors[c(3,1)], lty=c(1,1))
      
      for(r in repLoop[-1]){
        lines(0:maxIdx[m], c(0,pp[r,1:maxIdx[m]]), col=colors()[180+ r],lty=1)
        text(maxIdx[m]-2*r, pp[r,maxIdx[m]], r)
        lines(0:maxIdx[m], c(0,ppCC[r,1:maxIdx[m]]), col=colors()[180+r],lty=2)
        text(maxIdx[m]-2*r, ppCC[r,maxIdx[m]], r)
        if(thetaGene !=0 & thetaMethy!=0){
          lines(0:maxIdx[m], c(0,ppGene[r,1:maxIdx[m]]), col=colors()[180+r],lty=3)
          text(maxIdx[m]-2*r, ppGene[r,maxIdx[m]], r)
          #                    legend("bottomright", bty = "n", c("CC","BI","BIgene"), cex=1,  col=colors[c(3,1,6)], lty=c(1,1))
        }
      }
      
    }
    if (thetaGene==0)
      save(tp, tpCC, tpGene,tpCV,file=paste("~/BI/SAll/0_05_pp.RData", sep=""))
    if (thetaMethy==0)
      save(tp, tpCC, tpGene,tpCV,file=paste("~/BI/SAll/05_0_pp.RData", sep=""))
    if (thetaGene!=0 & thetaMethy!=0)
      save(tp, tpCC, tpGene,tpCV,file=paste("~/BI/SAll/02_02_pp.RData", sep=""))
    
  }
  dev.off()
  
  
  for(q in 1:length(qcutAll)){
    
    for(n in 1:length(all_n)){
      mse[n,1,q] <- signif( mean(yHatAll[n,,q],na.rm=TRUE), 3)
      mse[n,2,q] <- signif( sd(yHatAll[n,,q],na.rm=TRUE)/sqrt(length(repLoop)), 3)
      mseCC[n,1,q] <- signif( mean(yHatAll_CC[n,,q],na.rm=TRUE), 3)
      mseCC[n,2,q] <- signif( sd(yHatAll_CC[n,,q],na.rm=TRUE)/sqrt(length(repLoop)), 3)
      #   mseFull[n,1,q] <- signif( mean(yHatAll_Full[n,,q],na.rm=TRUE), 3)
      #   mseFull[n,2,q] <- signif( sd(yHatAll_Full[n,,q],na.rm=TRUE)/sqrt(length(repLoop)), 3)
      if(thetaGene !=0 & thetaMethy!=0){
        mseGene[n,1,q] <- signif( mean(yHatAll_Gene[n,,q],na.rm=TRUE), 3)
        mseGene[n,2,q] <- signif( sd(yHatAll_Gene[n,,q],na.rm=TRUE)/sqrt(length(repLoop)), 3)
      }
    }
  }
  save(mse, mseCC,mseFull,mseGene, yHatAll, yHatAll_CC, yHatAll_Full,yHatAll_Gene,all_n, decisions, file=paste(wd, "/mse.RData", sep=""))
  
  ##     print(mse)
  ##     print(mseCC)
  ##     print(mseFull)
  ##     print(mseGene)
  
  ##     pdf(paste(wd,"/MSEyHat.pdf",sep=""))
  ##     par(mfrow=c(1,1),mar=c(5.1,4.1,1.9,2.1),pin=c(4,4))
  ##     x <- 1:length(all_n)
  ##     barEpsilon <- 0.04
  
  ##                                         #plot MSE(yhat) and feature selections separatedly
  ##     for(q in 1:length(qcutAll)){    
  ##                                         #        ylimTmp <- max(mse[x,1,q],mseCC[x,1,q], mseFull[x,1,q]) + max(mse[x,2,q],mseCC[x,2,q], mseFull[x,2,q])
  ##         ylimTmp <- max(mse[x,1,q],mseCC[x,1,q]) + max(mse[x,2,q],mseCC[x,2,q]+10) 
  ##         ylimTmp <- ylimTmp*1.2
  
  ##                                         #       ylimTmpNeg <- min(mse[x,1,q],mseCC[x,1,q],mseFull[x,1,q])-max(mse[x,2,q],mseCC[x,2,q],mseFull[x,2,q])
  ##         ylimTmpNeg <- min(mse[x,1,q],mseCC[x,1,q])-max(mse[x,2,q],mseCC[x,2,q])
  ##         ylimTmpNeg<-min(ylimTmpNeg*1.1,-10)
  
  ##         plot(x, mse[x,1,q], ylim=c(ylimTmpNeg,ylimTmp),xlab="N",main=paste("FDR <",qcutAll[q]),
  ##              ylab="MSE",xaxt="n",col=colors[1],type="l")
  ##         axis(1,at=x,label=all_n[x],las=2,cex.axis=0.8)
  ##         segments(x,mse[x,1,q]-mse[x,2,q], x, mse[x,1,q]+mse[x,2,q],col=colors[1])
  ##         segments(x-barEpsilon, mse[x,1,q]-mse[x,2,q], x+barEpsilon, mse[x,1,q]-mse[x,2,q],col=colors[1] )
  ##         segments(x-barEpsilon, mse[x,1,q]+mse[x,2,q], x+barEpsilon, mse[x,1,q]+mse[x,2,q],col=colors[1] )
  
  ##         lines(x+barEpsilon, mseCC[x,1,q],col=colors[3],type="l")
  ##         segments(x+barEpsilon,mseCC[x,1,q]-mseCC[x,2,q], x+barEpsilon, mseCC[x,1,q]+mseCC[x,2,q],col=colors[3])
  ##         segments(x, mseCC[x,1,q]-mseCC[x,2,q], x+2*barEpsilon, mseCC[x,1,q]-mseCC[x,2,q],col=colors[3] )
  ##         segments(x, mseCC[x,1,q]+mseCC[x,2,q], x+2*barEpsilon, mseCC[x,1,q]+mseCC[x,2,q],col=colors[3] )
  
  ## #        lines(x+barEpsilon, mseFull[x,1,q],col=colors[4],type="l")
  ## #        segments(x+barEpsilon,mseFull[x,1,q]-mseFull[x,2,q], x+barEpsilon, mseFull[x,1,q]+mseFull[x,2,q],col=colors[4])
  ## #        segments(x, mseFull[x,1,q]-mseFull[x,2,q], x+2*barEpsilon, mseFull[x,1,q]-mseFull[x,2,q],col=colors[4] )
  ## #        segments(x, mseFull[x,1,q]+mseFull[x,2,q], x+2*barEpsilon, mseFull[x,1,q]+mseFull[x,2,q],col=colors[4] )
  
  ##         if(thetaGene !=0 & thetaMethy!=0){
  ##             lines(x+barEpsilon, mseGene[x,1,q],col=colors[6],type="l")
  ##             segments(x+barEpsilon,mseGene[x,1,q]-mseGene[x,2,q], x+barEpsilon, mseGene[x,1,q]+mseGene[x,2,q],col=colors[6])
  ##             segments(x, mseGene[x,1,q]-mseGene[x,2,q], x+2*barEpsilon, mseGene[x,1,q]-mseGene[x,2,q],col=colors[6] )
  ##             segments(x, mseGene[x,1,q]+mseGene[x,2,q], x+2*barEpsilon, mseGene[x,1,q]+mseGene[x,2,q],col=colors[6] )
  ##                                         #            legend("topright", bty = "n", methods, cex=0.8, col=colors[c(3,1,4,6)], lty=c(1,1,1,1))
  ##             legend("topright", bty = "n", methods[c(1,2,4)], cex=0.8, col=colors[c(3,1,6)], lty=c(1,1,1))
  
  ##         } else{
  ##             legend("topright", bty = "n", methods[1:3], cex=0.8, col=colors[c(3,1)], lty=c(1,1))
  ##         }
  ##                                         #        abline(epsilon^2,0,col=colors[2],lty=2) #mse of yhat should shrink to sigma^2
  
  
  ##     }
  
  ##     dev.off()
  
  print("Done for varying nSample")
  
  
}



