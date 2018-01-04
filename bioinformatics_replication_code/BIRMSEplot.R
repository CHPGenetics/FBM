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
qcutAll <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)

localTest <- 0
singleTest <- 0
multipleReplicates <- 0
fullRun <- 0


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

allthetaTitle <- paste( c(rep("Missing gene: ",1),rep("Missing methy: ",1), rep("Missing both: ",1)), c(50,50,20),rep("%",3),sep="")
colors<- c(colors()[26],  colors()[261], colors()[35], colors()[614],  colors()[621], colors()[121]) #blue, grey, red, green, yellow. light blue[6]

pdf("./SAll/asthma_RMSE_RMSE.pdf")
    par(mfrow=c(2,2), oma=c(2.5,3,0.5,0.5)+1, mar=c(2,1,1.5,1)+1)
#par(mfrow=c(2,2),mar=c(5.1,4.1,1.9,2.1),pin=c(4,4))
load("./SAll/mse_05_0.RData")
pType <- rep(4,10)
colType <- rep(1,10)
for(i in 1:10){
    if((decisions[i]==2 & yHatAll[i] < yHatAll_CC[i]) |(decisions[i]==1 & yHatAll[i] > yHatAll_CC[i]) )
        pType[i] <- 1
    if(yHatAll[i] > yHatAll_CC[i])
        colType[i] <- 3
}
theMax <- max(yHatAll, yHatAll_CC)
plot(yHatAll_CC, yHatAll, xlim=c(0,theMax), ylim=c(0,theMax), xlab="RMSE of CC", ylab="RMSE of BI", pch=pType,col=colors[colType], main=allthetaTitle[1])
lines(0:theMax, 0:theMax, col=colors[2],type="l")

mtext("RMSE (BI)", side=2, line=2.2, cex=1)
mtext("RMSE (CC)", side=1, line=2.2, cex=1)

load("./SAll/mse_0_05.RData")
pType <- rep(4,10)
colType <- rep(1,10)
for(i in 1:10){
    if((decisions[i]==2 & yHatAll[i] < yHatAll_CC[i]) |(decisions[i]==1 & yHatAll[i] > yHatAll_CC[i]) )
        pType[i] <- 1
    if(yHatAll[i] > yHatAll_CC[i])
        colType[i] <- 3
}
theMax <- max(yHatAll, yHatAll_CC)
plot(yHatAll_CC, yHatAll, xlim=c(0,theMax), ylim=c(0,theMax), xlab="RMSE of CC", ylab="RMSE of BI", pch=pType,col=colors[colType], main=allthetaTitle[2])
lines(0:theMax, 0:theMax, col=colors[2],type="l")

mtext("RMSE (BI)", side=2, line=2.2, cex=1)
mtext("RMSE (CC)", side=1, line=2.2, cex=1)


load("./SAll/mse_02_02.RData")
pType <- rep(4,10)
colType <- rep(1,10)
for(i in 1:10){
    if((decisions[i]!=1 & yHatAll[i] < yHatAll_CC[i]) |(decisions[i]!=2 & yHatAll[i] > yHatAll_CC[i]) )
        pType[i] <- 1
    if(yHatAll[i] > yHatAll_CC[i])
        colType[i] <- 3
}
theMax <- max(yHatAll, yHatAll_CC)
plot(yHatAll_CC, yHatAll, xlim=c(0,theMax), ylim=c(0,theMax), xlab="RMSE of CC", ylab="RMSE of BI", pch=pType,col=colors[colType], main=allthetaTitle[3])
lines(0:theMax, 0:theMax, col=colors[2],type="l")
mtext("RMSE (BI)", side=2, line=2.2, cex=1)
mtext("RMSE (CC)", side=1, line=2.2, cex=1)


pType <- rep(4,10)
colType <- rep(3,10)
for(i in 1:10){
    if((decisions[i]!=1 & yHatAll_Gene[i] < yHatAll_CC[i]) |(decisions[i]!=3 & yHatAll_Gene[i] > yHatAll_CC[i]) )
        pType[i] <- 1
    if(yHatAll_Gene[i] < yHatAll_CC[i])
        colType[i] <- 6
}
theMax <- max(yHatAll, yHatAll_CC)
plot(yHatAll_CC, yHatAll_Gene, xlim=c(0,theMax), ylim=c(0,theMax), xlab="RMSE of CC", ylab="RMSE of BIgene", pch=pType,col=colors[colType], main=allthetaTitle[3])
lines(0:theMax, 0:theMax, col=colors[2],type="l")
mtext("RMSE (BIgene)", side=2, line=2.2, cex=1)
mtext("RMSE (CC)", side=1, line=2.2, cex=1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(10, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottomright", c("CC","BI","BIgene") , xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n",col=colors[c(3,1,6)],pch=c(1,1,1), cex = 1.5)

dev.off()
