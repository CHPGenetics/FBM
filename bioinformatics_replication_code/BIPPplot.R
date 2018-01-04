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



pdf("~/BI/SAll/PP.pdf",width=8.6, height=3.5)
par(mfrow=c(1,3),mar=c(1.3,2,1,1)+0.5, oma=c(3.5,1,1,-0.3)+0.5)

colors<- c(colors()[26],  colors()[274], colors()[35], colors()[128],  colors()[621])

n=1
allthetaTitle <- paste( c(rep("Missing gene: ",1),rep("Missing methy: ",1), rep("Missing both: ",1)), c(50,50,20),rep("%",3),sep="")
maxIdx <- 100
load("~/BI/SAll/05_0_pp.RData")

m=1
ylimTmp <- 20
plot(0:maxIdx[m], c(0,tp[1:maxIdx[m]]), ylim=c(0, ylimTmp),xlab="Positive features",ylab="Overlapped featu",,col=colors[1],type="l", main=allthetaTitle[1])
                                        #  abline(h=1, col=colors[2],lty=2)
            curve(0.001*x^2, 0,maxIdx[m], add = TRUE, col = colors[2], lwd=0.7)

            lines(0:maxIdx[m], c(0,tpCC[1:maxIdx[m]]), col=colors[3],type="l")
            lines(0:maxIdx[m], c(0,tpCV[1:maxIdx[m]]), col=colors[5],type="l")
            
mtext("Positive features", side=1, line=2.2, cex=0.8)
mtext("Overlapped features", side=2, line=2.4, cex=0.8)


load("~/BI/SAll/0_05_pp.RData")
plot(0:maxIdx[m], c(0,tp[1:maxIdx[m]]), ylim=c(0, ylimTmp), xlab="Positive features",ylab="Overlapped featu",,col=colors[1],type="l", main=allthetaTitle[2])
                                        #  abline(h=1, col=colors[2],lty=2)
            curve(0.001*x^2, 0,maxIdx[m], add = TRUE, col = colors[2], lwd=0.7)

            lines(0:maxIdx[m], c(0,tpCC[1:maxIdx[m]]), col=colors[3],type="l")
            lines(0:maxIdx[m], c(0,tpCV[1:maxIdx[m]]), col=colors[5],type="l")

mtext("Positive features", side=1, line=2.2, cex=0.8)
mtext("Overlapped features", side=2, line=2.4, cex=0.8)

load("~/BI/SAll/02_02_pp.RData")
plot(0:maxIdx[m], c(0,tp[1:maxIdx[m]]), ylim=c(0, ylimTmp),xlab="Positive features",ylab="Overlapped featu",,col=colors[1],type="l", main=allthetaTitle[3])
                                        #  abline(h=1, col=colors[2],lty=2)
            curve(0.001*x^2, 0,maxIdx[m], add = TRUE, col = colors[2], lwd=0.7)

            lines(0:maxIdx[m], c(0,tpCC[1:maxIdx[m]]), col=colors[3],type="l")
            lines(0:maxIdx[m], c(0,tpCV[1:maxIdx[m]]), col=colors[5],type="l")

lines(0:maxIdx[m], c(0,tpGene[1:maxIdx[m]]), col=colors[4],type="l")
mtext("Positive features", side=1, line=2.2, cex=0.8)
mtext("Overlapped features", side=2, line=2.4, cex=0.8)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(10, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottomright", c("CC","BI","BIgene","CV"), xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n",col=colors[c(3,1,4,5)],lty=c(1,1,1,1), cex = 1.2)




    dev.off()
