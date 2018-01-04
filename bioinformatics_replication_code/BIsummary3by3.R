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
                                        #    all_n <- c(50,100,200,500)
    methods <- c("CC", "BI", "Full")
    colors<- c(colors()[26],  colors()[261], colors()[35], colors()[614],  colors()[621])
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

    pdf(paste(wd,"/BI/ROC3by3.pdf",sep=""))
    par(mfrow=c(3,3), oma=c(1.5,3,1,1)+0.5, mar=c(1,1,1,1)+0.5)
                                        #par(mfrow=c(3,3),mar=c(5.1,4.1,1.9,2.1),pin=c(4,4))

        allthetaName <- c("01_0","025_0","05_0", "0_01", "0_025", "0_05", "01_01", "02_02", "03_03")
    allthetaTitle <- paste( c(rep("Missing gene: ",3),rep("Missing methy: ",3), rep("Missing both: ",3)), c(10,25,50,10,25,50,10,20,30),rep("%",9),sep="")
                                        #    alltheta <- c("01_0","025_0","05_0", "0_015", "0_025", "0_05", "02_02", "025_025", "03_03")
        for(t in 1:length(allthetaName)){
            for (n in 1:3){ #50 100 200
                load(paste(wd,"/BI_", allthetaName[t],"/data_n",n,".RData",sep=""))

                auc <- signif(auc,3)
                aucCC <- signif(aucCC,3)
                aucFull <- signif(aucFull,3)
                
                if( (t-1)%%3 ==0 & ( (t)%%3==0 & (n-1)%%3==0 )){
                    plot(0:1,0:1,main=paste("nSample",all_n[n]),
                         xlab="1-Specificity",ylab="Sensitivity",xaxt="n",col=colors[2],type="l",lty=2)
                }else if( (t-1)%%3 ==0 & ( (t)%%3==0 & (n-1)%%3!=0 ) ) {
                    plot(0:1,0:1,main=paste("nSample",all_n[n]),
                         xlab="1-Specificity",ylab=NA,xaxt="n",col=colors[2],type="l",lty=2)
                }else if( (t-1)%%3 ==0 & ( (t)%%3!=0 & (n-1)%%3==0 ) ) {
                    plot(0:1,0:1,main=paste("nSample",all_n[n]),
                         xlab=NA,ylab="Sensitivity",xaxt="n",col=colors[2],type="l",lty=2)
                }else if( (t-1)%%3 ==0 & ( (t)%%3!=0 & (n-1)%%3!=0 ) ) {
                    plot(0:1,0:1,main=paste("nSample",all_n[n]),
                         xlab=NA,ylab=NA,xaxt="n",col=colors[2],type="l",lty=2)
                }             else if( (t-1)%%3 !=0 & ( (t)%%3==0 & (n-1)%%3==0 )){ ####divide here
                    plot(0:1,0:1,main=NA,
                         xlab="1-Specificity",ylab="Sensitivity",xaxt="n",col=colors[2],type="l",lty=2)
                }else if( (t-1)%%3 !=0 & ( (t)%%3==0 & (n-1)%%3!=0 ) ) {
                    plot(0:1,0:1,main=NA,
                         xlab="1-Specificity",ylab=NA,xaxt="n",col=colors[2],type="l",lty=2)
                }else if( (t-1)%%3 !=0 & ( (t)%%3!=0 & (n-1)%%3==0 ) ) {
                    plot(0:1,0:1,main=NA,
                         xlab=NA,ylab="Sensitivity",xaxt="n",col=colors[2],type="l",lty=2)
                }else if( (t-1)%%3 !=0 & ( (t)%%3!=0 & (n-1)%%3!=0 ) ) {
                    plot(0:1,0:1,main=NA,
                         xlab=NA,ylab=NA,xaxt="n",col=colors[2],type="l",lty=2)
                }

                if((n)%%3==0 ){
                                        #                axis(4)            
                    mtext(allthetaTitle[t], side=4, line=0.8)
                }

                if((n-1)%%3==0 ){
                                        #                axis(4)            
                    mtext("Sensitivity", side=2, line=2.2, cex=0.8)
                }
                if(t%%3==0 ){
                                        #                axis(4)            
                    mtext("1-Specificity", side=1, line=2.2, cex=0.8)
                }
                                        #       xlab="1-Specificity",ylab="Sensitivity",xaxt="n",sub=paste("nSample",all_n[n]),col=colors[2],type="l",lty=2)
                
                axis(1,at=c(0,0.25,0.5,0.75,1),label=c(0,0.25,0.5,0.75,1),cex.axis=1)

                polygon(c(rev(roc[,4]),roc[,4]), c(rev(roc[,1]), roc[,3]), col = alpha(colors[1],0.2), border=NA)
                polygon(c(rev(rocCC[,4]),rocCC[,4]), c(rev(rocCC[,1]), rocCC[,3]), col = alpha(colors[3],0.2), border=NA)
                polygon(c(rev(rocFull[,4]),rocFull[,4]), c(rev(rocFull[,1]), rocFull[,3]), col = alpha(colors[4],0.2), border=NA)

                lines(roc[,4], roc[,2], col=colors[1],type="l", lwd=0.7)
                lines(rocCC[,4], rocCC[,2], col=colors[3],type="l", lwd=0.7)
                lines(rocFull[,4], rocFull[,2], col=colors[4],type="l", lwd=0.7)

                points(roc[which(youden==max(youden))[1],4], roc[which(youden==max(youden))[1],2], pch=8, cex=1.5) #changed
                points(rocCC[which(youdenCC==max(youdenCC))[1],4], rocCC[which(youdenCC==max(youdenCC))[1],2], pch=8, cex=1.5)
                points(rocFull[which(youdenFull==max(youdenFull))[1],4], rocFull[which(youdenFull==max(youdenFull))[1],2], pch=8, cex=1.5)
                legend("bottomright", bty = "n", c(paste("CC(AUC=",aucCC,")",sep=""), paste("BI(AUC=",auc,")",sep=""), paste("Full(AUC=",aucFull,")",sep="")), cex=1,  col=colors[c(3,1,4)], lty=c(1,1,1))

                ##       plot(c(0,1),c(0,1)) #changed 
                ##       text(min(max(roc[which(youden==max(youden))[1],4],0.1),0.9), min(max(roc[which(youden==max(youden))[1],2]*1.1,0.1),0.9), labels = signif(max(youden),3), col=colors[1])
                ##       text(min(max(rocCC[which(youdenCC==max(youdenCC))[1],4],0.1),0.9),min(max(rocCC[which(youdenCC==max(youdenCC))[1],2]*1.1,0.1),0.9), labels = signif(max(youdenCC),3), col=colors[3])
                ## text(min(max(rocFull[which(youdenFull==max(youdenFull))[1],4],0.1),0.9),min(max(rocFull[which(youdenFull==max(youdenFull))[1],2]*1.1,0.1),0.9), labels = signif(max(youdenFull),3), col=colors[4])
                                        #            legend("bottomright", bty = "n", c(paste("CC(AUC=",aucCC,")",sep=""), paste("BI(AUC=",auc,")",sep=""), paste("Full(AUC=",aucFull,")",sep="")), cex=1,  col=colors[c(3,1,4)], lty=c(1,1,1))
            }
        }

        dev.off()

        barEpsilon=0.02

        pdf(paste(wd,"/BI/RMSE3by3.pdf",sep=""))
        par(mfrow=c(3,3))
        par(oma = c(4, 1.5, 1, 0)+0.1, mar=c(2,2,2,0)+0.2)
 #       par(mfrow=c(3,3), oma=c(5,4,0,0)+0.1, mar=c(0,0,1,1.5)+0.1)

                                        #    par(mfrow=c(3,3),mar=c(5.1,4.1,1.9,2.1),pin=c(4,4))
        i=7
        xZoom <- list()
        xZoom[[1]] <- seq(1:length(all_n))

        for(t in 1:length(allthetaName)){
            load(paste(wd,"/BI_", allthetaName[t],"/mse.RData",sep=""))
          
            x <- xZoom[[1]]
            ylimTmp <- max(sqrt(mseFull[i,x,1]),sqrt(mseCC[i,x,1]),sqrt(mse[i,x,1])) + max(sqrt(mseCC[i,x,2]),sqrt(mseFull[i,x,2]),sqrt(mse[i,x,2]))
            ylimTmp <- ylimTmp*1.2

            
            ylimTmpNeg <- min(sqrt(mseCC[i,x,1]),sqrt(mseFull[i,x,1]),sqrt(mse[i,x,1]))-max(sqrt(mseCC[i,x,2]),sqrt(mseFull[i,x,2]),sqrt(mse[i,x,2]))
            ylimTmpNeg<-ylimTmpNeg*1.1

            ylimAll <- c(rep(90,3),rep(100,3),rep(110,3))
                                        #            plot(x, mse[i,x,1], ylim=c(ylimTmpNeg,ylimTmp),xlab="N",main=parNames[i],
            plot(x, sqrt(mse[i,x,1]), ylim=c(0,100),xlab="N",main=allthetaTitle[t],
                 ylab="RMSE",xaxt="n",col=colors[1],type="l")
            axis(1,at=x,label=all_n[x],las=2,cex.axis=0.8)

            if((t-1)%%3==0 ){
                mtext(expression(paste("RMSE(",hat(Y),")",sep="")), side=2, line=2.2, cex=0.8)
            }
                if(t%in% c(7,8,9) ){
                    mtext("nSample", side=1, line=2.2, cex=0.8)
                }

            segments(x,sqrt(mse[i,x,1])-sqrt(mse[i,x,2]), x, sqrt(mse[i,x,1])+sqrt(mse[i,x,2]),col=colors[1])
            segments(x-barEpsilon, sqrt(mse[i,x,1])-sqrt(mse[i,x,2]), x+barEpsilon, sqrt(mse[i,x,1])-sqrt(mse[i,x,2]),col=colors[1] )
            segments(x-barEpsilon, sqrt(mse[i,x,1])+sqrt(mse[i,x,2]), x+barEpsilon, sqrt(mse[i,x,1])+sqrt(mse[i,x,2]),col=colors[1] )
            
            lines(x+barEpsilon, sqrt(mseCC[i,x,1]),col=colors[3],type="l")
            segments(x+barEpsilon,sqrt(mseCC[i,x,1])-sqrt(mseCC[i,x,2]), x+barEpsilon, sqrt(mseCC[i,x,1])+sqrt(mseCC[i,x,2]),col=colors[3])
            segments(x, sqrt(mseCC[i,x,1])-sqrt(mseCC[i,x,2]), x+2*barEpsilon, sqrt(mseCC[i,x,1])-sqrt(mseCC[i,x,2]),col=colors[3] )
            segments(x, sqrt(mseCC[i,x,1])+sqrt(mseCC[i,x,2]), x+2*barEpsilon, sqrt(mseCC[i,x,1])+sqrt(mseCC[i,x,2]),col=colors[3] )

            lines(x+barEpsilon, sqrt(mseFull[i,x,1]),col=colors[4],type="l")
            segments(x+barEpsilon,sqrt(mseFull[i,x,1])-sqrt(mseFull[i,x,2]), x+barEpsilon, sqrt(mseFull[i,x,1])+sqrt(mseFull[i,x,2]),col=colors[4])
            segments(x, sqrt(mseFull[i,x,1])-sqrt(mseFull[i,x,2]), x+2*barEpsilon, sqrt(mseFull[i,x,1])-sqrt(mseFull[i,x,2]),col=colors[4] )
            segments(x, sqrt(mseFull[i,x,1])+sqrt(mseFull[i,x,2]), x+2*barEpsilon, sqrt(mseFull[i,x,1])+sqrt(mseFull[i,x,2]),col=colors[4] )

            abline(0,0,col=colors[2],lty=2) #mse of yhat should shrink to sigma^2
                                        #            legend("topright", bty = "n", methods, cex=0.8, col=colors[c(3,1,4)], lty=c(1,1,1))
            

        }
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(10, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        legend("bottomright", methods, xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n",col=colors[c(3,1,4)],lty=c(1,1,1), cex = 1.5)
        
        dev.off()
        print("Done for varying nSample")

    }



