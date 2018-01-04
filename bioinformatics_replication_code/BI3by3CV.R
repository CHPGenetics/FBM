args = commandArgs(trailingOnly=TRUE)

require(MCMCpack)
library(doParallel)
library(foreach)
library(cluster)
library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(scales)
library(ggplot2)
                                        #Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
                                        #sourceCpp("~/BI/BI.cpp")

localTest <- 0
singleTest <- 0
multipleReplicates <- 0
fullRun <- 0
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

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
    methods <- c("CC", "BI", "BIgene")
    colors<- c(colors()[26],  colors()[261], colors()[35], colors()[614],  colors()[621], colors()[121])
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

                                        #par(mfrow=c(3,3),mar=c(5.1,4.1,1.9,2.1),pin=c(4,4))

                                        #        allthetaName <- c("01_0","025_0","05_0", "0_01", "0_025", "0_05", "01_01", "02_02", "03_03")
                                        #   allthetaTitle <- paste( c(rep("Missing gene: ",3),rep("Missing methy: ",3), rep("Missing both: ",3)), c(10,25,50,10,25,50,10,20,30),rep("%",9),sep="")

    allthetaName <- c("05_0","0_05", "025_025")
    allthetaTitle <- paste( c(rep("Missing gene: ",1),rep("Missing methy: ",1), rep("Missing both: ",1)), c(50,50,25),rep("%",3),sep="")

                                        #    allthetaName <- c("05_0")
                                        #    allthetaTitle <- paste( c(rep("Missing gene: ",1)), c(50),rep("%",1),sep="")

                                        #    alltheta <- c("01_0","025_0","05_0", "0_015", "0_025", "0_05", "02_02", "025_025", "03_03")

    barEpsilon=0.02

    pdf(paste(wd,"/BICV/Decision3by3.pdf",sep=""))
    par(mfrow=c(1,3))
    par(oma = c(4, 1.5, 1, 0)+0.1, mar=c(2,2,2,0)+0.2)
                                        #       par(mfrow=c(3,3), oma=c(5,4,0,0)+0.1, mar=c(0,0,1,1.5)+0.1)

                                        #    par(mfrow=c(3,3),mar=c(5.1,4.1,1.9,2.1),pin=c(4,4))
    i=7
    xZoom <- list()
    xZoom[[1]] <- seq(1:length(all_n))

    for(t in 1:length(allthetaName)){
        load(paste(wd,"/BICV_", allthetaName[t],"/mse.RData",sep=""))

        x <- xZoom[[1]]           

                                        #            plot(x, mse[i,x,1], ylim=c(ylimTmpNeg,ylimTmp),xlab="N",main=parNames[i],
        decisionTrueMode <- rep(0,4)
        for(xidx in x){
            print(decisionTrue[x[xidx],])
            print(Mode(decisionTrue[xidx,]))
            decisionTrueMode[xidx] <- Mode(decisionTrue[x[xidx],])
        }
        print(decisionTrueMode)
        print(x)
        print(decisionErr)
        print(decisionErrAdjust)
        print(newStat)
        plot(x, decisionErr,xlab="N", ylim = c(0,1), main=allthetaTitle[t],
             ylab="CV error rate",xaxt="n",col=(colors[c(3,1,6)])[decisionTrueMode],pch=decisionTrueMode-1, cex=2)
        axis(1,at=x,label=all_n[x],las=2,cex.axis=0.8)

        points(x, decisionErrAdjust, col=(colors[c(3,1,6)])[decisionTrueMode],pch=(decisionTrueMode+14), cex=1)

                                        #        lines(x, newStat[1,], col=colors()[300], type="l")

        abline(h=0.1, lty=2, col=colors()[330])
                                        #        if(sum(newStat[2,])!=0){
                                        #            lines(x, newStat[2,], col=colors()[320], type="l")
                                        #        }
        if((t-1)%%3==0 ){
            mtext("CV error rate", side=2, line=2.2, cex=0.8)
        }
        if(t%in% c(7,8,9) ){
            mtext("nSample", side=1, line=2.2, cex=0.8)
        }
    }       

    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(10, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottomright", c(methods, "RMSE statistics"), xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n",col=c(colors[c(3,1,6)], colors()[300]),pch=c(0,1,2,20), cex = 1.5)
    
    dev.off()
    print("Done for varying nSample")



    colors<- colors()[c(26,30,28,  505,36,35,  121,131,128, 261)]

    t=1
    load(paste(wd,"/BICV_", allthetaName[t],"/mse.RData",sep=""))

    pdf(paste(wd,"/BICV/RMSEplot.pdf",sep=""))
    par(mfrow=c(2,2), oma=c(2.5,3,1,1)+1, mar=c(1.5,1,1.5,1)+1)

    theMax <- sqrt(max(yHatAll[1,,1], yHatAll_CC[1,,1]))
    for(n in 1:2){
        pType <- rep(4,35)
        colType <- rep(1,35)
        for(i in 1:35){
            if( (decisionCV[n,i] == decisionTrue[n,i]))
                pType[i] <- 1
            if(yHatAll[n,i,1] > yHatAll_CC[n,i,1])
                colType[i] <- 4
        }
        thecolor <- colors[colType+n-1]
        if(n==1){
            plot(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll[n,,1]), xlim=c(15,theMax), ylim=c(15,theMax), pch=pType, col=thecolor, main =allthetaTitle[1])
            lines(0:theMax, 0:theMax, col=colors[10], type="l")
            mtext("RMSE (BI)", side=2, line=2.2, cex=1)
            mtext("RMSE (CC)", side=1, line=2.2, cex=1)
        }else{
            points(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll[n,,1]), pch=pType, col=thecolor )
        }
    }

    t=2
    load(paste(wd,"/BICV_", allthetaName[t],"/mse.RData",sep=""))

    theMax <- sqrt(max(yHatAll[1,,1], yHatAll_CC[1,,1]))
    for(n in 1:2){
        pType <- rep(4,35)
        colType <- rep(1,35)
        for(i in 1:35){
            if( (decisionCV[n,i] == decisionTrue[n,i]))
                pType[i] <- 1
            if(yHatAll[n,i,1] > yHatAll_CC[n,i,1])
                colType[i] <- 4
        }
        thecolor <- colors[colType+n-1]
        if(n==1){
            plot(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll[n,,1]), xlim=c(15,theMax), ylim=c(15,theMax), pch=pType, col=thecolor, main =allthetaTitle[2] )
            lines(0:theMax, 0:theMax, col=colors[10], type="l")
            mtext("RMSE (BI)", side=2, line=2.2, cex=1)
            mtext("RMSE (CC)", side=1, line=2.2, cex=1)
        }else{
            points(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll[n,,1]), pch=pType, col=thecolor )
        }
    }


    t=3
    load(paste(wd,"/BICV_", allthetaName[t],"/mse.RData",sep=""))

    theMax <- sqrt(max(yHatAll[1,,1], yHatAll_CC[1,,1]))
    for(n in 1:2){
        pType <- rep(4,35)
        colType <- rep(1,35)
        for(i in 1:35){
            if( (decisionCV[n,i] == decisionTrue[n,i]) | decisionTrue[n,i]==3)
                pType[i] <- 1
            if(yHatAll[n,i,1] > yHatAll_CC[n,i,1])
                colType[i] <- 4
        }
        thecolor <- colors[colType+n-1]
        if(n==1){
            plot(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll[n,,1]), xlim=c(15,theMax), ylim=c(15,theMax), pch=pType, col=thecolor , main =allthetaTitle[3])
            lines(0:theMax, 0:theMax, col=colors[10], type="l")
            mtext("RMSE (BI)", side=2, line=2.2, cex=1)
            mtext("RMSE (CC)", side=1, line=2.2, cex=1)
        }else{
            points(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll[n,,1]), pch=pType, col=thecolor )
        }
    }

    theMax <- sqrt(max(yHatAll_Gene[1,,1], yHatAll_CC[1,,1]))
    for(n in 1:2){
        ptype <- rep(4,35)
        colType <- rep(4,35)
        for(i in 1:35){
            if( (decisionCV[n,i] == decisionTrue[n,i]) | decisionTrue[n,i]==2 )
                pType[i] <- 1
            if(yHatAll_Gene[n,i,1] < yHatAll_CC[n,i,1])
                colType[i] <- 7
        }
        thecolor <- colors[colType+n-1]
        if(n==1){
            plot(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll_Gene[n,,1]), xlim=c(15,theMax), ylim=c(15,theMax), pch=pType, col=thecolor , main =allthetaTitle[4])
            lines(0:theMax, 0:theMax, col=colors[10], type="l")
            mtext("RMSE (BIgene)", side=2, line=2.2, cex=1)
            mtext("RMSE (CC)", side=1, line=2.2, cex=1)
        }else{
            points(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll_Gene[n,,1]), pch=pType, col=thecolor)
        }
    }

    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(10, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
                                        #    legend("bottomright", c("CC(50)","CC(100)","CC(200)","BI(50)","BI(100)","BI(200)","BIgene(50)","BIgene(100)","BIgene(200)") , xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n",col=colors[c(4:6,1:3,7:9)],pch=c(1,1,1), cex = 1)
    legend("bottomright", c("CC(50)","CC(100)","BI(50)","BI(100)","BIgene(50)","BIgene(100)") , xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n",col=colors[c(4:5,1:2,7:8)],pch=rep(1,6), cex = 1.3)

    dev.off()


    pdf(paste(wd,"/BICV/","n_RMSEplot.pdf",sep=""))

    for(n in 1:4){
        par(mfrow=c(2,2), oma=c(2.5,3,0.2,1)+1, mar=c(2,1,0.5,1)+1)

        t=1
        load(paste(wd,"/BICV_", allthetaName[t],"/mse.RData",sep=""))


        theMax <- sqrt(max(yHatAll[n,,1], yHatAll_CC[n,,1]))
        theMin <- sqrt(min(yHatAll[n,,1], yHatAll_CC[n,,1]))
        pType <- rep(4,35)
        colType <- rep(1,35)
        for(i in 1:35){
            if( (decisionCV[n,i] == decisionTrue[n,i]))
                pType[i] <- 1
            if(yHatAll[n,i,1] > yHatAll_CC[n,i,1])
                colType[i] <- 4
        }
        thecolor <- colors[colType]
        plot(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll[n,,1]), xlim=c(theMin,theMax), ylim=c(theMin,theMax), pch=pType, col=thecolor, main=allthetaTitle[1])
        lines(0:theMax, 0:theMax, col=colors[10], type="l")
        mtext("RMSE (BI)", side=2, line=2.2, cex=1)
        mtext("RMSE (CC)", side=1, line=2.2, cex=1)
        
        t=2
        load(paste(wd,"/BICV_", allthetaName[t],"/mse.RData",sep=""))

        theMax <- sqrt(max(yHatAll[n,,1], yHatAll_CC[n,,1]))
        theMin <- sqrt(min(yHatAll[n,,1], yHatAll_CC[n,,1]))
        pType <- rep(4,35)
        colType <- rep(1,35)
        for(i in 1:35){
            if( (decisionCV[n,i] == decisionTrue[n,i]))
                pType[i] <- 1
            if(yHatAll[n,i,1] > yHatAll_CC[n,i,1])
                colType[i] <- 4
        }
        thecolor <- colors[colType]
        plot(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll[n,,1]), xlim=c(theMin,theMax), ylim=c(theMin,theMax), pch=pType, col=thecolor, main =allthetaTitle[2] )
        lines(0:theMax, 0:theMax, col=colors[10], type="l")
        mtext("RMSE (BI)", side=2, line=2.2, cex=1)
        mtext("RMSE (CC)", side=1, line=2.2, cex=1)


        t=3
        load(paste(wd,"/BICV_", allthetaName[t],"/mse.RData",sep=""))

        theMax <- sqrt(max(yHatAll[n,,1], yHatAll_CC[n,,1]))
        theMin <- sqrt(min(yHatAll[n,,1], yHatAll_CC[n,,1]))
        pType <- rep(4,35)
        colType <- rep(1,35)
        for(i in 1:35){
            if( (decisionCV[n,i] == decisionTrue[n,i]) | decisionTrue[n,i]==3)
                pType[i] <- 1
            if(yHatAll[n,i,1] > yHatAll_CC[n,i,1])
                colType[i] <- 4
        }
        thecolor <- colors[colType]
        plot(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll[n,,1]), xlim=c(theMin,theMax), ylim=c(theMin,theMax), pch=pType, col=thecolor , main =allthetaTitle[3])
        lines(0:theMax, 0:theMax, col=colors[10], type="l")
        mtext("RMSE (BI)", side=2, line=2.2, cex=1)
        mtext("RMSE (CC)", side=1, line=2.2, cex=1)

        theMax <- sqrt(max(yHatAll_Gene[n,,1], yHatAll_CC[n,,1]))
        theMin <- sqrt(min(yHatAll_Gene[n,,1], yHatAll_CC[n,,1]))
        ptype <- rep(4,35)
        colType <- rep(4,35)
        for(i in 1:35){
            if( (decisionCV[n,i] == decisionTrue[n,i]) | decisionTrue[n,i]==2 )
                pType[i] <- 1
            if(yHatAll_Gene[n,i,1] < yHatAll_CC[n,i,1])
                colType[i] <- 7
        }
        thecolor <- colors[colType]
        plot(sqrt(yHatAll_CC[n,,1]), sqrt(yHatAll_Gene[n,,1]), xlim=c(theMin,theMax), ylim=c(theMin,theMax), pch=pType, col=thecolor , main =allthetaTitle[4])
        lines(0:theMax, 0:theMax, col=colors[10], type="l")
        mtext("RMSE (BIgene)", side=2, line=2.2, cex=1)
        mtext("RMSE (CC)", side=1, line=2.2, cex=1)

        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(10, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        legend("bottomright", c("CC","BI","BIgene") , xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n",col=colors[c(3,1,7)],pch=c(1,1,1), cex = 1)

    }
    dev.off()



    colors<- c(colors()[26],  colors()[261], colors()[35], colors()[614],  colors()[621], colors()[121])


    pdf(paste(wd,"/BICV/","vioplot.pdf",sep=""),width=5,height=3)

    par(mfrow=c(3,4), oma=c(2.5,3,0.2,1)+1, mar=c(2,1,0.5,1)+1)

    for(t in 1:2){
        load(paste(wd,"/BICV_", allthetaName[t],"/mse.RData",sep=""))
        for(n in 1:4){
            yHatTrue <- rep(0, 35)
            yHatCCtmp <- sqrt(yHatAll_CC[n,,1])
            yHattmp <- sqrt(yHatAll[n,,1])
            yHatCV <- yHatCCtmp

            for(i in 1:35){
                yHatTrue[i] <- min(yHatCCtmp[i], yHattmp[i])
                if(decisionCV[n,i]==2)
                    yHatCV[i] = yHattmp[i]
            }
            print(paste(mean(yHatTrue), mean(yHatCV), mean(yHatCCtmp)))
            yHat <- data.frame()
            yHat <- as.numeric(c(yHatTrue, yHatCV, yHatCCtmp, yHattmp))
            print(yHat[1:10])
            yHat <- cbind(yHat , c(rep("Truth",35),rep("CV",35),rep("CC",35),rep("BI",35)))
#            allcolor <- c(rep(colors[4],35),rep(colors[5],35),rep(colors[3],35),rep(colors[1],35))
            yHat <-  as.data.frame(yHat)
            colnames(yHat) <- c("RMSE", "methods")
            yHat$RMSE <- as.numeric(as.character(yHat$RMSE))
            yHat$methods <- as.factor(yHat$methods)

#            transform(yHat, RMSE = as.numeric(RMSE))

            print(mode(yHat$methods))
            print(yHat[1:10,])
            dp <- ggplot(yHat, aes(x=methods, y=RMSE, fill=methods)) +
                geom_violin(trim=FALSE)+
                geom_boxplot(width=0.1, fill="white")+
                labs(title=paste(allthetaTitle[t]," (n=",all_n[n],")",sep=""),x="methods", y = "RMSE")
            dp <- dp + scale_fill_brewer(palette="RdBu") + theme_minimal()
            print(dp)
   
#            vioplot(yHatTrue, yHatCV, yHatAll_CC[n,,1], yHatAll[n,,1], names=c("Truth","CV","CC","BI"), col=colors[c(4,5,3,1)])
        }
    }
    t=3
    load(paste(wd,"/BICV_", allthetaName[t],"/mse.RData",sep=""))
        for(n in 1:4){
            yHatTrue <- rep(0, 35)
            yHatCCtmp <- sqrt(yHatAll_CC[n,,1])
            yHattmp <- sqrt(yHatAll[n,,1])
            yHatGenetmp <- sqrt(yHatAll_Gene[n,,1])
            yHatCV <- yHatCCtmp

            for(i in 1:35){
                yHatTrue[i] <- min(yHatCCtmp[i], yHattmp[i], yHatGenetmp[i])
                if(decisionCV[n,i]==2)
                    yHatCV[i] = yHattmp[i]
                if(decisionCV[n,i]==3)
                    yHatCV[i] = yHatGenetmp[i]
            }
            print(paste(mean(yHatTrue), mean(yHatCV), mean(yHatCCtmp)))
            yHat <- data.frame()
            yHat <- as.numeric(c(yHatTrue, yHatCV, yHatCCtmp, yHattmp,yHatGenetmp))
            print(yHat[1:10])
            yHat <- cbind(yHat , c(rep("Truth",35),rep("CV",35),rep("CC",35),rep("BI",35),rep("BIgene",35)))
#            allcolor <- c(rep(colors[4],35),rep(colors[5],35),rep(colors[3],35),rep(colors[1],35))
            yHat <-  as.data.frame(yHat)
            colnames(yHat) <- c("RMSE", "methods")
            yHat$RMSE <- as.numeric(as.character(yHat$RMSE))
            yHat$methods <- as.factor(yHat$methods)

#            transform(yHat, RMSE = as.numeric(RMSE))

            print(mode(yHat$methods))
#            print(yHat[1:10,])
            dp <- ggplot(yHat, aes(x=methods, y=RMSE, fill=methods)) +
                geom_violin(trim=FALSE)+
                geom_boxplot(width=0.1, fill="white")+
                labs(title=paste(allthetaTitle[t]," (n=",all_n[n],")",sep=""),x="methods", y = "RMSE")
            dp <- dp + scale_fill_brewer(palette="RdBu") + theme_minimal()
            print(dp)
   
#            vioplot(yHatTrue, yHatCV, yHatAll_CC[n,,1], yHatAll[n,,1], names=c("Truth","CV","CC","BI"), col=colors[c(4,5,3,1)])
        }
        
#        vioplot(yHatTrue, yHatCV, yHatAll_CC[n,,1], yHatAll[n,,1],yHatAll_Gene[n,,1], names=c("Truth","CV","CC","BI","BIGene"), col=colors[c(4,5,3,1,6)])
    dev.off()
}

