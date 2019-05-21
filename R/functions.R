###########################################################################
####################  FUNCTIONS, MODIFIED 040917      #####################
###########################################################################


# Remedy with built-in martingale increment calculator

DuplicateRows <- function(frame,indices){
        if(is.null(dim(frame)))
                return(frame[rep(1:length(indices), diff(c(0,indices)))])
      return(frame[rep(1:length(indices), diff(c(0,indices))),])
}


boostFrame <- function(fFrame){
      intTimes <- sort(unique(fFrame$to))
      numTimes <- length(intTimes)
      # daInds <- unlist(by(fFrame$to,INDICES = fFrame$id,
      #                     FUN = function(info)unique(c(1,which(intTimes%in%info)))))
      
      daInds <- sapply(unique(fFrame$id),
                        function(ids)unique(c(1,which(intTimes%in%fFrame$to[fFrame$id==ids])) + (1-ids)*numTimes))
      
      
      return(DuplicateRows(fFrame,daInds))
}





refineFr <- function(fFrame){
        intTimes <- sort(unique(fFrame$to))
        
        thaIds <- unique(fFrame$id)
        fTable <- as.data.table(fFrame)
        setkey(fTable,"id","to")
        
        # Which columns to replicate:
        replicateColumns <- which(!(names(fTable) %in%c("from","to","to.state")))
        
        retIdList <- function(idds){
                
                # Finding possible event times and indices for individual idds
                # whichAddInds <- fTable[id==idds,from.state] != toStateVal
                toTimes <- fTable[id==idds,to]
                # intTimesLoc <- sort(unique(c(intTimes[intTimes<=max(toTimes)],toTimes)))
                toInds <- which(intTimes %in% toTimes)
                
                # Extending frame from individual idds
                thisIndTab <- DuplicateRows(fTable[id==idds,replicateColumns,with=F],toInds)
                # locTimes <- intTimes[1:nrow(thisIndTab)]
                
                # Setting the other time-variying columns
                thisIndTab$to <- intTimes
                thisIndTab$from <- c(0,intTimes[-length(intTimes)])
                # thisIndTab$to <- locTimes
                # thisIndTab$from <- c(0,locTimes[-length(locTimes)])
                
                thisIndTab$to.state <- rep(0,nrow(thisIndTab))
                thisIndTab$to.state[toInds] <- fTable[id==idds,to.state]
                
                return(thisIndTab)
        }
        toIndsPos <- lapply(thaIds,retIdList)
        frame <- do.call(rbind,toIndsPos)
        return(frame)
}






# 01.03.17 - attempt on final function for expanding data frame and calculating treatment weights, lanbda_tilde and b given
# endState <- 5
twCalc <- function(fFrame,fFit,lambda_tilde,b,endState){
        wst <- !(fFrame$from.state%in%endState)
        intTimes <- sort(unique(fFrame$to[wst]))
        
        cumB <- fFit$cum[,]
        
        # Extracting cumulative hazard estimates evaluated at the event times
        repVec <- sapply(intTimes,function(t)max(which(cumB[,1] <= t )))
        B <- cumB[repVec,]
        
        # B <- lapply(intTimes,function(t)cumB[max(which(cumB[,1] <=t)),])
        # B <- do.call(rbind,B)
        B[,1] <- intTimes
        
        fFrame$atRisk <- attr(fFit,"weights")
        
        # Refining the factual data frame
        frame <- refineCourse(fFrame,endState)
        # frame <- refineFr(fFrame)
        # frame <- frame[which(frame$from.state!=5)]
        
        idHaz <- function(idds){
                subfr <- frame[id == idds,]
                intinds <- match(subfr$to,intTimes)
                
                haz <- subfr$atRisk * ( diff(c(0,B[intinds,2])) + subfr$L*diff(c(0,B[intinds,3])))
                hazti <- subfr$atRisk * predict(lambda_tilde,intTimes[intinds])$y*diff(c(0,intTimes[intinds]))
                
                dN <- rep(0,nrow(subfr))
                dN[subfr$from.state%in%c(1,3) & subfr$to.state%in%c(2,4)] <- 1
                
                theta_b <- backDiffEst( subfr[which(dN==1),from] ,cumsum(haz),lambda_tilde,1,subfr[,to])
                
                cumPr <- 1 + (theta_b-1)*dN + haz - hazti
                cumPr <- c(1,cumPr[-length(cumPr)])
                
                return(cumprod(cumPr ))
        }
        Lam <- sapply(unique(frame$id),idHaz)
        
        frame$weights <- unlist(Lam)
        
        return(frame)
}


refineCourse <- function(fFrame,endState){
        fFrame <- fFrame[fFrame$from.state!=endState,]
        intTimes <- sort(unique(fFrame$to))

        thaIds <- unique(fFrame$id)
        fTable <- as.data.table(fFrame)
        setkey(fTable,"id","to")
        
        # Which columns to replicate:
        replicateColumns <- which(!(names(fTable) %in%c("from","to","to.state")))
        
        retIdList <- function(idds){
                
                # Finding possible event times and indices for individual idds
                # whichAddInds <- fTable[id==idds,from.state] != toStateVal
                toTimes <- fTable[id==idds,to]
                # intTimesLoc <- sort(unique(c(intTimes[intTimes<=max(toTimes)],toTimes)))
                toInds <- which(intTimes %in% toTimes)
                
                # Extending frame from individual idds
                thisIndTab <- DuplicateRows(fTable[id==idds,replicateColumns,with=F],toInds)
                locTimes <- intTimes[intTimes <= toTimes[length(toTimes)]]
                
                # Setting the other time-variying columns
                # thisIndTab$to <- intTimes
                # thisIndTab$from <- c(0,intTimes[-length(intTimes)])
                thisIndTab$to <- locTimes
                thisIndTab$from <- c(0,locTimes[-length(locTimes)])
        
                thisIndTab$to.state <- rep(0,nrow(thisIndTab))
                thisIndTab$to.state[toInds] <- fTable[id==idds,to.state]
                
                return(thisIndTab)
        }
        toIndsPos <- lapply(thaIds,retIdList)
        frame <- do.call(rbind,toIndsPos)
        return(frame)
}


#factualHaz <- cumsum(haz);cfProc<-lambda_tilde;times <- subfr$to;time <- subfr[which(dN==1),from]
backDiffEst <- function(time,factualHaz,cfProc,b,times){
        inds <- which(time >= times & times > time-b)
        if(length(inds)<=1)
                return(1)
        lambda_tilde_t <- predict(cfProc,time)$y
        DA <- factualHaz[max(inds)] - factualHaz[min(inds)]
        
        theta_b <- ifelse(DA==0,0,lambda_tilde_t/DA) * min(time,b)
        # convex modification for small times
        # if(time<b)
        #         theta_b <- (b-time)/b + (time/b)*theta_b
        return(theta_b)
}


# 24.02.17 - added cf-intensity input
RefTreatWeights2 <- function(fFrame,fFit,lambda_tilde,b){
      
      intTimes <- sort(unique(fFrame$to))
      
      cumB <- fFit$cum[,]
      covs <- all.vars(attr(fFit,"Formula")[[3]])
      AtRisk <- attr(fFit,"weights")

      # Extract estimate of covariates
      CovEst <- cumB[,colnames(cumB)%in%covs]
      
      fFrame$atRisk <- AtRisk
      
      thaIds <- unique(fFrame$id)
      numIds <- length(thaIds)
      numTimes <- length(intTimes)
      fTable <- as.data.table(fFrame)
      setkey(fTable,"id","to")
      
      dN_Indicator <- unique(fFrame$id[fFrame$to.state%in%c(2,4)])
      toIndsB <- which(intTimes %in% cumB[,1])
      if(toIndsB[length(toIndsB)] != length(intTimes))
            toIndsB <- c(toIndsB,length(intTimes))
      BRep <- which(!(colnames(cumB) %in% "time"))
      
      
      # Extending B at intTimes
      BEst <- DuplicateRows(cumB[,BRep],toIndsB)
      
      # dM2 <- fFit$residuals$dM[-1,]
      # dM2 <- apply(dM2,2,cumsum)
      # dMEst <- DuplicateRows(dM2,toIndsB)
      
      
      toIndsId <- rep(thaIds,each=numIds)
      
      replicateColumns <- which(!(names(fTable) %in%c("from","to","to.state")))
      
      retIdList <- function(idds){
            toInds <- which(intTimes %in% fTable[id==idds,to])
            thisIndTab <- DuplicateRows(fTable[id==idds,replicateColumns,with=F],toInds)
            #thisIndTab1 <- fTable[id==idds,replicateColumns,with=F][rep(1:length(toInds),diff(c(0,toInds)))]
            
            toIndsVec <- rep(0,numTimes)
            toIndsVec[toInds] <- fTable[id==idds,to.state]
            thisIndTab$to.state <- toIndsVec
            
            # matProduct <- BEst[,-1] * thisIndTab[,names(thisIndTab) %in% covs,with=F]
            
            # 28.09.16 ~ remedy
            matProduct <- diff(c(0,BEst[,-1])) * thisIndTab[,covs,with=F]*thisIndTab[,atRisk]
            matProduct <- cumsum(as.matrix(matProduct))
            
            if(!is.null(dim(matProduct)))
                  matProduct <- apply(matProduct,1,sum)
            
            # Calculating the martingale increments
            dN <- abs(diff(c(thisIndTab[,atRisk],0))) * (idds %in% dN_Indicator)
            # dM <- dN - (BEst[,1] + matProduct)*diff(c(0,intTimes)) * thisIndTab[,atRisk]
            LambdaHat <- cumsum(thisIndTab[,atRisk]* diff(c(0,BEst[,1]))) + matProduct
            dM <- dN - diff(c(0,LambdaHat)) * thisIndTab[,atRisk]
            
            
            # Estiamtion of weights - quickfix
            
            ##
            # aTiVec <- ifelse(thisIndTab$from.state%in%c(1,3),lati,0)
            # aVec <- ifelse(thisIndTab$from.state%in%c(1),la12,ifelse(thisIndTab$from.state%in%c(3),la34,0))
            # theta <- aTiVec/aVec
            # theta[is.na(theta)] <- 1
            
            # theta <-  lati / ( (1-thisIndTab$L)*la12 + thisIndTab$L*la34)
            
            # Estimation of weights - better fix, given cf-intensity
            
            # covProc <- thisIndTab[,covs,with=F]
            # atRiskProc <- thisIndTab[,atRisk]
            # dB0 <- diff(c(0,BEst[,1]))
            # X_dBj <-  covProc*diff(c(0,BEst[,-1]))
            # if(!is.null(X_dBj))
            #         X_dBj <- apply(X_dBj, 1, sum)
            # dA <- dB0 + X_dBj
            # A <- cumsum(atRiskProc*dA)
            
            DL <- diff(c(0,LambdaHat))
            l_tilde <- predict(lambda_tilde,c(0,intTimes[-length(intTimes)]))$y
            DL_tilde <- l_tilde * diff(c(0,intTimes))
            
            # estimating theta at the jump times
            theta_b <- sapply(intTimes[dN != 0],backDiffEst,factualHaz = LambdaHat,cfProc=lambda_tilde,b=b,times=intTimes)
            
            # Writing out the cumulative product
            Wi <- 1 + DL - DL_tilde
            if(any(dN != 0))
                Wi[dN != 0] <- Wi[dN != 0] + (theta_b - 1)
           
            thisIndTab$treatWeights <- cumprod(Wi)
            # thisIndTab <- thisIndTab[,!(names(thisIndTab) %in% "atRisk" ),with=F]
            return(thisIndTab)
      }
      
      toIndsPos <- lapply(thaIds,retIdList)
      frame <- do.call(rbind,toIndsPos)
      frame$to <- rep(intTimes,numIds)
      frame$from <- rep(c(0,intTimes)[-(numTimes+1)],numIds)
      frame <- subset(frame,select = !(names(frame) %in% c("covariates","atRisk","dM")))
      return(frame)
}




backDiff <- function(vec,times,t,b){
        relTimes <- which(times <= t & times >= t-b)
        if(length(relTimes) == 0)
                return(NA)
        if(min(relTimes)==1)
                Dvec <- vec[max(relTimes)] - vec[1]
        if(min(relTimes)!=1)
                Dvec <- vec[max(relTimes)] - vec[min(relTimes)-1]
        return(Dvec)
}


# toStateVal <- 5
RefTreatWeights3 <- function(fFrame,fFit,toStateVal){
      
      # Extracing the event times at interest, and the total event times
      intTimes <- sort(unique(fFrame$to[fFrame$to.state%in%toStateVal]))
      totIntTimes <- sort(unique(fFrame$to))
      
      # Extracting cum haz estimate, covariate column names
      cumB <- fFit$cum[,]
      covs <- all.vars(attr(fFit,"Formula")[[3]])
      AtRisk <- attr(fFit,"weights")
      
      fFrame$atRisk <- AtRisk
      
      thaIds <- unique(fFrame$id)
      fTable <- as.data.table(fFrame)
      setkey(fTable,"id","to")
      
      
      # Quick fix - is the individual treated?
      dN_Indicator <- unique(fFrame$id[fFrame$to.state%in%c(2,4)])
      
      
      # Finding estimate of cumulative hazard at all hypothetical event times
      toInds <- which(totIntTimes %in% cumB[,1])
      if(toInds[length(toInds)] != length(totIntTimes))
            toInds <- c(toInds,length(totIntTimes))
      
      BRep <- which(!(colnames(cumB) %in% "time"))
      
      BEst <- DuplicateRows(cumB[,BRep],toInds)
      BEst <- cbind(totIntTimes,BEst)
      
      #BEst <- cumB[,BRep][rep(1:length(toInds),diff(c(0,toInds))),]
      
      # Which columns to replicate:
      replicateColumns <- which(!(names(fTable) %in%c("from","to","to.state")))
      
      retIdList <- function(idds){
            
            # Finding possible event times and indices for individual idds
            whichAddInds <- fTable[id==idds,from.state] != toStateVal
            toTimes <- fTable[id==idds,to][whichAddInds]
            intTimesLoc <- sort(unique(c(intTimes[intTimes<=max(toTimes)],toTimes)))
            toInds <- which(intTimesLoc %in% toTimes)
            
            # Extending frame from individual idds
            thisIndTab <- DuplicateRows(fTable[id==idds,replicateColumns,with=F],toInds)
            
            # Extracting the cumulative hazard estimate for idds for its possible event times
            BEstLocInds <- which(totIntTimes %in% intTimesLoc)
            # matProduct <- BEst[BEstLocInds,-1] * thisIndTab[,names(thisIndTab) %in% covs,with=F]
            
            ##
            matProduct <- diff(c(0,BEst[,-1])) * thisIndTab[,covs,with=F]
            matProduct <- cumsum(thisIndTab[,atRisk]*matProduct)
            ##
            if(!is.null(dim(matProduct)))
                  matProduct <- apply(matProduct,1,sum)
            
            # Calculating the martingale increments
            dN <- abs(diff(c(1,thisIndTab[,atRisk]))) * (idds %in% dN_Indicator)
            #dM <- dN - (BEst[BEstLocInds,1] + matProduct)*diff(c(0,intTimesLoc)) * thisIndTab[,atRisk]
            dM <- dN - diff(c(0,BEst[BEstLocInds,1] + matProduct)) * thisIndTab[,atRisk]
            
            dM <- dN - thisIndTab[,atRisk]*(  (1-thisIndTab$L)*la12 + thisIndTab$L*la34  )*diff(c(0,intTimesLoc))
            
            # removing redundant columns
            thisIndTab <- thisIndTab[,!(names(thisIndTab)%in%"atRisk"),with=F]
            
            # Estiamtion of weights - quickfix
            # lati <- predict(cfTreatTransitionMatrix[[1]]$smoothspline[[2]],0)$y
            # la12 <- predict(fTreatTransitionMatrix[[1]]$smoothspline[[2]],0)$y
            # la34 <- predict(fTreatTransitionMatrix[[3]]$smoothspline[[4]],0)$y
            theta <-  lati / ( (1-thisIndTab$L)*la12 + thisIndTab$L*la34)
            thisIndTab$treatWeights <- cumprod(1 + (theta - 1)*dM)
            
            # Setting the other time-variying columns
            thisIndTab$to <- intTimesLoc
            thisIndTab$from <- c(0,intTimesLoc)[-(length(intTimesLoc)+1)]
            
            thisIndTab$to.state <- rep(0,nrow(thisIndTab))
            thisIndTab$to.state[toInds] <- fTable[id==idds,to.state][whichAddInds]
            
            
            return(thisIndTab)
      }
      toIndsPos <- lapply(thaIds,retIdList)
      frame <- do.call(rbind,toIndsPos)
      return(frame)
}


RefTreatWeights4 <- function(fFrame,fFit){
      
      intTimes <- sort(unique(fFrame$to))
      
      cumB <- fFit$cum[-1,]
      covs <- all.vars(attr(fFit,"Formula")[[3]])
      AtRisk <- attr(fFit,"weights")
      
      # Extract estimate of covariates
      CovEst <- cumB[,colnames(cumB)%in%covs]
      
      fFrame$atRisk <- AtRisk
      
      thaIds <- unique(fFrame$id)
      numIds <- length(thaIds)
      numTimes <- length(intTimes)
      fTable <- as.data.table(fFrame)
      setkey(fTable,"id","to")
      
      dN_Indicator <- unique(fFrame$id[fFrame$to.state%in%c(2,4)])
      toInds <- which(intTimes %in% cumB[,1])
      if(toInds[length(toInds)] != length(intTimes))
            toInds <- c(toInds,length(intTimes))
      BRep <- which(!(colnames(cumB) %in% "time"))
      
      
      # Extending B at intTimes
      BEst <- DuplicateRows(cumB[,BRep],toInds)
      
      toIndsId <- rep(thaIds,each=numIds)
      
      replicateColumns <- which(!(names(fTable) %in%c("from","to","to.state")))
      
      daInds <- unlist(by(fTable[,to],INDICES = fTable[,id],FUN = function(info)which(intTimes%in%info)))
      
      daFrame <- DuplicateRows(fTable[,replicateColumns,with=F],dadaInds)
      
      
      
      retIdList <- function(idds){
            toInds <- which(intTimes %in% fTable[id==idds,to])
            thisIndTab <- DuplicateRows(fTable[id==idds,replicateColumns,with=F],toInds)
            #thisIndTab1 <- fTable[id==idds,replicateColumns,with=F][rep(1:length(toInds),diff(c(0,toInds)))]
            
            toIndsVec <- rep(0,numTimes)
            toIndsVec[toInds] <- fTable[id==idds,to.state]
            thisIndTab$to.state <- toIndsVec
            
            matProduct <- BEst[,-1] * thisIndTab[,names(thisIndTab) %in% covs,with=F]
            if(!is.null(dim(matProduct)))
                  matProduct <- apply(matProduct,1,sum)
            
            # Calculating the martingale increments
            dN <- abs(diff(c(1,thisIndTab[,atRisk]))) * (idds %in% dN_Indicator)
            # dM <- dN - (BEst[,1] + matProduct)*diff(c(0,intTimes)) * thisIndTab[,atRisk]
            dM <- dN - diff(c(0,BEst[,1] + matProduct)) * thisIndTab[,atRisk]
            
            
            # Estiamtion of weights - quickfix
            theta <-  lati / ( (1-thisIndTab$L)*la12 + thisIndTab$L*la34)
            thisIndTab$weights <- cumprod(1 + (theta - 1)*dM)
            thisIndTab <- thisIndTab[,!(names(thisIndTab) %in% "atRisk" ),with=F]
            return(thisIndTab)
      }
      
      toIndsPos <- lapply(thaIds,retIdList)
      frame <- do.call(rbind,toIndsPos)
      frame$to <- rep(intTimes,numIds)
      frame$from <- rep(c(0,intTimes)[-(numTimes+1)],numIds)
      frame <- subset(frame,select = !(names(frame) %in% c("covariates","atRisk","dM")))
      return(frame)
}



RefTreatWeights5 <- function(fFrame,fFit){
      
      intTimes <- sort(unique(fFrame$to))
      
      cumB <- fFit$cum
      covs <- all.vars(attr(fFit,"Formula")[[3]])
      AtRisk <- attr(fFit,"weights")
      
      # Extract estimate of covariates
      CovEst <- cumB[,colnames(cumB)%in%covs]
      
      fFrame$atRisk <- AtRisk
      
      thaIds <- unique(fFrame$id)
      numIds <- length(thaIds)
      numTimes <- length(intTimes)
      fTable <- as.data.table(fFrame)
      setkey(fTable,"id","to")
      
      dN_Indicator <- unique(fFrame$id[fFrame$to.state%in%c(2,4)])
      toIndsB <- which(intTimes %in% cumB[,1])
      if(toIndsB[length(toIndsB)] != length(intTimes))
            toIndsB <- c(toIndsB,length(intTimes))
      BRep <- which(!(colnames(cumB) %in% "time"))
      
      
      # Extending B at intTimes
      BEst <- DuplicateRows(cumB[,BRep],toIndsB)
      
      dMRes <- fFit$residuals$dM
      dMRes <- apply(dMRes,2,cumsum)
      dMEst <- DuplicateRows(dMRes,toIndsB)
      dMEst <- apply(rbind(rep(0,ncol(dMEst)),dMEst),2,diff)
      
      
      
      toIndsId <- rep(thaIds,each=numIds)
      
      replicateColumns <- which(!(names(fTable) %in%c("from","to","to.state")))
      
      retIdList <- function(idds){
            toInds <- which(intTimes %in% fTable[id==idds,to])
            thisIndTab <- DuplicateRows(fTable[id==idds,replicateColumns,with=F],toInds)
            #thisIndTab1 <- fTable[id==idds,replicateColumns,with=F][rep(1:length(toInds),diff(c(0,toInds)))]
            
            toIndsVec <- rep(0,numTimes)
            toIndsVec[toInds] <- fTable[id==idds,to.state]
            thisIndTab$to.state <- toIndsVec
            
            matProduct <- BEst[,-1] * thisIndTab[,names(thisIndTab) %in% covs,with=F]
            if(is.null(ncol(BEst[,-1]))){
                  matProduct <- diff(c(BEst[,-1],0)) * thisIndTab[,names(thisIndTab) %in% covs,with=F]
                  matProduct <- cumsum(matProduct)
                  matProduct <- as.matrix(matProduct)
            }
            if(!is.null(dim(BEst[,-1]))){
                  matProduct <- apply(rbind(BEst[,-1],rep(0,ncol(BEst[,-1]))),2,diff) * thisIndTab[,names(thisIndTab) %in% covs,with=F]
                  matProduct <- apply(matProduct,2,cumsum)
                  matProduct <- apply(matProduct,1,sum)
            }
                  
            
            # Calculating the martingale increments
            dN <- abs(diff(c(thisIndTab[,atRisk],0))) * (idds %in% dN_Indicator)
            # dM <- dN - (BEst[,1] + matProduct)*diff(c(0,intTimes)) * thisIndTab[,atRisk]
            dB <- diff(c(BEst[,1] + matProduct,0)) * thisIndTab[,atRisk]
            dM <- dN - dB
            
            colLol <- rep(1:length(toInds),diff(c(0,which(intTimes %in% fFrame$to[fFrame$id==idds]))))
            dM2 <- sapply(1:nrow(thisIndTab),function(ipsy)dMEst[,fTable[,id]==idds][ipsy,colLol[ipsy]])
            # fTable[fTable[,id]==idds,]
            # toInds
            # thisIndTab[12,]
            #dM2 <- apply(dMEst[,fFrame$id==idds],1,sum)
            # dMEst[,fFrame$id==idds]
            cbind(dM,dM2,dN,dB,intTimes,thisIndTab[,c(2,4,5,6),with=F],diff(c(BEst[,1],0)))
            plot(intTimes[intTimes<15],dM[intTimes<15],type="l")
            lines(intTimes[intTimes<15],dM2[intTimes<15],lty=3,col="red")
            
            # Estiamtion of weights - quickfix
            theta <-  lati / ( (1-thisIndTab$L)*la12 + thisIndTab$L*la34)
            thisIndTab$weights <- cumprod(1 + (theta - 1)*dM)
            # thisIndTab <- thisIndTab[,!(names(thisIndTab) %in% "atRisk" ),with=F]
            return(thisIndTab)
      }
      
      toIndsPos <- lapply(thaIds,retIdList)
      frame <- do.call(rbind,toIndsPos)
      frame$to <- rep(intTimes,numIds)
      frame$from <- rep(c(0,intTimes)[-(numTimes+1)],numIds)
      frame <- subset(frame,select = !(names(frame) %in% c("covariates","atRisk","dM")))
      return(frame)
}






RefTreatWeights6 <- function(fFrame,fFit){
      
      intTimes <- sort(unique(fFrame$to))
      
      cumB <- fFit$cum[,]
      covs <- all.vars(attr(fFit,"Formula")[[3]])
      AtRisk <- attr(fFit,"weights")
      
      # Extract estimate of covariates
      CovEst <- cumB[,colnames(cumB)%in%covs]
      
      fFrame$atRisk <- AtRisk
      
      thaIds <- unique(fFrame$id)
      numIds <- length(thaIds)
      numTimes <- length(intTimes)
      fTable <- as.data.table(fFrame)
      setkey(fTable,"id","to")
      
      dN_Indicator <- unique(fFrame$id[fFrame$to.state%in%c(2,4)])
      toIndsB <- which(intTimes %in% cumB[,1])
      if(toIndsB[length(toIndsB)] != length(intTimes))
            toIndsB <- c(toIndsB,length(intTimes))
      BRep <- which(!(colnames(cumB) %in% "time"))
      
      # Extending B at intTimes
      BEst <- DuplicateRows(cumB[,BRep],toIndsB)
      BEstTilde <- DuplicateRows(fFitTilde$cum[,],toIndsB)
      BEstTilde <- BEstTilde[,2]

      toIndsId <- rep(thaIds,each=numIds)
      
      replicateColumns <- which(!(names(fTable) %in%c("from","to","to.state")))
      
      retIdList <- function(idds){
            toInds <- which(intTimes %in% fTable[id==idds,to])
            thisIndTab <- DuplicateRows(fTable[id==idds,replicateColumns,with=F],toInds)
            #thisIndTab1 <- fTable[id==idds,replicateColumns,with=F][rep(1:length(toInds),diff(c(0,toInds)))]
            
            toIndsVec <- rep(0,numTimes)
            toIndsVec[toInds] <- fTable[id==idds,to.state]
            thisIndTab$to.state <- toIndsVec
            
            # matProduct <- BEst[,-1] * thisIndTab[,names(thisIndTab) %in% covs,with=F]
            
            # 28.09.16 ~ remedy
            matProduct <- diff(c(0,BEst[,-1])) * thisIndTab[,names(thisIndTab) %in% covs,with=F]
            matProduct <- cumsum(as.matrix(matProduct))
            
            if(!is.null(dim(matProduct)))
                  matProduct <- apply(matProduct,1,sum)
            
            # Calculating the martingale increments
            dN <- abs(diff(c(thisIndTab[,atRisk],0))) * (idds %in% dN_Indicator)
            # dM <- dN - (BEst[,1] + matProduct)*diff(c(0,intTimes)) * thisIndTab[,atRisk]
            dLambda <- diff(c(0,BEst[,1] + matProduct)) * thisIndTab[,atRisk]
            Lambda <- cumsum(dLambda)
            ##
            # dLambda <- (la12*(thisIndTab$from.state==1) + la34 * (thisIndTab$from.state==3)) * diff(c(0,intTimes))
            # dLambdaTilde <- lati * diff(c(0,intTimes))
            ##
            
            dM <- dN - dLambda
            dLambdaTilde <- lati * diff(c(0,intTimes))*thisIndTab[,atRisk]
            #dLambdaTilde <- diff(c(0,BEstTilde))
            
            # Estiamtion of weights - quickfix
            
            ##
            # aTiVec <- ifelse(thisIndTab$from.state%in%c(1,3),lati,0)
            # aVec <- ifelse(thisIndTab$from.state%in%c(1),la12,ifelse(thisIndTab$from.state%in%c(3),la34,0))
            # thetaExact <- aTiVec/aVec
            # thetaExact[is.na(thetaExact)] <- 1
            
            cumP <- 1 + dLambda - dLambdaTilde
            
            if(any(dN!=0)){
                    #cfProc <- smooth.spline(0:10,rep(lati,11))
                    theta <- backDiffEst(intTimes[dN!=0],Lambda,cfProc,1,intTimes)
                    cumP[dN != 0] <- cumP[dN != 0] + (theta - 1)
            }
            # 
            # theta <- thetaExact
            # # 02.11.16 - backdiff for intensities, adaptive backstep for small times
            # b <- 1
            # La <- as.numeric(unlist(BEst[,1] + thisIndTab[,covs,with=F] * BEst[,2]))
            # dLa <- diff(c(0,La))
            # La <- cumsum(thisIndTab[,atRisk]*dLa)
            # #iftr <- sapply((fFit$cum[,1] >= min(thisIndTab$from[thisIndTab$from.state==3],na.rm=T)),function(lat)ifelse(lat,1,0))
            # #La <- fFit$cum[,2] +  iftr * fFit$cum[,3]
            # DLa <- sapply(intTimes,function(t)(t<b)*backDiff(La,fFit$cum[,1],t,t) + (t>=b)*backDiff(La,fFit$cum[,1],t,b))
            # theta <- aTiVec*c(intTimes[intTimes<b],rep(b,sum(intTimes>=b))) / DLa
            # theta[theta==Inf] <- 1
            # theta[is.na(theta)] <- 1
            # theta[! (thisIndTab$from.state%in%c(1,3))] <- 1
            # 
            # # plot(intTimes,thetaExact,type="l",xlim=c(0,5),ylim=c(0.5,4.5))
            # # lines(intTimes,theta,col="red")
            # 
            # # theta <-  lati / ( (1-thisIndTab$L)*la12 + thisIndTab$L*la34)
            # 
            # ##
            # 
            # countInt <- ifelse(cumsum(dN)==0,0,cumsum(dN) * (theta[dN==1] - 1))
            # dK <- (countInt + dLambda - dLambdaTilde) * thisIndTab$atRisk
            # #thisIndTab$weights <- 
            # # v1 <- cumprod(1 + (thetaExact - 1)*dM)
            # v2 <- cumprod(1 + dK)
            # # plot(intTimes,v1,type="l",xlim=c(0,6),ylim=c(0.5,1.7))
            # # lines(intTimes,v2,col="red")
            # thisIndTab$weights <- v2
            # # thisIndTab <- thisIndTab[,!(names(thisIndTab) %in% "atRisk" ),with=F]
            
            thisIndTab$weights <- cumprod(cumP)
            
            return(thisIndTab)
      }
      
      toIndsPos <- lapply(thaIds,retIdList)
      frame <- do.call(rbind,toIndsPos)
      frame$to <- rep(intTimes,numIds)
      frame$from <- rep(c(0,intTimes)[-(numTimes+1)],numIds)
      frame <- subset(frame,select = !(names(frame) %in% c("covariates","atRisk","dM")))
      return(frame)
}


refineFrame <- function(fFrame,toStateVal){
      # Extracing the event times at interest, and the total event times
      intTimes <- sort(unique(fFrame$to[fFrame$to.state%in%toStateVal]))
      totIntTimes <- sort(unique(fFrame$to))
      
      thaIds <- unique(fFrame$id)
      fTable <- as.data.table(fFrame)
      setkey(fTable,"id","to")

      # Which columns to replicate:
      replicateColumns <- which(!(names(fTable) %in%c("from","to","to.state")))
      
      retIdList <- function(idds){
            
            # Finding possible event times and indices for individual idds
            whichAddInds <- fTable[id==idds,from.state] != toStateVal
            toTimes <- fTable[id==idds,to][whichAddInds]
            intTimesLoc <- sort(unique(c(intTimes[intTimes<=max(toTimes)],toTimes)))
            toInds <- which(intTimesLoc %in% toTimes)
            
            # Extending frame from individual idds
            thisIndTab <- DuplicateRows(fTable[id==idds,replicateColumns,with=F],toInds)

            # Setting the other time-variying columns
            thisIndTab$to <- intTimesLoc
            thisIndTab$from <- c(0,intTimesLoc)[-(length(intTimesLoc)+1)]
            
            thisIndTab$to.state <- rep(0,nrow(thisIndTab))
            thisIndTab$to.state[toInds] <- fTable[id==idds,to.state][whichAddInds]
            
            return(thisIndTab)
      }
      toIndsPos <- lapply(thaIds,retIdList)
      frame <- do.call(rbind,toIndsPos)
      return(frame)
}






toEvents <- function(apEndTime=15,plotIt=F){
      tts <- seq(0,apEndTime,0.5)
      survEst <- sapply(tts,function(tim)
                        sum(fFrame$to[fFrame$from.state%in%c(1,3)&fFrame$to.state!=3]>tim))
      survEst <- survEst/length(unique(fFrame$id))
      if(plotIt){
            plot(tts,survEst,ylim=c(0,1),type = "l",ann=F)
            title("Proportion at risk")
            mtext(expression(t),side=1,line=2)
            mtext(expression(hat(S)[t]),side=2,line=3)
      }
      cutofTime <- min(tts[survEst < 0.05])
      return(cutofTime)
}

pltWts <- function(gFrame,twPlot=T,wPlot=T,WnPlot=F,n=100,ymax = 3){
      g <- gFrame$id
      ug <- unique(g)[1:n]
      tot <- gFrame$to; pTimes <- seq(0,cutofTime,length.out = 40);w <- gFrame$weights;tw <- gFrame$treatWeights
      wn <- gFrame$Wn
      ylm <- c(0,ymax) #added
      cols <- colorRampPalette(c("grey","red"));cll <- cols(10)[3];clol <- cols(10)[3]
      plot(NULL, xlim=c(min(tot),min(max(tot),2990)), 
           ylim=ylm,ann=F)
      for (v in ug) {
            J <- which(v == g)
            if(twPlot)
                  lines(tot[J],tw[J], type='l', col='#00000022')
            if(wPlot)
                  lines(tot[J],w[J], type='l', col=cll)
            if(WnPlot)
                  lines(tot[J],w[J], type='l', col=clol)
      }
      ywt <- sapply(pTimes,function(tim)mean(frame$treatWeights[frame$from<=tim & frame$to >= tim]))
      # ywt_approx <- sapply(pTimes,function(tim)mean(frame$weights[frame$from<=tim & frame$to >= tim]))
      ywt_approx2 <- sapply(pTimes,function(tim)mean(frame$Wn[frame$from<=tim & frame$to >= tim]))
      #ywt_disc <- sapply(pTimes,function(tim)mean(frame$discWeights[frame$from<=tim & frame$to >= tim]))
      lines(pTimes,ywt_approx2,col=clol,lwd=2)
      lines(pTimes,ywt,col='#000022',lwd=2)
      # lines(pTimes,ywt_approx,col=cols(10)[6],lwd=2)
      
}


pltWtsFun <- function(gFrame,aw,ww,onePlot=T,twoPlot=T,n=100,ymax = 3){
      g <- gFrame$id
      ug <- unique(g)[1:n]
      tot <- gFrame$to; pTimes <- seq(0,cutofTime,length.out = 40)
      ylm <- c(0,ymax)
      cols <- colorRampPalette(c("grey","red"));cll <- cols(10)[3]
      plot(NULL, xlim=c(min(tot),min(max(tot),2990)), 
           ylim=ylm,ann=F)
      for (v in ug) {
            J <- which(v == g)
            if(onePlot)
                  lines(tot[J],aw[J], type='l', col='#00000022')
            if(twoPlot)
                  lines(tot[J],ww[J], type='l', col=cll)
      }
      ywt <- sapply(pTimes,function(tim)mean(aw[gFrame$from<=tim & gFrame$to >= tim]))
      ywt_approx2 <- sapply(pTimes,function(tim)mean(ww[gFrame$from<=tim & gFrame$to >= tim]))
      lines(pTimes,ywt_approx2,col=cols(10)[6],lwd=2)
      lines(pTimes,ywt,col='#000022',lwd=2)
      
}



# xx <- data$to[ids]; ssObj <- loessFit
smDeriv <- function(xx,ssObj,epsilon=0.1){
  highOrderDeriv <- function(xVal){
    xGrid <- xVal + epsilon * (-4):4
    cdVec <- c(1/280,	-4/105,1/5,-4/5,0, 4/5, -1/5,	4/105, -1/280)
    yVals <- predict(ssObj,xGrid)
    if(class(ssObj)=="smooth.spline")
      yVals <- yVals$y
    # if(max(is.na(yVals))==1){
    #   naInd <- which(is.na(yVals))
    #   repInd <- ifelse(max(naInd) < 8,max(naInd)+1,min(naInd)-1)
    #   # if(repInd==0)
    #   #   return(highOrderDeriv(xVal - epsilon))
    #   yVals[naInd] <- yVals[repInd]
    # }
    return(yVals %*% cdVec/epsilon)
  }
  
  lowOrderDeriv <- function(xVal){
    xGrid <- xVal + epsilon * c(-1,0)
    yVals <- predict(ssObj,xGrid)
    return(yVals %*% c(-1,1)/epsilon)
  }
  yDeriv <- sapply(xx,highOrderDeriv)
  
  naAction <- function(ind,vec){
    iter <- 1
    while(T){
      subVec <- (ind-iter):(ind+iter)
      # remInd <- ifelse(min(subVec)<=0,1,ifelse(max(subVec)>=length(vec)+1,length(vec),NA))
      subVec <- subVec[which((subVec >=1)*(subVec <= length(vec)) > 0)]
      if(min(is.na(vec[subVec]))==0)
        return(mean(vec[subVec],na.rm=T))
      iter <- iter+1
    }
  }
  yDeriv[is.na(yDeriv)] <- sapply(which(is.na(yDeriv)),naAction,vec=yDeriv)
  return(yDeriv)
}



getWcomb2 <- function(frame){
  
  fullfFrame <- frame
  frame <- frame[order(frame$to),]
  frame <- RefineTimeScale(frame,startTimeName,stopTimeName,frame$to,endStatusName,NULL)
  frame$treatWeights <- 1
  ids <- which((frame$from.state%in%c(1,3)) * (frame$to.state %in% c(2,4)) > 0)
  
  frame$treatWeights[ids] <- sapply(frame$to[ids],function(t)lambdaTilde$coefficients %*% c(1,t,t^2,t^3)/lambda.fit$coefficients %*% c(1,2*t,3*t^2))
  LR <-  cumProdGroup(frame$treatWeights, frame[,idName])
  
  expInt <- sapply(frame$to,W.contComb)
  frame$treatWeights <- LR * expInt
  
  # timeIntervals <- seq(0,endTime,endTime/numIntervals)[-(numIntervals+1)]
  # frs <- sapply(timeIntervals,iptw,regFr=fullfFrame,
  #               startTimeName = startTimeName, stopTimeName=stopTimeName,endStatusName=endStatusName,
  #               endStatus=endStatus)
  # frs <- t(apply(frs,1,cumprod))
  
  # NB!
  # frame <- frame[order(frame$id),]
  
  # frame <- smoothWeights(frame,frs)
  
  # frame$discWeights <- 1
  # for(i in unique(frame$id)){
  #   whichIntervals <- apply( outer( frame[frame$id==i,][,names(frame[frame$id==i,])==stopTimeName], timeIntervals, ">"), 1, sum)
  #   frame$discWeights[frame$id==i] <- frs[i,whichIntervals]
  # }
  
  return(frame)
}



smoothWeights <- function(data,frs,timeIntervals){
  w <- w1 <- w2 <- w3 <- data$treatWeights
  
  for(z in unique(data$id)){
    
    # toTimes <- data$to[data$id==z]
    # smooth.spline
    # if(length(toTimes)>3){
    #   lts <- length(data$to[data$id==z])
    #   numKnots <- round(min(lts,max(lts/10,10)))
    #   optSpline <- smooth.spline(data$to[data$id==z],w[data$id==z],df=4)
    #   w[data$id == z] <- predict(optSpline,data$to[data$id==z])$y
    # }
    
    whichIntervals <- apply( outer( data[data$id==z,][,names(data[data$id==z,])==stopTimeName], timeIntervals, ">"), 1, sum)
    w3[data$id == z] <- frs[z,whichIntervals]
    
  }
  data$smoothWeights <- w
  data$discWeights <- w3
  return(data)
}


plotWeights <- function(w,t,g,n=100){
  ug <- unique(g)[1:n]
  ylm <- c(min(w[frame$id %in% ug]),max(w[frame$id %in% ug])) #added
  plot(NULL, xlim=c(min(t),min(max(t),2990)), 
       ylim=ylm, xlab='time', ylab='weights')
  for (v in ug) {
    J <- which(v == g)
    lines(t[J],w[J], type='s', col='#00000022')
  }
}

plotMean <- function(weights, startTime, stopTime, pTimes){
  n <- length(pTimes)
  y <- rep(0,n)
  for (j in 1:n){
    J <- (startTime <= pTimes[j]) & (stopTime >= pTimes[j])
    y[j] <- mean(weights[J],na.rm=T)
  }
  lines(pTimes,y,col='red')
}

plotWeightsMod <- function(frame,apEndTime,w,t,g,n=100){
  ug <- unique(g)[1:n]
  startTime <- frame$from;g <- frame$id; pTimes <- seq(0,apEndTime,.1)
  ylm <- c(max(min(w[frame$id %in% ug],na.rm=T),0),min(max(w[frame$id %in% ug],na.rm=T),3)) #added
  plot(NULL, xlim=c(min(t),min(max(t),2990)), 
       ylim=ylm,ann=F)
  for (v in ug) {
    J <- which(v == g)
    lines(t[J],w[J], type='s', col='#00000022')
  }
  n <- length(pTimes)
  ywt <- sapply(pTimes,function(tim)mean(frame$treatWeights[frame$from<=tim & frame$to >= tim]))
  ywt_approx <- sapply(pTimes,function(tim)mean(frame$weights[frame$from<=tim & frame$to >= tim]))
  ywt_disc <- sapply(pTimes,function(tim)mean(frame$discWeights[frame$from<=tim & frame$to >= tim]))
  # for (j in 1:n){
  #   J <- (startTime <= pTimes[j]) & (t >= pTimes[j])
  #   y[j] <- mean(w[J],na.rm=T)
  # }
  lines(pTimes,ywt,col='red',lty=2)
  lines(pTimes,ywt_approx,col='red')
  lines(pTimes,ywt_disc,col='green')
}



hazPlot <- function(){
  ymx <- max(c(predict(cfTreatTransitionMatrix[[1]]$smoothspline[[2]][[1]],xx)$y,
               predict(fTreatTransitionMatrix[[1]]$smoothspline[[2]][[1]],xx)$y,
               predict(fTreatTransitionMatrix[[3]]$smoothspline[[4]][[1]],xx)$y))
               
  plot(predict(cfTreatTransitionMatrix[[1]]$smoothspline[[2]][[1]],xx),type="l",ylim=c(0,ymx),ann=F,bty="n",lwd=2)
  lines(predict(fTreatTransitionMatrix[[1]]$smoothspline[[2]][[1]],xx),col="orange")
  lines(predict(fTreatTransitionMatrix[[3]]$smoothspline[[4]][[1]],xx),col="darkorange")
  legend("topleft",c(expression(lambda[12]),expression(lambda[34]),expression(tilde(lambda))),
         col=c("darkorange","orange",1),lwd=c(1,1,2),bty="n")
  title(expression(treatment~intensities))
}


transfer <- function(num, individual, integralMatrix,
                     times, endTime, transitionMatrix,initConditionMatr){
  jumpTimes <- c(0)
  currentTime <- 0
  # bl <- initConditionMatr[num,1]
  currentState <- initConditionMatr[num,1]
  thisTransitionMatrix <- integralMatrix[[currentState]]
  numberOfStates <- length(transitionMatrix[[currentState]]$neighbours)
  jumpPositions <- c(currentState)
  
  while(currentTime < endTime){
    nextEvent <- transitionFunction2(currentState,currentTime,
                                    transitionMatrix,thisTransitionMatrix,
                                    numberOfStates,times,endTime)
    currentState <- nextEvent$nextState
    thisTransitionMatrix <- integralMatrix[[currentState]]
    numberOfStates <- length(transitionMatrix[[currentState]]$neighbours)
    currentTime <- nextEvent$jumpTime
    jumpTimes <- c(jumpTimes,currentTime)
    jumpPositions <- c(jumpPositions,currentState)
  }
  # jumpTimes;jumpPositions
  # jumpTimes <- jumpTimes[-length(jumpTimes)]
  # jumpPositions <- jumpPositions[-length(jumpTimes)]
  
  

  # Commented out 10.06.2018 ~ needed for recurrent events
  # (handling absorbed states)
  # if(length(unique(jumpPositions)) == 1)jumpPositions <- currentState
  len <- length(jumpPositions)
  
  
  Delta <- diff(times)[1]
  # Adding noise to remove ties:
  # jumpTimes[-c(1,len)] <- jumpTimes[-c(1,len)] + runif(length(jumpTimes[-c(1,len)])) * 1/(100*length(times))
  
  jumpTimes[-c(1,len)] <- jumpTimes[-c(1,len)] + runif(length(jumpTimes[-c(1,len)])) * Delta/5
  
  # frame <- data.frame(id = rep(individual,len), from = c(0,jumpTimes[-length(jumpTimes)]),
  #                     to = c(jumpTimes[-1],endTime),from.state = jumpPositions,
  #                     to.state = c(jumpPositions[-1],tail(jumpPositions,1)),
  #                     u = rep(bl,len))
  
  frame <- data.frame(id = rep(individual,len-1), from = jumpTimes[-len],
                      to = jumpTimes[-1],from.state = jumpPositions[-len],
                      to.state = jumpPositions[-1])#,
                      #u = rep(bl,len-1))
  # frame
  return(frame)
}

getIntegralMatrix <- function(transitionMatrix){
      
      mx <- 0
      mx <- max(sapply(1:length(transitionMatrix),function(iis)max(transitionMatrix[[iis]]$neighbours,mx)))
      
      # integralMatrix is a list containing matrices, one for each baseline value.
      # For given baseline value, each row of the matrix contains the cumulative 
      # cause-specific hazards by increments as defined in the variable 'increments'.
      # The last row contains the 'total' cumulative hazard.
      integralMatrix <- vector(length=mx,"list")
      
      for(z in 1:mx){
            
            currentState <- z
            neighbours <- transitionMatrix[[currentState]]$neighbours
            baselineValues <- transitionMatrix[[currentState]]$baseline
            numBaselineValues <- length(baselineValues)
            numberOfStates <- length(neighbours)
            if(numberOfStates > 0){
                  
                  replist <- vector(length=numBaselineValues,"list")
                  
                  for(w in 1:numBaselineValues)replist[[w]] <- matrix(0,numberOfStates+1,maxNumberOfTimesteps)
                  
                  #integralMatrix[[z]] <- list(matrix(0,numberOfStates+1,maxNumberOfTimesteps),
                  #                            matrix(0,numberOfStates+1,maxNumberOfTimesteps))
                  
                  for(k in 1:numBaselineValues){
                        for(i in 1:numberOfStates){
                              yy <- predict(transitionMatrix[[currentState]]$smoothspline[[neighbours[i]]][[k]],times)
                              replist[[k]][i,] <- cumsum(yy$y)
                        }
                  }
                  for(k in 1:numBaselineValues){
                        replist[[k]] <- replist[[k]] * Delta
                        replist[[k]][numberOfStates+1,] = apply(replist[[k]],2,sum)
                  }
                  
                  integralMatrix[[z]] <- replist
                  
            }
            
      }
      return(integralMatrix)
      
}


getIntegralMatrix2 <- function(transitionMatrix,times,Delta){
      
      maxNumberOfTimesteps <- length(times)
      mx <- max(sapply(1:length(transitionMatrix),function(iis)max(transitionMatrix[[iis]]$neighbours,0)))
      
      # integralMatrix is a list containing matrices, one for each baseline value.
      # For given baseline value, each row of the matrix contains the cumulative 
      # cause-specific hazards by increments as defined in the variable 'increments'.
      # The last row contains the 'total' cumulative hazard.
      integralMatrix <- vector(length=mx,"list")
      
      for(z in 1:mx){
            
            currentState <- z
            neighbours <- transitionMatrix[[currentState]]$neighbours
            # baselineValues <- transitionMatrix[[currentState]]$baseline
            # numBaselineValues <- length(baselineValues)
            numberOfStates <- length(neighbours)
            if(numberOfStates > 0){
                  
                  # replist <- vector(length=numBaselineValues,"list")
                  
                  # for(w in 1:numBaselineValues)replist[[w]] <- matrix(0,numberOfStates+1,maxNumberOfTimesteps)
                  replist <- matrix(0,numberOfStates+1,maxNumberOfTimesteps)
                  
                  # for(k in 1:numBaselineValues){
                        for(i in 1:numberOfStates){
                              neighb <- ifelse(length(neighbours)==1,neighbours,neighbours[i])
                              intensty_pred <- transitionMatrix[[currentState]]$smoothspline[[neighb]][[1]]#[[k]]
                              if(class(intensty_pred)!="smooth.spline")
                                      intensty_pred <- transitionMatrix[[currentState]]$smoothspline[[neighb]]
                              yy <- predict(intensty_pred,times)$y
                              # replist[[k]][i,] <- yy
                              replist[i,] <- yy
                        }
                  # }
                  # for(k in 1:numBaselineValues){
                        # replist[[k]][numberOfStates+1,] <- apply(replist[[k]],2,sum) * Delta
                        # replist[[k]][numberOfStates+1,][1] <- 0
                        # replist[[k]][numberOfStates+1,] <- cumsum(replist[[k]][numberOfStates+1,])
                        replist[numberOfStates+1,] <- apply(replist,2,sum) * Delta
                        replist[numberOfStates+1,][1] <- 0
                        replist[numberOfStates+1,] <- cumsum(replist[numberOfStates+1,])
                  # }
                  
                  integralMatrix[[z]] <- replist
                  
            }
            
      }
      return(integralMatrix)
      
}


generateFrame <- function(n, endTime, transitionMatrix,
                          startingState, baselineProbabilities=NULL){
  
  # Number of intervals:
  increments <- min(2e2,100*n)
  Delta <- endTime/increments
  times <- seq(0,endTime,length.out=increments)
  
  integralMatrix <- getIntegralMatrix2(transitionMatrix,times,Delta)
  
  # baselineValues <- sapply(startingState,function(s)transitionMatrix[[s]]$baseline)
  
  # Obtaining n realisations of the process(es).
  # if(length(baselineValues)==1)baselineValues <- c(baselineValues,baselineValues);baselineProbabilities <- c(1,1)/2
  # if(length(startingState)==1)startingState <- c(startingState,startingState)
  # baselineVec <- sample(baselineValues,n,replace=T,prob=baselineProbabilities)
  startingStates <- sample(startingState,n,replace=T)
  # initConditionMatr <- matrix(c(baselineVec,startingStates),ncol=2)
  initConditionMatr <- as.matrix(startingStates,ncol=1)
  
  obj <- lapply(1:n,transfer,individual=0,
                integralMatrix=integralMatrix,times=times,
                endTime=endTime,transitionMatrix=transitionMatrix,
                initConditionMatr=initConditionMatr)
  frame <- do.call("rbind",obj)
  
  # # Dealing with ties
  # duplicates <- duplicated(frame$to)
  # duplicates[frame$to == endTime] <- F
  # numDuplicates <- sum(duplicates)
  # if(numDuplicates > 0){
  #   dupNoise <- 1e-4*runif(numDuplicates)
  #   frame$to[duplicates] <- frame$to[duplicates] + dupNoise
  #   frame$from[duplicates] <- frame$from[duplicates] + dupNoise
  # }
  
  # Adding id
  # counter <- 1
  # frame$id[1] <- 1
  # for(i in 2:(dim(frame)[1])){
  #   if(frame$to[i-1] >= frame$to[i])
  #     counter <- counter + 1
  #   frame$id[i] <- counter
  # }
  
  
  # frame$id <- 1 + cumsum(1*(diff(c(0,frame$to))<0))
  
  frame$id <- ifelse(frame$from == 0,1,0)
  frame$id <- cumsum(frame$id)
  return(frame)
}

# This function samples jump times and states.
transitionFunction <- function(baseline, currentState, currentTime,
                               transitionMatrix, thisTransitionMatrix, 
                               numberOfStates, times, endTime){
  # Handling cases where there are no neighbours or no baseline is provided.
  neighbours <- transitionMatrix[[currentState]]$neighbours
  if(is.null(neighbours)){
    return(list(jumpTime=endTime, nextState=currentState))
  }
  if(!(baseline %in% transitionMatrix[[currentState]]$baseline)){
    stop('invalid baseline')
  }
  
  thisTransitionMatrix <- thisTransitionMatrix[[baseline]]
  # Drawing the jump time using the inversion method. If the jump time is larger than
  # the study time (endTime), we indicate that the data is right censored by just
  # returning the current state (currentState) and endTime. Otherwise, we randomly
  # draw the next state based ratios of the cumulative hazards.
  U <- runif(1)
  valids <- which(thisTransitionMatrix[numberOfStates+1,] + log(U) >= 0)
  if(length(valids) == 0){
    return(list(jumpTime=endTime, nextState=currentState))
  }
  jumpTimePos <- min(valids)
  transitionProbs <- rep(0,length = numberOfStates)
  for(i in 1:numberOfStates){
    transitionProbs[i] <- thisTransitionMatrix[i,jumpTimePos]/
      thisTransitionMatrix[numberOfStates+1,jumpTimePos]
  }
  
  V <- runif(1)
  nextState <- neighbours[min(which(V < cumsum(transitionProbs)))]
  jumpTime <- currentTime + times[jumpTimePos]
  if(jumpTime==currentTime)
    jumpTime <- jumpTime + 1e-2*runif(1)
  return(list(jumpTime=jumpTime, nextState=nextState))
}


# This function samples jump times and states.


# transitionFunction2 <- function(baseline, currentState, currentTime,
#                                 transitionMatrix, thisTransitionMatrix, 
#                                 numberOfStates, times, endTime){
transitionFunction2 <- function(currentState, currentTime,
                                transitionMatrix, thisTransitionMatrix, 
                                numberOfStates, times, endTime){
      neighbours <- transitionMatrix[[currentState]]$neighbours
      
      if(is.null(neighbours)){
            return(list(jumpTime=endTime, nextState=currentState))
      }
      # if(!(baseline %in% transitionMatrix[[currentState]]$baseline)){
      #       stop('invalid baseline')
      # }
      
      
      # thisTransitionMatrix <- thisTransitionMatrix[[1]]
      
      # Drawing the jump time using the inversion method. If the jump time is larger than
      # the study time (endTime), we indicate that the data is right censored by just
      # returning the current state (currentState) and endTime. Otherwise, we randomly
      # draw the next state based ratios of the cumulative hazards.
      
      U <- runif(1)
      
      startTimePos <- min(which(currentTime <= times))
      adjustSum <- thisTransitionMatrix[numberOfStates+1,startTimePos]
      
      thisTransitionMatrix[numberOfStates+1,] <- thisTransitionMatrix[numberOfStates+1,] - adjustSum
      
      
      valids <- which(thisTransitionMatrix[numberOfStates+1,] + log(U) >= 0)
      if(length(valids) == 0){
            return(list(jumpTime=endTime, nextState=currentState))
      }
      jumpTimePos <- min(valids)
      # jumpTimePos <- jumpTimePos + startTimePos - 1
      
      transitionProbs <- thisTransitionMatrix[,jumpTimePos]
      transitionProbs <- transitionProbs[-length(transitionProbs)]
      transitionProbs <- transitionProbs/sum(transitionProbs)
      
      V <- runif(1)
      
      nextState <- neighbours[min(which(V < cumsum(transitionProbs)))]
      jumpTime <- times[jumpTimePos]
      return(list(jumpTime=jumpTime, nextState=nextState))
}




# This function samples jump times and states.
transitionFunction3 <- function(baseline, currentState, currentTime,
                               transitionMatrix, thisTransitionMatrix, 
                               numberOfStates, times, endTime){
      # Handling cases where there are no neighbours or no baseline is provided.
      neighbours <- transitionMatrix[[currentState]]$neighbours
      if(is.null(neighbours)){
            return(list(jumpTime=endTime, nextState=currentState))
      }
      if(!(baseline %in% transitionMatrix[[currentState]]$baseline)){
            stop('invalid baseline')
      }
      
      thisTransitionMatrix <- thisTransitionMatrix[[baseline]]
      # Drawing the jump time using the inversion method. If the jump time is larger than
      # the study time (endTime), we indicate that the data is right censored by just
      # returning the current state (currentState) and endTime. Otherwise, we randomly
      # draw the next state based ratios of the cumulative hazards.
      
      jumpTime <- currentTime
      while(jumpTime <= currentTime){
            U <- runif(1)
            
            valids <- which(thisTransitionMatrix[numberOfStates+1,] + log(U) >= 0)
            if(length(valids) == 0){
                  return(list(jumpTime=endTime, nextState=currentState))
            }
            jumpTimePos <- min(valids)
            jumpTime <- times[jumpTimePos]

      }
      
      
      # startTimePos <- min(which(currentTime <= times))
      # adjustSum <- thisTransitionMatrix[numberOfStates+1,startTimePos]
      # 
      # thisTransitionMatrix[numberOfStates+1,] <- thisTransitionMatrix[numberOfStates+1,] - adjustSum

      transitionProbs <- thisTransitionMatrix[,jumpTimePos]
      transitionProbs <- transitionProbs[-length(transitionProbs)]
      transitionProbs <- transitionProbs/sum(transitionProbs)
      
      
      V <- runif(1)
      nextState <- neighbours[min(which(V < cumsum(transitionProbs)))]
      # jumpTime <- currentTime + times[jumpTimePos]
      # clumsy tie-dealing
      # if(jumpTime==currentTime)
      #       jumpTime <- jumpTime + 1e-2*runif(1)
      return(list(jumpTime=jumpTime, nextState=nextState))
}


# # The celebrated MakeCensWeights-function.
# MakeCensWeights2 <- function(fFit, cfFit, data, startTimeName, stopTimeName,
#                              endStatusName, intStatusName, idName){
#   
#   # Extracting and merging the residuals for the f and cf model.
#   fRes <- aalenRes(fFit)
#   cfRes <- aalenRes(cfFit)
#   names(fRes)[1] <- 'fRes.dM'
#   names(cfRes)[1] <- 'cfRes.dM'
#   B <- merge(fRes, cfRes, by=c('id','time'), all=T, sort=T)
#   
#   colnames(B)[which(names(B) == 'time')] <- stopTimeName
#   colnames(B)[which(names(B) == 'id')] <- paste('temp',idName,sep='.')
#   data[paste('temp',idName,sep='.')] <-  attr(fFit,'id')
#   
#   # Test parameters for insertion into RefineTimeScale:
#   # A <- data; intTimes <- sort(unique(B[,stopTimeName]))
#   
#   data <- RefineTimeScale(data, startTimeName, stopTimeName,
#                           sort(unique(B[,stopTimeName])), endStatusName, 
#                           intStatusName)
#   
#   data <- merge(data, B, by=c(paste('temp', idName, sep='.'), stopTimeName),
#                 all.x=T, sort=T)
#   data$fRes.dM[is.na(data$fRes.dM)] <- 0
#   data$cfRes.dM[is.na(data$cfRes.dM)] <- 0
#   LR <-  cumProdGroup(1 + data$cfRes.dM - data$fRes.dM, data[,idName])
#   data$weights <- rightShiftGroup(LR, data[,idName],
#                                   initVal=rep(1,ncol(fFit$residuals$dM)))
#   
#   data <- data[-which(names(data) %in% c(paste('temp', idName, sep='.'),
#                                          'fRes.dM', 'cfRes.dM'))]
#   return(data)
# }


cf.est <- function(frame){
  stateind <- which((frame$from.state = 1)*(frame$to.state == 2) > 0)
  ord <- order(frame$to[stateind])
  stateind <- stateind[ord]
  tim <- frame$to[stateind]
  n <- length(stateind)
  alpha <- rep(0,n)
  w <- frame$weights
  if(is.null(w))w <- rep(1,dim(frame)[1])
  #getEst <- function(Ti){
  #  eventInd <- which((frame$from < Ti)*(frame$to >= Ti) > 0)
  #  beta_Ti <- sum(w[eventInd])
  #  pos <- stateind[frame$to[stateind]==Ti]
  #  alpha_i <- w[pos]/beta_Ti
  #  return(alpha_i)
  #}
  #alpha <- sapply(tim,getEst)
  for(i in 1:n){
    Ti <- tim[i]
    eventInd <- which((frame$from < Ti)*(frame$to >= Ti) > 0)
    beta_Ti <- sum(w[eventInd])
    if(beta_Ti != 0)
      alpha[i] <- w[stateind[i]]/beta_Ti
  }
  return(data.frame(t=tim,A = cumsum(alpha)))
}

  myConstWeights <- function(frame,rates){
  times <- sort(frame$to)
  frame <- frame[order(frame$to),]
  ind <- which(frame$to.state == 3)
  
  timeFun <- function(t){
    index <- times >=t
    W <- exp(t*rates[frame$u[index],1])
    dA <- exp(t*rates[frame$u[which(frame$to==t)],1])/sum(W)
    return(dA)
  }
  dA <- sapply(times[ind],timeFun)
  frame$A <- rep(0,dim(frame)[1])
  frame$A[ind] <- dA
  frame$A <- cumsum(frame$A)
  frame$weights <- exp(frame$to * rates[frame$u,1])
  return(frame)
}


playFun <- function(rates,num,tMax){
  u <- sample(1:2,num,replace=T,prob=c(0.5,0.5))
  n1 <- sum(u==1);n2 <- num-n1
  times1 <- rexp(n1,sum(rates[1,]))
  times2 <- rexp(n2,sum(rates[2,]))
  
  times <- to.state <- rep(0,num)
  times[u==1] <- times1
  times[u==2] <- times2
  switchInd <- tMax < times
  times[switchInd] <- tMax
  
  if(sum(rates[,1])==0){
    to.state <- rep(3,num)
    to.state[switchInd] <- 1
    frame <- data.frame(id = 1:num, from = rep(0,num),to = times,
                        from.state = rep(1,num),to.state = to.state,u = u)
    return(frame)
  }
  
  to.state1 <- sample(2:3,n1,prob = rates[1,]/sum(rates[1,]),replace=T)
  to.state2 <- sample(2:3,n2,prob = rates[2,]/sum(rates[2,]),replace=T)
  
  to.state[u==1] <- to.state1
  to.state[u==2] <- to.state2
  to.state[switchInd] <- 1
  
  frame <- data.frame(id = 1:num,from = rep(0,num),to = times,
                      from.state = rep(1,num),to.state = to.state,u = u)
  return(frame)
}

storeSet <- function(num, smoothSplines, baselineProbabilities){
  source("transitionMatrix.R")
  source("functions.R")
  endTime <- 10
  startingState <- 1
  fTransitionMatrix[[1]]$smoothspline <- smoothSplines
  cfTransitionMatrix[[1]]$smoothspline[[3]] <- smoothSplines[[3]]
  
  fFrame <- generateFrame(num, endTime, fTransitionMatrix,
                          startingState, baselineProbabilities)
  cfFrame <- generateFrame(num, endTime, cfTransitionMatrix,
                           startingState, baselineProbabilities)
  fFrame <- fFrame[fFrame$from.state == 1,]
  cfFrame <- cfFrame[cfFrame$from.state == 1,]
  
  fFit <- aalen(Surv(from,to,to.state==2) ~ 1 + u, data = fFrame, residuals = T,n.sim=0)
  cfFit <- aalen(Surv(from,to,to.state==2) ~ 1, data = fFrame, residuals = T,n.sim=0)
  
  fFrameNew <- MakeCensWeights2(fFit,cfFit,fFrame,"from","to","to.state",
                                NULL,"id")
  
  weightedNA <- survfit(Surv(from,to,to.state==3)~1,
                        data=fFrameNew,weights=fFrameNew$weights)
  weightedNA <- data.frame(weightedNA$time,
                           cumsum(c(0,weightedNA$n.event[-1]/weightedNA$n.risk[-1])))
  
  fNA <- survfit(Surv(from,to,to.state==3)~1,data=fFrame)
  fNA <- data.frame(fNA$time,
                    cumsum(c(0,fNA$n.event[-1]/fNA$n.risk[-1])))
  cfNA <- survfit(Surv(from,to,to.state==3)~1,data=cfFrame)
  cfNA <- data.frame(cfNA$time,
                     cumsum(c(0,cfNA$n.event[-1]/cfNA$n.risk[-1])))
  
  names(weightedNA)[1] <- names(fNA)[1] <- names(cfNA)[1] <-"t"
  names(weightedNA)[2] <- "weighted factual"
  names(fNA)[2] <- "factual"
  names(cfNA)[2] <- "counterfactual"
  
  #fr <- merge(weightedNA,fNA,by.x="t",all=T)
  fFrameNew <- 0
  storeframes <- list(fFrame,cfFrame,fFrameNew,fNA,cfNA,weightedNA)
  #}
  
  return(storeframes)
}

plotExamples <- function(results){
  #load("storeframes.RData")
  n <- length(results)
  tMax <- 10
  for(i in 1:n){
    frames <- results[[i]]
    fNA <- frames[[4]];cfNA <- frames[[5]];weightedNA <- frames[[6]]
    
    ymax <- max(tail(fNA[,2],1),tail(cfNA[,2],1),tail(weightedNA[,2],1))
    
    name <- paste("case",i,".pdf",sep="")
    pdf(file=name)
    plot(cfNA,type="s",xlim=c(0,tMax),
         ylim=c(0,ymax),bty="n",ann=F,col="blue")
    par(new=T)
    plot(weightedNA,type="s",ann=F,bty="n",
         xlim=c(0,tMax),ylim=c(0,ymax),col="red")
    par(new=T)
    plot(fNA,type="s",ann=F,
         xlim=c(0,tMax),ylim=c(0,ymax))
    legend("topleft",c("counterfactual","weighted factual","factual"),
           col=c("blue","red",1),lty=1,bty="n")
    dev.off()
  }
}

comparePlot <- function(num, intervalLengths,n.sim=1){
  #source("smoothsplineslist.R")
  source(paste(getwd(),"/R/smoothsplineslist.R",sep=""))
  startingState <- 1;endTime<-10;baselineProbabilities<-c(.5,.5)
  numCases <- length(smoothSplinesList)
  for(k in 1:numCases){
    smoothSplines <- smoothSplinesList[[k]]
    fTransitionMatrix[[1]]$smoothspline <- smoothSplines
    cfTransitionMatrix[[1]]$smoothspline[[3]] <- smoothSplines[[3]]
    endTime <- 10
    baselineProbabilities <- c(.5,.5)
    print("Generating data.frames...")
    fFrame <- generateFrame(num, endTime, fTransitionMatrix,
                            startingState, baselineProbabilities)
    fFrame <- fFrame[fFrame$from.state == 1,]; fFrame <- fFrame[order(fFrame$to),]
    
    cfFrame <- generateFrame(num, endTime, cfTransitionMatrix,
                             startingState, baselineProbabilities)
    cfFrame <- cfFrame[cfFrame$from.state == 1,]
    
    print("Calculating continuous-time weights and estimates...")
    fFit <- aalen(Surv(from,to,to.state==2) ~ 1 + u, data = fFrame, residuals = T,n.sim=0)
    cfFit <- aalen(Surv(from,to,to.state==2) ~ 1, data = fFrame, residuals = T,n.sim=0)
    
    # naObj <- naEstimate(fFit,cfFit,fFrame,"id",Surv(from,to,to.state==3)~1,n.sim)
    naObj <- naEstimate(fFit,cfFit,fFrame,"id",Surv(from,to,to.state==3)~1,4,n.sim,"continuous")
    
    # fFrameNew <- MakeCensWeights(fFit,cfFit,fFrame,"from","to","to.state",
    #                              NULL,"id")
    
    # survObj <- survfit(Surv(from,to,to.state==3)~1,data=fFrameNew,weights=fFrameNew$weights)
    # survObj <- data.frame(survObj$time,cumsum(c(0,survObj$n.event[-1]/survObj$n.risk[-1])))
    na.three <- survfit(Surv(from,to,to.state==3)~1,data=fFrame)
    na.three <- data.frame(na.three$time,cumsum(c(0,na.three$n.event[-1]/na.three$n.risk[-1])))
    compareFit <- survfit(Surv(from,to,to.state==3)~1,data=cfFrame)
    compareFit <- data.frame(compareFit$time,cumsum(c(0,compareFit$n.event[-1]/compareFit$n.risk[-1])))
    
    #print("Looping over the discretisations, finding discrete weights.")
    for(j in 1:length(intervalLengths)){
      timeIntervals <- seq(0,endTime,length.out=intervalLengths[j]+1)
      timeIntervals <- timeIntervals[-length(timeIntervals)]
      len <- length(timeIntervals)
      delta <- endTime/len
      weightMatrix <- sapply(timeIntervals, getProbs, frame=fFrame,
                             len=len, endTime=endTime)
      weightMatrix <- t(apply(weightMatrix,1,cumprod))
      
      ## added
      numIntervals<-intervalLengths[j]; regressionFrame<-fFrame; startTimeName<-"from"
      stopTimeName<-"to"; endStatusName<-"to.state"; endStatus<-3; freqTabl<-table(1:dim(fFrame)[1])
      bootTransitionTimes <- stopTimeTransitions <- unique(fFrame$to[fFrame$to.state==3])
      disc <- naDiscreteBoot(numIntervals, endTime, regressionFrame, startTimeName,
                             stopTimeName, endStatusName, endStatus, freqTabl,
                             bootTransitionTimes, stopTimeTransitions)
      ##
      
      # dA <- rep(0,dim(fFrame)[1])
      # for(i in 1:dim(fFrame)[1]){
      #   if(fFrame$to.state[i]==3){
      #     Ti <- fFrame$to[i]
      #     intervalPos <- which((Ti < c(timeIntervals[-1],endTime))*(Ti >=timeIntervals) > 0)
      #     W <- sum(weightMatrix[fFrame$to >= Ti,intervalPos])
      #     dA[i] <- weightMatrix[i,intervalPos]/W
      #   }
      # }
      # A <- cumsum(dA)
      # unique(A);as.vector(c(0,disc))
      
      ymax <- max(tail(compareFit[,2],1), tail(na.three[,2],1),A)
      
      #name <- paste("comparison-",k,"-",intervalLengths[j],".pdf",sep="")
      #pdf(file=name)
      # plot(sort(fFrame$to),A,type="s",xlim=c(0,endTime),
      #      ylim=c(0,ymax),col="green",ann=F,bty="n")
      plot(as.numeric(names(disc)),disc,type="s",xlim=c(0,endTime),
                 ylim=c(0,ymax),col="green",ann=F,bty="n")
      #par(new=T)
      plot(compareFit,type="s",xlim=c(0,tMax),
           ylim=c(0,ymax),bty="n",ann=F,col="blue")
      par(new=T)
      #plot(survObj,type="s",
      #     xlim=c(0,tMax),ylim=c(0,ymax),col="red",bty="n",ann=F)
      #par(new=T)
      plot(na.three,type="s",ann=F,
           xlim=c(0,tMax),ylim=c(0,ymax))
      par(new=T)
      plot(naObj$time,naObj$contCurve,type="s",ann=F,col="darkred",
           xlim=c(0,tMax),ylim=c(0,ymax))
      par(new=T)
      plot(naObj$time,naObj$contCurveCI[1,],type="s",
           xlim=c(0,tMax), ylim=c(0,ymax),ann=F,bty="n",col="darkred",lty=2)
      par(new=T)
      plot(naObj$time,naObj$contCurveCI[2,],type="s",
           xlim=c(0,tMax), ylim=c(0,ymax),ann=F,bty="n",col="darkred",lty=2)
      legend("topleft",c(paste("IPCW, timepoints:",intervalLengths[j]),"counterfactual","factual","weighted factual","weighted factual 95 percent CI"),
             col=c("green","blue",1,"darkred","darkred"),lty=c(1,1,1,1,2),bty="n")
      #dev.off()
    }
  }
  
}

# IPCW
getProbs <- function(time, frame, len, endTime){
  v <- rep(1,dim(frame)[1])
  if(sum(time<=frame$to)==0)return(v)
  
  delta <- endTime/(len-1)
  
  geqind <- which(time<=frame$to)
  regressionFrame <- frame[geqind,]
  leind <- which(delta + time > regressionFrame$to)
  
  regressionFrame$nc <- rep(0,dim(regressionFrame)[1])
  regressionFrame$nc[leind] <- 1*(regressionFrame$to.state[leind] %in% c(1,2))
  
  marginalMod <- glm(formula = nc ~ 1, data = regressionFrame, family = "binomial")
  conditionalMod <- glm(formula = nc ~ 1 + u, data = regressionFrame, family = "binomial")
  
  marginalProbabilities <- 1 - predict(marginalMod,type="response")
  conditionalProbabilities <- 1 - predict(conditionalMod,type="response")
  
  v[geqind] <- marginalProbabilities/conditionalProbabilities
  return(v)
}


lastStore <- function(){
  num <- 2e3
  n.sim <- 1e3
  intervals <- c(4,8,16,32)
  
  startingState <- 1
  endTime <- 10
  baselineProbabilities<-c(.5,.5)
  numCases <- length(smoothSplinesList)
  modelFormula <- Surv(from,to,to.state==3)~1
  mkList <- vector(mode = "list",length = numCases)
  for(k in 1:numCases){
    smoothSplines <- smoothSplinesList[[k]]
    fTransitionMatrix[[1]]$smoothspline <- smoothSplines
    cfTransitionMatrix[[1]]$smoothspline[[3]] <- smoothSplines[[3]]
    print("Generating data.frames...")
    fFrame <- generateFrame(num, endTime, fTransitionMatrix,
                            startingState, baselineProbabilities)
    fFrame <- fFrame[fFrame$from.state == 1,]
    # NB! Should not depend upon ordering of fFrame!!
    #fFrame <- fFrame[order(fFrame$to),]
    
    cfFrame <- generateFrame(num, endTime, cfTransitionMatrix,
                             startingState, baselineProbabilities)
    cfFrame <- cfFrame[cfFrame$from.state == 1,]
    
    print("Calculating continuous-time weights and estimates...")
    fFit <- aalen(Surv(from,to,to.state==2) ~ 1 + u, data = fFrame, residuals = T,n.sim=0)
    cfFit <- aalen(Surv(from,to,to.state==2) ~ 1, data = fFrame, residuals = T,n.sim=0)
    
    data <- fFrame
    
    naObj <- naEstimate(fFit,cfFit,data,"id",modelFormula,intervals,n.sim)
    
    na.three <- survfit(Surv(from,to,to.state==3)~1,data=fFrame)
    na.three <- data.frame(na.three$time,cumsum(c(0,na.three$n.event[-1]/na.three$n.risk[-1])))
    compareFit <- survfit(Surv(from,to,to.state==3)~1,data=cfFrame)
    compareFit <- data.frame(compareFit$time,cumsum(c(0,compareFit$n.event[-1]/compareFit$n.risk[-1])))
    
    naObj$factual <- na.three
    naObj$counterfactual <- compareFit
    
    mkList[[k]] <- naObj
  }
  save(mkList, file = "mkList.RData")
}
