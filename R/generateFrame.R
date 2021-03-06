#' Survival data generator
#'
#' @description Generates event history data
#'
#' @param n Total number of indiivduals
#' @param endTime End of follow-up time
#' @param transitionMatrix from and to states, with transition intensities
#' @param startingState starting state
#'
#' @author Pål Christie Ryalen \email{p.c.ryalen@@medisin.uio.no}
#'
#' @return \code{data.frame} with survival data
#'
#'
#' @examples
#'
#' # generating multistate data with 5 states
#'
#' n <- 100
#'
#'
#' a12 <- 0.8;a13 <- 0.5;a15 <- 1
#' a21 <- 1.2;a24 <- 0.7;a25 <- 1.2
#' a34 <- 0.8;a35 <- 1.4;a43 <- 0.6;a45 <- 1.5
#'
#' endTime <- 10
#' xx <- seq(0,endTime,length.out = 20)
#'
#' spl12 <- smooth.spline(xx,seq(a12,a12,length.out=length(xx)))
#' spl13 <- smooth.spline(xx,seq(a13,a13,length.out=length(xx)))
#' spl15 <- smooth.spline(xx,seq(a15,a15,length.out=length(xx)))
#' spl21 <- smooth.spline(xx,seq(a21,a21,length.out=length(xx)))
#' spl24 <- smooth.spline(xx,seq(a24,a24,length.out=length(xx)))
#' spl25 <- smooth.spline(xx,seq(a25,a25,length.out=length(xx)))
#' spl34 <- smooth.spline(xx,seq(a34,a34,length.out=length(xx)))
#' spl35 <- smooth.spline(xx,seq(a35,a35,length.out=length(xx)))
#' spl43 <- smooth.spline(xx,seq(a43,a43,length.out=length(xx)))
#' spl45 <- smooth.spline(xx,seq(a45,a45,length.out=length(xx)))
#'
#'
#'
#' transitionMatrix <- vector("list",length = 5)
#' transitionMatrix[[1]]$smoothspline[[1]] <- list(NULL)
#' transitionMatrix[[1]]$smoothspline[[2]] <- list(spl12)
#' transitionMatrix[[1]]$smoothspline[[3]] <- list(spl13)
#' transitionMatrix[[1]]$smoothspline[[4]] <- list(NULL)
#' transitionMatrix[[1]]$smoothspline[[5]] <- list(spl15)
#' transitionMatrix[[1]]$neighbours <- c(2,3,5)
#' transitionMatrix[[2]]$smoothspline[[1]] <- list(spl21)
#' transitionMatrix[[2]]$smoothspline[[2]] <- list(NULL)
#' transitionMatrix[[2]]$smoothspline[[3]] <- list(NULL)
#' transitionMatrix[[2]]$smoothspline[[4]] <- list(spl24)
#' transitionMatrix[[2]]$smoothspline[[5]] <- list(spl25)
#' transitionMatrix[[2]]$neighbours <- c(1,4,5)
#' transitionMatrix[[3]]$smoothspline[[1]] <- list(NULL)
#' transitionMatrix[[3]]$smoothspline[[2]] <- list(NULL)
#' transitionMatrix[[3]]$smoothspline[[3]] <- list(NULL)
#' transitionMatrix[[3]]$smoothspline[[4]] <- list(spl34)
#' transitionMatrix[[3]]$smoothspline[[5]] <- list(spl35)
#' transitionMatrix[[3]]$neighbours <- c(4,5)
#' transitionMatrix[[4]]$smoothspline[[1]] <- list(NULL)
#' transitionMatrix[[4]]$smoothspline[[2]] <- list(NULL)
#' transitionMatrix[[4]]$smoothspline[[3]] <- list(spl43)
#' transitionMatrix[[4]]$smoothspline[[4]] <- list(NULL)
#' transitionMatrix[[4]]$smoothspline[[5]] <- list(spl45)
#' transitionMatrix[[4]]$neighbours <- c(3,5)
#' transitionMatrix[[5]]$smoothspline[[1]] <- list(NULL)
#' transitionMatrix[[5]]$smoothspline[[2]] <- list(NULL)
#' transitionMatrix[[5]]$smoothspline[[3]] <- list(NULL)
#' transitionMatrix[[5]]$smoothspline[[4]] <- list(NULL)
#' transitionMatrix[[5]]$smoothspline[[5]] <- list(NULL)
#' transitionMatrix[[5]]$neighbours <- NULL
#'
#'
#' # setting starting state 1
#' dfr <- generateFrame(n,endTime,transitionMatrix,1)
#'
#' # 5 is a terminating state
#' dfr <- dfr[!(dfr$from.state == 5),]
#'
#' head(dfr)
#'
#' @export
#'
#' @importFrom stats predict runif

generateFrame <- function(n, endTime, transitionMatrix,startingState){

  # Number of intervals:
  increments <- min(2e2,100*n)
  Delta <- endTime/increments
  times <- seq(0,endTime,length.out=increments)

  integralMatrix <- getIntegralMatrix(transitionMatrix,times,Delta)

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
