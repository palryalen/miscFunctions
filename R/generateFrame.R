#' Survival data generator
#'
#' @description
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
#'
#'
#' @export

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
