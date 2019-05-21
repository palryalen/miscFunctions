getIntegralMatrix <- function(transitionMatrix,times,Delta){

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
