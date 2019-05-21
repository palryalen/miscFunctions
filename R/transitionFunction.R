transitionFunction <- function(currentState, currentTime,
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
