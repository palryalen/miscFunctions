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
    nextEvent <- transitionFunction(currentState,currentTime,
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
