#Finding Hardest possible n sized test set script
FindHardest <- function(Distance,size,number){
  #Distance is the TF distance matrix
  #size is the size of the test set you want to build
  #first minmax normalize the Distance matrix
  #number : number of testsets generated
  testInd <- matrix(nrow = number, ncol = size)
  trainInd <- matrix(nrow = number,ncol=(nrow(Distance)-size))
  dismn <- min(disMat[Distance != 0])
  dismx <- max(Distance)
  disrng <- dismx - dismn
  disMatNorm <- (Distance-dismn)/disrng
  rownames(disMatNorm) <- c(1:nrow(Distance))
  HarmDisHold <- matrix(nrow=number,ncol = size)
  #
  for (num in 1:number){
    NewDistance <- disMatNorm
    
    for (i in 1:size){
      minDis <- numeric(nrow(NewDistance))
      for (j in 1:length(minDis)){
        #getting the (nrow(Distance)-size)th entry of the sorted distance values for each condition
        minDis[j] <- sort(NewDistance[j,],decreasing = T)[(nrow(Distance)-size)]
      }
      #getting the (num)th entry of the sorted minDis as testInd
      testInd[num,i] <- as.numeric(rownames(NewDistance)[sort(minDis,decreasing = T,index.return = T)$ix[num]])
      #indexKeeper[which.max(minDis):nrow(Distance)] = indexKeeper[which.max(minDis):nrow(Distance)] +1
      #removing the chosen one from the distance matrix
      NewDistance <- NewDistance[-c(sort(minDis,decreasing = T,index.return = T)$ix[num]),-c(sort(minDis,decreasing = T,index.return = T)$ix[num])]
      
    }
    trainInd[num,] <- which(! 1:ncol(Distance) %in% testInd[num,])
    dishol <- matrix(nrow = size, ncol = (ncol(Distance) - size))
    
    for (eac in 1: length(testInd[num,] )){
      dishol[eac,1:length(trainInd[num,])] <- disMatNorm[testInd[num,eac],trainInd[num,]]
      HarmDisHold[num,eac] <- harmMean(dishol[eac,])
    }
  }
  results <- list()
  results[[1]] <- testInd
  results[[2]] <- HarmDisHold
  return(results)
  }
##############################
disMat <- rdist(TFexpTran)
fff2 <- FindHardest(disMat,100,531)
