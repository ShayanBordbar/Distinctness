#Finding Hardest possible n sized test set script
FindHardest <- function(Distance,size,number){
  #Distance is the TF distance matrix
  #size is the size of the test set you want to build
  #first minmax normalize the Distance matrix
  #number : number of testsets produced
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
par(mfrow = c(1,1),mar = c(4,5,4,4))
#boxplot.matrix(fff[[2]],use.cols = F,ylim = c(0.2,1))
NewDifGrad <- matrix(nrow = nrow(fff2[[2]]), ncol = ncol(fff2[[2]]))
for (i in 1:nrow(fff2[[2]])){
  NewDifGrad[(nrow(fff2[[2]])-i +1),] <- fff2[[2]][i,]
}
boxplot.matrix(NewDifGrad,use.cols = F,ylim = c(0,1))
NewDifGradSam <- matrix(nrow = 50, ncol = 100)
for(i in 1:50){
  NewDifGradSam[i,] <- NewDifGrad[(10*i)+31,]
}
boxplot.matrix(NewDifGradSam,use.cols = F,ylim = c(0,1),ylab = "hardness")
abline(h = 0,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 0.1,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 0.2,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 0.3,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 0.4,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 0.5,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 0.6,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 0.7,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 0.8,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 0.9,col = "lightgray", lty = "dotted",lwd = par("lwd"))
abline(h = 1,col = "lightgray", lty = "dotted",lwd = par("lwd"))
# grid(nx = 0, ny = 10, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
boxplot.matrix(AstrixTestInput2HarmNotAvg,ylim = c(0,1))
abline(h = 0.2,col=2)
abline(h = 0.3,col=2)
abline(h = 0.4,col=2)
abline(h = 0.5,col=2)
abline(h = 0.6,col=2)
abline(h = 0.7,col=2)
abline(h = 0.8,col=2)

 minDis <- numeric(631)
for (i in 1:631){
     minDis[i] <- sort(disMatNorm[i,],decreasing = T)[501]
   }
 hist(minDis)
 minDisSorted <- sort(minDis, decreasing = T, index.return = T)$ix
 plot(1:631,minDis[minDisSorted] )
 testInd <- minDisSorted[532:631]
 trainInd <- minDisSorted[1:531]
 size = 100
 dishol <- matrix(nrow = size, ncol = (631 - size))
 HarmDisHold <- matrix(nrow=1,ncol = size)
 for (eac in 1: length(testInd )){
   dishol[eac,1:length(trainInd)] <- disMatNorm[testInd[eac],trainInd]
   HarmDisHold[1,eac] <- harmMean(dishol[eac,])
 }

 HardHiuSample <- matrix(nrow = 50, ncol = 100)
   
 for(i in 1:50){
   HardHiuSample[i,] <- fff2[[1]][(10*i)+31,]
 }
 HardHiuSample <- t(HardHiuSample)
 #############
 ffPANCAN <- FindHardest(disMat,200,664)
 NewDifGrad <- matrix(nrow = nrow(ffPANCAN[[2]]), ncol = ncol(ffPANCAN[[2]]))
 for (i in 1:nrow(ffPANCAN[[2]])){
   NewDifGrad[(nrow(ffPANCAN[[2]])-i +1),] <- ffPANCAN[[2]][i,]
 }
 boxplot.matrix(NewDifGrad,use.cols = F,ylim = c(0,0.8))
 abline ( h = seq(0.2,0.8,0.1),col = 2)
 