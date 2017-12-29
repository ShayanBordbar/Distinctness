#PanCan Distinctness calculations:
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#calculating hardness of the  cluster based partitionings
FoldHardnessCalculator <- function(CLusterList, DistanceMatrix,no.folds,HarmDisRandom){
  #clusterList is a list of partitionings
  #HarmDisRandom is a vector of harmonic distances for random partitioning to be used in the comparison plots or Put Zero
  par(mfrow = c(3,4),mar = c(0.2, 0.2, 0.2, 0.2) )
  dismn <- min(DistanceMatrix[DistanceMatrix != 0])
  dismx <- max(DistanceMatrix)
  disrng <- dismx - dismn
  no.Clus <- length(CLusterList)
  #initializing the matrix that each row is one partitioning
  #/each column is a condition indexed based on the actual index of the conditons
  HarmDisHoldSorted <- matrix(nrow = no.Clus,ncol = ncol(DistanceMatrix))
  #indexed based on clusters
  HarmDisHold <- matrix(nrow = no.Clus,ncol = ncol(DistanceMatrix))
  HarmDisHoldPerfold <- matrix(nrow = no.Clus,ncol = no.folds)
  #for each clustering enters the cluster 
  #and computes the harmonic distance for all conditions based on their distance to their training points
  for (clus.cnt in 1:no.Clus){
    #
    dishol <- matrix(nrow = ncol(DistanceMatrix), ncol= (ncol(DistanceMatrix) - min(CLusterList[[clus.cnt]]$size)))
    rownames(dishol) <- numeric(ncol(DistanceMatrix))
    indcount <- 0
    #for each fold get the test and training index
    for (eaf in 1:no.folds){
      testInd <- which(CLusterList[[clus.cnt]]$cluster %in% eaf)
      trainInd <- which(! CLusterList[[clus.cnt]]$cluster %in% eaf)
      #for each test point gather its distance to its training points
      #dishold each row is the distance of all training points of each point
      for (eac in 1: length(testInd )){
        dishol[indcount + eac,1:length(trainInd)] <- DistanceMatrix[testInd[eac],trainInd]
        rownames(dishol)[indcount + eac] <- testInd[eac]
      }
      indcount <- indcount + length(testInd )
    }
    #minmax norm the train test distance matrix
    disholNorm <- (dishol - dismn)/disrng
    #calculating harmonic mean of distance of each test condition from all of its training conditions
    HarmDis <- numeric(ncol(DistanceMatrix)) # vector containing the harmonic mean of the distance of each condition from the its training set
    #computes the harmonic distance for each point but the order is not as the first index
    for (i in 1: length(HarmDis)){
      HarmDis[i] <- 1/mean(1/disholNorm[i,],na.rm = T)
    }
    #return the index system to the actual first one[The current function returns the clustered indexing change if you want otherwise] 
    IndConv <- sort(x =rownames(dishol),decreasing = F,index.return = T)$ix
    HarmDisHoldSorted[clus.cnt,] <- HarmDis[IndConv]
    HarmDisHold[clus.cnt,] <- HarmDis
    tIndcnt <- 0
    for (harmPf in 1:no.folds){
      HarmDisHoldPerfold[clus.cnt,harmPf] <- mean(HarmDis[(tIndcnt + 1):(tIndcnt + CLusterList[[clus.cnt]]$size[harmPf])],na.rm = T)
      tIndcnt <- tIndcnt + CLusterList[[clus.cnt]]$size[harmPf]
    }
    if(! HarmDisRandom == 0){
      #plotting the harmonic distances comparing to a random fold choice##[NEEDS TO BE SORTED --> [NO/Yes] [Draw border/Don't draw border]
      plotHarmDistance(RandomHarmDis = HarmDisRandom,ClusHarmDis = HarmDis,isSort = F,CLuster = CLusterList[[clus.cnt]])
    }
  }#end of each partitioning
  #boxplot harmonic distance for all the partitionings 
  par(mfrow = c(1,1),mar = c(2, 2, 2, 2) )
  boxplot.matrix(HarmDisHold ,use.cols = F,axes = T,oma = c(5, 5, 5, 5), xlab = "Partitioning no.", main = "Harmonic Distance")
  #calculating the mean harmonic distance in order to find the cluster with maximum mean
  meanHarmhol <- numeric(nrow(HarmDisHold))
  for (i in 1: nrow(HarmDisHold)){
    meanHarmhol[i] <- mean(HarmDisHold[i,])
  }
  resultss <- list()
  resultss[[1]] <- HarmDisHold
  resultss[[2]] <- HarmDisHoldPerfold
  
  return(resultss)
}
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#Function to plot two vectors of harmonic distances in comparison 
plotHarmDistance <- function(RandomHarmDis,ClusHarmDis,isSort,CLuster){
  #isSort is either True(if the vector has been sorted back to initial index ) or False [otherwise]
  plot(RandomHarmDis, type ="l", xaxt = "n", yaxt = "n",col = "blue", ylim = c(0,1))
  lines(ClusHarmDis, col = "red")
  if(! isSort) {
    bord <- 0
    for (i in 1:length(CLuster$size)){
      bord = bord +CLuster$size[i]
      abline(v= (bord),col=3)
    }
  }
}
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
DistanceClustered <- FoldHardnessCalculator(CLusterList = CLchosen,DistanceMatrix = disMat,no.folds = 6,HarmDisRandom = 0)
DistanceClustered2 <- FoldHardnessCalculator(CLusterList = ClList2,DistanceMatrix = disMat,no.folds = 6,HarmDisRandom = 0)

DistanceIden <- matrix(nrow = 100,ncol = 100)
for (i in 1:100){
  for (j in i:100){
    DistanceIden[i,j] <- identical(sort(DistanceClustered2[i,]),sort(DistanceClustered2[j,]))
  }
}

#calculate the hardness of each cluster, to use for correlation calculations
clusterHardness2 <- FoldHardnessCalculator(CLchosen2, disMat, 6,0)
boxplot.matrix(clusterHardness2[[2]],use.cols = F)
ClusterPfoldHardness <- matrix(t(clusterHardness2[[2]]),nrow = 1,ncol=60)
#construct the hardness matrix for 500 genes: basically 500 copies of the above matrix
ClusterPfoldHardness500genes <- matrix(nrow = 500,ncol = 60)
for (i in 1:500){
  ClusterPfoldHardness500genes[i,]<- ClusterPfoldHardness
}

#Calculate the hardness of random partitionings
RandomHardnessCalculator <- function(indices,borders,DistanceMatrix){
  #indices is a list containing the numeric index of the conditions as each entry of the list
  #borders is a list containing the index of the border between folds of each partitioning as each entry of the list, each entry of the list has (no.folds +1) length
  dismn <- min(DistanceMatrix[DistanceMatrix != 0])
  dismx <- max(DistanceMatrix)
  disrng <- dismx - dismn
  no.Clus <- length(indices)
  no.folds <- length(borders[[1]])-1
  #initializing the matrix that each row is one partitioning
  #/each column is a condition indexed based on the actual index of the conditons
  HarmDisHoldSorted <- matrix(nrow = no.Clus,ncol = ncol(DistanceMatrix))
  #indexed based on clusters
  HarmDisHold <- matrix(nrow = no.Clus,ncol = ncol(DistanceMatrix))
  HarmDisHoldPerfold <- matrix(nrow = no.Clus,ncol = no.folds)
  #for each clustering enters the cluster 
  #and computes the harmonic distance for all conditions based on their distance to their training points
  for (clus.cnt in 1:no.Clus){
    #
    dishol <- matrix(nrow = ncol(DistanceMatrix), ncol= (ncol(DistanceMatrix) - min((borders[[clus.cnt]][3] - borders[[clus.cnt]][2]),(borders[[clus.cnt]][2] - borders[[clus.cnt]][1]) )))
    rownames(dishol) <- numeric(ncol(DistanceMatrix))
    indcount <- 0
    #for each fold get the test and training index
    for (eaf in 1:no.folds){
      testInd <- indices[[clus.cnt]][(borders[[clus.cnt]][eaf]+1):borders[[clus.cnt]][eaf+1]]
      trainInd <- setdiff(x = c(1:ncol(DistanceMatrix)),y = testInd)
      #for each test point gather its distance to its training points
      #dishold each row is the distance of all training points of each point
      for (eac in 1: length(testInd )){
        dishol[indcount + eac,1:length(trainInd)] <- DistanceMatrix[testInd[eac],trainInd]
        rownames(dishol)[indcount + eac] <- testInd[eac]
      }
      indcount <- indcount + length(testInd )
    }
    #minmax norm the train test distance matrix
    disholNorm <- (dishol - dismn)/disrng
    #calculating harmonic mean of distance of each test condition from all of its training conditions
    HarmDis <- numeric(ncol(DistanceMatrix)) # vector containing the harmonic mean of the distance of each condition from the its training set
    #computes the harmonic distance for each point but the order is not as the first index
    for (i in 1: length(HarmDis)){
      HarmDis[i] <- 1/mean(1/disholNorm[i,],na.rm = T)
    }
    #return the index system to the actual first one[The current function returns the clustered indexing change if you want otherwise] 
    IndConv <- sort(x =rownames(dishol),decreasing = F,index.return = T)$ix
    HarmDisHoldSorted[clus.cnt,] <- HarmDis[IndConv]
    HarmDisHold[clus.cnt,] <- HarmDis
    
    for (harmPf in 1:no.folds){
      HarmDisHoldPerfold[clus.cnt,harmPf] <- mean(HarmDis[(borders[[clus.cnt]][harmPf]+1):borders[[clus.cnt]][harmPf+1]],na.rm = T)
    }
#     if(! HarmDisRandom == 0){
#       #plotting the harmonic distances comparing to a random fold choice##[NEEDS TO BE SORTED --> [NO/Yes] [Draw border/Don't draw border]
#       plotHarmDistance(RandomHarmDis = HarmDisRandom,ClusHarmDis = HarmDis,isSort = F,CLuster = CLusterList[[clus.cnt]])
#     }
  }#end of each partitioning
  #boxplot harmonic distance for all the partitionings 
  par(mfrow = c(1,1),mar = c(2, 2, 2, 2) )
  boxplot.matrix(HarmDisHold ,use.cols = F,axes = T,oma = c(5, 5, 5, 5), xlab = "Partitioning no.", main = "Harmonic Distance")
  #calculating the mean harmonic distance in order to find the cluster with maximum mean
  meanHarmhol <- numeric(nrow(HarmDisHold))
  for (i in 1: nrow(HarmDisHold)){
    meanHarmhol[i] <- mean(HarmDisHold[i,])
  }
  resultss <- list()
  resultss[[1]] <- HarmDisHold
  resultss[[2]] <- HarmDisHoldPerfold
  return(resultss)
}
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
######Calculating the random 
RandomHardness2 <- RandomHardnessCalculator(indices = indR,borders = bordR,DistanceMatrix = disMat)
RandomPfoldHardness <- matrix(t(RandomHardness2[[2]]),nrow = 1,ncol=60)
#construct the hardness matrix for 500 genes: basically 500 copies of the above matrix
RandomPfoldHardness500genes <- matrix(nrow = 500,ncol = 60)
for (i in 1:500){
  RandomPfoldHardness500genes[i,]<- RandomPfoldHardness
}
#####This is the matrix to be used for correlation calculations of performance and hardness
HardnessClusRand <- cbind(ClusterPfoldHardness500genes,RandomPfoldHardness500genes)
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#compute the hardness of the validation set
#creating the distance matrix 
TFexpValidNormSamp <- GENETFValidExpNormalizedSampled[478:1716,]
TFexpTrTesVal <- cbind(SampledTFExpDataImputed,TFexpValidNormSamp)
library(fields)
TFexpTrTesValTran <- t(TFexpTrTesVal)
disMatTrTesVal <- rdist(TFexpTrTesValTran)
boxplot.matrix(disMatTrTesVal)
IndexValidation <- list()
IndexValidation[[1]] <- c(1:936)
bordersValidation <- list()
bordersValidation[[1]] <- c(0,864,936)
######
#####Second validation set hardness calculations
disMatProstateVal <- rdist(t(TFmatValidationProstateFinal)) 
disMatProstateVal2mat <- rdist(t(TFmatValidationProstateFinalMatlab))
IndexValidation2 <- list()
IndexValidation2[[1]] <- c(1:1064)
bordersValidation2 <- list()
bordersValidation2[[1]] <- c(0,200,1064)
ValidationSetHardnessPros3 <- RandomHardnessCalculator(indices = IndexValidation2,borders = bordersValidation2,DistanceMatrix = disMatProstateVal2mat)
#ValidationSetHardnessMatPros <- matrix(t(ValidationSetHardness[[1]]),nrow = 72,ncol=1)
ValidationSetHardnessPros[[1]][1:200] # is the hardness of each validation condition
#####
ValidationSetHardness <- RandomHardnessCalculator(indices = IndexValidation,borders = bordersValidation,DistanceMatrix = disMatTrTesVal)
ValidationSetHardnessMat <- matrix(t(ValidationSetHardness[[1]]),nrow = 72,ncol=1)
######Hardness of the validation set is 0.4640503 given 864 conditions as the training set
#####this value is close to some clusters in clustering based approach
#calculate the varinace of TF expression accross conditions old vs prostate
OldTFVariance2 <- matrix(nrow = 864,ncol=1)
ProstTFVariance2 <- matrix(nrow = 200,ncol=1)
for(i in 201:1064){
  OldTFVariance2[i-200,1] <- var(TFmatValidationProstateFinal[,i])
  #ProstTFVariance2[i,1] <- var(TFmatValidationProstateFinal[,i])
}
boxplot.matrix(cbind(OldTFVariance,ProstTFVariance))
ProstTFVarianceB <- matrix(nrow = 477,ncol=1)
for(i in 1:477){
  ProstTFVarianceB[i,1]<- var(dataProstateNum477[i,SampledCondsProstate])
}

##################################################################################################


