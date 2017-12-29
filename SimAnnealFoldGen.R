#functions: harmMean, InitialDisRugDif, InitialDisRugDifONEGENE, RugDisDifChange, RugDisDifChangeONEGENE, simAnn, simAnnONEGENE, accProb
##In this script I want to implemet simulated annealing in order to creat cross validation 
#folds (for monte-carlo cross-validation) with different difficulty levels.
##first I need to define the function which computes the difficulty of a fold.
##I assume that distance between conditions and the raw ruggedness of each gene condition is known.
###UPDATE: RUGGEDNESS won't be consisdered for calculating distinctness
###The main function to be used for constructing folds is: simAnn
############################################################################################
#function to compute the harmonic mean of the distance for the test point
#given min-max normalized distance vector
#calculate the harmonic mean given a vector
harmMean <- function(vec){
  HarmS <- 0
  for (i in 1:length(vec)){
    #if (vec[i] == 0){
    #  stop("Division by zero")
   # }
    HarmS<- HarmS + 1/vec[i]
  }
  HarMM <- i/HarmS
  return(HarMM)
}
############################################################################################
############################################################################################
InitialDis <- function(TFexp, sampleNumber,initCond){
  library(fields)
  #TFexp is the tf expression matrix: each row one TF, each column one condition
  #rawRugg is raw ruggedness for each gene, condition pair
  #initCond are the initial conditions chosen randomly
  #compute the harmonic distance of the test set
  #getting the tf distance matrix
  TFdist <- rdist(t(TFexp))
  #minmax normalizing tf distance matrix
  dismn <- min(TFdist[TFdist != 0])
  dismx <- max(TFdist)
  disrng <- dismx - dismn
  TFdistNorm <- (TFdist-dismn)/disrng
  #creating matrix to hold distance of each test condition
  dishol <- matrix(nrow = length(initCond), ncol= (ncol(TFexp) - length(initCond)))
  testInd <- which(1:ncol(TFexp) %in% initCond)
  trainInd <- which(! 1:ncol(TFexp) %in% initCond)
  HarmDisHold <- matrix(nrow=1,ncol = length(initCond))
  for (eac in 1: length(testInd )){
    dishol[eac,1:length(trainInd)] <- TFdistNorm[testInd[eac],trainInd]
    HarmDisHold[1,eac] <- harmMean(dishol[eac,])
  }
  colnames(HarmDisHold) <- testInd
  #compute the initial ruggedness of the test set
  #minmax normalize ruggedness per gene
#   rugmnGene = numeric(nrow(rawRugg))
#   rugmxGene = numeric(nrow(rawRugg))
#   rugrngGene = numeric(nrow(rawRugg))
#   for (i in 1: nrow(rawRugg)){
#     rugmnGene[i] <- min(rawRugg[i,])
#     rugmxGene[i] <- max(rawRugg[i,])
#     rugrngGene[i] <- rugmxGene[i] - rugmnGene[i]
#   }
#   RawRugNorm <- (rawRugg - rugmnGene)/rugrngGene
  
  testInd <- which(1:ncol(TFexp) %in% initCond)
  trainInd <- which(! 1:ncol(TFexp) %in% initCond)
#   RugMat <- matrix(nrow = nrow(rawRugg), ncol = length(testInd) )#matrix that hold the ruggedness of prediction of each gene in each test condition
#   weightedRug <- matrix(nrow = nrow(rawRugg), ncol = length(testInd) )#matrix that hold the weightedRuggedness[not average weighted] of prediction of each gene in each test condition
#   weightOfRug <-matrix(nrow = nrow(rawRugg), ncol = length(testInd) )#matrix that hold the weight of ruggedness of prediction of each gene in each test condition
#   #compute the weighted average ruggedness for each gene in each test point
#   for (test.gene in 1:nrow(rawRugg)){
#     for (test.con in 1: length(testInd)){
#       har <- 0
#       wei <- 0
#       for (train.con in 1:length(trainInd)){
#         x <- (1 - TFdistNorm[testInd[test.con],trainInd[train.con]])* (RawRugNorm[test.gene,trainInd[train.con]])
#         har <- har + x
#         wei <- wei + (1 - TFdistNorm[testInd[test.con],trainInd[train.con]])
#       }
#       RugMat[test.gene,test.con] <- har/wei
#       weightedRug[test.gene,test.con] <- har
#       weightOfRug[test.gene,test.con] <- wei
#     }
#   }
#   colnames(RugMat)<-testInd
  #calculate the hardness using distance and ruggedness
#   DifficultyTest <- matrix(nrow = nrow(RugMat), ncol = ncol(RugMat) )
#   for (i in 1: nrow(RugMat)){
#     for (j in 1: ncol(RugMat)){
#       DifficultyTest[i,j] <- HarmDisHold[1,j] * RugMat[i,j]
#     }
#   }
  Results <- list()
  Results[[1]] <- testInd #test index
  Results[[2]] <- trainInd # training index
  Results[[3]] <- HarmDisHold
#   Results[[4]] <- RugMat
#   Results[[5]] <- weightedRug
#   Results[[6]] <- weightOfRug
#   Results[[7]] <- DifficultyTest
  Results[[4]] <- TFdistNorm
  # Results[[9]] <- RawRugNorm
  return(Results)
}
############################################################################################
############################################################################################
DisChange <- function(testInd,subtit,HarmDis,TFDistNorm){
  #testInd is a vector containing the index of the test set conditions before change
  #subtit is a vector of length 2, first entry is the current 
  #test condition which is going to be changed with second entry
  #HarmDis is the harmonic distance matrix(1*lenghth(testInd))
  #Ruggedness is the ruggedness matrix (500*lenghth(testInd)) of the test set
  #Difficulty is the difficulty matrix (500*lenghth(testInd)) of the test set
  #TFDistNorm is the normalized euclidean distance matrix (631*631) of TF expression vectors 
  #rawRuggedness is the normalized raw ruggedness matrix (500 * 631)
  NewTestInd <- testInd
  #substitute the old and new index
  NewTestInd[which(testInd == subtit[1])] <- subtit[2]
  NewTrainInd <- which(! 1:ncol(TFDistNorm) %in% NewTestInd)
  NewHarmDis <- HarmDis
  colnames(NewHarmDis) <- NewTestInd
  #change the harmonic distance measure of each of the test points using the substitution
  for (i in 1:length(NewTestInd)){
    NewHarmDis[1,i] <- (length(NewTrainInd)/((length(NewTrainInd)/NewHarmDis[1,i]) + (1/TFDistNorm[NewTestInd[i],subtit[1]]) - (1/TFDistNorm[NewTestInd[i],subtit[2]])) )
  }
  #measuring the harmonic distance of the new test point
  NewHarmDis[1,which(as.numeric(colnames(NewHarmDis))== subtit[2])] <- harmMean(TFDistNorm[subtit[2],NewTrainInd])
  #####Measuring the new ruggedness
#   NewRugged <- Ruggedness
#   NewWightedRug <- weightedRug
#   NewWightOfRug <- weightOfRug
#   colnames(NewRugged) <- NewTestInd
#   for (geneNum in 1:nrow(NewRugged)){
#     for (test.con in 1:ncol(NewRugged)){
#       NewWightedRug[geneNum,test.con] <- weightedRug[geneNum,test.con] - ((1-TFDistNorm[test.con,subtit[2]]) * rawRuggednessNorm[geneNum,subtit[2]]) + ((1-TFDistNorm[test.con,subtit[1]]) * rawRuggednessNorm[geneNum,subtit[1]])
#       NewWightOfRug[geneNum,test.con] <- weightOfRug[geneNum,test.con] - (1-TFDistNorm[test.con,subtit[2]]) + (1-TFDistNorm[test.con,subtit[1]])
#       NewRugged[geneNum,test.con] <- NewWightedRug[geneNum,test.con] / NewWightOfRug[geneNum,test.con]
#     }
#     har <- 0
#     wei <- 0
#     for (train.con in 1:length(NewTrainInd)){
#       x <- (1 - TFdistNorm[subtit[2],NewTrainInd[train.con]])* (rawRuggednessNorm[geneNum,NewTrainInd[train.con]])
#       har <- har + x
#       wei <- wei + (1 - TFdistNorm[subtit[2],NewTrainInd[train.con]])
#     }
#     NewWightedRug[geneNum,which(as.numeric(colnames(NewRugged))== subtit[2])] <- har
#     NewWightOfRug[geneNum,which(as.numeric(colnames(NewRugged))== subtit[2])] <- wei
#     NewRugged[geneNum,which(as.numeric(colnames(NewRugged))== subtit[2])] <- har/wei
#   }
#   #####Measuring the new difficulty
#   NewDif <- Difficulty
#   ChangeOfDif <- matrix(nrow = nrow(NewDif), ncol = ncol(NewDif))
#   for (i in 1:nrow(Difficulty)){
#     NewDif[i,] <- NewRugged[i,] * NewHarmDis[1,]
#   }
#   ChangeOfDif <- NewDif - Difficulty
  Results <- list()
  Results[[1]] <- NewTestInd #test index
  Results[[2]] <- NewTrainInd # training index
  Results[[3]] <- NewHarmDis
#   Results[[4]] <- NewRugged
#   Results[[5]] <- NewWightedRug
#   Results[[6]] <- NewWightOfRug
#   Results[[7]] <- NewDif
#   Results[[8]] <- ChangeOfDif
  return(Results)
}
#################################################################################################################
#################################################################################################################
#function that performs the simulated annealing
simAnn <- function(itertemp, TempMin, alPha,sampleNumber,TFexpres,initConds){
  #itertemp is the number of iterations in each temperature
  #TempMin is the minimum temperature
  #alpha is the temp change coefficient
  #TFexpres is the TF expression matrix (236 * 631)
  #rawRugged is the Raw ruggedness matrix (500*631)
  #sampleNumber is the number of conditions in the test set
  #sample the conditions for initial condistions
  #first choose a random set of samplenNmber conditions as test conditions and compute the 
  #optoption is the option in which the optimization is base on: one of the following: "HarmDis","RuggedRandom","HarmRuggedRandom"
  #initConds is a vector of conditions to start with or zero if you want to start randomly
  if (initConds == 0){
    initCondit <- sample(x = 1:631 , size = sampleNumber,replace = F)
  }else{
    initCondit <- initConds
  }
  initResults <- InitialDis(TFexp = TFexpres,initCond= initCondit)
  Tmp <- 1
  minTemp <-TempMin
  ALPHA <- alPha
  randGene <- sample(1:500,size = 1,replace = F)
  substit <- numeric(2)
  substit[1] <- sample(x = initResults[[1]] , size = 1,replace = F)
  substit[2] <- sample(x = initResults[[2]] , size = 1,replace = F)
  FirstChange <- DisChange(testInd = initResults[[1]],subtit = substit,HarmDis = initResults[[3]], TFDistNorm = initResults[[4]] )
  change <- list()
  change[[1]] <- FirstChange
  cnt <- 2
  while (Tmp > minTemp){
    for (tempitr in 1:itertemp){
      print(paste("iteration", tempitr ,"/",itertemp,"in temperature",Tmp,sep = " "))
      substit[1] <- sample(x = change[[cnt-1]][[1]] , size = 1,replace = F)
      substit[2] <- sample(x = change[[cnt-1]][[2]] , size = 1,replace = F)
      change[[cnt]] <- DisChange(testInd = change[[cnt-1]][[1]],subtit = substit,HarmDis = change[[cnt-1]][[3]], TFDistNorm = initResults[[4]])
      #acceptProb <- exp(1)^((mean(colMeans(t(change[[cnt-1]][[8]]),na.rm = T),na.rm = T))/Tmp)
      #acceptProb <- exp(1)^(((mean(change[[cnt]][[3]],na.rm = T) - mean(change[[cnt-1]][[3]],na.rm = T))*100)/Tmp)
      acceptProb <- accProbDis(CurSol = change[[cnt]],PrevSol = change[[cnt-1]],temprature = Tmp)
      if(acceptProb > runif(1)){
        print("YES",quote = F)
        cnt = cnt + 1
        print(paste("Current solution number is:", cnt,sep = " "))
      }
    }
    Tmp = Tmp * ALPHA
  }
  #print(paste("Random gene index chosen to optimize ruggedness is:",randGene,sep = " "))
  return(change)
}
#############################################################################################################
#############################################################################################################
####Acceptance probability function
accProb <- function(CurSol,PrevSol,temprature,Option,randGen){
  #CurSol is the solution after change
  #PrevSol is the solution before change
  #Option is the option to use for calculating probability: it can be "HarmDis","RuggedRandom","HarmRuggedRandom"
  if (Option == "HarmDis"){
    acceptProb <- exp(1)^(((mean(CurSol[[3]],na.rm = T) - mean(PrevSol[[3]],na.rm = T))*100)/temprature)
  }else if (Option == "RuggedRandom"){
    acceptProb <- exp(1)^(((mean(CurSol[[4]],na.rm = T) - mean(PrevSol[[4]],na.rm = T))*100)/temprature)
  }else if (Option == "HarmRuggedRandom"){
    acceptProb1 <- exp(1)^(((mean(CurSol[[3]],na.rm = T) - mean(PrevSol[[3]],na.rm = T))*100)/temprature)
    acceptProb2 <- exp(1)^(((mean(CurSol[[4]][randGen,],na.rm = T) - mean(PrevSol[[4]][randGen,],na.rm = T))*100)/temprature)
    acceptProb <- acceptProb1 * acceptProb2
  }
  return(acceptProb)
}
###########
accProbDis <- function(CurSol,PrevSol,temprature){
  #CurSol is the solution after change
  #PrevSol is the solution before change
  #Option is the option to use for calculating probability: it can be "HarmDis","RuggedRandom","HarmRuggedRandom"
    acceptProb <- exp(1)^(((mean(CurSol[[3]],na.rm = T) - mean(PrevSol[[3]],na.rm = T))*100)/temprature)
  return(acceptProb)
}
########
