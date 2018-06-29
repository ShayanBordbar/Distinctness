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
    HarmS<- HarmS + 1/vec[i]
  }
  HarMM <- i/HarmS
  return(HarMM)
}
############################################################################################
############################################################################################
InitialDis <- function(TFexp, sampleNumber,initCond){
  library(fields)
  # TFexp is the tf expression matrix: each row one TF, each column one condition
  # initCond are the initial conditions chosen randomly
  # compute the harmonic distance of the test set
  # getting the tf distance matrix
  TFdist <- rdist(t(TFexp))
  # minmax normalizing tf distance matrix
  dismn <- min(TFdist[TFdist != 0])
  dismx <- max(TFdist)
  disrng <- dismx - dismn
  TFdistNorm <- (TFdist - dismn) / disrng
  # creating matrix to hold distance of each test condition
  dishol <- matrix(nrow = length(initCond),
                   ncol= (ncol(TFexp) - length(initCond)))
  testInd  <- which(1:ncol(TFexp) %in% initCond)
  trainInd <- which(! 1:ncol(TFexp) %in% initCond)
  HarmDisHold <- matrix(nrow=1, ncol = length(initCond))
  for (eac in 1: length(testInd )){
    dishol[eac, 1:length(trainInd)] <- TFdistNorm[testInd[eac], trainInd]
    HarmDisHold[1, eac] <- harmMean(dishol[eac, ])
  }
  colnames(HarmDisHold) <- testInd
  
  testInd <- which(1:ncol(TFexp) %in% initCond)
  trainInd <- which(! 1:ncol(TFexp) %in% initCond)

  Results <- list()
  Results[[1]] <- testInd # test index
  Results[[2]] <- trainInd # training index
  Results[[3]] <- HarmDisHold
  Results[[4]] <- TFdistNorm
  return(Results)
}
############################################################################################
############################################################################################
DisChange <- function(testInd, subtit, HarmDis, TFDistNorm){
  # testInd is a vector containing the index of the test set conditions before change
  # subtit is a vector of length 2, first entry is the current 
  # test condition which is going to be substituted with the second entry
  # HarmDis is the harmonic distance matrix(1*lenghth(testInd))
  # TFDistNorm is the normalized euclidean distance matrix (631*631) of TF expression vectors 
  NewTestInd <- testInd
  # substitute the old and new index
  NewTestInd[which(testInd == subtit[1])] <- subtit[2]
  NewTrainInd <- which(! 1:ncol(TFDistNorm) %in% NewTestInd)
  NewHarmDis <- HarmDis
  colnames(NewHarmDis) <- NewTestInd
  # change the harmonic distance measure of each of the test points using the substitution
  for (i in 1:length(NewTestInd)){
    NewHarmDis[1, i] <- (length(NewTrainInd)/((length(NewTrainInd)/NewHarmDis[1, i])
                                              + (1/TFDistNorm[NewTestInd[i], subtit[1]])
                                              - (1/TFDistNorm[NewTestInd[i],subtit[2]])))
  }
  # measuring the harmonic distance of the new test point
  NewHarmDis[1, which(as.numeric(colnames(NewHarmDis))== subtit[2])] <- harmMean(TFDistNorm[subtit[2],NewTrainInd])

  Results <- list()
  Results[[1]] <- NewTestInd # test index
  Results[[2]] <- NewTrainInd # training index
  Results[[3]] <- NewHarmDis
  return(Results)
}
#################################################################################################################
#################################################################################################################
#function that performs the simulated annealing
simAnn <- function(itertemp, TempMin, alPha, sampleNumber, TFexpres, initConds){
  # itertemp is the number of iterations in each temperature
  # TempMin is the minimum temperature
  # alpha is the temp change coefficient
  # sampleNumber is the number of conditions in the test set
  # TFexpres is the TF expression matrix (236 * 631)
  # initConds is a vector of conditions to start with or zero if you want to start randomly

  # sample the conditions for initial conditions
  if (initConds == 0){
    initCondit <- sample(x = 1:631 , size = sampleNumber, replace = F)
  }else{
    initCondit <- initConds
  }
  initResults <- InitialDis(TFexp = TFexpres, initCond= initCondit)
  Tmp <- 1
  minTemp <-TempMin
  ALPHA <- alPha
  substit <- numeric(2)
  substit[1] <- sample(x = initResults[[1]], size = 1, replace = F)
  substit[2] <- sample(x = initResults[[2]], size = 1, replace = F)
  FirstChange <- DisChange(testInd = initResults[[1]],
                           subtit = substit,
                           HarmDis = initResults[[3]],
                           TFDistNorm = initResults[[4]] )
  change <- list()
  change[[1]] <- FirstChange
  cnt <- 2
  while (Tmp > minTemp){
    for (tempitr in 1:itertemp){
      print(paste("iteration", tempitr , "/", 
                  itertemp, "in temperature", Tmp, sep = " "))
      substit[1] <- sample(x = change[[cnt-1]][[1]], size = 1, replace = F)
      substit[2] <- sample(x = change[[cnt-1]][[2]], size = 1, replace = F)
      change[[cnt]] <- DisChange(testInd = change[[cnt-1]][[1]],
                                 subtit = substit,
                                 HarmDis = change[[cnt-1]][[3]],
                                 TFDistNorm = initResults[[4]])
      acceptProb <- accProbDis(CurSol = change[[cnt]],
                               PrevSol = change[[cnt-1]],
                               temprature = Tmp)
      if(acceptProb > runif(1)){
        print("YES", quote = F)
        cnt = cnt + 1
        print(paste("Current solution number is:", cnt, sep = " "))
      }
    }
    Tmp = Tmp * ALPHA
  }
  return(change)
}
#############################################################################################################
#############################################################################################################
####Acceptance probability function

accProbDis <- function(CurSol, PrevSol, temprature){
  # CurSol is the solution after change
  # PrevSol is the solution before change
  # Option is the option to use for calculating probability: it can be "HarmDis","RuggedRandom","HarmRuggedRandom"
    acceptProb <- exp(1)^(((mean(CurSol[[3]], na.rm = T) -
                            mean(PrevSol[[3]], na.rm = T)) * 100)/temprature)
  return(acceptProb)
}
########
