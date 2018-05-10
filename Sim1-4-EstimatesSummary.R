# Setting the current directory and loading/sourcing files
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/Estimates")
sigma2hatBayesPLS <- get(load("sigma2hatBayesPLS_full.RData"))
sigma2hatPCR <- get(load("sigma2hatPCR_full.RData"))
sigma2hatPLSRkrylov <- get(load("sigma2hatPLSRkrylov_red.RData"))
sigma2hatPLSRnaive <- get(load("sigma2hatPLSRnaive_full.RData"))
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")

# ------------------------ #
#                          #
#   Summary of estimates   #
#                          #
# ------------------------ #

# Design points (dp)/simrel levels
# n = 50 / 15
# p = 25
# m = 3
# q = 25
# relpos = 1-2-3 / 3-5-7
# gamma = 0.9 / 0.2
# R2 = 0.7 / 0.2
# k = 8

# dplevels
nvec <- rep(c(rep(50,4), rep(15,4)),2)
p <- 25
relposmat <- rbind(rep(c(rep(1,2), rep(3,2)),4), rep(c(rep(2,2), rep(5,2)),4), rep(c(rep(3,2), rep(7,2)),4))
gammavec <- rep(c(0.9, 0.2),8)
R2vec <- c(rep(0.7,8), rep(0.2,8))
sigma2vec <- 1 - R2vec

# Number of components to include in the model
k <- 8

methods <- c("BayesPLS", "PCR", "PLSRkrylov", "PLSRnaive")

# Calculating means, standardized means, MSE, RMSE, standard error and biases of the estimates
Measures <- c("Means", "STMeans", "MSE", "RMSE", "SE", "Bias")
estSummary <- vector("list", 16)
names(estSummary) <- paste("dp", 1:16, sep="")
for(i in 1:16){
  sigma2 <- sigma2vec[i]
  
  dp <- array(rep(0,(length(methods)*length(Measures)*k)), dim=c(length(methods),k,length(Measures)))
  dimnames(dp) <- list(methods, paste("Comp ", 1:k, sep=""), Measures)
  
  for(j in 1:k){
    # Means
    dp[1,j,1] <- mean(sigma2hatBayesPLS[[i]][,,j], na.rm=TRUE)
    dp[2,j,1] <- mean(sigma2hatPCR[[i]][,,j], na.rm=TRUE)
    dp[3,j,1] <- mean(sigma2hatPLSRkrylov[[i]][,,j], na.rm=TRUE)
    dp[4,j,1] <- mean(sigma2hatPLSRnaive[[i]][,,j], na.rm=TRUE)
    
    # Standardized means
    dp[,,2] <- dp[,,1]/sigma2
    
    # Mean square error (MSE)
    dp[1,j,3] <- mean((sigma2hatBayesPLS[[i]][,,j] - sigma2)^2, na.rm=TRUE)
    dp[2,j,3] <- mean((sigma2hatPCR[[i]][,,j] - sigma2)^2, na.rm=TRUE)
    dp[3,j,3] <- mean((sigma2hatPLSRkrylov[[i]][,,j] - sigma2)^2, na.rm=TRUE)
    dp[4,j,3] <- mean((sigma2hatPLSRnaive[[i]][,,j] - sigma2)^2, na.rm=TRUE)
    
    # Root mean square error (RMSE)
    dp[1,j,4] <- sqrt(mean((sigma2hatBayesPLS[[i]][,,j] - sigma2)^2, na.rm=TRUE))
    dp[2,j,4] <- sqrt(mean((sigma2hatPCR[[i]][,,j] - sigma2)^2, na.rm=TRUE))
    dp[3,j,4] <- sqrt(mean((sigma2hatPLSRkrylov[[i]][,,j] - sigma2)^2, na.rm=TRUE))
    dp[4,j,4] <- sqrt(mean((sigma2hatPLSRnaive[[i]][,,j] - sigma2)^2, na.rm=TRUE))
    
    # Standard error (SE)
    dp[1,j,5] <- sqrt(mean((sigma2hatBayesPLS[[i]][,,j] - dp[1,j,1])^2, na.rm=TRUE))
    dp[2,j,5] <- sqrt(mean((sigma2hatPCR[[i]][,,j] - dp[2,j,1])^2, na.rm=TRUE))
    dp[3,j,5] <- sqrt(mean((sigma2hatPLSRkrylov[[i]][,,j] - dp[3,j,1])^2, na.rm=TRUE))
    dp[4,j,5] <- sqrt(mean((sigma2hatPLSRnaive[[i]][,,j] - dp[4,j,1])^2, na.rm=TRUE))
    
    # Average bias
    dp[,,6] <- dp[,,1] - sigma2
  }
  estSummary[[i]] <- dp
}

save(estSummary, file="EstimatesSummary.RData")


