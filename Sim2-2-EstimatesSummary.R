# Setting the current directory and loading/sourcing files
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/Estimates_sim2")
sigma2hatPCR <- get(load("sigma2hatPCR_sim2.RData"))
sigma2hatPLSRnaive <- get(load("sigma2hatPLSRnaive_sim2.RData"))
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")

# ------------------------ #
#                          #
#   Summary of estimates   #
#          sim 2           #
#                          #
# ------------------------ #

# Design points (dp)/simrel levels
# n = 50
# p = 25
# m = 3
# q = 25
# relpos = 1-2-3 / 3-5-7
# gamma = 0.9 / 0.2
# R2 = 0.7 / 0.2
# k = 25

# dplevels (8 dp's, n is held fixed)
n <- 50
p <- 25
m <- 3
q <- 25
relposmat <- rbind(rep(c(rep(1,2), rep(3,2)),2), rep(c(rep(2,2), rep(5,2)),2), rep(c(rep(3,2), rep(7,2)),2))
gammavec <- rep(c(0.9, 0.2),4)
R2vec <- c(rep(0.7,4), rep(0.2,4))
sigma2vec <- 1 - R2vec

# Number of components to include in the model
k <- p

methods <- c("PCR", "PLSRnaive")

# Calculating means, MSE, variance, and standardized means of the estimates
Measures <- c("Means", "STMeans", "RMSE", "SD", "Bias")
estSummary <- vector("list", 8)
names(estSummary) <- paste("dp", c(seq(from=1, to=4), seq(from=9, to=12)), sep="")
for(i in 1:8){
  sigma2 <- sigma2vec[i]
  
  dp <- array(rep(0,(length(methods)*length(Measures)*k)), dim=c(length(methods),k,length(Measures)))
  dimnames(dp) <- list(methods, paste("Comp ", 1:k, sep=""), Measures)
  
  for(j in 1:k){
    # Means
    dp[1,j,1] <- mean(sigma2hatPCR[[i]][,,j], na.rm=TRUE)
    dp[2,j,1] <- mean(sigma2hatPLSRnaive[[i]][,,j], na.rm=TRUE)
    
    # Standardized means
    dp[,,2] <- dp[,,1]/sigma2
    
    # Root mean square error (RMSE)
    dp[1,j,3] <- sqrt(mean((sigma2hatPCR[[i]][,,j] - sigma2)^2, na.rm=TRUE))
    dp[2,j,3] <- sqrt(mean((sigma2hatPLSRnaive[[i]][,,j] - sigma2)^2, na.rm=TRUE))
    
    # Standard error (SD)
    dp[1,j,4] <- sqrt(mean((sigma2hatPCR[[i]][,,j] - dp[1,j,1])^2, na.rm=TRUE))
    dp[2,j,4] <- sqrt(mean((sigma2hatPLSRnaive[[i]][,,j] - dp[2,j,1])^2, na.rm=TRUE))
    
    # Average bias
    dp[,,5] <- dp[,,1] - sigma2
  }
  estSummary[[i]] <- dp
}

save(estSummary, file="EstimatesSummary_sim2.RData")


