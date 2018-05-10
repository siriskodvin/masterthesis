# Setting the current directory and loading/sourcing files
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")
library(simrel)
library(plsdof)
source("Sim1-2a-EstimatesFunctions.R")

# --------------------------------------- #
#                                         #
#              Simulation 2               #
#   PCR and PLSRnaive sigma^2 estimates   #
#          for k = 1, 2, ... , p          #
#                                         #
# --------------------------------------- #

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

seeds <- c(-49, 1, 13, 9, 2957, -4002, 486)
r <- 3  # number of replicates

# Number of components to include in the model
k <- p


# ------------- #
# PCR estimates #
# ------------- #
# Loop for all 8 dp's
sigma2hatPCR <- vector("list", 8)

for(i in 1:8){
  relpos <- relposmat[,i]
  gamma <- gammavec[i]
  R2 <- R2vec[i]
  
  sigma2hatPCRdp <- sigma2hatPCRfunc(n = n,
                                     p = p,
                                     m = m,
                                     q = q,
                                     relpos = relpos,
                                     gamma = gamma,
                                     R2 = R2,
                                     r = r,
                                     seeds = seeds,
                                     k = k)
  
  sigma2hatPCR[[i]] <- sigma2hatPCRdp
  if(i %in% seq(from=1, to=4)){
    names(sigma2hatPCR)[[i]] <- paste("dp", i, sep="")
  }else{
    names(sigma2hatPCR)[[i]] <- paste("dp", i+4, sep="")
  }
}

save(sigma2hatPCR, file="sigma2hatPCR_sim2.RData")


# --------------------- #
#  PLSRnaive estimates  #
# --------------------- #
# Loop for all 8 dp's
sigma2hatPLSRnaive <- vector("list", 8)

for(i in 1:8){
  relpos <- relposmat[,i]
  gamma <- gammavec[i]
  R2 <- R2vec[i]
  
  sigma2hatPLSRdp <- sigma2hatPLSRfunc(n = n,
                                       p = p,
                                       m = m,
                                       q = q,
                                       relpos = relpos,
                                       gamma = gamma,
                                       R2 = R2,
                                       r = r,
                                       seeds = seeds,
                                       k = k)
  
  sigma2hatPLSRnaive[[i]] <- sigma2hatPLSRdp[[1]]
  
  if(i %in% seq(from=1, to=4)){
    names(sigma2hatPLSRnaive)[[i]] <- paste("dp", i, sep="")
  }else{
    names(sigma2hatPLSRnaive)[[i]] <- paste("dp", i+4, sep="")
  }
}

save(sigma2hatPLSRnaive, file="sigma2hatPLSRnaive_sim2.RData")




