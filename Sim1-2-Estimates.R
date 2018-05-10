# Setting the current directory and loading/sourcing files
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")
library(simrel)
library(plsdof)
library(BayesPLS)
source("Sim1-2a-EstimatesFunctions.R")

# --------------------- #
#                       #
#   sigma^2 estimates   #
#                       #
# --------------------- #

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
m <- 3
q <- 25
relposmat <- rbind(rep(c(rep(1,2), rep(3,2)),4), rep(c(rep(2,2), rep(5,2)),4), rep(c(rep(3,2), rep(7,2)),4))
gammavec <- rep(c(0.9, 0.2),8)
R2vec <- c(rep(0.7,8), rep(0.2,8))

seeds <- c(-49, 1, 13, 9, 2957, -4002, 486)
r <- 3  # number of replicates

# Number of components to include in the model
k <- 8


# ------------- #
# PCR estimates #
# ------------- #
# Loop for all 16 dp's
sigma2hatPCR <- vector("list", 16)

for(i in 1:16){
  n <- nvec[i]
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
  names(sigma2hatPCR)[[i]] <- paste("dp", i, sep="")
}

save(sigma2hatPCR, file="sigma2hatPCR_full.RData")


# ---------------- #
#  PLSR estimates  #
# Krylov and naive #
#       DoF        #
# ---------------- #
# Loop for all 16 dp's
sigma2hatPLSRnaive <- vector("list", 16)
sigma2hatPLSRkrylov <- vector("list", 16)
DoFkrylov <- vector("list", 16)

for(i in 1:16){
  n <- nvec[i]
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
  sigma2hatPLSRkrylov[[i]] <- sigma2hatPLSRdp[[2]]
  DoFkrylov[[i]] <- sigma2hatPLSRdp[[3]]
  
  names(sigma2hatPLSRnaive)[[i]] <- paste("dp", i, sep="")
  names(sigma2hatPLSRkrylov)[[i]] <- paste("dp", i, sep="")
  names(DoFkrylov)[[i]] <- paste("dp", i, sep="")  
}

save(sigma2hatPLSRnaive, file="sigma2hatPLSRnaive_full.RData")
save(sigma2hatPLSRkrylov, file="sigma2hatPLSRkrylov_full.RData")
save(DoFkrylov, file="DoFkrylov.RData")


# ------------------ #
# BayesPLS estimates #
# ------------------ #
# Obtaining the estimates manually for each dp-comp
setwd("C:/Users/siris/Desktop/Estimation/BayesPLSdp14/8 comp/Estimation")

dp <- 14
n <- nvec[dp]
relpos <- relposmat[,dp]
gamma <- gammavec[dp]
R2 <- R2vec[dp]
k <- 8

sigma2hatBayesdp14comp8 <- matrix(rep(0, (r * length(seeds))), nrow=r, ncol=length(seeds))

# File names
combnames <- NULL
for(a in 1:length(seeds)){
  seed <- NULL
  for(b in 1:r){
    seed <- c(seed, paste("dp", dp, " comp ", k, " seed ", a, " r", b, ".Rdata", sep=""))
  }
  combnames <- cbind(combnames, seed)
}

for(i in 1:length(seeds)){
  for(j in 1:r){
    load(combnames[j,i])
    est <- estimate.BayesPLS(Bayes1, start=5000, stop=15000)
    sigma2hatBayesdp14comp8[j,i] <- est$sigma.sq
  }
}

load("dp14 comp 8 seed 4 r1.RData")
est <- estimate.BayesPLS(Bayes1, start=5000, stop=25000)
sigma2hatBayesdp14comp8[1,4] <- est$sigma.sq
sigma2hatBayesdp14comp8

colnames(sigma2hatBayesdp14comp8) <- paste("Seed", seq(1:(length(seeds))))
rownames(sigma2hatBayesdp14comp8) <- paste("r", seq(1:r), sep="")

setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")
save(sigma2hatBayesdp14comp8, file="sigma2hatBayesdp14comp8.RData")

# Organizing the estimates into an array for each dp
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/Estimates")
load("sigma2hatBayesdp14comp1.RData")
load("sigma2hatBayesdp14comp2.RData")
load("sigma2hatBayesdp14comp3.RData")
load("sigma2hatBayesdp14comp4.RData")
load("sigma2hatBayesdp14comp5.RData")
load("sigma2hatBayesdp14comp6.RData")
load("sigma2hatBayesdp14comp7.RData")
load("sigma2hatBayesdp14comp8.RData")
sigma2hatBayesdp14 <- sapply(list(sigma2hatBayesdp14comp1,
                                  sigma2hatBayesdp14comp2,
                                  sigma2hatBayesdp14comp3,
                                  sigma2hatBayesdp14comp4,
                                  sigma2hatBayesdp14comp5,
                                  sigma2hatBayesdp14comp6,
                                  sigma2hatBayesdp14comp7,
                                  sigma2hatBayesdp14comp8), identity, simplify="array")
dimnames(sigma2hatBayesdp14)[[3]] <- paste("comp", 1:8)
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")
save(sigma2hatBayesdp14, file="sigma2hatBayesdp14.RData")

# Collecting all arrays into a list (now containing all estimates)
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/Estimates")
sigma2hatBayesPLS <- vector("list", 16)
for(i in 1:16){
  filename <- paste("sigma2hatBayesdp", i, ".RData", sep="")
  sigma2hatBayes[[i]] <- get(load(filename))
  names(sigma2hatBayes)[[i]] <- paste("dp", i, sep="")
}
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/")
save(sigma2hatBayesPLS, file="sigma2hatBayesPLS_full.RData")

