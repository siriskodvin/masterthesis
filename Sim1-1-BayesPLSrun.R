library(simrel)
library(BayesPLS)

# -------------------- #
#                      #
#   BayesPLS objects   #
#                      #
# -------------------- #

# Design points (dp)/simrel levels
# n = 50 / 15
# p = 25
# m = 3
# q = 25
# relpos = 1-2-3 / 3-5-7
# gamma = 0.9 / 0.2
# R2 = 0.7 / 0.2
# k = 5/8

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

kmax <- 8

for(dp in 1:length(nvec)){           # Changing dp
  n <- nvec[dp]
  relpos <- relposmat[,dp]
  gamma <- gammavec[dp]
  R2 <- R2vec[dp]
  
  for(k in 1:kmax){       # Changing the number of components
    if(k == 1){
      burnin <- 5000
    }else if(k == 2){
      if(gamma == 0.9){
        burnin <- 5000
      }else if(gamma == 0.2){
        burnin <- 1000
      }
    }else{
      burnin <- 1000
    }
    
    for(g in 1:2){        # Changing between not approx/approx
      if(g == 1){
        # File names
        combnames <- NULL
        for(a in 1:length(seeds)){
          seed <- NULL
          for(b in 1:r){
            seed <- c(seed, paste("dp", dp, " comp ", k, " seed ", a, " r", b, " run 1", sep=""))
          }
          combnames <- rbind(combnames, seed)
        }
        approx=FALSE
      }else if(g == 2){
        combnames <- NULL
        for(a in 1:length(seeds)){
          seed <- NULL
          for(b in 1:r){
            seed <- c(seed, paste("dp", dp, " comp ", k, " seed ", a, " r", b, " run 1 approx", sep=""))
          }
          combnames <- rbind(combnames, seed)
        }
        approx=TRUE
      }
      
      # Running BayesPLS
      # setwd("C:/Users/siris/Desktop/BayesPLS")
      for(i in 1:length(seeds)){
        set.seed(seeds[i])
        sim <- simrel(n = (r * n), # simulating all replicates at once
                      p = p,
                      m = m,
                      q = q,
                      relpos = relpos,
                      gamma = gamma,
                      R2 = R2)
        for(j in 1:r){
        print(paste("dp", dp))
        print(paste("comp", k))
        print(paste("seed", i))
        print(paste("r =", j))
        
        X <- sim$X[((j-1)*n+1):(j*n), 1:p]
        Y <- sim$Y[((j-1)*n+1):(j*n)]
        
        Bayes1 <- BayesPLS(Y = Y,
                           X = X,
                           ncomp = k,
                           totiter=150000,
                           init.method = "PLS",
                           burnin=burnin,
                           dotrace=TRUE,
                           plotint=5000,
                           thin=10,
                           approx=approx,
                           appit=100,
                           approp=0.95,
                           compreduce = FALSE)
        
        savePlot(filename=combnames[i,j], type="png", device=dev.cur())
        save(Bayes1, file=paste(combnames[i,j], ".RData", sep=""))
        
        dev.off()
        }
      }
    }
  }
}

