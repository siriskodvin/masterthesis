library(simrel)
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")

# --------------------------- #
#                             #
#      Property plots of      #
#   examples of simulations   #
#        from each dp         #
#                             #
# --------------------------- #

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

# Number of components to include in the plots
k <- 15

# Generating 4 pdf's each with 4 (x 2) plots
for(i in 1:4){
  pdftitle <- paste("Propertyplot", i, ".pdf", sep="")
  pdf(pdftitle)
  
  frame <- matrix(c(1,2,3,4,5,6,7,8), nrow=4, ncol=2, byrow=TRUE)
  layout(mat=frame, heights=rep(0.25, 4))
  par(oma=c(0,0,2.5,0))
  
  # par(mar=c(1,1,1,1))
  # plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  # legend(x="top", inset=0, legend="Blablablubb   ", bty="n", cex=1.5, horiz=TRUE)
  
  for(j in 1:4){
    dp <- (i - 1)*4 + j
    
    n <- nvec[dp]
    relpos <- relposmat[,dp]
    gamma <- gammavec[dp]
    R2 <- R2vec[dp]
    
    set.seed(seeds[1])
    sim <- simrel(n=(r*n), p=p, m=m, q=q, relpos=relpos, gamma=gamma, R2=R2)
    
    X <- sim$X[1:n, 1:p]  # Following the same seed-r philosophy as before,
    Y <- sim$Y[1:n]       # to retrieve identical datasets
    
    lambda <- sim$lambda
    covZY <- sim$Sigma[2:(p+1),1]
    
    # Scaling the true covariances between Z and Y
    covabs <- abs(covZY)
    maxcov <- max(covabs)
    covsc <- covabs/maxcov
    
    # Estimating the eigenvalues and eigenvectors
    Sxx <- (t(X) %*% X)/(n - 1)
    Ehat <- eigen(Sxx)$vectors
    lambdahat <- eigen(Sxx)$values
    
    # Scaling the estimated eigenvalues
    lambdahatsc <- lambdahat/lambdahat[1]
    
    # Computing and scaling the estimated covariances between Z and Y
    Z <- X %*% Ehat
    covhat <- cov(Z, Y)
    covhatabs <- abs(covhat)
    maxcovhat <- max(covhatabs)
    covhatsc <- covhatabs/maxcovhat
    
    # Plotting regions
    ymin <- -0.02
    ymax <- 1.1
    
    # Plots
    par(mar=c(2.5,2,1,2))
    barplot(lambda[1:k], space=26, col="black", ylim=c(ymin, ymax), yaxs="i", xlab="Component no.")
    # mtext("Eigenvalues", side=2, line=2, las=3)
    # mtext("Covariances", side=4, las=3)
    par(new=TRUE)
    plot(covsc[1:k], pch=21, col="dodgerblue", bg="dodgerblue", ylim=c(ymin, ymax), yaxs="i", xlab="", ylab="", cex=1.4)
    legend("topright", legend=paste("dp", dp, "  ", sep=""), cex=1.2)
    par(mar=c(2.5,3,1,0.5))
    barplot(lambdahatsc[1:k], space=26, col="black", ylim=c(ymin, ymax), yaxs="i", xlab="Component no.")
    # mtext("Estimated eigenvalues", side=2, line=2, las=3)
    # mtext("Estimated covariances", side=4, las=3)
    par(new=TRUE)
    plot(covhatsc[1:k], pch=21, col="dodgerblue", bg="dodgerblue", ylim=c(ymin, ymax), yaxs="i", yaxt="n", xlab="", ylab="", cex=1.4)
    legend("topright", legend=paste("dp", dp, "  ", sep=""), cex=1.2)
    mtext(text="True (left) and estimated (right) eigenvalues and covariances of components",
          line=0.4, outer=TRUE, cex=1.1)
  }
  dev.off()
}


