# ------------------------ #
#                          #
#  PCR estimates function  #
#                          #
# ------------------------ #
sigma2hatPCRfunc <- function(n, p, m, q, relpos, gamma, R2, r, seeds, k){
  sigma2hatPCR <- array(rep(0,(r*length(seeds)*k)), dim=c(r, length(seeds), k))
  dimnames(sigma2hatPCR)[1] <- list(paste("r", 1:r, sep=""))
  dimnames(sigma2hatPCR)[2] <- list(paste("Seed", 1:(length(seeds))))
  dimnames(sigma2hatPCR)[3] <- list(paste("comp", 1:k))
  
  for(a in 1:k){
    for(b in 1:length(seeds)){
      set.seed(seeds[b])
      sim <- simrel(n=(r*n), p=p, m=m, q=q, relpos=relpos, gamma=gamma, R2=R2)
      
      for(j in 1:r){
        X <- sim$X[((j-1)*n+1):(j*n), 1:p]
        Y <- sim$Y[((j-1)*n+1):(j*n)]
        
        # Center the data
        Xc <- scale(X, scale=FALSE)
        Yc <- scale(Y, scale=FALSE)
        
        # Eigenmatrix of XtX
        E <- eigen((t(Xc) %*% Xc)/(n - 1))$vectors
        Ek <- E[,1:a]
        
        # Principal component scores
        Z <- Xc %*% Ek
        
        # Estimated coefficients
        alfahat <- solve(t(Z) %*% Z) %*% t(Z) %*% Yc
        betahatPCR <- Ek %*% alfahat
        
        sigma2hat <- sum((Yc - (Xc %*% betahatPCR))^2)/(n - (a + 1))

        sigma2hatPCR[j, b, a] <- sigma2hat
      }
    }
  }
  return(sigma2hatPCR)
}

# -------------------------------------------- #
#                                              #
#  PLSR (naive and Krämer) estimates function  #
#              (also returns DoF)              #
#                                              #
# -------------------------------------------- #
sigma2hatPLSRfunc <- function(n, p, m, q, relpos, gamma, R2, r, seeds, k){
  sigma2hatPLSRnaive <- array(rep(0,(r*length(seeds)*k)), dim=c(r, length(seeds), k))
  dimnames(sigma2hatPLSRnaive)[1] <- list(paste("r", 1:r, sep=""))
  dimnames(sigma2hatPLSRnaive)[2] <- list(paste("Seed", 1:(length(seeds))))
  dimnames(sigma2hatPLSRnaive)[3] <- list(paste("comp", 1:k))
  sigma2hatPLSRkrylov <- array(rep(0,(r*length(seeds)*k)), dim=c(r, length(seeds), k))
  dimnames(sigma2hatPLSRkrylov)[1] <- list(paste("r", 1:r, sep=""))
  dimnames(sigma2hatPLSRkrylov)[2] <- list(paste("Seed", 1:(length(seeds))))
  dimnames(sigma2hatPLSRkrylov)[3] <- list(paste("comp", 1:k))
  DoFkrylov <- array(rep(0,(r*length(seeds)*k)), dim=c(r, length(seeds), k))
  dimnames(DoFkrylov)[1] <- list(paste("r", 1:r, sep=""))
  dimnames(DoFkrylov)[2] <- list(paste("Seed", 1:(length(seeds))))
  dimnames(DoFkrylov)[3] <- list(paste("comp", 1:k))
  
  for(b in 1:length(seeds)){
    set.seed(seeds[b])
    sim <- simrel(n=(r*n), p=p, m=m, q=q, relpos=relpos, gamma=gamma, R2=R2)
    
    for(j in 1:r){
      X <- sim$X[((j-1)*n+1):(j*n), 1:p]
      Y <- sim$Y[((j-1)*n+1):(j*n)]
      
      plsmodel <- linear.pls.fit(X, Y, m=k)
      K <- X %*% t(X)   # kernel matrix
      dofobj <- pls.dof(plsmodel, n, Y, K, m=k, DoF.max=min(ncol(X)+1,nrow(X)-1))
      
      for(a in 1:k){
        sigma2hatPLSRnaive[j, b, a] <- (plsmodel$sigmahat[a+1])^2
        sigma2hatPLSRkrylov[j, b, a] <- (dofobj$sigmahat[a])^2
        DoFkrylov[j, b, a] <- dofobj$DoF[a]
      }
    }
  }
  return(list(sigma2hatPLSRnaive, sigma2hatPLSRkrylov, DoFkrylov))
}