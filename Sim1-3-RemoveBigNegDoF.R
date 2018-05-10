# Setting the current directory and loading/sourcing files
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/Estimates")
load("sigma2hatPLSRkrylov_full.RData")
load("DoFkrylov.RData")
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")

DoFkrylovCount <- vector("list", 16)
for(dp in 1:16){
  if(dp %in% seq(from=1, to=4) | dp %in% seq(from=9, to=12)){
    bound <- 26
  }else{
    bound <- 14
  }
  upperboundDoF <- NULL
  negativeDoF <- NULL
  for(k in 1:8){
    big <- 0
    neg <- 0
    for(s in 1:7){
      for(r in 1:3){
        if(DoFkrylov[[dp]][r,s,k] == bound){
          sigma2hatPLSRkrylov[[dp]][r,s,k] <- NA
          big <- big + 1
        }
        if(DoFkrylov[[dp]][r,s,k] < 0){
          sigma2hatPLSRkrylov[[dp]][r,s,k] <- NA
          neg <- neg + 1
        }
      }
    }
    upperboundDoF <- c(upperboundDoF, big)
    negativeDoF <- c(negativeDoF, neg)
  }
  DoFkrylovCount[[dp]] <- rbind(negativeDoF, upperboundDoF)
  names(DoFkrylovCount)[[dp]] <- paste("dp", dp, sep="")
}

save(sigma2hatPLSRkrylov, file="sigma2hatPLSRkrylov_red.RData")
save(DoFkrylovCount, file="DoFkrylovCount.RData")





