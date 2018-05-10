# ------------------------ #
#                          #
#   Plotting the summary   #
#     of the estimates     #
#                          #
# -----------------------  #

# For details, se the function call file "Plots.R".

SimulationPlots <- function(estSummary, dps, measNum, colors){
  # Information about the datasets (dp's)
  nvec <- rep(c(rep(50,4), rep(15,4)),2)
  relposmat <- rbind(rep(c(rep(1,2), rep(3,2)),4), rep(c(rep(2,2), rep(5,2)),4), rep(c(rep(3,2), rep(7,2)),4))
  gammavec <- rep(c(0.9, 0.2),8)
  R2vec <- c(rep(0.7,8), rep(0.2,8))
  sigma2 <- 1 - R2vec
  
  # Labels
  methods <- c("BayesPLS", "PCR", "PLSRkrylov", "PLSRnaive")
  measures <- c("Means", "STMeans", "MSE", "RMSE", "SE", "Bias")
  ylabMean <- expression(bar(hat(sigma)^2))
  ylabSTMean <- expression(paste(bar(hat(sigma)^2), "/", sigma^2))
  ylabMSE <- expression(paste("MSE(", hat(sigma)^2, ")"))
  ylabRMSE <- expression(paste("RMSE(", hat(sigma)^2, ")"))
  ylabSD <- expression(paste(hat("SE"), "(", hat(sigma)^2, ")"))
  ylabBias <- expression(paste(hat("Bias"), "(", hat(sigma)^2, ")"))
  ylabs <- c(ylabMean, ylabSTMean, ylabMSE, ylabRMSE, ylabSD, ylabBias)
  mainMean <- "Average of estimates"
  mainSTmean <- "Standardized average of estimates"
  mainMSE <- "Mean square error of estimates"
  mainRMSE <- "Root mean square error of estimates"
  mainSD <- "Standard error of estimates"
  mainBias <- "Bias of estimates"
  mains <- c(mainMean, mainSTmean, mainMSE, mainRMSE, mainSD, mainBias)
  
  # Determining plotting regions
  ymin <- min(estSummary[[dps[1]]][,,measNum], estSummary[[dps[2]]][,,measNum],
              estSummary[[dps[3]]][,,measNum], estSummary[[dps[4]]][,,measNum])
  ymax <- max(estSummary[[dps[1]]][,,measNum], estSummary[[dps[2]]][,,measNum],
              estSummary[[dps[3]]][,,measNum], estSummary[[dps[4]]][,,measNum])
  sigmasubs <- c(sigma2[dps[1]], sigma2[dps[2]], sigma2[dps[3]], sigma2[dps[4]])
  
  # Initializing a pdf
  pdfname = paste("Simulation", dps[1], "-", dps[2], "-", dps[3], "-", dps[4],
                  measures[measNum], ".pdf", sep="")
  pdf(pdfname, height=7.5)
  
  frame <- matrix(c(1,2,3,4,5,5), nrow=3, ncol=2, byrow=TRUE)
  layout(mat=frame, heights=c(0.46,0.46,0.08))
  par(oma=c(0,0,2.5,0))
  
  # Plotting
  for(j in 1:4){
    par(mar=c(5,5,4,2))
    dp <- dps[j]
    str1 <- paste("n = ", nvec[dp], ", p = 25, ", sep="")
    str2 <- paste("relpos = ", relposmat[1,dp], "-", relposmat[2,dp],
                  "-", relposmat[3,dp], ",", sep="")
    plot(estSummary[[dp]][1,,measNum], ylim=c(ymin, ymax), type="n",
         ylab=ylabs[measNum], xlab="number of components",
         main=bquote(atop(.(str1) ~ R^2 == .(R2vec[dp]), .(str2) ~ gamma == .(gammavec[dp]))))
    grid(col="gray")
    if(measNum == 1){
      abline(h=sigmasubs[j], col="gray26", lty=2, cex=0.8)
    }else if(measNum == 2){
      abline(h=1, col="gray26", lty=2, cex=0.8)
    }else if(measNum == 5){
      abline(h=0, col="gray26", lty=2, cex=0.8)
    }
    par(new=TRUE)
    plot(estSummary[[dp]][1,,measNum], ylim=c(ymin, ymax), type="l", col=colors[1], ylab="", xlab="")
    par(new=TRUE)
    plot(estSummary[[dp]][2,,measNum], ylim=c(ymin, ymax), type="l", col=colors[2], ylab="", xlab="")
    par(new=TRUE)
    plot(estSummary[[dp]][3,,measNum], ylim=c(ymin, ymax), type="l", col=colors[3], ylab="", xlab="")
    par(new=TRUE)
    plot(estSummary[[dp]][4,,measNum], ylim=c(ymin, ymax), type="l", col=colors[4], ylab="", xlab="")
  }
  mtext(mains[measNum], line=0.4, outer=TRUE, cex=1.1)
  
  par(mar=c(1,1,1,1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x="top", inset=0, legend=methods, col=colors, lwd=1.7, cex=0.9, horiz=TRUE)
  
  dev.off()
}


