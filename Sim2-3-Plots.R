# Setting the current directory and loading/sourcing files
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/Estimates_sim2")
load("EstimatesSummary_sim2.RData")
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")


# ------------------------ #
#                          #
#   Plotting the summary   #
#     of the estimates     #
#       from sim. 2        #
#                          #
# -----------------------  #

# nvec <- rep(c(rep(50,4), rep(15,4)),2)
relposmat <- rbind(rep(c(rep(1,2), rep(3,2)),4), rep(c(rep(2,2), rep(5,2)),4), rep(c(rep(3,2), rep(7,2)),4))
gammavec <- rep(c(0.9, 0.2),8)
R2vec <- c(rep(0.7,8), rep(0.2,8))
sigma2 <- 1 - R2vec

# Labels
methods <- c("PCR", "PLSRnaive")
colors <- c("red", "blue")


# -------------- #
# Plotting means #
# -------------- #

# Determining plotting regions for dp1-4
ymin1 <- min(estSummary[[1]][,,1], estSummary[[2]][,,1],
             estSummary[[3]][,,1], estSummary[[4]][,,1])
ymax1 <- max(estSummary[[1]][,,1], estSummary[[2]][,,1],
             estSummary[[3]][,,1], estSummary[[4]][,,1])
            
# Determining plotting regions for dp9-12
ymin2 <- min(estSummary[[5]][,,1], estSummary[[6]][,,1],
             estSummary[[7]][,,1], estSummary[[8]][,,1])
ymax2 <- max(estSummary[[5]][,,1], estSummary[[6]][,,1],
             estSummary[[7]][,,1], estSummary[[8]][,,1])

# Initializing a pdf
pdfname = "PCR-PLSRbias.pdf"
pdf(pdfname, height=7.5)

frame <- matrix(c(1,2,3,4,5,6,7,8,9,9), nrow=5, ncol=2, byrow=TRUE)
layout(mat=frame, heights=c(0.233,0.233,0.233,0.233,0.068))

# Plotting
for(dp in 1:4){
  par(mar=c(4,7,3,4))
  str1 <- paste("relpos = ", relposmat[1,dp], "-", relposmat[2,dp],
                "-", relposmat[3,dp], ",", sep="")
  str2 <- ", "
  plot(estSummary[[dp]][1,,1], ylim=c(ymin1, ymax1), type="n",
       ylab=expression(bar(hat(sigma)^2)), xlab="number of components",
       main=bquote(.(str1) ~ gamma == .(gammavec[dp]) ~ .(str2) ~ R^2 == .(R2vec[dp])))
  grid(col="gray")
  abline(h=sigma2[dp], col="gray26", lty=2, cex=0.8)
  par(new=TRUE)
  plot(estSummary[[dp]][1,,1], ylim=c(ymin1, ymax1), type="l", col=colors[1], ylab="", xlab="")
  par(new=TRUE)
  plot(estSummary[[dp]][2,,1], ylim=c(ymin1, ymax1), type="l", col=colors[2], ylab="", xlab="")
}

for(dp in 9:12){
  par(mar=c(4,7,3,4))
  str1 <- paste("relpos = ", relposmat[1,dp], "-", relposmat[2,dp],
                "-", relposmat[3,dp], ",", sep="")
  str2 <- ", "
  plot(estSummary[[dp-4]][1,,1], ylim=c(ymin2, ymax2), type="n",
       ylab=expression(bar(hat(sigma)^2)), xlab="number of components",
       main=bquote(.(str1) ~ gamma == .(gammavec[dp]) ~ .(str2) ~ R^2 == .(R2vec[dp])))
  grid(col="gray")
  abline(h=sigma2[dp], col="gray26", lty=2, cex=0.8)
  par(new=TRUE)
  plot(estSummary[[dp-4]][1,,1], ylim=c(ymin2, ymax2), type="l", col=colors[1], ylab="", xlab="")
  par(new=TRUE)
  plot(estSummary[[dp-4]][2,,1], ylim=c(ymin2, ymax2), type="l", col=colors[2], ylab="", xlab="")
}
par(mar=c(1,1,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="top", inset=0, legend=methods, col=colors, lwd=1.7, cex=0.9, horiz=TRUE)

dev.off()

# ------------- #
# Plotting RMSE #
# ------------- #

# Determining plotting regions for dp1-4
ymin1 <- min(estSummary[[1]][,,3], estSummary[[2]][,,3],
             estSummary[[3]][,,3], estSummary[[4]][,,3])
ymax1 <- max(estSummary[[1]][,,3], estSummary[[2]][,,3],
             estSummary[[3]][,,3], estSummary[[4]][,,3])

# Determining plotting regions for dp9-12
ymin2 <- min(estSummary[[5]][,,3], estSummary[[6]][,,3],
             estSummary[[7]][,,3], estSummary[[8]][,,3])
ymax2 <- max(estSummary[[5]][,,3], estSummary[[6]][,,3],
             estSummary[[7]][,,3], estSummary[[8]][,,3])

# Initializing a pdf (plotting means)
pdfname = "PCR-PLSRrmse.pdf"
pdf(pdfname, height=7.5)

frame <- matrix(c(1,2,3,4,5,6,7,8,9,9), nrow=5, ncol=2, byrow=TRUE)
layout(mat=frame, heights=c(0.233,0.233,0.233,0.233,0.068))

# Plotting
for(dp in 1:4){
  par(mar=c(4,7,3,4))
  str1 <- paste("relpos = ", relposmat[1,dp], "-", relposmat[2,dp],
                "-", relposmat[3,dp], ",", sep="")
  str2 <- ", "
  plot(estSummary[[dp]][1,,1], ylim=c(ymin1, ymax1), type="n",
       ylab=expression(paste("RMSE(", hat(sigma)^2, ")")), xlab="number of components",
       main=bquote(.(str1) ~ gamma == .(gammavec[dp]) ~ .(str2) ~ R^2 == .(R2vec[dp])))
  grid(col="gray")
  par(new=TRUE)
  plot(estSummary[[dp]][1,,3], ylim=c(ymin1, ymax1), type="l", col=colors[1], ylab="", xlab="")
  par(new=TRUE)
  plot(estSummary[[dp]][2,,3], ylim=c(ymin1, ymax1), type="l", col=colors[2], ylab="", xlab="")
}

for(dp in 9:12){
  par(mar=c(4,7,3,4))
  str1 <- paste("relpos = ", relposmat[1,dp], "-", relposmat[2,dp],
                "-", relposmat[3,dp], ",", sep="")
  str2 <- ", "
  plot(estSummary[[dp-4]][1,,1], ylim=c(ymin2, ymax2), type="n",
       ylab=expression(paste("RMSE(", hat(sigma)^2, ")")), xlab="number of components",
       main=bquote(.(str1) ~ gamma == .(gammavec[dp]) ~ .(str2) ~ R^2 == .(R2vec[dp])))
  grid(col="gray")
  par(new=TRUE)
  plot(estSummary[[dp-4]][1,,3], ylim=c(ymin2, ymax2), type="l", col=colors[1], ylab="", xlab="")
  par(new=TRUE)
  plot(estSummary[[dp-4]][2,,3], ylim=c(ymin2, ymax2), type="l", col=colors[2], ylab="", xlab="")
}
par(mar=c(1,1,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="top", inset=0, legend=methods, col=colors, lwd=1.7, cex=0.9, horiz=TRUE)

dev.off()





