# Setting the current directory and loading/sourcing files
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/Estimates")
sigma2hatBayesPLS <- get(load("sigma2hatBayesPLS_full.RData"))
sigma2hatPCR <- get(load("sigma2hatPCR_full.RData"))
sigma2hatPLSRnaive <- get(load("sigma2hatPLSRnaive_full.RData"))
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")


# ---------------------------- #
#                              #
#      STACKING THE DATA       #
#   using squared error as Y   #
#                              #
# ---------------------------- #

# Varied factors
Y <- NULL
Method <- NULL
n <- NULL
relpos <- NULL
gamma <- NULL
R2 <- NULL
Comp <- NULL
Seed <- NULL
r <- NULL

methodLevels <- c("BayesPLS", "PCR", "PLSRnaive")
nLevels <- rep(c(rep(50,4), rep(15,4)),2)
relposLevels <- rep(c(rep("1,2,3", 2), rep("3,5,7", 2)), 4)
gammaLevels <- rep(c(0.9, 0.2),8)
R2Levels <- c(rep(0.7,8), rep(0.2,8))
kmax <-8

for(i in 1:length(methodLevels)){
  varname <- paste("sigma2hat", methodLevels[i], sep="")
  methodtemp <- methodLevels[i]
  
  for(dp in 1:16){
    ntemp <- nLevels[dp]
    relpostemp <- relposLevels[dp]
    gammatemp <- gammaLevels[dp]
    R2temp <- R2Levels[dp]
    
    for(k in 1:kmax){
      for(s in 1:7){
        Ytemp <- NULL
        
        for(rs in 1:3){
          Ytempnew <- get(varname)[[dp]][rs,s,k]
          
          if(!is.na(Ytempnew)){
            Ytempnew <- (Ytempnew - (1 - R2temp))^2  # Squared error of the estimate
            Ytemp <- c(Ytemp, Ytempnew)
            r <- c(r, rs)
          }
        }
        j <- length(Ytemp)
        
        Y <- c(Y, Ytemp)
        Method <- c(Method, rep(methodtemp, j))
        n <- c(n, rep(ntemp, j))
        relpos <- c(relpos, rep(relpostemp, j))
        gamma <- c(gamma, rep(gammatemp, j))
        R2 <- c(R2, rep(R2temp, j))
        Comp <- c(Comp, rep(k, j))
        Seed <- c(Seed, rep(s, j))
      }
    }
  }
}

sigma2hatStacked <- data.frame(Y, Method, n, relpos, gamma, R2, Comp, Seed, r)

sigma2hatStacked$n <- as.factor(sigma2hatStacked$n)
sigma2hatStacked$gamma <- as.factor(sigma2hatStacked$gamma)
sigma2hatStacked$R2 <- as.factor(sigma2hatStacked$R2)
sigma2hatStacked$Comp <- as.factor(sigma2hatStacked$Comp)

save(sigma2hatStacked, file="sigma2hatStacked_withoutKrylov.RData")


# ----------------------------------- #
#                                     #
#                ANOVA                #
#         Mixed effects model         #
#        r nested within seed         #
#                                     #
# ----------------------------------- #

setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/Estimates")
sigma2hatStacked <- get(load("sigma2hatStacked_withoutKrylov.RData"))
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")
library(lme4)


# ------------ #
#  Full model  #
# ------------ #
options(contrasts = c("contr.sum","contr.poly"))
lmermodel3 <- lmer(Y ~ (Method + Comp + n + relpos + gamma + R2)^6
                  + (1|Seed/r), data=sigma2hatStacked)

save(lmermodel3, file="lmermodel3.RData")

anmodel3 <- anova(lmermodel3, data=sigma2hatStacked)

# Degrees of freedom
dfR <- sum(anmodel3$Df)
dfT <- dim(sigma2hatStacked)[1] - 1
dfE <- dfT - dfR

options(scipen=999)

anmodel3$`Sum Sq` <- round(anmodel3$`Sum Sq`, digits=4)
anmodel3$`Mean Sq` <- round(anmodel3$`Mean Sq`, digits=4)
anmodel3$`F value` <- round(anmodel3$`F value`, digits=4)

`p value` <- 1 - pf(anmodel3$`F value`, anmodel3$Df, dfE)
`p value` <- round(`p value`, digits=4)
anmodel3 <- cbind(anmodel3, `p value`)

save(anmodel3, file="anovamodel3.RData")


# --------------- #
#  Reduced model  #
# --------------- #
options(contrasts = c("contr.sum","contr.poly"))
lmermodel4 <- lmer(Y ~ Method + Comp + n + relpos + gamma + R2 + Method:Comp + Method:n
                   + Method:relpos + Method:gamma + Method:R2 + Comp:n + Comp:relpos
                   + Comp:gamma + Comp:R2 + n:gamma + n:R2 + relpos:gamma + relpos:R2 + gamma:R2
                   + Method:Comp:n + Method:Comp:relpos + Method:Comp:gamma + Method:Comp:R2
                   + Method:n:gamma + Method:R2:n + Method:relpos:gamma + Method:relpos:R2
                   + Method:gamma:R2 + Comp:gamma:n + Comp:R2:n + Comp:relpos:gamma + Comp:relpos:R2
                   + Comp:gamma:R2 + n:gamma:R2 + relpos:gamma:R2 + Method:Comp:relpos:R2
                   + Method:Comp:gamma:R2 + Method:gamma:R2:n + Comp:relpos:gamma:R2
                   + (1|Seed/r), data=sigma2hatStacked)

lmermodel4summary <- summary(lmermodel4)
options(max.print=1200)
lmermodel4summary

save(lmermodel4, file="lmermodel4.RData")
save(lmermodel4summary, file="lmermodel4summary.RData")

anmodel4 <- anova(lmermodel4, data=sigma2hatStacked)

# Degrees of freedom
dfR <- sum(anmodel4$Df)
dfT <- dim(sigma2hatStacked)[1] - 1
dfE <- dfT - dfR

options(scipen=999)

anmodel4$`Sum Sq` <- round(anmodel4$`Sum Sq`, digits=4)
anmodel4$`Mean Sq` <- round(anmodel4$`Mean Sq`, digits=4)
anmodel4$`F value` <- round(anmodel4$`F value`, digits=4)

`p value` <- 1 - pf(anmodel4$`F value`, anmodel4$Df, dfE)
`p value` <- round(`p value`, digits=7)
anmodel4 <- cbind(anmodel4, `p value`)

save(anmodel4, file="anovamodel4.RData")


# ------------------------------- #
#                                 #
#      Plotting the effects       #
#   of the analysis of variance   #
#                                 #
# ------------------------------- #

lmermodel <- lmermodel4
library(effects, pos=21)


# --------------- #
#   Main effects  #
# --------------- #
eff <- Effect("R2", lmermodel)

pdf("eff_R2.pdf", height=5)
plot(eff, main=NULL, ylab="Estimation error")
dev.off()


# ------------------------ #
#     Grid of plots of     #
#  remaining main effects  #
# ------------------------ #
library(lattice)
library(grid)
library(gridExtra)
library(ggplot2)

latticeGrob <- function(p){ 
  grob(p=p, cl="lattice") 
} 

eff1 <- Effect("n", lmermodel)
eff2 <- Effect("relpos", lmermodel)
eff3 <- Effect("gamma", lmermodel)
eff4 <- Effect("R2", lmermodel)

p1 <- plot(eff1, main=NULL, ylab="Estimation error")
p2 <- plot(eff2, main=NULL, ylab="Estimation error")
p3 <- plot(eff3, main=NULL, ylab="Estimation error")
p4 <- plot(eff4, main=NULL, ylab="Estimation error")

layout <- matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)

pdf("effmain_n_relpos_gamma_R2.pdf", height=6)
grid.arrange(latticeGrob(p1),
             latticeGrob(p2),
             latticeGrob(p3),
             latticeGrob(p4),
             layout_matrix=layout)
dev.off()


# -------------------------- #
#   Remaining higher order   #
#  interaction effect plots  #
# -------------------------- #

whichEffects <- c(5,6)  # 1 = Method, 2 = Comp, 3 = n, 4 = relpos, 5 = gamma, 6 = R2

# Factors and number of levels
factors <- c("Method", "Comp", "n", "relpos", "gamma", "R2")
factorLevels <- c(3, 8, 2, 2, 2, 2)

nr <- 1
nc <- 2

eff <- Effect(factors[whichEffects], lmermodel)
pdftitle <- paste("eff_", paste(factors[whichEffects], collapse='_'), ".pdf", sep="")

pdf(pdftitle, height=4)
plot(eff, layout=c(nc,nr), main=NULL, ylab="Estimation error")
dev.off()


