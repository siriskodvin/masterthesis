# Setting the current directory and loading/sourcing files
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim/Estimates")
estSummary <- get(load("EstimatesSummary.RData"))
setwd("C:/Users/siris/OneDrive/NMBU/Masteroppgave/R kode/Simulation/Simulation 13 Sigma2hatsim")
source("Sim1-5a-PlotsFunction.R")

# In this simulation 16 different datapoints (dp's) are simulated. Four key properties of the data are
# varied between a high and a low level:
# n = 50 - 15
# relpos = 1,2,3 - 3,5,7
# gamma = 0.9 - 0.2
# R2 = 0.7 - 0.2
# Below is an overview of the levels of the factors of each dp.

#        ------------------------------
#        | n  |  relpos | gamma | R2  |
# |------|----|---------|-------|-----|
# | dp1  | 50 |  1,2,3  |  0.9  | 0.7 |
# |------|----|---------|-------|-----|
# | dp2  | 50 |  1,2,3  |  0.2  | 0.7 |
# |------|----|---------|-------|-----|
# | dp3  | 50 |  3,5,7  |  0.9  | 0.7 |
# |------|----|---------|-------|-----|
# | dp4  | 50 |  3,5,7  |  0.2  | 0.7 |
# |------|----|---------|-------|-----|
# | dp5  | 15 |  1,2,3  |  0.9  | 0.7 |
# |------|----|---------|-------|-----|
# | dp6  | 15 |  1,2,3  |  0.2  | 0.7 |
# |------|----|---------|-------|-----|
# | dp7  | 15 |  3,5,7  |  0.9  | 0.7 |
# |------|----|---------|-------|-----|
# | dp8  | 15 |  3,5,7  |  0.2  | 0.7 |
# |------|----|---------|-------|-----|
# | dp9  | 50 |  1,2,3  |  0.9  | 0.2 |
# |------|----|---------|-------|-----|
# | dp10 | 50 |  1,2,3  |  0.2  | 0.2 |
# |------|----|---------|-------|-----|
# | dp11 | 50 |  3,5,7  |  0.9  | 0.2 |
# |------|----|---------|-------|-----|
# | dp12 | 50 |  3,5,7  |  0.2  | 0.2 |
# |------|----|---------|-------|-----|
# | dp13 | 15 |  1,2,3  |  0.9  | 0.2 |
# |------|----|---------|-------|-----|
# | dp14 | 15 |  1,2,3  |  0.2  | 0.2 |
# |------|----|---------|-------|-----|
# | dp15 | 15 |  3,5,7  |  0.9  | 0.2 |
# |------|----|---------|-------|-----|
# | dp16 | 15 |  3,5,7  |  0.2  | 0.2 |
# |------|----|---------|-------|-----|

# Function call
# The function creates a pdf with 4 plots.
SimulationPlots(estSummary = estSummary,   # Loaded data
                dps = c(7,8,15,16),          # dp numbers (see table above)
                measNum = 4,               # Which measure to plot: # 1 = Means, 2 = Standardized means,
                                                                    # 3 = MSE, 4 = RMSE, 5 = SE, 6 = Bias
                colors = c("magenta", "red", "green", "blue"))

