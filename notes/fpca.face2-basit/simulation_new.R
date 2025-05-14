library(refund)
library(splines)
library(mgcv)
library(Matrix)

source("GeneData.R")
source("fpca.face2.R")
source("simulation_functions.R")

n <- rep(c(100, 200, 500), 6)
J <- rep(rep(c(100, 1000, 10000), each = 3), 2)
design <- rep(c("regular", "irregular"), each = 9)


## Equally Spaced Data

sim_final_results(
  Nel = 100, K = 3, n = n, J = J, sigma = 0.5, equally = T, design = design,
  level = 0.5, parallel = T
)

## Unequally Spaced Data

sim_final_results(
  Nel = 100, K = 3, n = n, J = J, sigma = 0.5, equally = F, design = design,
  level = 0.5, parallel = T
)