################################################################################
######                           Data generation                          ######
################################################################################

# Setting working directory (where functions and Simulation.seed are saved)
setwd("") 

#---------------------- Needed Packages and Functions --------------------------
library(foreach)
library(doParallel)

source("0_CostumeFunctions_DataGeneration.R")

#-------------------- Same parameter in all scenarios  -------------------------
# Number of used cores for Simulations
num.cl <- 48

# Number of iterations for each sub-scenario
n_sim <- 10000 

# Seed for simulation
load("Simulation.seed.RData")
#seed = sample(x = 1:1000000000, size = it, replace = FALSE)
#save(seed, file = "Simulation.seed.RData")

# Path for saving generated data:
path_data <- "" # e.g. ".../Simulations/Data/"


#------------------------- Standard Scenario 1 ---------------------------------
#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- c(0.2, 0.4, 0.6) 
accrual.time <- 24 # months --> 2 years
fu.time      <- NA # will set to 2*med.C later


# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, 
                          med.C, HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,8],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", "r", 
                          "med.C", "HR.var", "HR", "n_sim", "n.T", "n.C")

# Sorting parameter matrix after the sample size
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters[,10],decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters, file = paste0(path_data, "Scen1_parameters.rds"))


#### Data generation
# Starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival")) %dopar% {
                     # Generating data "n_sim" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                        data.scen.i[[i]] <- 
                           DataEXP(
                              n.T            = para$n.T, 
                              n.C            = para$n.C,
                              median.control = para$med.C,
                              accrual.time   = para$accrual.time,
                              fu.time        = para$fu.time,
                              HR             = para$HR, 
                              HR.var         = para$HR.var, 
                              cens.rate      = para$cens.rate,
                              index.seed     = seed[i])
                     }
                     
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data, 
                                           "Scen1_Data", 
                                           "..beta.",        para$beta,
                                           ".accrual.time.", para$accrual.time,
                                           ".fu.time.",      para$fu.time,
                                           ".cens.rate.",    para$cens.rate,
                                           ".r.",            para$r,
                                           ".med.C.",        para$med.C,
                                           ".HR.var.",       para$HR.var,
                                           ".HR.",           para$HR,
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)

#------------------------------ Scenario 2 -------------------------------------
#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- c(0.8, 0.9, 1.1, 1.2)
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- c(0.2, 0.4, 0.6) 
accrual.time <- 24 # in months --> 2 years
fu.time      <- NA # will set to 2*med.C later

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, med.C, 
                          HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,8],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", "r", 
                          "med.C", "HR.var", "HR", "n_sim", "n.T", "n.C")

# Sorting parameter matrix after the sample size
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters[,10],decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path_data, "Scen2_parameters.rds"))


#### Data generation
# Starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival")) %dopar% {
                     # Generating data "n_sim" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                        data.scen.i[[i]] <- 
                           DataEXP(
                              n.T            = para$n.T, 
                              n.C            = para$n.C,
                              median.control = para$med.C,
                              accrual.time   = para$accrual.time,
                              fu.time        = para$fu.time,
                              HR             = para$HR, 
                              HR.var         = para$HR.var, 
                              cens.rate      = para$cens.rate,
                              index.seed     = seed[i])
                     }
                     
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data,
                                           "Scen2_Data",
                                           "..beta.",        para$beta,
                                           ".accrual.time.", para$accrual.time,
                                           ".fu.time.",      para$fu.time,
                                           ".cens.rate.",    para$cens.rate,
                                           ".r.",            para$r,
                                           ".med.C.",        para$med.C,
                                           ".HR.var.",       para$HR.var,
                                           ".HR.",           para$HR,
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)


#------------------------- Scenario 3a: WEIBULL --------------------------------
#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- c(0.2, 0.4, 0.6) 
shape.C.T    <- c(0.5, 1.5)   
accrual.time <- 24 # in months --> 2 years
fu.time      <- NA # will set to 2*med.C later


# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, med.C, 
                          shape.C.T, HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,9],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", "r", 
                          "med.C", "shape.C.T", "HR.var", "HR", "n_sim", "n.T", 
                          "n.C")

# Sorting parameter matrix after the sample size
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters[,11],decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path_data, "Scen3aWEIB_parameters.rds"))


#### Data generation
# Starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival")) %dopar% {
                     # Generating data "n_sim" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                        data.scen.i[[i]] <- 
                           DataWEIB(
                              n.T            = para$n.T, 
                              n.C            = para$n.C,
                              median.control = para$med.C,
                              accrual.time   = para$accrual.time,
                              fu.time        = para$fu.time,
                              HR             = para$HR, 
                              HR.var         = para$HR.var, 
                              cens.rate      = para$cens.rate,
                              shape.C.T      = para$shape.C.T,
                              index.seed     = seed[i])
                     }
                     
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data, 
                                           "Scen3aWEIB_Data",
                                           "..beta.",     para$beta,
                                           ".accrual.time.", para$accrual.time,
                                           ".fu.time.",      para$fu.time,
                                           ".cens.rate.", para$cens.rate,
                                           ".r.",         para$r,
                                           ".med.C.",     para$med.C,
                                           ".shape.C.T.", para$shape.C.T,
                                           ".HR.var.",    para$HR.var,
                                           ".HR.",        para$HR,
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)


#------------------------- Scenario 3b: Gompertz -------------------------------
#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- c(0.2, 0.4, 0.6) 
shape.C.T    <- c(-0.2, 0.2)   
accrual.time <- 24 # in months --> 2 years
fu.time      <- NA # will set to 2*med.C later


# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, med.C, 
                          shape.C.T, HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,9],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", "r", 
                          "med.C", "shape.C.T",  "HR.var", "HR", "n_sim", "n.T", 
                          "n.C")

# Sorting parameter matrix after the sample size
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters[,11],decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path_data, "Scen3bGOMP_parameters.rds"))


#### Data generation
# Starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival", "flexsurv")) %dopar% {
                     # Generating data "it" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                        data.scen.i[[i]] <- 
                           DataGomp(
                              n.T            = para$n.T, 
                              n.C            = para$n.C,
                              median.control = para$med.C,
                              accrual.time   = para$accrual.time,
                              fu.time        = para$fu.time,
                              HR             = para$HR, 
                              HR.var         = para$HR.var, 
                              cens.rate      = para$cens.rate,
                              shape.C.T      = para$shape.C.T,
                              index.seed     = seed[i])
                     }
                     
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data, 
                                           "Scen3bGOMP_Data",
                                           "..beta.",        para$beta,
                                           ".accrual.time.", para$accrual.time,
                                           ".fu.time.",      para$fu.time,
                                           ".cens.rate.",    para$cens.rate,
                                           ".r.",            para$r,
                                           ".med.C.",        para$med.C,
                                           ".shape.C.T.",    para$shape.C.T,
                                           ".HR.var.",       para$HR.var,
                                           ".HR.",           para$HR,
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 4 -------------------------------------
# treat.effect.start is set to 1/3*med.C due to ESMOs disadvantage of using 
# "gain"(med.T-med.C)

#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- c(0.2, 0.4, 0.6) 
accrual.time <- 24 # in months --> 2 years
fu.time      <- NA # will set to 2*med.C later

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, med.C, 
                          HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,8],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", 
                          "r", "med.C", "HR.var", "HR", "n_sim", "n.T", "n.C")

# Sorting parameter matrix after the sample size 
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters[,10],decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path_data, "Scen4_parameters.rds"))

#### Data generation
# starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl)
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival")) %dopar% {
                     # Generating data "it" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                        data.scen.i[[i]] <- 
                           DataEXPNonProp(
                              n.T            = para$n.T, 
                              n.C            = para$n.C,
                              median.control = para$med.C,
                              accrual.time   = para$accrual.time,
                              fu.time        = para$fu.time,
                              HR             = para$HR, 
                              HR.var         = para$HR.var, 
                              effect.start.T = 1/3 * para$med.C,
                              cens.rate      = para$cens.rate,
                              index.seed     = seed[i])
                     }
                     
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data,
                                           "Scen4_Data",
                                           "..beta.",        para$beta,
                                           ".accrual.time.", para$accrual.time,
                                           ".fu.time.",      para$fu.time,
                                           ".cens.rate.",    para$cens.rate,
                                           ".r.",            para$r,
                                           ".med.C.",        para$med.C,
                                           ".HR.var.",       para$HR.var,
                                           ".HR.",           para$HR,
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)
