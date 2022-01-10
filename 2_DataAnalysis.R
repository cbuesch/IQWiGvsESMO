################################################################################
#########            Analysis of the generated data                    #########                 
################################################################################

# Setting working directory (where functions and Simulation.seed are saved)
setwd("") 

#---------------------- Needed Packages and Functions --------------------------
library(tidyverse)
library(foreach)
library(doParallel)

source("0_CostumeFunctions_Analysis.R")

#-------------------- Same parameter in all scenarios  -------------------------
# Number of used cores for analysis
num.cl <- 4

# Loading / Saving paths
path_data <- "" # e.g. ".../Simulations/Data/"
path_ana  <- "" # e.g. ".../Simulations/Ana/"


#------------------------- Standard Scenario 1 ---------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen1_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results_Scen1 <- foreach(para = iter(parameters, by='row'), 
                        .packages = "survival", .combine = rbind) %dopar% {
   # Load data of sub-scenarios with "n_sim" times replicates
   df.loaded <- readRDS(paste0(path_data,
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
   
   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, 
                                                  ncol = dim(ana.data.scen.1)[2],
                                                  nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)
   
   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Perform summary of the analyzed data and return it
   return(
      data.frame(
         Scenario = "Scen_1",
         para,
         analysis.summary(df.raw = ana.data.scen.i, 
                          n_sim = para$n_sim)
      )
   )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 2 -------------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen2_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results_Scen2 <- foreach(para = iter(parameters, by='row'), 
                  .packages = "survival", .combine = rbind) %dopar% {
   # load data of sub-scenarios with "n_times" times replicates
   df.loaded <- readRDS(paste0(path_data,
                               "Scen2_Data",
                               "..beta.",     para$beta,
                               ".accrual.time.", para$accrual.time,
                               ".fu.time.",      para$fu.time,
                               ".cens.rate.", para$cens.rate,
                               ".r.",         para$r,
                               ".med.C.",     para$med.C,
                               ".HR.var.",    para$HR.var,
                               ".HR.",        para$HR,
                               ".rds"))

   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, 
                                                  ncol = dim(ana.data.scen.1)[2], 
                                                  nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)

   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Perform summary of the analyzed data and return it
   return(
      data.frame(
         Scenario = "Scen_2",
         para,
         analysis.summary(df.raw = ana.data.scen.i, 
                          n_sim = para$n_sim)
      )
   )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)


#------------------------- Scenario 3a: WEIBULL --------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen3aWEIB_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl)
registerDoParallel(cl)

# Starting parallel computing
results_Scen3a <- foreach(para = iter(parameters, by='row'),
                        .packages = "survival", .combine = rbind) %dopar% {
   # load data of sub-scenarios with "n_sim" times replicates
   df.loaded <- readRDS(paste0(path_data,
                               "Scen3aWEIB_Data",
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

   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, ncol = dim(ana.data.scen.1)[2], nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)

   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Perform summary of the analyzed data and return it
   return(
      data.frame(
         Scenario = "Scen_3aWEIB",
         para,
         analysis.summary(df.raw = ana.data.scen.i, 
                          n_sim = para$n_sim)
      )
   )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)


#------------------------- Scenario 3b: GOMPERTZ -------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen3bGOMP_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# setting up parallel computing
cl <- makeCluster(num.cl)
registerDoParallel(cl)

# Starting parallel computing
results_Scen3b <- foreach(para = iter(parameters, by='row'), 
                         .packages = "survival", .combine = rbind) %dopar% {
   # load data of sub-scenarios with "n_sim" times replicates
   df.loaded <- readRDS(paste0(path_data,
                               "Scen3bGOMP_Data",
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

   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, ncol = dim(ana.data.scen.1)[2], nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)

   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Perform summary of the analyzed data and return it
   return(
      data.frame(
         Scenario = "Scen_3bGOMP",
         para,
         analysis.summary(df.raw = ana.data.scen.i, 
                          n_sim = para$n_sim)
      )
   )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 4 -------------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen4_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results_Scen4 <- foreach(para = iter(parameters, by='row'), 
                  .packages = "survival", .combine = rbind) %dopar% {
   # load data of sub-scenarios with "n_sim" times replicates
   df.loaded <- readRDS(paste0(path_data,
                               "Scen4_Data",
                               "..beta.",     para$beta,
                               ".accrual.time.", para$accrual.time,
                               ".fu.time.",      para$fu.time,
                               ".cens.rate.", para$cens.rate,
                               ".r.",         para$r,
                               ".med.C.",     para$med.C,
                               ".HR.var.",    para$HR.var,
                               ".HR.",        para$HR,
                               ".rds"))

   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, ncol = dim(ana.data.scen.1)[2], nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)

   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Perform summary of the analyzed data and return it
   return(
      data.frame(
         Scenario = "Scen_4",
         para,
         analysis.summary(df.raw = ana.data.scen.i, 
                          n_sim = para$n_sim)
      )
   )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)


#-------------- Combining of every analyzed data and saving --------------------
# Combining of every analyzed data
results <- rbind(results_Scen1 %>% mutate(shape.C.T = NA) %>% 
                    relocate(names(results_Scen1)[1:7], shape.C.T, everything()),
                 results_Scen2 %>% mutate(shape.C.T = NA) %>% 
                    relocate(names(results_Scen2)[1:7], shape.C.T, everything()),
                 results_Scen3a, results_Scen3b, 
                 results_Scen4 %>% mutate(shape.C.T = NA) %>% 
                    relocate(names(results_Scen4)[1:7], shape.C.T, everything())
            )
# Saving
saveRDS(results, file = paste0(path_ana, "Results.rds"))
