################################################################################
######                           ROC curves                               ######
################################################################################
library(tidyverse)
library(ggpubr)
library(survival)

## Thresholds for assumed true added benefit
delta_deserved <- c(0.5, 0.6, 0.7, 0.8)

## Different thresholds for different threshold methods
## (Point estimate, upper CI, lower CI)
HR_thresh_methods <- seq(0.2,1, by = 0.01)

#------------------------------ Standard Scenario 1 -------------------------------------
## Working directory path to generated data sets
path_ROC  <- "" # e.g. ".../Simulations/Data/"

## Loading analysis functions
source("0_CostumeFunctions_Analysis.R")

## Loading parameters and order them
parameters  <- readRDS(paste0(path_ROC, "Scen1_parameters.rds"))
parameters_ordered <- parameters[order(parameters$beta, parameters$accrual.time,
                                       parameters$fu.time, parameters$cens.rate,
                                       parameters$r, parameters$med.C,
                                       parameters$HR.var, parameters$HR),]

#--------------------------------- ROC Generation ---------------------------------------
##### Preparations
number_HR_simulated_in_each_scenario <- length(unique(parameters_ordered$HR))
number_scenarios <-
  dim(parameters_ordered)[1]/number_HR_simulated_in_each_scenario

# Selecting sub-scenarios used for ROC generation
k_it=c(6,8)
# 6: Power 90%, censoring rate 20%, med_C = 12 months
# 8: Power 90%, censoring rate 60%, med_C = 12 months

# Needed function
elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

# Preparing / Setting up data frames for plotting
df_plot <- data.frame(
  k_it           = rep(k_it,
                       each = length(delta_deserved)*length(HR_thresh_methods)),
  Scen           = rep(c("Pow90_censrate20_medC12", "Pow90_censrate60_medC12"),
                       each = length(delta_deserved)*length(HR_thresh_methods)),
  HR_trueBenefit = rep(delta_deserved, each = length(HR_thresh_methods),
                       times = length(k_it)),
  Threshold      = rep(HR_thresh_methods,
                       times = length(delta_deserved)*length(k_it)),
  TPR.LL         = NA,
  FPR.LL         = NA,
  TPR.PE         = NA,
  FPR.PE         = NA,
  TPR.UL         = NA,
  FPR.UL         = NA
) %>%
  mutate(
    HR_thresh_methods_plot = ifelse(!rep(c(rep(FALSE, 5), TRUE),
                                         length(Threshold)/6), 0, Threshold),
    text_together = factor(paste0(ifelse(k_it == 6, "20'%'~censoring~rate,~",
                                         "60'%'~censoring~rate,~"),
                                  "delta[deserved]==",HR_trueBenefit))
  )

Methods <- data.frame(
  k_it           = rep(k_it, each = length(delta_deserved)*4),
  Scen           = rep(c("Pow90_censrate20_medC12", "Pow90_censrate60_medC12"),
                       each = length(delta_deserved)),
  Method         = rep(c("IQWiG", "Mod.IQWiG", "ESMO", "ESMO.RB"),
                       times = length(k_it)),
  Method_label   = rep(c("IQWiG[RR]", "IQWiG[HR]", "ESMO", "ESMO[RB]"),
                       times = length(k_it)),
  HR_trueBenefit = rep(delta_deserved, each = 4),
  TPR            = NA,
  FPR            = NA
)

##### For Loop to calculate FPR and TPR
for(k in k_it){
  parameters.k <-parameters_ordered[
    (k*number_HR_simulated_in_each_scenario - number_HR_simulated_in_each_scenario + 1) :
      (k*number_HR_simulated_in_each_scenario),]


  #### Setting up variables
  FN.IQWiG <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                   nrow = 1, dimnames = list(NULL, delta_deserved)))
  TP.IQWiG <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                   nrow = 1, dimnames = list(NULL, delta_deserved)))
  TN.IQWiG <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                   nrow = 1, dimnames = list(NULL, delta_deserved)))
  FP.IQWiG <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                   nrow = 1, dimnames = list(NULL, delta_deserved)))

  FN.Mod.IQWiG <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                       nrow = 1, dimnames = list(NULL, delta_deserved)))
  TP.Mod.IQWiG <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                       nrow = 1, dimnames = list(NULL, delta_deserved)))
  TN.Mod.IQWiG <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                       nrow = 1, dimnames = list(NULL, delta_deserved)))
  FP.Mod.IQWiG <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                       nrow = 1, dimnames = list(NULL, delta_deserved)))

  FN.ESMO <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                  nrow = 1, dimnames = list(NULL, delta_deserved)))
  TP.ESMO <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                  nrow = 1, dimnames = list(NULL, delta_deserved)))
  TN.ESMO <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                  nrow = 1, dimnames = list(NULL, delta_deserved)))
  FP.ESMO <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                  nrow = 1, dimnames = list(NULL, delta_deserved)))

  FN.ESMO.RB <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                     nrow = 1, dimnames = list(NULL, delta_deserved)))
  TP.ESMO.RB <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                     nrow = 1, dimnames = list(NULL, delta_deserved)))
  TN.ESMO.RB <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                     nrow = 1, dimnames = list(NULL, delta_deserved)))
  FP.ESMO.RB <- as.data.frame(matrix(0,ncol = length(delta_deserved),
                                     nrow = 1, dimnames = list(NULL, delta_deserved)))


  results_help <- list(
    `0.5` = as.data.frame(
      matrix(0, nrow = length(HR_thresh_methods), ncol = 13,
             dimnames = list(NULL, c("Threshold",
                                     "FN.LL", "TP.LL", "TN.LL", "FP.LL",
                                     "FN.PE", "TP.PE", "TN.PE", "FP.PE",
                                     "FN.UL", "TP.UL", "TN.UL", "FP.UL")))
    ) %>%
      mutate(
        Threshold = HR_thresh_methods
      ),
    `0.6` = as.data.frame(
      matrix(0, nrow = length(HR_thresh_methods), ncol = 13,
             dimnames = list(NULL, c("Threshold",
                                     "FN.LL", "TP.LL", "TN.LL", "FP.LL",
                                     "FN.PE", "TP.PE", "TN.PE", "FP.PE",
                                     "FN.UL", "TP.UL", "TN.UL", "FP.UL")))
    ) %>%
      mutate(
        Threshold = HR_thresh_methods
      ),
    `0.7` = as.data.frame(
      matrix(0, nrow = length(HR_thresh_methods), ncol = 13,
             dimnames = list(NULL, c("Threshold",
                                     "FN.LL", "TP.LL", "TN.LL", "FP.LL",
                                     "FN.PE", "TP.PE", "TN.PE", "FP.PE",
                                     "FN.UL", "TP.UL", "TN.UL", "FP.UL")))
    ) %>%
      mutate(
        Threshold = HR_thresh_methods
      ),
    `0.8` = as.data.frame(
      matrix(0, nrow = length(HR_thresh_methods), ncol = 13,
             dimnames = list(NULL, c("Threshold",
                                     "FN.LL", "TP.LL", "TN.LL", "FP.LL",
                                     "FN.PE", "TP.PE", "TN.PE", "FP.PE",
                                     "FN.UL", "TP.UL", "TN.UL", "FP.UL")))
    ) %>%
      mutate(
        Threshold = HR_thresh_methods
      )
  )

  #### Starting loop for calculating FN, TP, TN and FP
  for (j in 1:dim(parameters.k)[1]){
    para <- parameters.k[j,]
    ### Load data of sub-scenarios with "n_sim" times replicates
    df.loaded <- readRDS(paste0(path_ROC,
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
    ### Analyze data sets
    df.raw.1         <- analysis.assessment.methods(df = df.loaded[[1]])
    df.raw           <- data.frame(matrix(NA,
                                          ncol = dim(df.raw.1)[2],
                                          nrow = para$n_sim))
    df.raw[1,]       <- df.raw.1
    colnames(df.raw) <- names(df.raw.1)
    for (i in 2:para$n_sim) {
      df.raw[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
    }

    ## selecting only significant trials for further analysis
    index.sig <- which(df.raw$sig %in% c(1))
    df.raw.sig <- df.raw[index.sig,]

    ## Calculation of maximal score (yes/no) using only ESMO.RB
    df.raw.sig$ESMO.RB.max <-
      factor(ifelse(df.raw.sig$median.C<=12,
                    as.numeric(df.raw.sig$HR.CI.low <= 0.65),
                    as.numeric(df.raw.sig$HR.CI.low <= 0.7)),
             levels = c(0,1))

    ### FN, TP, TN and FP for different methods for true benefit
    ### "designHR<=delta_deserved vs designHR > delta_deserved"
    for (i.thresh in 1:length(delta_deserved)) {
      true.added.benefit <- (para$HR*para$HR.var) <= delta_deserved[i.thresh]

      ## IQWiG
      FN.IQWiG[i.thresh] <- FN.IQWiG[i.thresh] +
        sum(true.added.benefit==1 & df.raw.sig$IQWiG!=6)
      TP.IQWiG[i.thresh] <- TP.IQWiG[i.thresh] +
        sum(true.added.benefit==1 & df.raw.sig$IQWiG==6)
      TN.IQWiG[i.thresh] <- TN.IQWiG[i.thresh] +
        sum(true.added.benefit==0 & df.raw.sig$IQWiG!=6)
      FP.IQWiG[i.thresh] <- FP.IQWiG[i.thresh] +
        sum(true.added.benefit==0 & df.raw.sig$IQWiG==6)

      ## Mod.IQWiG
      FN.Mod.IQWiG[i.thresh] <- FN.Mod.IQWiG[i.thresh] +
        sum(true.added.benefit==1 & df.raw.sig$Mod.IQWiG!=6)
      TP.Mod.IQWiG[i.thresh] <- TP.Mod.IQWiG[i.thresh] +
        sum(true.added.benefit==1 & df.raw.sig$Mod.IQWiG==6)
      TN.Mod.IQWiG[i.thresh] <- TN.Mod.IQWiG[i.thresh] +
        sum(true.added.benefit==0 & df.raw.sig$Mod.IQWiG!=6)
      FP.Mod.IQWiG[i.thresh] <- FP.Mod.IQWiG[i.thresh] +
        sum(true.added.benefit==0 & df.raw.sig$Mod.IQWiG==6)

      ## ESMO
      FN.ESMO[i.thresh] <- FN.ESMO[i.thresh] +
        sum(true.added.benefit==1 & df.raw.sig$ESMO.OS!=4)
      TP.ESMO[i.thresh] <- TP.ESMO[i.thresh] +
        sum(true.added.benefit==1 & df.raw.sig$ESMO.OS==4)
      TN.ESMO[i.thresh] <- TN.ESMO[i.thresh] +
        sum(true.added.benefit==0 & df.raw.sig$ESMO.OS!=4)
      FP.ESMO[i.thresh] <- FP.ESMO[i.thresh] +
        sum(true.added.benefit==0 & df.raw.sig$ESMO.OS==4)

      ## ESMO.RB
      FN.ESMO.RB[i.thresh] <- FN.ESMO.RB[i.thresh] +
        sum(true.added.benefit==1 & df.raw.sig$ESMO.RB.OS.max==0)
      TP.ESMO.RB[i.thresh] <- TP.ESMO.RB[i.thresh] +
        sum(true.added.benefit==1 & df.raw.sig$ESMO.RB.OS.max==1)
      TN.ESMO.RB[i.thresh] <- TN.ESMO.RB[i.thresh] +
        sum(true.added.benefit==0 & df.raw.sig$ESMO.RB.OS.max==0)
      FP.ESMO.RB[i.thresh] <- FP.ESMO.RB[i.thresh] +
        sum(true.added.benefit==0 & df.raw.sig$ESMO.RB.OS.max==1)

      ## Different threshold methods
      for (i in 1:length(HR_thresh_methods)) {
        # Lower CI threshold
        results_help[[i.thresh]]$FN.LL[i] <- results_help[[i.thresh]]$FN.LL[i] +
          sum(true.added.benefit==1 & df.raw.sig$HR.CI.low >  results_help[[i.thresh]]$Threshold[i])
        results_help[[i.thresh]]$TP.LL[i] <- results_help[[i.thresh]]$TP.LL[i] +
          sum(true.added.benefit==1 & df.raw.sig$HR.CI.low <= results_help[[i.thresh]]$Threshold[i])
        results_help[[i.thresh]]$TN.LL[i] <- results_help[[i.thresh]]$TN.LL[i] +
          sum(true.added.benefit==0 & df.raw.sig$HR.CI.low >  results_help[[i.thresh]]$Threshold[i])
        results_help[[i.thresh]]$FP.LL[i] <- results_help[[i.thresh]]$FP.LL[i] +
          sum(true.added.benefit==0 & df.raw.sig$HR.CI.low <= results_help[[i.thresh]]$Threshold[i])

        # Point estimate threshold
        results_help[[i.thresh]]$FN.PE[i] <- results_help[[i.thresh]]$FN.PE[i] +
          sum(true.added.benefit==1 & df.raw.sig$HR.point >  results_help[[i.thresh]]$Threshold[i])
        results_help[[i.thresh]]$TP.PE[i] <- results_help[[i.thresh]]$TP.PE[i] +
          sum(true.added.benefit==1 & df.raw.sig$HR.point <= results_help[[i.thresh]]$Threshold[i])
        results_help[[i.thresh]]$TN.PE[i] <- results_help[[i.thresh]]$TN.PE[i] +
          sum(true.added.benefit==0 & df.raw.sig$HR.point >  results_help[[i.thresh]]$Threshold[i])
        results_help[[i.thresh]]$FP.PE[i] <- results_help[[i.thresh]]$FP.PE[i] +
          sum(true.added.benefit==0 & df.raw.sig$HR.point <= results_help[[i.thresh]]$Threshold[i])

        # Upper CI threshold
        results_help[[i.thresh]]$FN.UL[i] <- results_help[[i.thresh]]$FN.UL[i] +
          sum(true.added.benefit==1 & df.raw.sig$HR.CI.up >= results_help[[i.thresh]]$Threshold[i])
        results_help[[i.thresh]]$TP.UL[i] <- results_help[[i.thresh]]$TP.UL[i] +
          sum(true.added.benefit==1 & df.raw.sig$HR.CI.up <  results_help[[i.thresh]]$Threshold[i])
        results_help[[i.thresh]]$TN.UL[i] <- results_help[[i.thresh]]$TN.UL[i] +
          sum(true.added.benefit==0 & df.raw.sig$HR.CI.up >= results_help[[i.thresh]]$Threshold[i])
        results_help[[i.thresh]]$FP.UL[i] <- results_help[[i.thresh]]$FP.UL[i] +
          sum(true.added.benefit==0 & df.raw.sig$HR.CI.up <  results_help[[i.thresh]]$Threshold[i])
      }
    }
  }

  ### calculating of FPR and TPR
  ## If one of the denominators is zero, no TRP/TNR can be calculated
  ## --> setting value to: NA
  # IQWiG
  Methods[which(Methods$k_it==k & Methods$Method=="IQWiG"),"TPR"] <-
    as.vector(unlist(ifelse((TP.IQWiG + FN.IQWiG)==0,
                            NA, TP.IQWiG/(TP.IQWiG+FN.IQWiG))))
  Methods[which(Methods$k_it==k & Methods$Method=="IQWiG"),"FPR"] <-
    as.vector(unlist(ifelse((TN.IQWiG + FP.IQWiG)==0, NA,
                            1- TN.IQWiG/(TN.IQWiG+FP.IQWiG))))

  # Mod.IQWiG
  Methods[which(Methods$k_it==k & Methods$Method=="Mod.IQWiG"),"TPR"] <-
    as.vector(unlist(ifelse((TP.Mod.IQWiG + FN.Mod.IQWiG)==0,
                            NA, TP.Mod.IQWiG/(TP.Mod.IQWiG+FN.Mod.IQWiG))))
  Methods[which(Methods$k_it==k & Methods$Method=="Mod.IQWiG"),"FPR"] <-
    as.vector(unlist(ifelse((TN.Mod.IQWiG + FP.Mod.IQWiG)==0,
                            NA, 1- TN.Mod.IQWiG/(TN.Mod.IQWiG+FP.Mod.IQWiG))))

  # ESMO.RB
  Methods[which(Methods$k_it==k & Methods$Method=="ESMO.RB"),"TPR"] <-
    as.vector(unlist(ifelse((TP.ESMO.RB + FN.ESMO.RB)==0,
                            NA,  TP.ESMO.RB/(TP.ESMO.RB+FN.ESMO.RB))))
  Methods[which(Methods$k_it==k & Methods$Method=="ESMO.RB"),"FPR"] <-
    as.vector(unlist(ifelse((TN.ESMO.RB + FP.ESMO.RB)==0,
                            NA, 1- TN.ESMO.RB/(TN.ESMO.RB+FP.ESMO.RB))))

  # ESMO
  Methods[which(Methods$k_it==k & Methods$Method=="ESMO"),"TPR"] <-
    as.vector(unlist(ifelse((TP.ESMO + FN.ESMO)==0,
                            NA, TP.ESMO/(TP.ESMO+FN.ESMO))))
  Methods[which(Methods$k_it==k & Methods$Method=="ESMO"),"FPR"] <-
    as.vector(unlist(ifelse((TN.ESMO + FP.ESMO)==0,
                            NA, 1- TN.ESMO/(TN.ESMO+FP.ESMO))))

  # Different threshold methods
  for (i.thresh in 1:length(delta_deserved)) {
    results_help[[i.thresh]]$TPR.LL <-
      ifelse((results_help[[i.thresh]]$TP.LL+results_help[[i.thresh]]$FN.LL)==0,
             NA,
             results_help[[i.thresh]]$TP.LL/(results_help[[i.thresh]]$TP.LL+results_help[[i.thresh]]$FN.LL))
    results_help[[i.thresh]]$FPR.LL <-
      ifelse((results_help[[i.thresh]]$TN.LL+results_help[[i.thresh]]$FP.LL)==0,
             NA,
             1 - results_help[[i.thresh]]$TN.LL/(results_help[[i.thresh]]$TN.LL+results_help[[i.thresh]]$FP.LL))

    results_help[[i.thresh]]$TPR.PE <-
      ifelse((results_help[[i.thresh]]$TP.PE+results_help[[i.thresh]]$FN.PE)==0,
             NA,
             results_help[[i.thresh]]$TP.PE/(results_help[[i.thresh]]$TP.PE+results_help[[i.thresh]]$FN.PE))
    results_help[[i.thresh]]$FPR.PE <-
      ifelse((results_help[[i.thresh]]$TN.PE+results_help[[i.thresh]]$FP.PE)==0,
             NA,
             1 - results_help[[i.thresh]]$TN.PE/(results_help[[i.thresh]]$TN.PE+results_help[[i.thresh]]$FP.PE))

    results_help[[i.thresh]]$TPR.UL <-
      ifelse((results_help[[i.thresh]]$TP.UL+results_help[[i.thresh]]$FN.UL)==0,
             NA,
             results_help[[i.thresh]]$TP.UL/(results_help[[i.thresh]]$TP.UL+results_help[[i.thresh]]$FN.UL))
    results_help[[i.thresh]]$FPR.UL <-
      ifelse((results_help[[i.thresh]]$TN.UL+results_help[[i.thresh]]$FP.UL)==0,
             NA,
             1 - results_help[[i.thresh]]$TN.UL/(results_help[[i.thresh]]$TN.UL+results_help[[i.thresh]]$FP.UL))



    df_plot[which(df_plot$k_it==k & df_plot$HR_trueBenefit == delta_deserved[i.thresh]),
            c("Threshold", "TPR.LL", "FPR.LL", "TPR.PE", "FPR.PE", "TPR.UL", "FPR.UL")] <-
      results_help[[i.thresh]] %>%
      dplyr::select(Threshold, TPR.LL, FPR.LL, TPR.PE, FPR.PE, TPR.UL, FPR.UL)
  }
}





##### Plotting of calculated data
size_points <- 0.00005
size_points_Method <- 0.4
size_line <- 0.2
size_thresholds <- 1


## Final changes to data sets
df_plot %>% mutate(
  # Only plot every 6th threshold, so that the plot is not overfilled
  HR_thresh_methods_plot = ifelse(!rep(c(rep(FALSE, 5), TRUE),
                                       length(Threshold)/6), 0, Threshold),
  text_together = factor(paste0(ifelse(k_it == 6, "'20%'~censoring~rate*','~",
                                       "'60%'~censoring~rate*','~"),
                                "delta[deserved]==",HR_trueBenefit))
) -> df_plot

Methods %>% mutate(
  text_together = factor(paste0(ifelse(k_it == 6, "'20%'~censoring~rate*','~",
                                       "'60%'~censoring~rate*','~"),
                                "delta[deserved]==",HR_trueBenefit))
) -> Methods


## ROC plot
df_plot %>%
  ggplot() +
  # LL thresholds
  geom_line(aes(x = FPR.LL,
                y = TPR.LL,
                color = "LL thresholds"), size = size_line) +
  geom_point(data = subset(df_plot,
                           elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
             aes(x = FPR.LL, y = TPR.LL,  color = "LL thresholds"), size = size_points) +
  geom_text(data = subset(df_plot, elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
            aes(x = FPR.LL - 0.05,
                y = TPR.LL + 0.03,
                label = Threshold, color = "LL thresholds"),
            check_overlap = TRUE, size = size_thresholds) +
  # PE thresholds
  geom_line(aes(x = FPR.PE,
                y = TPR.PE,
                color = "PE thresholds"), size = size_line) +
  geom_point(data = subset(df_plot,
                           elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
             aes(x = FPR.PE, y = TPR.PE, color = "PE thresholds"), size = size_points) +
  geom_text(data = subset(df_plot, elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
            aes(x = FPR.PE + 0.07,
                y = TPR.PE - 0.01,
                label = Threshold, color = "PE thresholds"),
            check_overlap = TRUE, size = size_thresholds) +
  # UL thresholds
  geom_line(aes(x = FPR.UL,
                y = TPR.UL, color = "UL thresholds"), size = size_line) +
  geom_point(data = subset(df_plot,
                           elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
             aes(x = FPR.UL, y = TPR.UL, color = "UL thresholds"), size = size_points) +
  geom_text(data = subset(df_plot, elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
            aes(x = FPR.UL + 0.08,
                y = TPR.UL - 0.02,
                label = Threshold, color = "UL thresholds"),
            check_overlap = TRUE, size = size_thresholds) +
  # Methods (IQWiG, Mod.IQWiG, ESMO, ESMO.RB)
  geom_point(data = Methods,
             aes(x = FPR, y = TPR, color = factor(Method),
                 fill = factor(Method), shape = factor(Method)),
             size = size_points_Method) +
  # Legend
  scale_color_manual(
    name = "Maximal score rule using",
    values = c("LL thresholds" = "blue", "PE thresholds" = "orange",
               "UL thresholds" = "black", "IQWiG" = "green4", "Mod.IQWiG" = "green",
               "ESMO" = "blueviolet", "ESMO.RB" = "cyan"),
    labels = c(bquote(HR^"-" ~ "thresholds"), "HR  thresholds",
               bquote(HR^"+" ~ "thresholds"), expression(IQWiG[RR]),
               expression(IQWiG[HR]), "ESMO", expression(ESMO[RB]))
  ) +
  scale_shape_manual(
    name = "Method",
    values = c("IQWiG" = 24, "Mod.IQWiG" = 24, "ESMO" = 24, "ESMO.RB" = 24),
    labels = c(expression(IQWiG[RR]), expression(IQWiG[HR]), "ESMO", expression(ESMO[RB]))
  ) +
  scale_fill_manual(
    name = "Method",
    values = c("IQWiG" = "green4", "Mod.IQWiG" = "green", "ESMO" = "blueviolet", "ESMO.RB" = "cyan"),
    labels = c(expression(IQWiG[RR]), expression(IQWiG[HR]), "ESMO", expression(ESMO[RB]))
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = c(16,16,16,17,17,17,17), size = 0.3)),
    shape = guide_legend(override.aes = list(size = 0.4,
                                             fill = c("IQWiG" = "green4",
                                                      "Mod.IQWiG" = "green",
                                                      "ESMO" = "blueviolet",
                                                      "ESMO.RB" = "cyan")))
  ) +
  # Labels
  labs(title = "",
       subtitle = "",
       x = "FPR",
       y = "TPR",
       color = "") +
  # Changing range of x-axis and y-axis
  scale_y_continuous(limits=c(0, 1.06), breaks = seq(0,1, by=0.2)) +
  scale_x_continuous(limits=c(-0.06, 1), breaks = seq(0,1, by=0.2)) +
  # Theme
  theme_bw() +
  theme(
    text = element_text(size = 5,
                        family = "serif" # TT Times New Roman
                        #family = "sans" # TT Arial
    ),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position = "bottom",
    axis.ticks = element_line(size=0.2),
    axis.ticks.length = unit(0.05, "cm"),
    strip.background = element_rect(colour="black",size=0.4), # change line width of label of facet_wrap
    strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")) # change label width of facet_wrap
  ) +
  facet_wrap(text_together~., labeller = label_parsed, nrow = 2) -> Figure3
Figure3



## Better legend of ROC plot
df_plot %>%
  ggplot() +
  # LL thresholds
  geom_line(aes(x = FPR.LL,
                y = TPR.LL,
                color = "LL thresholds"), size = size_line) +
  geom_point(data = subset(df_plot,
                           elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
             aes(x = FPR.LL, y = TPR.LL,  color = "LL thresholds"), size = size_points) +
  geom_text(data = subset(df_plot, elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
            aes(x = FPR.LL - 0.05,
                y = TPR.LL + 0.03,
                label = Threshold, color = "LL thresholds"),
            check_overlap = TRUE, size = size_thresholds) +
  # PE thresholds
  geom_line(aes(x = FPR.PE,
                y = TPR.PE,
                color = "PE thresholds"), size = size_line) +
  geom_point(data = subset(df_plot,
                           elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
             aes(x = FPR.PE, y = TPR.PE, color = "PE thresholds"), size = size_points) +
  geom_text(data = subset(df_plot, elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
            aes(x = FPR.PE + 0.07,
                y = TPR.PE - 0.01,
                label = Threshold, color = "PE thresholds"),
            check_overlap = TRUE, size = size_thresholds) +
  # UL thresholds
  geom_line(aes(x = FPR.UL,
                y = TPR.UL, color = "UL thresholds"), size = size_line) +
  geom_point(data = subset(df_plot,
                           elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
             aes(x = FPR.UL, y = TPR.UL, color = "UL thresholds"), size = size_points) +
  geom_text(data = subset(df_plot, elementwise.all.equal(Threshold, HR_thresh_methods_plot)),
            aes(x = FPR.UL + 0.08,
                y = TPR.UL - 0.02,
                label = Threshold, color = "UL thresholds"),
            check_overlap = TRUE, size = size_thresholds) +
  # Methods
  geom_point(data = Methods,
             aes(x = FPR, y = TPR, color = factor(Method), shape = factor(Method)),
             size = size_points_Method) +
  # Legend
  scale_color_manual(
    name = "Maximal score rule using",
    values = c("LL thresholds" = "blue", "PE thresholds" = "orange",
               "UL thresholds" = "black"),
    labels = c(bquote(HR^"-" ~ "thresholds"), "HR  thresholds",
               bquote(HR^"+" ~ "thresholds"))
  ) +
  scale_shape_manual(
    name = "Method",
    values = c("IQWiG" = 17, "Mod.IQWiG" = 17, "ESMO" = 17, "ESMO.RB" = 17),
    labels = c(expression(IQWiG), expression("Mod-IQWiG"[HR]~"scores"),
               "ESMO", expression(ESMO[RB]))
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = c(16,16,16), size = 0.4)),
    shape = guide_legend(override.aes = list(size = 0.5,
                                             color = c("IQWiG" = "green4",
                                                       "Mod.IQWiG" = "green",
                                                       "ESMO" = "blueviolet",
                                                       "ESMO.RB" = "cyan")))
  ) +
  # Labels
  labs(title = "",
       subtitle = "",
       x = "FPR",
       y = "TPR",
       color = "") +
  scale_y_continuous(limits=c(0, 1.06), breaks = seq(0,1, by=0.2)) +
  scale_x_continuous(limits=c(-0.06, 1), breaks = seq(0,1, by=0.2)) +
  theme_bw() +
  theme(
    text = element_text(size = 5,
                        family = "serif" # TT Times New Roman
                        #family = "sans" # TT Arial
    ),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position = "none",
    legend.text.align = 0, # align legend to the left
    axis.ticks = element_line(size=0.2),
    axis.ticks.length = unit(0.05, "cm"),
    strip.background = element_rect(colour="black",size=0.4), # change line width of label of facet_wrap
    strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")) # change label width of facet_wrap
  ) +
  facet_wrap(text_together~., labeller = label_parsed, nrow = 2) -> Paper1_ROC_legend

legend <- ggpubr::get_legend(Paper1_ROC_legend + theme(legend.position = "right") +
                       guides(guide_legend(nrow = 1)))
ggpubr::as_ggplot(legend)



####### for respond to Reviewer 1, major 6:
df_plot %>% filter(Scen=="Pow90_censrate20_medC12" & HR_trueBenefit %in% c(0.5,0.7) & as.character(df_plot$Threshold)=="0.43")

