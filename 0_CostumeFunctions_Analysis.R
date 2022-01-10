################################################################################
####        Costume functions needed for analysis of generated data         ####
################################################################################

#--------------------- Analysis of generated data sets -------------------------
analysis.assessment.methods <- function(df){
  #### Analysis of generated data sets
  ## Input:
  # df: Generated data set
  
  ## Output: 
  # HR.point:    Estimated hazard ratio point estimate
  # HR.CI.low:   Estimated lower limit of the 95% confidence interval of the 
  #              hazard ratio
  # HR.CI.up:    Estimated upper limit of the 95% confidence interval of the 
  #              hazard ratio
  # median.C:    Median survival time of control group
  # median.T:    Median survival time of treatment group 
  # median.gain: Difference of median.T - median.C
  # surv.gain:   Survival rate difference between both treatment groups 
  #              (time point depending on observed median.C)
  # IQWiG:       IQWiG's additional benefit assessment method
  # Mod.IQWiG:   Modified IQWiG's additional benefit assessment method
  # ESMO:        ESMO's additional benefit assessment method
  # sig:         Significant difference (in the correct direction) between 
  #              treatment groups present
  # If data set does not show a significant difference between treatment and 
  # control group using the log-rank test, "sig" is set to 0 (not significant)
  # and all other outputs are set to NA.
  
  ## Information for data set construction
  # Data.TandC: Survival / Censoring time
  # Data.status: Indication of censoring and event/death 
  #              (0: censored, 1: event/death)
  # Data.Z: Treatment allocation / group 
  
  #### Function start
  ### Perform Cox-Regression
  cox.sum <- summary(coxph(Surv(Data.TandC, Data.status != 0) ~ Data.Z, 
                           data = df))
  ### Checking if Treatment is significant in the correct direction
  if(cox.sum$waldtest[[3]]>0.05 | cox.sum$conf.int[3]>1){
    ## Not significant, hence no further analysis is performed
    Ana.out <- data.frame(HR.point = NA, HR.CI.low = NA, HR.CI.up = NA,
                          median.C = NA, median.T = NA, median.gain = NA,
                          surv.gain = NA,
                          IQWiG = NA, Mod.IQWiG = NA, ESMO = NA,
                          sig = 0)
  }else{
    ## Significant, hence further analysis is performed
    
    ## Calculation of HR, upper and lower CI
    HR.point  <- cox.sum$conf.int[1]
    HR.CI.low <- cox.sum$conf.int[3]
    HR.CI.up  <- cox.sum$conf.int[4]
    
    ## Calculation of median survival times
    # Due to the scenarios with 60% censoring rate some of the survival curves
    # of the simulated data did not fall underneath 50% and, therefore, ESMO
    # score could not be calculated. To solve this issue a conservative 
    # approach was implemented using the last observed maximal event or
    # censoring time present. 
    survtest <- summary(
      survfit(Surv(Data.TandC, Data.status != 0) ~ Data.Z, data = df)
    )$table
    # Treatment group (median.T)
    if(is.na(survtest[2,7])){
      median.T <- max(df$Data.TandC[df$Data.Z=="Treatment"])
    }else{
      median.T <- survtest[2,7]
    }
    # Control group (median.C)
    if(is.na(survtest[1,7])){
      median.C <- max(df$Data.TandC[df$Data.Z=="Control"])
    }else{
      median.C <- survtest[1,7]
    }
    # median.gain
    median.gain <- abs(median.T - median.C)
    
    ## Calculation of med.C (decides which part of the ESMO scale 
    ## is used, see Table 1)
    if (median.C <=12){
      med.C = "<=12"
    }else if (median.C > 12 && median.C <= 24) {
      med.C = "12.24"
    } else if (median.C > 24) {
      med.C = ">24"
    }
    
    ## Survival gain
    # Some of the survival curves of the simulated data did not reach 2, 3 or 5
    # year survival and hence the ESMO score could not be calculated. To solve 
    # this issue we extended the survival curve, meaning that either the last
    # censoring time was carried forward or a survival rate of 0% was used.
    fit.C <- survfit(Surv(Data.TandC, Data.status != 0) ~ 1, 
                     data = df[which(df$Data.Z=="Control"),])
    fit.T <- survfit(Surv(Data.TandC, Data.status != 0) ~ 1, 
                     data = df[which(df$Data.Z=="Treatment"),])
    
    surv.gain <-
      switch (as.character(med.C),
              "<=12" = {
                surv.C = summary(fit.C, times=2*12, extend=TRUE)$surv 
                surv.T = summary(fit.T, times=2*12, extend=TRUE)$surv
                if(is.null(surv.C) && is.null(surv.T)){
                  0
                }else if(is.null(surv.C)){
                  abs(0 - surv.T)
                }else if(is.null(surv.T)){
                  abs(surv.C - 0)
                }else{
                  abs(surv.C - surv.T)
                }
              },
              "12.24" = {
                surv.C = summary(fit.C, times=3*12, extend=TRUE)$surv
                surv.T = summary(fit.T, times=3*12, extend=TRUE)$surv
                if(is.null(surv.C) && is.null(surv.T)){
                  0
                }else if(is.null(surv.C)){
                  abs(0 - surv.T)
                }else if(is.null(surv.T)){
                  abs(surv.C - 0)
                }else{
                  abs(surv.C - surv.T)
                }
              },
              ">24" = {
                surv.C = summary(fit.C, times=5*12, extend=TRUE)$surv
                surv.T = summary(fit.T, times=5*12, extend=TRUE)$surv
                if(is.null(surv.C) && is.null(surv.T)){
                  0
                }else if(is.null(surv.C)){
                  abs(0 - surv.T)
                }else if(is.null(surv.T)){
                  abs(surv.C - 0)
                }else{
                  abs(surv.C - surv.T)
                }
              })
    
    ## Additional benefit methods
    # IQWiG
    IQWiG <- (HR.CI.up<0.85)*6 + 
      (HR.CI.up<0.95 & HR.CI.up>=0.85)*5 +
      (HR.CI.up<1 & HR.CI.up>=0.95)*4
    
    # Modified IQWiG (Mod.IQWiG)
    Mod.IQWiG <- (HR.CI.up<0.7908765)*6 + 
      (HR.CI.up<0.9286669 & HR.CI.up>=0.7908765)*5 + 
      (HR.CI.up<1 & HR.CI.up>=0.9286669)*4
    
    # ESMO
    ESMO <-
      switch (as.character(med.C),
              "<=12" = {
                if(surv.gain>=0.1){
                  4
                }else{
                  (HR.CI.low>0.7 || median.gain<1.5)*1 +
                    (  (HR.CI.low<=0.65 && median.gain<2 && median.gain>=1.5) ||
                         (HR.CI.low<0.7 && HR.CI.low>0.65 && median.gain>=1.5)  )*2 +
                    (HR.CI.low<=0.65 && median.gain<3 && median.gain>=2)*3 +
                    (HR.CI.low<=0.65 && median.gain>=3)*4
                }
              },
              "12.24" = {
                if(surv.gain>=0.1){
                  4
                }else{
                  (HR.CI.low>0.75 || median.gain<1.5)*1 +
                    ( (HR.CI.low<=0.7 && median.gain<3 && median.gain>=1.5) ||
                        (HR.CI.low<0.75 && HR.CI.low>0.7 && median.gain>=1.5)  )*2 +
                    (HR.CI.low<=0.7 && median.gain<5 && median.gain>=3)*3 +
                    (HR.CI.low<=0.7 && median.gain>=5)*4
                }
              },
              ">24" = {
                if(surv.gain>=0.1){
                  4
                }else{
                  (HR.CI.low>0.75 || median.gain<4)*1 +
                    ( (HR.CI.low<=0.7 && median.gain<6 && median.gain>4) ||
                        (HR.CI.low<0.75 && HR.CI.low>0.7 && median.gain>=4)  )*2 +
                    (HR.CI.low<=0.7 && median.gain<9 && median.gain>=6)*3 +
                    (HR.CI.low<=0.7 && median.gain>=9)*4
                }
              })
    
    ## Analyzed data
    Ana.out <- data.frame(HR.point, HR.CI.low, HR.CI.up,
                          median.C, median.T, median.gain,
                          surv.gain,
                          IQWiG, Mod.IQWiG, ESMO,
                          sig = 1)
  }
  ### Return analyzed data 
  return(Ana.out)
}


#---------------------- Summary of analyzed data sets --------------------------
analysis.summary <- function(df.raw, n_sim){
  ### Summary of analyzed generated data sets
  ## Input:
  # df.raw: All n_sim time analysed data sets of one specific sub-scenario using the function analysis.assessment.methods()
  # n_sim: Amount of generated data set of each sub-scenario
  
  #### Function start
  ## Power
  pow <- sum(df.raw$sig == 1)/n_sim # all significant trials
  
  ## Selecting only significant trials for further analysis
  index.sig <- which(df.raw$sig ==1)
  df.raw.sig <- df.raw[index.sig,]
  
  ## Calculating the mean of important continous measures
  means <- colMeans(df.raw.sig[,c(1,2,3,4,5,6,7)])
  mean.HRpoint      <- means[[1]]
  mean.HRCIlow      <- means[[2]]
  mean.HRCIup       <- means[[3]]
  mean.MedianC      <- means[[4]]
  mean.MedianT      <- means[[5]]
  mean.MedianGain   <- means[[6]]
  mean.SurvGain     <- means[[7]]
  
  ## Maximal observed HR (point estimate) of significant trials
  max.HRpoint <- max(df.raw.sig$HR.point)
  
  ## IQWiG
  # Percentage maximal IQWiG score
  df.raw.sig$IQWiG.max <- df.raw.sig$IQWiG==6
  perc.IQWiG.max       <- sum(df.raw.sig$IQWiG.max==1)/dim(df.raw.sig)[1]
  
  # Number of maximal and not maximal IQWiG score
  n.IQWiG.max     <- sum(df.raw.sig$IQWiG.max==1)
  n.IQWiG.not.max <- sum(df.raw.sig$IQWiG.max!=1)
  
  # Calculation of maximal HR and corresponding CI for the maximum IQWiG scores
  if(sum(df.raw.sig$IQWiG.max) != 0){
    max.HRpoint.IQWiG <-
      max(df.raw.sig[which(df.raw.sig$IQWiG.max==1),]$HR.point)
    max.HRCIup.IQWiG  <-
      max(df.raw.sig[which(df.raw.sig$IQWiG.max==1),]$HR.CI.up)
    max.HRCIlow.IQWiG <-
      max(df.raw.sig[which(df.raw.sig$IQWiG.max==1),]$HR.CI.low)
  }else{
    # If no maximal score was achieved, setting to "NA"
    max.HRpoint.IQWiG <- NA
    max.HRCIup.IQWiG  <- NA
    max.HRCIlow.IQWiG <- NA
  }
  
  ## Modified IQWiG (Mod.IQWiG)
  # Percentage of maximal Mod.IQWiG score
  df.raw.sig$Mod.IQWiG.max <- df.raw.sig$Mod.IQWiG==6
  perc.Mod.IQWiG.max       <- sum(df.raw.sig$Mod.IQWiG.max==1)/dim(df.raw.sig)[1]
  
  # Number of maximal and not maximal Mod.IQWiG score
  n.Mod.IQWiG.max     <- sum(df.raw.sig$Mod.IQWiG==1)
  n.Mod.IQWiG.not.max <- sum(df.raw.sig$Mod.IQWiG!=1)
  
  # Calculation of maximal HR and corresponding CI for the maximum
  # Mod.IQWiG scores
  if(sum(df.raw.sig$Mod.IQWiG.max) != 0){
    max.HRpoint.Mod.IQWiG <-
      max(df.raw.sig[which(df.raw.sig$Mod.IQWiG.max==1),]$HR.point)
    max.HRCIup.Mod.IQWiG  <-
      max(df.raw.sig[which(df.raw.sig$Mod.IQWiG.max==1),]$HR.CI.up)
    max.HRCIlow.Mod.IQWiG <-
      max(df.raw.sig[which(df.raw.sig$Mod.IQWiG.max==1),]$HR.CI.low)
  }else{
    # If no maximal score was achieved, setting to "NA"
    max.HRpoint.Mod.IQWiG <- NA
    max.HRCIup.Mod.IQWiG  <- NA
    max.HRCIlow.Mod.IQWiG <- NA
  }
  
  ## ESMO
  # Percentage maximal ESMO score
  df.raw.sig$ESMO.max <- as.numeric(df.raw.sig$ESMO==4)
  perc.ESMO.max       <- sum(df.raw.sig$ESMO==4)/dim(df.raw.sig)[1]
  
  # Number of maximal and not maximal ESMO score
  n.perc.ESMO.max     <-  sum(df.raw.sig$ESMO==4)
  n.perc.ESMO.not.max <-  sum(df.raw.sig$ESMO!=4)
  
  # Calculation of maximal HR and corresponding CI for the maximum ESMO scores
  if(sum(df.raw.sig$ESMO.max) != 0){
    max.HRpoint.ESMO <-
      max(df.raw.sig[which(df.raw.sig$ESMO.max==1),]$HR.point)
    max.HRCIup.ESMO  <-
      max(df.raw.sig[which(df.raw.sig$ESMO.max==1),]$HR.CI.up)
    max.HRCIlow.ESMO <-
      max(df.raw.sig[which(df.raw.sig$ESMO.max==1),]$HR.CI.low)
  }else{
    # If no maximal score was achieved, setting to "NA"
    max.HRpoint.ESMO <- NA
    max.HRCIup.ESMO  <- NA
    max.HRCIlow.ESMO <- NA
  }
  
  ## ESMO's relative benefit rule (ESMO.RB)
  # Percentage maximal ESMO.RB score
  df.raw.sig$ESMO.RB.max <-
    ifelse(df.raw.sig$median.C<=12,
           as.numeric(df.raw.sig$HR.CI.low <= 0.65),
           as.numeric(df.raw.sig$HR.CI.low <= 0.7))
  
  perc.ESMO.RB.max <- sum(df.raw.sig$ESMO.RB.max==1)/dim(df.raw.sig)[1]
  
  # Number of maximal and not maximal maximal ESMO.RB score
  n.perc.ESMO.RB.max     <-  sum(df.raw.sig$ESMO.RB.max==1)
  n.perc.ESMO.RB.not.max <-  sum(df.raw.sig$ESMO.RB.max!=1)
  
  # Calculation of maximal HR and corresponding CI for the maximum
  # ESMO.RB scores
  if(sum(df.raw.sig$ESMO.RB.max) != 0){
    max.HRpoint.ESMO.RB <-
      max(df.raw.sig[which(df.raw.sig$ESMO.RB.max==1),]$HR.point)
    max.HRCIup.ESMO.RB  <-
      max(df.raw.sig[which(df.raw.sig$ESMO.RB.max==1),]$HR.CI.up)
    max.HRCIlow.ESMO.RB <-
      max(df.raw.sig[which(df.raw.sig$ESMO.RB.max==1),]$HR.CI.low)
  }else{
    # If no maximal score was achieved, setting to "NA"
    max.HRpoint.ESMO.RB <- NA
    max.HRCIup.ESMO.RB  <- NA
    max.HRCIlow.ESMO.RB <- NA
  }
  
  ## Percentages of maximal scores / significant trials compared to number of
  ## iterations (n_sim)
  perc.Mod.IQWiG.max.n_sim <- sum(df.raw.sig$Mod.IQWiG.max==1)/n_sim
  perc.IQWiG.max.n_sim     <- sum(df.raw.sig$IQWiG.max==1)/n_sim
  perc.ESMO.max.n_sim      <- sum(df.raw.sig$ESMO.max==1)/n_sim
  perc.ESMO.RB.max.n_sim   <- sum(df.raw.sig$ESMO.RB.max==1)/n_sim
  
  ## Correlation between ESMO / IQWiG, ESMO / Mod.IQWiG and Mod.IQWiG / IQWiG
  if(dim(table(as.numeric(df.raw.sig$ESMO))) == 1 |
     length(levels(factor(df.raw.sig$IQWiG))) == 1){
    # Setting correlation to NA if only the same score was given by
    # ESMO or IQWiG
    Cor.IQWiG.ESMO <- NA
  }else{
    Cor.IQWiG.ESMO  <- cor(df.raw.sig$IQWiG, df.raw.sig$ESMO,
                           method = "spearman")
  }
  
  if(dim(table(as.numeric(df.raw.sig$ESMO))) == 1 |
     length(levels(factor(df.raw.sig$Mod.IQWiG))) == 1){
    # Setting correlation to NA if only the same score was given by
    # ESMO or Mod.IQWiG
    Cor.Mod.IQWiG.ESMO <- NA
  }else{
    Cor.Mod.IQWiG.ESMO  <- cor(df.raw.sig$Mod.IQWiG, df.raw.sig$ESMO,
                               method = "spearman")
  }
  
  if(length(levels(factor(df.raw.sig$IQWiG))) == 1 | 
     length(levels(factor(df.raw.sig$Mod.IQWiG))) == 1 ){
    # Setting correlation to NA if only the same score was given by
    # IQWiG or Mod.IQWiG
    Cor.Mod.IQWiG.IQWiG <- NA
  }else{
    Cor.Mod.IQWiG.IQWiG <- cor(df.raw.sig$Mod.IQWiG, df.raw.sig$IQWiG,
                               method = "spearman")
  }
  
  
  ## Combining everything to edited and analyzed results
  results <- data.frame(
    pow,
    perc.Mod.IQWiG.max.n_sim, perc.IQWiG.max.n_sim, perc.ESMO.max.n_sim,
    perc.ESMO.RB.max.n_sim,
    n.Mod.IQWiG.max, n.Mod.IQWiG.not.max, n.IQWiG.max, n.IQWiG.not.max,
    n.perc.ESMO.RB.max, n.perc.ESMO.RB.not.max, n.perc.ESMO.max, 
    n.perc.ESMO.not.max,
    perc.IQWiG.max, max.HRpoint.IQWiG, max.HRCIup.IQWiG, max.HRCIlow.IQWiG,
    perc.Mod.IQWiG.max, max.HRpoint.Mod.IQWiG, max.HRCIup.Mod.IQWiG, 
    max.HRCIlow.Mod.IQWiG,
    perc.ESMO.max, max.HRpoint.ESMO, max.HRCIup.ESMO, max.HRCIlow.ESMO,
    perc.ESMO.RB.max, max.HRpoint.ESMO.RB, max.HRCIup.ESMO.RB, 
    max.HRCIlow.ESMO.RB,
    Cor.IQWiG.ESMO, Cor.Mod.IQWiG.ESMO, Cor.Mod.IQWiG.IQWiG,
    mean.HRpoint, mean.HRCIlow, mean.HRCIup, mean.MedianC, mean.MedianT,
    mean.MedianGain, mean.SurvGain,
    max.HRpoint
  )
  
  ## Return edited and analyzed data
  return(results)
}
