################################################################################
####           Costume functions needed for data generation                 ####
################################################################################

#-------------------------- Sample Size Calculation ----------------------------
sample.size.func <- function(alpha, beta, n.proportion, HR, median.control,
                             fu.time, accrual.time, cens.rate) {
# Calculation of sample size for survival data (Log-Rank Test) via Schoenfeld
# approach
# Parameters:
  # alpha:          type-I-error
  # beta:           type-II-error
  # n.proportion:   sample size proportion of the two intervention groups: n.T/n.C
  # HR:             hazard ratio of treatment against control
  # median.control: median survival time of the control group
  # cens.rate:      censoring rate

## Calculating lambda parameters for exponential distribution:
  lambda.C <- log(2)/median.control
  lambda.T <- HR*log(2)/median.control

## Shoenefeld approach:
  # events
  d.raw <- (1+n.proportion)^2/n.proportion * (qnorm(1-alpha/2) + qnorm(1-beta))^2/(log(HR)^2)
  # rounding events up so it can be divided by n.proportion
  d <- ceiling(d.raw)-1
  response <- 1
  while (response!=0) {
    d  <- d + 1
    response <- d %% (n.proportion+1)
  }

  # Probability of an event for combination of administrative censoring and
  # specific censoring rate
  P.D <- 1 - 1/(6*(1+n.proportion)) * (
    exp(-lambda.C*fu.time) +
    n.proportion*exp(-lambda.T*fu.time) +
    4*( exp(-lambda.C*(accrual.time/2 + fu.time)) +
          n.proportion*exp(-lambda.T*(accrual.time/2 + fu.time))
      ) +
    exp(-lambda.C*(accrual.time+fu.time)) +
    n.proportion*exp(-lambda.T*(accrual.time+fu.time))
  )

  # raw sample size
    # if probability of an event using administrative censoring is smaller than
    # (1-cens.rate), than use only P.D for calculation of n.raw, otherwise use
    # only cens.rate.
    # In other words: If censoring rate of administrative censoring is larger
    #                 than from specific censoring, than use only cens.rate for
    #                 calculation of n.raw, otherwise use only P.D.
  if(P.D < (1-cens.rate)){
    n.raw = d/P.D
  }else{
    n.raw = d/(1-cens.rate)
  }

  # Rounding sample size up so it can be divided by n.proportion
  n.schoen <- ceiling(n.raw)-1
  response <- 1
  while (response!=0) {
    n.schoen <- n.schoen + 1
    response <- n.schoen %% (n.proportion+1)
  }

  return(list(n.T = n.proportion*n.schoen/(n.proportion+1),
              n.C = n.schoen/(n.proportion+1)))
}

#----------- Data Generation of exponential distributed failure times ----------
DataEXP <- function(n.T, n.C, median.control, accrual.time, fu.time, HR, HR.var,
                    cens.rate, index.seed) {
## Generating of survival data with exponential distributed failure times
## with censoring (administrative and specific censoring)
## Parameters:
  # n.T:            number of patients in treatment group
  # n.C:            number of patients in control group
  # median.control: median survival time of the control group
  # accrual.time:   accrual time
  # fu.time:        follow-up time
  # HR:             hazard ratio of treatment against control
  # HR.var:         variation of the true HR and design HR.
  #                     1 --> design.HR = true.HR
  # cens.rate:      censoring rate
  # index.seed:     seed for the simulation

## Setting seed
  set.seed(index.seed)

## Event times
  mytimes.C <- rexp(n.C, rate = log(2)/median.control)
  mytimes.T <- rexp(n.T, rate = (HR*HR.var)*log(2)/median.control)
  # Combining event times of both groups
  mytimes   <- c(mytimes.C, mytimes.T)


## Censoring with combination of administrative censoring and specific censoring
## rate
  # Administrative censoring
  cens.times.AC <- runif(n.C+n.T, min = 0, max = accrual.time) + fu.time

  # Help data frame
  Data.help <- data.frame(id = 1:(n.C + n.T), T = mytimes,
                          Z = factor(c(rep("Control",n.C), rep("Treatment",n.T))),
                          C = cens.times.AC)
  Data.help$status <- factor(as.numeric(Data.help$T <= Data.help$C))
  Data.help$TandC  <- pmin(Data.help$T,Data.help$C)

  # Specific censoring rate of the remaining events "mytimes.AC"
    mytimes.AC <- Data.help$T[which(Data.help$status==1)]
    # calculating needed censoring rate to receive overall censoring rate of
    # cens.rate
    current.CensoringRate <- 1 # of Data.help$status==1 --> always 100%
    needed.CensoringRate  <-
      (cens.rate * dim(Data.help)[1] -
         sum(Data.help$status==0)*current.CensoringRate) /
      sum(Data.help$status==1)
    # if administrative censoring rate is already larger than wanted censoring
    # rate of cens.rate, no specific censoring rate is added
    if (needed.CensoringRate > 0){
      cens.times.SC <-
        sapply(mytimes.AC,
               function(mytimes.AC, needed.CensoringRate){
                 rexp(1, rate = -log(1-needed.CensoringRate)/mytimes.AC)
               },
               needed.CensoringRate)
    }else{
      # if administrative censoring rate is already larger than wanted censoring
      # rate of cens.rate, no specific censoring rate is added
      # --> Setting censoring added censoring rates to Infinity
      # --> event rates are always smaller than specific censoring rates
      cens.times.SC <- rep(Inf, length(mytimes.AC))
    }

## Final data frame
  Data <- Data.help
  Data$C[which(Data$status==1)] <- cens.times.SC
  Data$status <- factor(as.numeric(Data$T <= Data$C))
  Data$TandC  <- pmin(Data$T,Data$C)

  Data.out <- data.frame(Data$TandC, Data$status, Data$Z)

## Returning data
return(Data.out)
}

#------------- Data Generation of Weibull distributed failure times ------------
DataWEIB <- function(n.T, n.C, median.control, accrual.time, fu.time, HR, HR.var,
                    cens.rate, shape.C.T, index.seed) {
## Generating of survival data with Weibull distributed failure times
## with censoring (administrative and specific censoring)
## Parameters:
  # n.T:            number of patients in treatment group
  # n.C:            number of patients in control group
  # median.control: median survival time of the control group
  # accrual.time:   accrual time
  # fu.time:        follow-up time
  # HR:             hazard ratio of treatment against control
  # HR.var:         variation of the true HR and design HR.
  #                     1 --> design.HR = true.HR
  # cens.rate:      censoring rate, censoring exponential distributed
  # shape.C.T:      shape parameter for control and treatment group
  # index.seed:     seed for the simulation

## Setting seed
  set.seed(index.seed)

## Event times
  mytimes.C <- rweibull(
    n.C, shape = shape.C.T,
    scale = median.control/((log(2))^(1/shape.C.T)) )
  mytimes.T <- rweibull(
    n.C, shape = shape.C.T,
    scale = (median.control^shape.C.T / (HR*HR.var*log(2)))^(1/shape.C.T) )
  # Combining event times of both groups
  mytimes <- c(mytimes.C, mytimes.T)

## Censoring with combination of administrative censoring and specific censoring
## rate
  # Administrative censoring
  cens.times.AC <- runif(n.C+n.T, min = 0, max = accrual.time) + fu.time

  # Help Data Frame
  Data.help <- data.frame(id = 1:(n.C + n.T), T = mytimes,
                          Z = factor(c(rep("Control",n.C), rep("Treatment",n.T))),
                          C = cens.times.AC)
  Data.help$status <- factor(as.numeric(Data.help$T <= Data.help$C))
  Data.help$TandC  <- pmin(Data.help$T,Data.help$C)

  # Specific censoring rate of the remaining events "mytimes.AC"
    mytimes.AC    <- Data.help$T[which(Data.help$status==1)]
    # calculating needed censoring rate to receive overall censoring rate of
    # cens.rate
    current.CensoringRate <- 1 # of Data.help$status==1 --> always 100%
    needed.CensoringRate  <-
      (cens.rate * dim(Data.help)[1] -
         sum(Data.help$status==0)*current.CensoringRate) /
      sum(Data.help$status==1)
    # if administrative censoring rate is already larger than wanted censoring
    # rate of cens.rate, no specific censoring rate is added
    if (needed.CensoringRate > 0){
      cens.times.SC <-
        sapply(mytimes.AC,
               function(mytimes.AC, needed.CensoringRate){
                 rexp(1, rate = -log(1-needed.CensoringRate)/mytimes.AC)
               },
               needed.CensoringRate)
    }else{
      # if administrative censoring rate is already larger than wanted censoring
      # rate of cens.rate, no specific censoring rate is added
      # --> Setting censoring added censoring rates to Infinity
      # --> event rates are always smaller than specific censoring rates
      cens.times.SC <- rep(Inf, length(mytimes.AC))
    }

## Final data frame
  Data <- Data.help
  Data$C[which(Data$status==1)] <- cens.times.SC
  Data$status <- factor(as.numeric(Data$T <= Data$C))
  Data$TandC  <- pmin(Data$T,Data$C)

  Data.out <- data.frame(Data$TandC, Data$status, Data$Z)

## Returning data
return(Data.out)
}

#------------ Data Generation of Gompertz distributed failure times ------------
DataGomp <- function(n.T, n.C, median.control, accrual.time, fu.time, HR, HR.var,
                     cens.rate, shape.C.T, index.seed) {
## Generating of survival data with Gompertz distributed failure times
## with censoring (administrative and specific censoring)
## Parameters:
  # n.T:            number of patients in treatment group
  # n.C:            number of patients in control group
  # median.control: median survival time of the control group
  # accrual.time:   accrual time
  # fu.time:        follow-up time
  # HR:             hazard ratio of treatment against control
  # HR.var:         variation of the true HR and design HR.
  #                     0 --> design.HR = true.HR
  # cens.rate:      censoring rate, censoring exponential distributed
  # shape.C.T:      shape parameter for control and treatment group
  # index.seed:     seed for the simulation

## Setting seed
  set.seed(index.seed)

## Event times
  mytimes.C <-
    rgompertz(n.C, shape = shape.C.T,
              rate = shape.C.T*log(2) / (exp(median.control*shape.C.T)-1) )
  mytimes.T <-
    rgompertz(n.C, shape = shape.C.T,
              rate = (HR*HR.var)*shape.C.T*log(2) / (exp(median.control*shape.C.T)-1) )
  # Combining event times of both groups
  mytimes   <- c(mytimes.C, mytimes.T)

  ## Censoring with combination of administrative censoring and specific
  ## censoring rate
  # Administrative censoring
  cens.times.AC <- runif(n.C+n.T, min = 0, max = accrual.time) + fu.time

  # Help Data Frame
  Data.help <- data.frame(id = 1:(n.C + n.T), T = mytimes,
                          Z = factor(c(rep("Control",n.C), rep("Treatment",n.T))),
                          C = cens.times.AC)
  Data.help$status <- factor(as.numeric(Data.help$T <= Data.help$C))
  Data.help$TandC  <- pmin(Data.help$T,Data.help$C)

  # Specific censoring rate of the remaining events "mytimes.AC"
      mytimes.AC <- Data.help$T[which(Data.help$status==1)]
      # calculating needed censoring rate to receive overall censoring rate of
      # cens.rate
      current.CensoringRate <- 1 # of Data.help$status==1 --> always 100%
      needed.CensoringRate  <-
        (cens.rate * dim(Data.help)[1] -
           sum(Data.help$status==0)*current.CensoringRate) /
        sum(Data.help$status==1)
      # if administrative censoring rate is already larger than wanted censoring
      # rate of cens.rate, no specific censoring rate is added
      if (needed.CensoringRate > 0){
        cens.times.SC <-
          sapply(mytimes.AC,
                 function(mytimes.AC, needed.CensoringRate){
                   rexp(1, rate = -log(1-needed.CensoringRate)/mytimes.AC)
                 },
                 needed.CensoringRate)
      }else{
        # if administrative censoring rate is already larger than wanted
        # censoring rate of cens.rate, no specific censoring rate is added
        # --> Setting censoring added censoring rates to Infinity
        # --> event rates are always smaller than specific censoring rates
        cens.times.SC <- rep(Inf, length(mytimes.AC))
      }

## Final data frame
  Data <- Data.help
  Data$C[which(Data$status==1)] <- cens.times.SC
  Data$status <- factor(as.numeric(Data$T <= Data$C))
  Data$TandC  <- pmin(Data$T,Data$C)

  Data.out = data.frame(Data$TandC, Data$status, Data$Z)

## Returning data
  return(Data.out)
}

#---------- Data Generation of exponential distributed failure times -----------
#-------------------------- (non-proportional hazards) -------------------------
DataEXPNonProp = function(n.T, n.C, median.control, accrual.time, fu.time, HR,
                          HR.var, effect.start.T, cens.rate, index.seed) {
## Generating of survival data with exponential distributed failure times
## with censoring (administrative and specific censoring) and non-proportional
## hazards (late treatment effect)
## Parameters:
  # n.T:            number of patients in treatment group
  # n.C:            number of patients in control group
  # median.control: median survival time of the control group
  # accrual.time:   accrual time
  # fu.time:        follow-up time
  # HR:             hazard ratio of treatment against control
  # HR.var:         variation of the true HR and design HR.
  #                     0 --> design.HR = true.HR
  # effect.start.T: start of treatment effect of treatment group
  #                 (late treatment effect)
  # cens.rate:      censoring rate, censoring exponential distributed
  # index.seed:     seed for the simulation

## Setting seed
  set.seed(index.seed)

## Event times
  # Control group
  mytimes.C <- rexp(n.C, rate = log(2)/median.control)
  # Treatment group using inversion method
  # Inverse of probability distribution F
  Inv.T <- function(y, lambda.C, lambda.T, treat.start){
    return(
      ifelse(y <= 1-exp(-lambda.C*treat.start),
             -log(1-y)/lambda.C,
             -(log(1-y) + lambda.C*treat.start)/lambda.T + treat.start)
    )
  }
  U <- runif(n.T, min = 0, max = 1)
  mytimes.T <- Inv.T(y = U,
                     lambda.C = log(2)/median.control,
                     lambda.T = (HR*HR.var)*log(2)/median.control,
                     treat.start = effect.start.T)
  # Combining event times of both groups
  mytimes   <- c(mytimes.C, mytimes.T)



## Censoring with combination of administrative censoring and specific censoring
## rate
  # Administrative censoring
  cens.times.AC <- runif(n.C+n.T, min = 0, max = accrual.time) + fu.time

  # Help data frame
  Data.help <- data.frame(id = 1:(n.C + n.T), T = mytimes,
                          Z = factor(c(rep("Control",n.C), rep("Treatment",n.T))),
                          C = cens.times.AC)
  Data.help$status <- factor(as.numeric(Data.help$T <= Data.help$C))
  Data.help$TandC  <- pmin(Data.help$T,Data.help$C)

  # Specific censoring rate of the remaining events "mytimes.AC"
    mytimes.AC    <- Data.help$T[which(Data.help$status==1)]
    # calculating needed censoring rate to recieve overall censoring rate of
    # cens.rate
    current.CensoringRate <- 1 # of Data.help$status==1 --> always 100%
    needed.CensoringRate  <-
      (cens.rate * dim(Data.help)[1] -
         sum(Data.help$status==0)*current.CensoringRate) /
      sum(Data.help$status==1)
    # if administrative censoring rate is already larger than wanted censoring
    # rate of cens.rate, no specific censoring rate is added
    if (needed.CensoringRate > 0){
      cens.times.SC <-
        sapply(mytimes.AC,
               function(mytimes.AC, needed.CensoringRate){
                 rexp(1, rate = -log(1-needed.CensoringRate)/mytimes.AC)
               },
               needed.CensoringRate)
    }else{
      # if administrative censoring rate is already larger than wanted censoring
      # rate of cens.rate, no specific censoring rate is added
      # --> Setting censoring added censoring rates to Infinity
      # --> event rates are always smaller than specific censoring rates
      cens.times.SC <- rep(Inf, length(mytimes.AC))
    }

## Final data frame
  Data <- Data.help
  Data$C[which(Data$status==1)] <- cens.times.SC
  Data$status <- factor(as.numeric(Data$T <= Data$C))
  Data$TandC  <- pmin(Data$T,Data$C)

  Data.out <- data.frame(Data$TandC, Data$status, Data$Z)

## Returning data
  return(Data.out)
}




