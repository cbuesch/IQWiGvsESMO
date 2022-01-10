################################################################################
######                Further illustrations and information               ######
################################################################################

#------------------------------ Preparations ------------------------------------------
## Packages
library(survival)
library(tidyverse)
library(foreach)
library(doParallel)

## Number of used cores for Simulations
num.cl <- 45

## Data loading
results_final <- readRDS("Results.rds")

#--------------------------- Figure 2 with extra row ----------------------------------
# Figure 2 with an additional figure which has the same layout but displays the absolute
# difference.
results_final %>%
   filter(Scenario %in%  c("Scen_1") & cens.rate %in% c(0.2, 0.6) & beta == 0.1) %>%
   dplyr::select(beta, HR, med.C, Scenario, shape.C.T, HR.var, cens.rate,
                 pow, perc.IQWiG.max, perc.Mod.IQWiG.max, perc.ESMO.max,
                 perc.ESMO.RB.max
   ) %>%
   gather("methods", "value", -beta, -HR, -med.C , -shape.C.T, -Scenario, -HR.var,
          -pow, -cens.rate) %>%
   mutate(
      methods = factor(methods, levels(factor(methods))[c(1,2,3,4)],
                       labels = c("ESMO", "ESMO.RB", "Mod.IQWIG", "IQWIG.HR")),

      med.C.text = paste0("Median Control: ", med.C),
      pow.text   = paste0(100*(1-as.numeric(as.character(beta))),"% Power"),
      text_together = factor(paste0(pow.text, ", ", med.C.text),
                             levels = c("90% Power, Median Control: 6",
                                        "90% Power, Median Control: 12",
                                        "90% Power, Median Control: 18",
                                        "90% Power, Median Control: 24",
                                        "90% Power, Median Control: 30"
                             ),
                             labels = c("90% Power, Median Control: 6",
                                        "90% Power, Median Control: 12",
                                        "90% Power, Median Control: 18",
                                        "90% Power, Median Control: 24",
                                        "90% Power, Median Control: 30"
                             )
      ),
      cens.rate = factor(cens.rate, levels = c(0.1, 0.2, 0.4, 0.6),
                         labels = c("10%", "20%", "40%", "60%")),
      value = value*100,
      barlabel = paste0(round(value,2), "%"),
      HR = as.numeric(as.character(HR))
   ) %>%
   ggplot() +
   geom_line(aes(x=HR, y=value, colour = methods, linetype = cens.rate,
                 group = interaction(methods, cens.rate)), stat="identity",
             size = 0.4) +
   geom_point(aes(x=HR, y=value, colour = methods, shape = cens.rate,
                  group = interaction(methods, cens.rate)),
              size = 0.9) +
   # Setting tick intervals on x-axis
   scale_x_continuous(breaks = seq(0.3, 0.9, by = 0.06),
                      minor_breaks = seq(0.3, 0.9, by = 0.02)) +
   #Legend and changing colors
   scale_colour_manual(name = "Percentages of",
                       values = c("blueviolet", "cyan", "green4", "green"),
                       labels = c("maximal ESMO scores",
                                  expression(maximal~ESMO[RB]~scores),
                                  expression(maximal~IQWIG~scores),
                                  expression("maximal"~"Mod-IQWIG"[HR]~"scores"))) +
   #Labels
   labs(title = "",
        subtitle = "",
        linetype = "Censoring rate",
        shape = "Censoring rate",
        x = "designHR",
        y = "Percentages of maximal scores") +
   scale_y_continuous(limits=c(0, 101), breaks = seq(0,100,by=20)) +
   theme_bw() +
   theme(
      text = element_text(#size = 5,
         family = "serif" # TT Times New Roman
         #family = "sans" # TT Arial
      ),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      #legend.position = "none",
      legend.text.align = 0, # align legend to the left
      #axis.ticks = element_line(size=0.2),
      #axis.ticks.length = unit(0.05, "cm"),
      strip.background = element_rect(colour="black",size=0.4), # change line width of label of facet_wrap
      #strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")) # change label width of facet_wrap
      #axis.ticks.x = element_blank()
   ) +
   facet_wrap(text_together~., nrow = 1) -> Figure2


results_final %>%
   filter(Scenario %in%  c("Scen_1") & cens.rate %in% c(0.2, 0.6) & beta == 0.1) %>%
   mutate(
      Diff_IQWiG     = perc.IQWiG.max - perc.ESMO.max,
      Diff_Mod.IQWiG = perc.Mod.IQWiG.max - perc.ESMO.max,
      Diff_ESMO.RB   = perc.ESMO.RB.max - perc.ESMO.max,
   ) %>%
   dplyr::select(beta, HR, med.C, Scenario, shape.C.T, HR.var, cens.rate,
                 pow, Diff_IQWiG, Diff_Mod.IQWiG, Diff_ESMO.RB
   ) %>%
   gather("methods", "value", -beta, -HR, -med.C , -shape.C.T, -Scenario, -HR.var,
          -pow, -cens.rate) %>%
   mutate(
      methods = factor(methods),
      med.C.text = paste0("Median Control: ", med.C),
      pow.text   = paste0(100*(1-as.numeric(as.character(beta))),"% Power"),
      text_together = factor(paste0(pow.text, ", ", med.C.text),
                             levels = c("90% Power, Median Control: 6",
                                        "90% Power, Median Control: 12",
                                        "90% Power, Median Control: 18",
                                        "90% Power, Median Control: 24",
                                        "90% Power, Median Control: 30"
                             ),
                             labels = c("90% Power, Median Control: 6",
                                        "90% Power, Median Control: 12",
                                        "90% Power, Median Control: 18",
                                        "90% Power, Median Control: 24",
                                        "90% Power, Median Control: 30"
                             )
      ),
      cens.rate = factor(cens.rate, levels = c(0.2, 0.6),
                         labels = c("20%", "60%")),
      value = value*100,
      barlabel = paste0(round(value,2), "%"),
      HR = as.numeric(as.character(HR))
   ) %>%
   ggplot() +
   geom_line(aes(x=HR, y=value, colour = methods, linetype = cens.rate,
                 group = interaction(methods, cens.rate)), stat="identity",
             size = 0.4) +
   geom_point(aes(x=HR, y=value, colour = methods, shape = cens.rate,
                  group = interaction(methods, cens.rate)),
              size = 0.9) +
   # Setting tick intervals on x-axis
   scale_x_continuous(breaks = seq(0.3, 0.9, by = 0.06),
                      minor_breaks = seq(0.3, 0.9, by = 0.02)) +
   #Legend and changing colours
   scale_colour_manual(name = "Difference of maximal scoring\npercentages between",
                       values = c("cyan", "green4", "green"),
                       labels = c(expression(ESMO~and~ESMO[RB]),
                                  expression(ESMO~and~IQWIG),
                                  expression(ESMO~and~"Mod-IQWIG"[HR]))) +
   #Labels
   labs(title = "",
        subtitle = "",
        linetype = "Censoring rate",
        shape = "Censoring rate",
        x = "designHR",
        y = "Difference of maximal scoring percentages") +
   scale_y_continuous(limits=c(-90, 50), breaks = seq(-100,50,by=20)) +
   theme_bw() +
   theme(
      text = element_text(#size = 5,
         family = "serif" # TT Times New Roman
         #family = "sans" # TT Arial
      ),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      #legend.position = "none",
      legend.text.align = 0, # align legend to the left
      #axis.ticks = element_line(size=0.2),
      #axis.ticks.length = unit(0.05, "cm"),
      strip.background = element_rect(colour="black",size=0.4), # change line width of label of facet_wrap
      #strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")) # change label width of facet_wrap
      #axis.ticks.x = element_blank()
   ) +
   facet_wrap(text_together~., nrow = 1) -> Figure2_Diff

ggarrange(Figure2, Figure2_Diff, nrow = 2) -> Figure2_comb
Figure2_comb



#---------------------- Bias due to censoring mechanism -------------------------------
##### Generate data of Standard Scenario 1 with no sample size calculation (N=400)
n_sim = 10000

# Setting working directory (where functions and Simulation.seed are saved)
setwd("")

# Seed for simulation
load("Simulation.seed.RData")
# Functions for data generation
source("0_CostumeFunctions_DataGeneration.R")

# setting parameters
med.C          <- c(6,12,18,24,30) # in months
HR             <- seq(0.3, 0.9, by = 0.02)
HR.var         <- 1
beta           <- 0.1
alpha          <- 0.05
r              <- 1
cens.rate      <- c(0.2, 0.4, 0.6)
accrual.time   <- 24 # in months --> 2 years
fu.time        <- NA # will set to 2*med.C later


# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, med.C, HR.var, HR,
                          n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Setting sample size for each sub-scenario and adding it to the parameter matrix
parameters$n.T <- 200
parameters$n.C <- 200

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", "r",
                          "med.C", "HR.var", "HR", "n_sim", "n.T", "n.C")

#### Data generation and analysis
start.time.1 <- Sys.time()

# setting up parallel computing
cl <- makeCluster(num.cl)
registerDoParallel(cl)

# starting parallel computing
Data.out.allforBias <- NULL
results_CensMechBias <- foreach(para = iter(parameters, by='row'),
                               .combine = rbind,
                               .packages = c("tidyverse", "survival")) %dopar% {
   data.scen.i <- rep(list(NA), para$n_sim)
   summary.i   <- data.frame(matrix(NA, ncol = 9, nrow = para$n_sim))
   colnames(summary.i) <- c("HR.point", "HR.CI.low", "HR.CI.up",
                            "Prop.Death.Treatment", "sig",
                            "Cens.T.admin.i", "Cens.C.admin.i", "Cens.T.specific.i", "Cens.C.specific.i")
   for (i in 1:para$n_sim) {
      # Generating Data "n_sim" times for each sub-scenario
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

      # Analysis of the times simulated trial
      cox.sum.i.help <- summary(coxph(Surv(Data.TandC, Data.status != 0) ~ Data.Z,
                                      data = data.scen.i[[i]]))
      Prop.Death.Treatmenti <-
         sum(data.scen.i[[i]]$Data.status[which(data.scen.i[[i]]$Data.Z=="Treatment")]==1)/
         sum(data.scen.i[[i]]$Data.Z=="Treatment")

      Cens.T.admin.i    <- sum(data.scen.i[[i]]$Data.statusDetailed[
         which(data.scen.i[[i]]$Data.Z=="Treatment")]=="Admin")/sum(data.scen.i[[i]]$Data.Z=="Treatment")
      Cens.C.admin.i    <- sum(data.scen.i[[i]]$Data.statusDetailed[
         which(data.scen.i[[i]]$Data.Z=="Control")]=="Admin")/sum(data.scen.i[[i]]$Data.Z=="Control")
      Cens.T.specific.i <- sum(data.scen.i[[i]]$Data.statusDetailed[
         which(data.scen.i[[i]]$Data.Z=="Treatment")]=="Specific")/sum(data.scen.i[[i]]$Data.Z=="Treatment")
      Cens.C.specific.i <- sum(data.scen.i[[i]]$Data.statusDetailed[
         which(data.scen.i[[i]]$Data.Z=="Control")]=="Specific")/sum(data.scen.i[[i]]$Data.Z=="Control")


      summary.i[i,] <- c(cox.sum.i.help$conf.int[c(1,3,4)], Prop.Death.Treatmenti,
                         ifelse(cox.sum.i.help$waldtest[[3]]>0.05 | cox.sum.i.help$conf.int[3]>1, 0,
                                ifelse(sum(data.scen.i[[i]]$Data.status[data.scen.i[[i]]$Data.Z=="Treatment"]==0)==0 ||
                                          sum(data.scen.i[[i]]$Data.status[data.scen.i$Data.Z=="Control"]==0)==0,
                                       2,1)),
                         Cens.T.admin.i, Cens.C.admin.i, Cens.T.specific.i, Cens.C.specific.i)
   }

   # Combining results
   Data.out.allforBias <- Data.out.allforBias %>%
      rbind(
         tibble(
            beta         = para$beta,
            accrual.time = para$accrual.time,
            fu.time      = para$fu.time,
            cens.rate    = para$cens.rate,
            r            = para$r,
            med.C        = para$med.C,
            HR.var       = para$HR.var,
            HR           = para$HR,
            trueHR       = para$HR.var*para$HR,
            Diff_HRestimate_trueHR = summary.i$HR.point - para$HR.var*para$HR
         )
      )

   Data.out <- tibble(
      beta         = para$beta,
      accrual.time = para$accrual.time,
      fu.time      = para$fu.time,
      cens.rate    = para$cens.rate,
      r            = para$r,
      med.C        = para$med.C,
      HR.var       = para$HR.var,
      HR           = para$HR,
      trueHR       = para$HR.var*para$HR,

      mean.HR.pointALL  = mean(summary.i$HR.point),
      mean.HR.CI.lowALL = mean(summary.i$HR.CI.low),
      mean.HR.CI.upALL  = mean(summary.i$HR.CI.up),

      mean.HR.pointSIG  = mean(
         ifelse(summary.i$sig %in% c(1,2), summary.i$HR.point, NA),
         na.rm = TRUE),
      mean.HR.CI.lowSIG = mean(
         ifelse(summary.i$sig %in% c(1,2), summary.i$HR.CI.low, NA),
         na.rm = TRUE),
      mean.HR.CI.upSIG  = mean(
         ifelse(summary.i$sig %in% c(1,2), summary.i$HR.CI.up, NA),
         na.rm = TRUE),

      mean.Cens.T.admin.i = mean(summary.i$Cens.T.admin.i, na.rm = TRUE),
      mean.Cens.C.admin.i = mean(summary.i$Cens.C.admin.i, na.rm = TRUE),
      mean.Cens.T.specific.i = mean(summary.i$Cens.T.specific.i, na.rm = TRUE),
      mean.Cens.C.specific.i = mean(summary.i$Cens.C.specific.i, na.rm = TRUE),

      mean.Prop.Death.TreatmentALL = mean(summary.i$Prop.Death.Treatment),
      mean.Prop.Death.TreatmentSIG = mean(
         ifelse(summary.i$sig %in% c(1,2), summary.i$Prop.Death.Treatment, NA),
         na.rm = TRUE)
   )
   return(list(Data.out, Data.out.allforBias))
}

stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)

## Data processing
results_Scen1_CensMechBias_allforBias_SameN <- bind_rows(results_CensMechBias[,2])
results_CensMechBias            <- bind_rows(results_CensMechBias[,1])


##### Illustrating results
results_Scen1_CensMechBias_allforBias_SameN %>%
   filter(beta==0.1, cens.rate %in% c(0.2,0.4,0.6)) %>%
   group_by(cens.rate, med.C) %>%
   summarize(
      beta         = first(beta),
      accrual.time = first(accrual.time),
      fu.time      = first(fu.time),
      cens.rate    = first(cens.rate),
      r            = first(r),
      med.C        = first(med.C),
      HR.var       = first(HR.var),
      HR           = first(HR),
      trueHR       = first(trueHR),

      n            = n(),
      Sum_Diff     = sum(Diff_HRestimate_trueHR),
      Sum_Diff_quad = sum(Diff_HRestimate_trueHR^2)
   ) %>%
   ungroup() %>%
   mutate(
      Bias    = Sum_Diff/n,
      SD_Bias = sqrt(1/(n*(n-1)) * Sum_Diff_quad)
   ) -> df_table_Bias

# Overall bias
results_Scen1_CensMechBias_allforBias_SameN %>%
   filter(beta==0.1, cens.rate %in% c(0.2,0.4,0.6)) -> help

1/(dim(help)[1]) * sum(help$Diff_HRestimate_trueHR)
1/(dim(help)[1]*(dim(help)[1]-1)) * sum(help$Diff_HRestimate_trueHR^2)

# Plot
results_Scen1_CensMechBias_allforBias_SameN %>%
   filter(beta==0.1, cens.rate %in% c(0.2,0.4,0.6)) %>%
   group_by(cens.rate, med.C, trueHR) %>%
   summarize(
      beta         = first(beta),
      accrual.time = first(accrual.time),
      fu.time      = first(fu.time),
      cens.rate    = first(cens.rate),
      r            = first(r),
      med.C        = first(med.C),
      HR.var       = first(HR.var),
      HR           = first(HR),
      trueHR       = first(trueHR),

      n            = n(),
      Sum_Diff     = sum(Diff_HRestimate_trueHR),
      Sum_Diff_quad = sum(Diff_HRestimate_trueHR^2)
   ) %>%
   ungroup() %>%
   mutate(
      Bias    = Sum_Diff/n,
      SD_Bias = sqrt(1/(n*(n-1)) * Sum_Diff_quad)
   ) -> df_plot_Bias

df_plot_Bias %>%
   mutate(
      med.C.text = factor(paste0("Median Control: ", med.C),
                          levels = c("Median Control: 6", "Median Control: 12",
                                     "Median Control: 18", "Median Control: 24",
                                     "Median Control: 30"),
                          labels = c("Median Control: 6", "Median Control: 12",
                                     "Median Control: 18", "Median Control: 24",
                                     "Median Control: 30"))
   ) %>%
   ggplot(aes(x = trueHR, y=Bias, col=factor(cens.rate))) +
   geom_line() +
   scale_x_continuous(breaks = seq(0.3, 0.9, by = 0.06),
                      minor_breaks = seq(0.3, 0.9, by = 0.02)) +
   labs(y = "HR Bias") + guides(col=guide_legend(title="Censoring rate")) +
   theme(
      text = element_text(#size = 5,
         family = "serif" # TT Times New Roman
         #family = "sans" # TT Arial
      ),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      #legend.position = "none",
      #legend.text.align = 0, # align legend to the left
      #axis.ticks = element_line(size=0.2),
      #axis.ticks.length = unit(0.05, "cm"),
      strip.background = element_rect(colour="black",size=0.4), # change line width of label of facet_wrap
      #strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")) # change label width of facet_wrap
      #axis.ticks.x = element_blank()
   ) +
   facet_wrap(~med.C.text) -> p_Bias
p_Bias


#-------------------------- Maturity of the Data --------------------------------------
##### Analysis data of Standard Scenario 1 to calculate the proportion of death in the treatment group
path_data <- "" # e.g. ".../Simulations/Data/"

## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen1_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl)
registerDoParallel(cl)

# Starting parallel computing
Results_Scen1_PropDeathTreatment <- foreach(para = iter(parameters, by='row'),
                              .packages = c("survival", "tidyverse"),
                              .combine = rbind,
                              .multicombine = TRUE) %dopar% {
   # Load data of sub-scenarios with "n_sim" times replicates
   df <- readRDS(paste0(path_data,
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
   cox.sum.1 <- summary(coxph(Surv(Data.TandC, Data.status != 0) ~ Data.Z, data = df[[1]]))
   Prop.Death.Treatment1 <- sum(df[[1]]$Data.status[which(df[[1]]$Data.Z=="Treatment")]==1)/sum(df[[1]]$Data.Z=="Treatment")

   cox.sum.i     <- data.frame(matrix(NA, ncol = 5, nrow = para$it))
   cox.sum.i[1,] <- c(cox.sum.1$conf.int[c(1,3,4)], Prop.Death.Treatment1,
                      ifelse(cox.sum.1$waldtest[[3]]>0.05 | cox.sum.1$conf.int[3]>1, 0,
                             ifelse(sum(df[[1]]$Data.status[df[[1]]$Data.Z=="Treatment"]==0)==0 ||
                                       sum(df[[1]]$Data.status[df$Data.Z=="Control"]==0)==0,
                                    2,1)))
   colnames(cox.sum.i) <- c("HR.point", "HR.CI.low", "HR.CI.up", "Prop.Death.Treatment", "sig")

   for (i in 2:para$it) {
      cox.sum.i.help <- summary(coxph(Surv(Data.TandC, Data.status != 0) ~ Data.Z, data = df[[i]]))
      Prop.Death.Treatmenti <- sum(df[[i]]$Data.status[which(df[[i]]$Data.Z=="Treatment")]==1)/sum(df[[i]]$Data.Z=="Treatment")
      cox.sum.i[i,] <- c(cox.sum.i.help$conf.int[c(1,3,4)], Prop.Death.Treatmenti,
                         ifelse(cox.sum.i.help$waldtest[[3]]>0.05 | cox.sum.i.help$conf.int[3]>1, 0,
                                ifelse(sum(df[[i]]$Data.status[df[[i]]$Data.Z=="Treatment"]==0)==0 ||
                                          sum(df[[i]]$Data.status[df$Data.Z=="Control"]==0)==0,
                                       2,1)))
   }

   Data.out <- tibble(
      beta         = para$beta,
      accrual.time = para$accrual.time,
      fu.time      = para$fu.time,
      cens.rate    = para$cens.rate,
      r            = para$r,
      med.C        = para$med.C,
      HR.var       = para$HR.var,
      HR           = para$HR,
      trueHR       = para$HR.var*para$HR,

      mean.Prop.Death.TreatmentALL = mean(cox.sum.i$Prop.Death.Treatment),
      mean.Prop.Death.TreatmentSIG = mean(
         ifelse(cox.sum.i$sig %in% c(1,2), cox.sum.i$Prop.Death.Treatment, NA),
         na.rm = TRUE)
      )
   return(Data.out)
}
stopCluster(cl)

(time.1 <- Sys.time() - start.time.1)


##### Illustrating results
results_final %>% filter(Scenario == "Scen_1") %>%
   merge(Results_Scen1_PropDeathTreatment,
         by.x = c("beta", "cens.rate", "fu.time", "accrual.time", "med.C", "HR"),
         by.y = c("beta", "cens.rate", "fu.time", "accrual.time", "med.C", "HR")) %>%
   filter(cens.rate %in% c(0.2,0.6) & beta == 0.1) %>%
   mutate(cens.rate.label = paste0(100*as.numeric(as.character(cens.rate)), "%"),

          med.C.text = paste0("Median Control: ", med.C),
          pow.text = paste0(100*(1-as.numeric(as.character(beta))),"% Power"),
          text_together = factor(paste0(pow.text, ", ", med.C.text),
                                 levels = c("90% Power, Median Control: 6",
                                            "90% Power, Median Control: 12",
                                            "90% Power, Median Control: 18",
                                            "90% Power, Median Control: 24",
                                            "90% Power, Median Control: 30"
                                 ),
                                 labels = c("90% Power, Median Control: 6",
                                            "90% Power, Median Control: 12",
                                            "90% Power, Median Control: 18",
                                            "90% Power, Median Control: 24",
                                            "90% Power, Median Control: 30"
                                 )
          ),
          HR = as.numeric(as.character(HR)),
          trueHR_plot = ifelse(!rep(c(rep(FALSE, 1), TRUE), length(trueHR)/2), 0, trueHR)
   ) -> help


elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})
help %>%
   ggplot(aes(x = mean.Prop.Death.TreatmentSIG, y = Cor.IQWiG.ESMO,
              group = interaction(cens.rate.label),
              linetype = cens.rate.label, shape = cens.rate.label)) +
   geom_point() + geom_line() +
   geom_text(data = subset(help, elementwise.all.equal(trueHR, trueHR_plot)),
             aes(x = mean.Prop.Death.TreatmentSIG - 0.03,
                 y = Cor.IQWiG.ESMO + 0.02,
                 label = trueHR),#, color = "UL thresholds"),
             check_overlap = TRUE, size = 2.5) +
   ylim(c(0,1)) +
   scale_x_continuous(limits = c(0.2,0.8), breaks = seq(0, 1, by = 0.1),
                      minor_breaks = seq(0, 1, by = 0.05)) +
   labs(y = "Spearman correlation",
        x = "Proportion of deaths observed in the intervention arm") +
   guides(linetype=guide_legend(title="Censoring rate"),
          shape=guide_legend(title="Censoring rate")) +
   theme(
      text = element_text(#size = 5,
         family = "serif" # TT Times New Roman
         #family = "sans" # TT Arial
      ),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      #legend.position = "none",
      #legend.text.align = 0, # align legend to the left
      #axis.ticks = element_line(size=0.2),
      #axis.ticks.length = unit(0.05, "cm"),
      strip.background = element_rect(colour="black",size=0.4), # change line width of label of facet_wrap
      #strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")) # change label width of facet_wrap
      #axis.ticks.x = element_blank()
   ) +
   facet_wrap(text_together~., nrow=2, scales = "fixed") -> p_MaturityOfData
p_MaturityOfData


