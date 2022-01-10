################################################################################
######                  Visualization of analyzed data                    ######
################################################################################

# Setting working directory (where final results are saved)
setwd("")

#----------------- Needed Packages and further setting ups ---------------------
library(tidyverse)
library(ggpubr)
library(survival)
library(flexsurv)

# Loading results of analyze
results <- readRDS("Results.rds")


#------------ Figure 1: Spearman correlation Standard Scenario 1 ---------------
results %>% filter(Scenario == "Scen_1" &
                      cens.rate %in% c(0.2, 0.6) & beta == 0.1) %>%
   arrange(beta, cens.rate, med.C, HR) %>%
   mutate(cens.rate.label = paste0(100*as.numeric(as.character(cens.rate)), "%"),
          med.C.text = paste0("Median Control: ", med.C),
          pow.text = paste0(100*(1-as.numeric(as.character(beta))),"% Power"),
          text_together = factor(paste0(pow.text, ", ", med.C.text),
                                 levels = c("90% Power, Median Control: 6",
                                            "90% Power, Median Control: 12",
                                            "90% Power, Median Control: 18",
                                            "90% Power, Median Control: 24",
                                            "90% Power, Median Control: 30"),
                                 labels = c("90% Power, Median Control: 6",
                                            "90% Power, Median Control: 12",
                                            "90% Power, Median Control: 18",
                                            "90% Power, Median Control: 24",
                                            "90% Power, Median Control: 30")
          )
   ) %>%
   ggplot(aes(x = HR, y = Cor.IQWiG.ESMO,
              group = interaction(cens.rate.label),
              linetype = cens.rate.label, shape = cens.rate.label)) +
   geom_point(size = 0.2) + geom_line(size = 0.2) +
   labs(title = "", x = "designHR", y = "Spearman correlation",
        colour = "Power", linetype = "Censoring rate", shape = "Censoring rate") +
   # Setting tick intervals on x-axis
   scale_x_continuous(breaks = seq(0.3, 0.9, by = 0.06),
                      minor_breaks = seq(0.3, 0.9, by = 0.02)) +
   # Setting tick intervals on y-axis
   scale_y_continuous(breaks = c(0,0.25,0.5,0.75),
                      minor_breaks = c(0,0.125,0.25,0.375,0.5,0.625,0.75,0.875)) +
   facet_wrap(text_together~., nrow=1, scales = "fixed") +
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
      #axis.ticks.x = element_blank()
   ) -> Figure1
Figure1


#------------------ Figure 2: Line chart Standard Scenario 1 -------------------
results %>%
   filter(Scenario %in%  c("Scen_1") & cens.rate %in% c(0.2, 0.6) & beta == 0.1) %>%
   select(beta, HR, med.C, Scenario, shape.C.T, HR.var, cens.rate, pow,
          perc.IQWiG.max, perc.Mod.IQWiG.max, perc.ESMO.max, perc.ESMO.RB.max
   ) %>%
   gather("methods", "value",
          -beta, -HR, -med.C , -shape.C.T, -Scenario, -HR.var, -pow, -cens.rate) %>%
   mutate(
      med.C.text = paste0("Median Control: ", med.C),
      pow.text   = paste0(100*(1-as.numeric(as.character(beta))),"% Power"),
      text_together = factor(paste0(pow.text, ", ", med.C.text),
                             levels = c("90% Power, Median Control: 6",
                                        "90% Power, Median Control: 12",
                                        "90% Power, Median Control: 18",
                                        "90% Power, Median Control: 24",
                                        "90% Power, Median Control: 30"),
                             labels = c("90% Power, Median Control: 6",
                                        "90% Power, Median Control: 12",
                                        "90% Power, Median Control: 18",
                                        "90% Power, Median Control: 24",
                                        "90% Power, Median Control: 30")
      ),
      cens.rate = factor(cens.rate, levels = c(0.2, 0.6), labels = c("20%", "60%")),
      methods = factor(methods, levels(factor(methods))[c(1,2,3,4)],
                       labels = c("ESMO", "ESMO.RB", "IQWIG.RR", "IQWIG.HR")),
      value = value*100
   ) %>%
   ggplot() +
   geom_line(aes(x=HR, y=value, colour = methods, linetype = cens.rate,
                 group = interaction(methods, cens.rate)), stat="identity",
             size = 0.2) +
   geom_point(aes(x=HR, y=value, colour = methods, shape = cens.rate,
                  group = interaction(methods, cens.rate)),
              size = 0.2) +
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
   labs(title = "", subtitle = "", linetype = "Censoring rate",
        shape = "Censoring rate", x = "designHR",
        y = "Percentages of maximal scores") +
   scale_y_continuous(limits=c(0, 101), breaks = seq(0, 100,by=20)) +
   facet_wrap(text_together~., nrow = 1) +
   theme_bw() +
   theme(
      text = element_text(size = 5,
                          family = "serif" # TT Times New Roman
                          #family = "sans" # TT Arial
      ),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.position = "bottom",
      legend.text.align = 0, # align legend to the left
      axis.ticks = element_line(size=0.2),
      axis.ticks.length = unit(0.05, "cm"),
      strip.background = element_rect(colour="black",size=0.4), # change line width of label of facet_wrap
      strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")) # change label width of facet_wrap
   ) -> Figure2
Figure2


#--------------------- Figure 3: ROC Standard Scenario 1 -----------------------
# This plot is generated in an separate R-Code file: "3_VisualizationsOfResults_ROC.R"


#--------- Figure 4: Bar chart of Scenario 2 (incorrect assumed HR) ------------
results %>%
   filter(Scenario %in%  c("Scen_1", "Scen_2") & cens.rate %in% c(0.2,0.6) &
             beta == 0.1 & med.C == 12 & HR %in% c(0.3, 0.5, 0.7, 0.8, 0.9) &
             HR.var %in% c(0.9, 1, 1.1, 1.2)) %>%
   select(beta, HR.var, HR, med.C, cens.rate, pow,
          perc.IQWiG.max, perc.Mod.IQWiG.max, perc.ESMO.max, perc.ESMO.RB.max) %>%
   gather("methods", "value", -beta, -HR.var, -HR, -med.C, -pow, -cens.rate) %>%
   mutate(
      Methods = factor(methods, levels(factor(methods))[c(1,2,3,4)],
                       labels = c("ESMO", "ESMO.RB", "IQWIG.RR", "IQWIG.HR")),
      Methods_with_cens = factor(
         paste0(methods, cens.rate),
         levels(factor(paste0(methods, cens.rate))),
         labels = c("ESMO (20% cens. rate)",     "ESMO (60% cens. rate)",
                    "ESMO.RB (20% cens. rate)",  "ESMO.RB (60% cens. rate)",
                    "IQWIG.RR (20% cens. rate)", "IQWIG.RR (60% cens. rate)",
                    "IQWIG.HR (20% cens. rate)", "IQWIG.HR (60% cens. rate)")),
      label = factor(
         paste0("HR[var]==", HR.var, "~(designHR==", HR,"*','~trueHR==", HR*HR.var, ")"),
         levels = c(
            "HR[var]==0.9~(designHR==0.3*','~trueHR==0.27)", "HR[var]==0.9~(designHR==0.5*','~trueHR==0.45)",
            "HR[var]==0.9~(designHR==0.7*','~trueHR==0.63)", "HR[var]==0.9~(designHR==0.8*','~trueHR==0.72)",
            "HR[var]==0.9~(designHR==0.9*','~trueHR==0.81)",
            "HR[var]==1~(designHR==0.3*','~trueHR==0.3)",    "HR[var]==1~(designHR==0.5*','~trueHR==0.5)",
            "HR[var]==1~(designHR==0.7*','~trueHR==0.7)",    "HR[var]==1~(designHR==0.8*','~trueHR==0.8)",
            "HR[var]==1~(designHR==0.9*','~trueHR==0.9)",
            "HR[var]==1.1~(designHR==0.3*','~trueHR==0.33)", "HR[var]==1.1~(designHR==0.5*','~trueHR==0.55)",
            "HR[var]==1.1~(designHR==0.7*','~trueHR==0.77)", "HR[var]==1.1~(designHR==0.8*','~trueHR==0.88)",
            "HR[var]==1.1~(designHR==0.9*','~trueHR==0.99)",
            "HR[var]==1.2~(designHR==0.3*','~trueHR==0.36)", "HR[var]==1.2~(designHR==0.5*','~trueHR==0.6)",
            "HR[var]==1.2~(designHR==0.7*','~trueHR==0.84)", "HR[var]==1.2~(designHR==0.8*','~trueHR==0.96)",
            "HR[var]==1.2~(designHR==0.9*','~trueHR==1.08)"
         ),
         labels = c(
            "HR[var]==0.9~(designHR==0.3*','~trueHR==0.27)", "HR[var]==0.9~(designHR==0.5*','~trueHR==0.45)",
            "HR[var]==0.9~(designHR==0.7*','~trueHR==0.63)", "HR[var]==0.9~(designHR==0.8*','~trueHR==0.72)",
            "HR[var]==0.9~(designHR==0.9*','~trueHR==0.81)",
            "HR[var]==1~(designHR==0.3*','~trueHR==0.3)",    "HR[var]==1~(designHR==0.5*','~trueHR==0.5)",
            "HR[var]==1~(designHR==0.7*','~trueHR==0.7)",    "HR[var]==1~(designHR==0.8*','~trueHR==0.8)",
            "HR[var]==1~(designHR==0.9*','~trueHR==0.9)",
            "HR[var]==1.1~(designHR==0.3*','~trueHR==0.33)", "HR[var]==1.1~(designHR==0.5*','~trueHR==0.55)",
            "HR[var]==1.1~(designHR==0.7*','~trueHR==0.77)", "HR[var]==1.1~(designHR==0.8*','~trueHR==0.88)",
            "HR[var]==1.1~(designHR==0.9*','~trueHR==0.99)",
            "HR[var]==1.2~(designHR==0.3*','~trueHR==0.36)", "HR[var]==1.2~(designHR==0.5*','~trueHR==0.6)",
            "HR[var]==1.2~(designHR==0.7*','~trueHR==0.84)", "HR[var]==1.2~(designHR==0.8*','~trueHR==0.96)",
            "HR[var]==1.2~(designHR==0.9*','~trueHR==1.08)"
         )
      ),
      OverUnderPowered_label =
         ifelse(HR.var < 1, "Overpowered",
                ifelse(HR.var>1, "Underpowered", "Planned~power~achieved")
         ),
      value = ifelse(is.na(value), 0, value*100),
      barlabel = paste0(round(value,1), "%")
   ) %>%
   ggplot(aes(x=Methods_with_cens, y=value, fill = Methods)) +
   geom_bar(stat="identity") +
   geom_text(aes(x=Methods_with_cens, y=0, label=barlabel), angle = 90, nudge_y = 25, size=1, fontface = "bold") +
   scale_fill_manual(
      name = "Method",
      values =c("blueviolet", "cyan", "green4", "green"),
      labels = c("ESMO",
                 expression(ESMO[RB]),
                 expression(IQWIG),
                 expression("Mod-IQWIG"[HR]~"scores"))
   ) +
   # x-axis ticks
   scale_x_discrete(labels=c("20%","60%","20%","60%","20%","60%","20%","60%")) +
   #Labels
   labs(title = "",
        x = "Censoring rate",
        y = "Proportion of maximal scores") +
   scale_y_continuous(limits=c(0, 101), breaks = c(0, 30, 60, 90)) +
   facet_wrap(label~OverUnderPowered_label, nrow = 4, labeller = label_parsed) +
   theme_bw() +
   theme(
      text = element_text(size = 4,
                          family = "serif" # TT Times New Roman
                          #family = "sans" # TT Arial
      ),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.position = "bottom",
      legend.text.align = 0, # align legend to the left
      axis.ticks = element_line(size=0.2),
      axis.ticks.length = unit(0.05, "cm"),
      strip.background = element_rect(colour="black",size=0.4), # change line width of label of facet_wrap
      strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")) # change label width of facet_wrap
   ) -> Figure4
Figure4


#------------------ Figure 5: Line chart Scenario 3 and 4 ----------------------
#----------------- (different failure time distributions) ----------------------
results %>%
   filter(Scenario %in%  c("Scen_1", "Scen_3aWEIB", "Scen_3bGOMP", "Scen_4") &
             med.C == 12 & cens.rate %in% c(0.2,0.6) &
             shape.C.T %in% c(NA, 0.5, 1.5, -0.2, 0.2) & beta == 0.1
   ) %>%
   select(beta, HR, med.C, Scenario, shape.C.T, HR.var, cens.rate, pow,
          perc.IQWiG.max, perc.Mod.IQWiG.max, perc.ESMO.max, perc.ESMO.RB.max
   ) %>%
   gather("methods", "value",
          -beta, -HR, -med.C , -shape.C.T, -Scenario, -HR.var, -pow, -cens.rate) %>%
   mutate(
      methods = factor(methods, levels(factor(methods))[c(1,2,3,4)],
                       labels = c("ESMO", "ESMO.RB", "IQWIG.RR", "IQWIG.HR")),
      Scen.text = factor(
         ifelse(Scenario == "Scen_1",
                paste0("Exp. prop. haz."),
                ifelse(Scenario == "Scen_4",
                       paste0("Piece-wise exp. (non. prop. haz.)"),
                       ifelse(Scenario == "Scen_3aWEIB",
                              paste0("Weib. prop. haz. (shape=", shape.C.T ,")"),
                              paste0("Gomp. prop. haz. (shape=", shape.C.T ,")")
                       )
                )
         ),
         labels = c("Exp. prop. haz.",
                    "Weib. prop. haz. (shape=0.5)", "Weib. prop. haz. (shape=1.5)",
                    "Gomp. prop. haz. (shape=-0.1)", "Gomp. prop. haz. (shape=0.1)",
                    "Gomp. prop. haz. (shape=-0.2)", "Gomp. prop. haz. (shape=0.2)",
                    "Piece-wise exp. (non. prop. haz.)"),
         levels = c("Exp. prop. haz.",
                    "Weib. prop. haz. (shape=0.5)", "Weib. prop. haz. (shape=1.5)",
                    "Gomp. prop. haz. (shape=-0.1)", "Gomp. prop. haz. (shape=0.1)",
                    "Gomp. prop. haz. (shape=-0.2)", "Gomp. prop. haz. (shape=0.2)",
                    "Piece-wise exp. (non. prop. haz.)")
      ),
      cens.rate = factor(cens.rate, levels = c(0.2, 0.6), labels = c("20%", "60%")),
      med.C.text = factor(med.C, levels = levels(med.C),
                          labels = paste0("med.C=", levels(med.C))),
      pow.text = factor(beta, levels = levels(beta),
                        labels = paste0(100*(1-as.numeric(levels(beta))),"% Power")),
      value = value*100
   ) %>%
   ggplot() +
   geom_line(aes(x=HR, y=value, colour = methods, linetype = cens.rate, group = interaction(methods, cens.rate)), stat="identity", size = 0.2) +
   geom_point(aes(x=HR, y=value, colour = methods, shape = cens.rate, group = interaction(methods, cens.rate)), size = 0.2) +
   # Setting tick intervals on x-axis
   scale_x_continuous(breaks = seq(0.3, 0.9, by = 0.06), minor_breaks = seq(0.3, 0.9, by = 0.02)) +
   #Legend and changing colours
   scale_colour_manual(name = "Percentages of",
                       values = c("blueviolet", "cyan", "green4", "green"),
                       labels = c("maximal ESMO scores",
                                  expression(maximal~ESMO[RB]~scores),
                                  expression(maximal~IQWIG~scores),
                                  expression("maximal"~"Mod-IQWIG"[HR]~"scores"))) +
   #Labels
   labs(title = "", subtitle = "", linetype = "Censoring rate",
        shape = "Censoring rate", x = "designHR",
        y = "Percentages of maximal scores") +
   scale_y_continuous(limits=c(0, 101), breaks = seq(0,100,by=20)) +
   facet_wrap(pow.text~Scen.text, nrow = 1, scales = "fixed") +
   theme_bw() +
   theme(
      text = element_text(size = 4,
                          family = "serif" # TT Times New Roman
                          #family = "sans" # TT Arial
      ),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.position = "bottom",
      legend.text.align = 0, # align legend to the left
      axis.ticks = element_line(size=0.2),
      axis.ticks.length = unit(0.05, "cm"),
      strip.background = element_rect(colour="black",size=0.3), # change line width of label of facet_wrap
      strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm")) # change label width of facet_wrap
   ) -> Figure5
Figure5


#----------- APPENDIX: Hazard functions for exponential, Weibull ---------------
#---------------- and Gompertz distribution for Appendix -----------------------
# Functions needed for calculation of hazard function
Hazard_Weib_C <- function(shape.C.T.WEIB, time_points, median.control){
   return(
      shape.C.T.WEIB *
         median.control/((log(2))^(1/shape.C.T.WEIB)) *
         (time_points*median.control/((log(2))^(1/shape.C.T.WEIB)))^(shape.C.T.WEIB-1)
   )
}
Hazard_Weib_T <- function(shape.C.T.WEIB, time_points, median.control, HR, HR.var){
   return(
      shape.C.T.WEIB *
         (median.control^shape.C.T.WEIB / (HR*HR.var*log(2)))^(1/shape.C.T.WEIB) *
         (time_points*(median.control^shape.C.T.WEIB / (HR*HR.var*log(2)))^(1/shape.C.T.WEIB))^(shape.C.T.WEIB-1)
   )
}
Hazard_Gomp_C <- function(shape.C.T.GOMP, time_points, median.control){
   return(
      hgompertz(
         x=time_points,
         shape = shape.C.T.GOMP,
         rate = shape.C.T.GOMP*log(2) / (exp(median.control*shape.C.T.GOMP)-1)
      )
   )
}
Hazard_Gomp_T <- function(shape.C.T.GOMP, time_points, median.control, HR, HR.var){
   return(
      hgompertz(
         x=time_points,
         shape = shape.C.T.GOMP,
         rate = (HR*HR.var)*shape.C.T.GOMP*log(2) / (exp(median.control*shape.C.T.GOMP)-1)
      )
   )
}

# Assumed parameters
HR <- 0.9
HR.var <- 1
median.control <- 6
time_points <- seq(1,
                   40,
                   by = 0.2)

# Data frame for plot
df_plot <- data.frame(
   t            = rep(time_points, times=2),
   grp          = c(rep("Control", length(time_points)),
                    rep("Treatment", length(time_points))),
   exp          = c(rep(log(2)/median.control, length(time_points)),
                    rep((HR*HR.var)*log(2)/median.control, length(time_points))),
   weib_shape0_5        = c(Hazard_Weib_C(shape.C.T.WEIB = 0.5,
                                          time_points = time_points,
                                          median.control = median.control),
                            Hazard_Weib_T(shape.C.T.WEIB = 0.5,
                                          time_points = time_points,
                                          median.control = median.control,
                                          HR = HR, HR.var = HR.var)),
   weib_shape1_5        = c(Hazard_Weib_C(shape.C.T.WEIB = 1.5,
                                          time_points = time_points,
                                          median.control = median.control),
                            Hazard_Weib_T(shape.C.T.WEIB = 1.5,
                                          time_points = time_points,
                                          median.control = median.control,
                                          HR = HR, HR.var = HR.var)),
   gomp_shape0_2        = c(Hazard_Gomp_C(shape.C.T.GOMP = 0.2,
                                          time_points = time_points,
                                          median.control = median.control),
                            Hazard_Gomp_T(shape.C.T.GOMP = 0.2,
                                          time_points = time_points,
                                          median.control = median.control,
                                          HR = HR, HR.var = HR.var)),
   gomp_shape_Minus_0_2 = c(Hazard_Gomp_C(shape.C.T.GOMP = -0.2,
                                          time_points = time_points,
                                          median.control = median.control),
                            Hazard_Gomp_T(shape.C.T.GOMP = -0.2,
                                          time_points = time_points,
                                          median.control = median.control,
                                          HR = HR, HR.var = HR.var))
) %>%
   pivot_longer(cols = c("exp", "weib_shape0_5", "weib_shape1_5",
                         "gomp_shape0_2", "gomp_shape_Minus_0_2"),
                names_to  = "dist",
                values_to = "hazard") %>%
   mutate(
      dist = factor(dist,
                    levels = c("exp", "weib_shape0_5", "weib_shape1_5",
                               "gomp_shape0_2", "gomp_shape_Minus_0_2"),
                    labels = c("Exponential",
                               "Weibull\n(shape = 0.5)", "Weibull\n(shape = 1.5)",
                               "Gompertz\n(shape = 0.2)", "Gompertz\n(shape = -0.2)"))
   )

# Plots
df_plot %>%
   filter(
      dist %in% c("Exponential", "Weibull\n(shape = 1.5)",
                  "Gompertz\n(shape = 0.2)"
      ) & grp == "Treatment"
   ) %>%
   ggplot(aes(x=t, y=hazard, col = dist)) +
   geom_line(size=0.4) +
   labs(x = "Time (in months)", y = "Hazard function") +
   scale_color_manual("Failure time\ndistribution",
                      values = c("black", "red", "blue")) +
   theme_bw() +
   theme(
      text = element_text(size = 8,
                          family = "serif" # TT Times New Roman
                          #family = "sans" # TT Arial
      ),
      legend.position = "bottom",
      axis.ticks = element_line(size=0.4),
      axis.ticks.length = unit(0.05, "cm")
   ) -> hazard_increasing


df_plot %>%
   filter(
      dist %in% c("Exponential", "Weibull\n(shape = 0.5)",
                  "Gompertz\n(shape = -0.2)"
      ) & grp == "Treatment"
   ) %>%
   ggplot(aes(x=t, y=hazard, col = dist)) +
   geom_line(size=0.4) +
   labs(x = "Time (in months)", y = "Hazard function") +
   scale_color_manual("", values = c("black", "orangered4", "royalblue4")) +
   theme_bw() +
   theme(
      text = element_text(size = 8,
                          family = "serif" # TT Times New Roman
                          #family = "sans" # TT Arial
      ),
      legend.position = "bottom",
      axis.ticks = element_line(size=0.4),
      axis.ticks.length = unit(0.05, "cm")
   ) -> hazard_decreasing

ggarrange(hazard_increasing, hazard_decreasing) -> AppendixFigure1_hazards
AppendixFigure1_hazards


#----------- APPENDIX: Survival function of late treatment effect --------------
# Functions:
Survival_Treatment <- function(x, l_c, l_t, c){
   # x:   time points
   # l_c: exponential parameter of control group (and treatment group until c)
   # l_t: exponential parameter of treatment group from c onwoards
   # c:   starting point of treatment effect in treatment group
   return(
      ifelse(x <= c,
             exp(-l_c*x),
             exp(-l_c*c) * exp(-l_t*(x-c)))
   )
}

# Assumed parameters
HR <- 0.7
HR.var <- 1
median.control <- 12
start_time <- 1/3*median.control
time_points <- seq(1,
                   50,
                   by = 0.2)
df_plot_Late_treatment <- data.frame(
   t            = rep(time_points, times=2),
   grp          = c(rep("Control: exponential failure times",
                        length(time_points)),
                    rep("Treatment: piece-wise exponential failure times",
                        length(time_points))),
   non.prop.exp = c(1-pexp(q=time_points, rate = log(2)/median.control),
                    Survival_Treatment(x=time_points,
                                       l_c = log(2)/median.control,
                                       l_t = (HR*HR.var)*log(2)/median.control,
                                       c=start_time))
) %>%
   mutate(non.prop.exp_100 = non.prop.exp*100)

# Data frame for plot
df_plot_Late_treatment %>%
   ggplot(aes(x=t, y=non.prop.exp_100, col = grp)) +
   geom_line(size=0.4) +
   labs(x = "Time (in months)", y = "Survival function") +
   scale_y_continuous(breaks = seq(0,100, by=20)) +
   scale_color_manual("Group", values = c("black", "red", "blue")) +
   theme_bw() +
   theme(
      text = element_text(size = 8,
                          family = "serif" # TT Times New Roman
                          #family = "sans" # TT Arial
      ),
      legend.position = "bottom",
      axis.ticks = element_line(size=0.4),
      axis.ticks.length = unit(0.05, "cm")
   ) -> AppendixFigure2_Survival_late_treatment
AppendixFigure2_Survival_late_treatment
