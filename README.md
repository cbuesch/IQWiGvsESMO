# IQWiGvsESMO
Further information (R-Code, ADEMP structure of simulations) of Paper "A comprehensive comparison of additional benefit assessment methods applied by IQWiG and ESMO for time-to-event endpoints after significant phase III trials – A simulation study".
In this github repository the paper and associated appendix including the ADEMP structure of our simulations can be found. Furthermore, the folder "FiguresAndTable" include all Figures and Tables of our paper as well as some plots illustrating further results of our simulations. And last but not least, the R-Code of our simulations are shared in the folder "RPrograms". Please read the R-Code instruction below before performing it!

## R-Code: 
Please keep in mind that the following R-packages should be installed: tidyverse, ggpubr, survival, flexsurv, foreach and doParallel. The provided programs are using the doParallel package for parallel computing in windows systems. If you have a unix-like system you need to install doMC instead of doParallel and replace 
```r
cl <- makeCluster(num.cl)
registerDoParallel(cl)
```
by
```r
registerDoMC(num.cl)
```
Furthermore, you need to remove all lines containing 
```r
stopCluster(cl)
```
In addition, the running time of the programs (especially the data generation and data analysis) is very long. We performed it with parallel computing using 48 cores and still needed several days. Therefore, please use many cores for the parallel computing or reduce the number of iterations (n<sub>sim</sub>).

In the following several R-Code-Scripts are explained, which can be found in the folder "RPrograms":
- Costume functions ("0_CostumeFunctions_Analysis.R" and "0_CostumeFunctions_DataGeneration.R"):
  These two scripts contain all needed functions for the data generation and analysis. Hence, they need to be loaded before conducting the simulations. 
  
- Data Generation ("1_DataGeneration.R"):
  This script conducts the data generation for all four scenarios. At the beginning the working directory needs to be set, where the file "Simulation.seed" and the script "0_CostumeFunctions_DataGeneration.R" are saved. Furthermore, as mentioned above please keep the number of cores (num.cl) / number of iterations (n_sim) in mind so that the running time is not too long. At last the folder path for saving of the generated data should be supplied. 
  
- Data Analysis ("2_DataAnalysis.R"):
  This script conducts the data analysis of the generated data for all four scenarios. At the beginning the working directory needs to be set, where the script "0_CostumeFunctions_DataGeneration.R" is saved. Furthermore, as mentioned above please keep the number of cores (num.cl) in mind so that the running time is not too long. At last the folder path for saving of the analysed data (path_ana) and the folder path for the generated data (path_data) should be supplied. 

- Visualizations ("3_VisualizationsOfResults.R", "3_VisualizationsOfResults_ROC.R"):
  The script "3_VisualizationsOfResults.R" creates the Figures 1, 2, 4 and 5 of the paper. In addition, it provides the Figures 1 and 2 of the Appendix. At the beginning you need to set your working directory to the folder where the final results of the script "2_DataAnalysis.R" is saved. 
  Figure 3 is created by script "3_VisualizationsOfResults_ROC.R" because further calculations of the True Positive Rate (TPR) and False Positive Rate (FPR) is needed and hence the generated data sets need to be loaded again. Therefore, at the beginning the working directory needs to be set (path_ROC), where the results of the script "1_DataGeneration.R" is saved. 

- Further illustrations and information ("4_FurtherIllustrationsAndInformation.R"):
    - Figure 2 with additional row: 
      Supplement Figure 2 with an additional row displaying the absolute difference of maximal scoring percentages between the methods.
    - Bias due to censoring mechanism: 
      Our censoring mechanism is partly depending on the event times (see Appendix) and consequently introduces bias to the HR estimation. To investigate the amount of introduced bias, we performed an additional simulation. Therefore, we simulated data of the Standard Scenario with n<sub>sim</sub>= 10,000 (iterations) and without sample size calculation, i.e. with a fixed sample size of 400 (200 per group) for each sub-scenario. The provided script performs this simulation and illustrates the results, which can also be found under the name "BiasDueToCensoringMechanism.tiff" in folder "FiguresAndTable".        
      As the results show, our censoring mechanism for specific censoring rates introduces bias for the HR estimation. Nevertheless, since we implemented a combination of administrative censoring (not dependent on the event time) and specific censoring rate (dependent on the event time), this bias is reduced, leading only to a slightly underestimated HR (overall bias of Standard Scenario 1: -0.01139, Monte Carlo SE of bias: 0.00000000160). Furthermore, this introduced bias is not affecting the method comparison to a substantial degree because it is affecting all compared methods equally. 
      In addition, below the table shows the bias and Monte Carlo SE of bias for different sub-scenarios of Standard Scenario 1 (without sample size calculation (N=400; 200 per group) and n<sub>sim</sub> = 10,000 (iterations), non-significant studies were still included in this analysis): 
      
      
      | censoring rate | med<sub>C</sub> | HR bias  | monte carlo SE of bias |
      |----------------|-----------------|----------|------------------------|
      | 0.2            | 6               | -0.00107 | 0.000129               |
      | 0.2            | 12              | 0.00110  | 0.000130               |
      | 0.2            | 18              | 0.00181  | 0.000131               |
      | 0.2            | 24              | 0.00198  | 0.000132               |
      | 0.2            | 30              | 0.00217  | 0.000133               |
      | 0.4            | 6               | -0.0158  | 0.000149               |
      | 0.4            | 12              | -0.0104  | 0.000147               |
      | 0.4            | 18              | -0.00748 | 0.000146               |
      | 0.4            | 24              | -0.00578 | 0.000146               |
      | 0.4            | 30              | -0.00467 | 0.000146               |
      | 0.6            | 6               | -0.0350  | 0.000186               |
      | 0.6            | 12              | -0.0285  | 0.000183               |
      | 0.6            | 18              | -0.0250  | 0.000181               |
      | 0.6            | 24              | -0.0229  | 0.000181               |
      | 0.6            | 30              | -0.0214  | 0.000180               |
      

## Further information
- Further Figures:
  In folder "FiguresAndTable" other versions of the figures in the paper can be found illustrating additional results for 80% power, 40% censoring and HR<sub>var</sub> of 0.8 (factor for deviance between designHR and trueHR):
    - Figure1_with80PowerAnd40Censoring.tiff
    - Figure2_with80PowerAnd40Censoring.tiff
    - Figure4_withHRvar08.tiff
    - Figure5_with80PowerAnd40Censoring.tiff
- In addition, in folder "FiguresAndTable" further figures giving a greater understanding of the ESMO method are provided:
    - AddFig_ESMO_Max_DueToSurvGain.tiff
    - AddFig_StdErrKaplanMeier.tiff
  
- Further interpretation of Figure 1 and Figure 2:
  In Figure 1 it can be seen that ESMO’s and IQWiG’s method are most highly correlated when the effect size is “moderate” (or ~0.80). Figure 2 then shows that for low designHR values, the lack of correlation between the two methods is due to ESMO’s method having a higher proportion of maximal scores (almost all being maximal). This is also illustrated by the table below, showing a cross table between ESMO’s dual rule and IQWiG’s method of Standard Scenario 1 for a designHRs of 0.36 with 90% power, censoring rate of 20%, and med<sub>C</sub> equal to 6 months (first column of Figure 2). It can clearly be seen that ESMO’s dual rule almost exclusively awards maximal scores for large treatment effects. IQWiG’s method on the other hand assigns studies to the whole range of categories. Hence, in this case the Spearman correlation between the methods are very small (0.0407). 


    | ESMO/IQWiG | Minor | Considerable | Major |
    |------------|-------|--------------|-------|
    | 1          | 0     | 0            | 0     |
    | 2          | 1     | 0            | 0     |
    | 3          | 1     | 0            | 0     |
    | 4          | 315   | 992          | 7767  |

  In cases of “moderate” designHRs, the Spearman correlation is the largest between the two methods due to the fact that both methods assign scores over the complete range of categories (e.g. see table below with designHRs = 0.78, 90% power, 20% censoring rate and med<sub>C</sub> = 6 months and hence a Spearman correlation of 0.7500). For the maximal category, this can also be seen in Figure 2, showing similar proportions in the area of designHRs where the Spearman correlation is the largest.

    | ESMO/IQWiG | Minor | Considerable | Major |
    |------------|-------|--------------|-------|
    | 1          | 1558  | 1457         | 31    |
    | 2          | 1     | 2892         | 217   |
    | 3          | 0     | 794          | 1143  |
    | 4          | 1     | 184          | 871   |

  In cases of large designHRs (small treatment effect), IQWiG’s method appears to have a higher proportion of maximal scores (see Figure 2). This again can be seen in the table below (Standard Scenario 1 for designHR of 0.84 with 90% power, censoring rate of 20% and med<sub>C</sub> equal to 6 months) showing that ESMO almost never awards the maximal score anymore (0.0549%). IQWiG’s method appears to assign a higher proportion of maximal scores, but overall the probability for a maximal score is small (5.3920%). 
Again, IQWiG’s method awards the whole range of categories and almost always awards the lowest category. Hence, the Spearman correlation between the methods are again smaller (0.4105 for designHR of 0.84 and 0.2118/ for designHR of 0.86). 

    | ESMO/IQWiG | Minor | Considerable | Major |
    |------------|-------|--------------|-------|
    | 1          | 2598  | 5652         | 41    |
    | 2          | 0     | 364          | 415   |
    | 3          | 0     | 0            | 31    |
    | 4          | 0     | 1            | 4     |

  Overall, when the designHR is very small / large (effect size is strong / weak) and hence the classification seems to be the simplest, the methods disagree the most. As shown above, the reason for this behaviour is due to the fact that IQWiG’s method always awards the whole range of possible categories. Note that this does not mean that IQWiG always awards the same proportion for every category (minor, considerable, major) over each scenario. It just means that IQWiG, over the n<sub>sim</sub> iterations, at least sometimes awards every possible category. And since ESMO’s dual rule, in extreme cases, awards only one category to almost every iteration, the Spearman correlation is small. Due to space limitations, we did not include these cross tables and explanations into our manuscript but added it here.
  
  
- Interpretation of Figure AddFig_ESMO_Max_DueToSurvGain.tiff:
We implemented the complete ESMO method; hence the ESMO method can achieve a maximal score if the survival rate increases by 10% or more at 2, 3 or 5 years (depending on the occurred median survival time in the control group; see table 1). This figure illustrates the proportion of maximal ESMO scores due to the survival gain rule, dual rule or both rules (number of maximal ESMO scores due to survival gain rule (OR dual rule OR Both) divided by number of maximal ESMO scores) of our simulations for Standard Scenario 1, showing that most of the times a maximal ESMO score was assigned by its dual rule as well as its survival gain rule. Only in scenarios with small treatment effects (designHR close to 1) the reason that ESMO assigns a maximal score is solely due to its survival gain rule (green line at 100%). The number of assigned maximal scores in these sub-scenarios, however, is small. For example, in the sub-scenario with med<sub>C</sub> = 12, censoring proportion of 20% and designHR=0.86 (top row, second panel from the left), ESMO’s method assigns 73 trials a maximal score; in sub-scenario with med<sub>C</sub> = 12, censoring proportion of 40% and designHR=0.88 (middle row, second panel from the left) ESMOs method assigns 3 trials a maximal score. Some values/points are missing from the figure, because in these scenarios no maximal scores were assigned by the ESMO dual rule. Overall, it is safe to say that in our simulation all parts of the ESMO methods contribute to the ESMO assessment.

- Interpretation of Figure AddFig_StdErrKaplanMeier.tiff:
To examine if the estimated 2/3/5-year survival rates for ESMO's method are reliable, we took a look at the standard error of the Kaplan-Meier estimates at 2/3/5 years of our Standard Scenario 1 (see figure AddFig_StdErrKaplanMeier.tiff). Due to our sample size calculation resulting in larger sample sizes with smaller treatment effects, the standard error is getting smaller with smaller treatment effects as well. Furthermore, in the scenario with med<sub>C</sub> = 6 and large treatment effects, a clear group difference for the standard error can be seen, which can be interpreted as unreliable survival rate estimation. Nevertheless, as the figure AddFig_ESMO_Max_DueToSurvGain.tiff illustrates, in these extreme scenarios with such a large treatment effect, both parts of ESMOs method do assign a maximal score and, hence, the survival gain part is not even needed. Reason for this group difference is the fact that the survival curve of the control group is more often decreasing to 0% before the 2-year mark, which leads to a standard error of 0% and thus reduces the mean. This difference, however, is getting smaller with larger underlying med<sub>C</sub> and censoring rates because less survival curves are dropping to 0% before the key milestone 2, 3 or 5 year. 
Since we implemented the complete ESMO method and, in our opinion, the 2/3/5-year survival rates are reliable. 

