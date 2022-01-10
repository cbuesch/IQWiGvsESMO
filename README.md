# IQWiGvsESMO
Further information (R-Code, ADEMP structur of simulations) of Paper "A comprehensive comparison of additional benefit assessment methods applied by IQWiG and ESMO for time-to-event endpoints after significant phase III trials – A simulation study".
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
In addition, the running time of the programs (escpecially the data generation and data analysis) is very long. We performed it with parallel computing using 48 cores and still needed several days. Therefore, please use many cores for the parallel computing or reduce the number of iterations (n<sub>sim</sub>).

In the following several R-Code-Scripts are explained, which can be found in the folder "RPrograms":
- Costume functions ("0_CostumeFunctions_Analysis.R" and "0_CostumeFunctions_DataGeneration.R"):
  These two scripts conatin all needed functions for the data generation and analysis. Hence, they need to be loaded before conducting the simulations. 
- Data Generation ("1_DataGeneration.R"):
  
- Data Analysis ("2_DataAnalysis.R"):
- Visualizations ("3_VisualizationsOfResults.R", "3_VisualizationsOfResults_ROC.R", "4_FurtherIllustrationsAndInformation.R")

## Further information
- Figure 2 with additional row: TODO

- Bias due to cesnoring mechanism: TODO



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
