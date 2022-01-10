# IQWiGvsESMO
Further information (R-Code, ADEMP structur of simulations) of Paper "A comprehensive comparison of additional benefit assessment methods applied by IQWiG and ESMO for time-to-event endpoints after significant phase III trials – A simulation study"


## R-Code: 
- Instruction: ff


## Appenmdix

## Further information
- Further interpretation of Figure 1 and Figure 2:
We agree with the interpretation of Figure 1: ESMO’s and IQWiG’s method are most highly correlated when the effect size is “moderate” (or ~0.80). Figure 2 then shows that for low designHR values, the lack of correlation between the two methods is due to ESMO’s method having a higher proportion of maximal scores (almost all being maximal). This is also illustrated by the table below, showing a cross table between ESMO’s dual rule and IQWiG’s method of Standard Scenario 1 for different designHRs with 90% power, censoring rate of 20%, and medC equal to 6 months (first column of Figure 2). It can clearly be seen that ESMO’s dual rule almost exclusively awards maximal scores for large treatment effects. IQWiG’s method on the other hand assigns studies to the whole range of categories. Hence, the Spearman correlation between the methods are very small in these cases. 
| ESMO/IQWiG | Minor | Considerable | Major |
|------------|-------|--------------|-------|
| 1          | 0     | 0            | 0     |
| 2          | 1     | 0            | 0     |
| 3          | 1     | 0            | 0     |
| 4          | 315   | 992          | 7767  |


- 
