# SchistoDrivers

 ### Power analysis for identifying drivers of schistosomiasis praziquantel treatment failure

 The script PowerAnalysis.R estimates power across a range of model parameter assumptions and sample sizes, detailed in comments are in the R script. Results are output as CSV to the results directory and plotted to schisto_power.pdf.

 The aim of the power analysis is to estimate power to detect drivers of praziquantel treatment failure in individuals infected with schistosomiasis. This is a simulation-based power analysis, where the study data is simulated and analysed multiple times under varying study design scenarios, and power is estimated as the proportion of simulated analyses that achieve the desired outcome (detecting a true driver of treatment failure). The association between the outcome (treatment failure) and each driver is estimated and tested in a multivariable GLM. For this analysis, power is defined as the proportion of drivers that are significantly associated with the outcome, averaged across 10 simulated data analyses per scenario.

 The following assumptions are made:
- 10 drivers are associated with the outcome, of which 4 are binary and 6 continuous. The prevalences of the binary drivers are: 0.2, 0.2, 0.2, 0.5, representing drivers such as co-infection (HIV, soil-transmitted helminths, malaria) and hybrid/resistance presence).
 - The drivers are correlated with each other, with a common correlation coefficient of