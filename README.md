## Sample size for Aim 2a: identifying drivers of schistosomiasis praziquantel treatment failure

 ### Methods

 The aim of this power analysis is to estimate power to detect drivers of praziquantel treatment failure in individuals infected with schistosomiasis. This is a simulation-based power analysis, where the study data is simulated and analysed multiple times under varying study design scenarios, and power is estimated as the proportion of simulated analyses that achieve the desired outcome (detecting a true driver of treatment failure). The association between the outcome (treatment failure) and each driver is estimated and tested in a multivariable GLM. For this analysis, power is defined as the proportion of drivers that are significantly associated with the outcome, averaged across 1000 simulated data analyses per scenario.

 The following assumptions are made:
- 10 drivers are associated with the outcome, of which 3 are binary and 7 continuous. The prevalences of the binary drivers are: 0.5, 0.2, 0.2, representing  malaria, soil-transmitted helminths, hybrid/resistance presence .
 - The drivers are correlated with each other, with a common correlation coefficient of 0.25. We don’t know what the true correlation is among drivers, but moderate correlations are likely and neglecting them will give optimistic power estimates.
 - In order to control inflation of the number of false positive results due to multiple testing of 10 drivers, the significance threshold of 0.05 is Bonferroni-adjusted to 0.005, i.e. a driver is significant if P < 0.005.

 We explore the effect on power of varying the following study design choices/assumptions:
 - Sample size (number of infected and treated individuals): 600, 1200, 1800, 2400, 3600, 4800, 6000, 7200, 8400, 9600.
 - The proportion of these that fail to clear: 0.1, 0.2, 0.3.
 - The strength of association between each driver and failure to clear, defined as an odds ratio for binary drivers and as an odds ratio per standard deviation unit for continuous drivers: 1.25, 1.5, 1.75, 2. To give a sense of what these effect sizes mean:
   - If a binary driver has an odds ratio of 1.5, then if 5% of people fail to clear in the absence of a driver, that proportion will be 7.3% among those who are exposed to the driver. If the prevalence of failure to clear is 25% without the driver, it will be 33% with the driver.
   - If we raise the odds ratio from 1.5 to 2, then the prevalences of failure to clear in the exposed population will be 9.5% relative to 5% in the unexposed population, and 40% relative to 25% in the unexposed population.

 Full details are provided in the script [PowerAnalysis.R](https://github.com/pcdjohnson/SchistoDrivers/blob/main/PowerAnalysis.R). Results are output as CSV to the [results](https://github.com/pcdjohnson/SchistoDrivers/tree/main/results) directory and plotted to [schisto_power2a.png](https://github.com/pcdjohnson/SchistoDrivers/blob/main/schisto_power2a.png).

 ### Results
 ![Power2aCurve](schisto_power2a.png) 


## Sample size for Aim 2b-i: identifying individual-level drivers of schistosomiasis infection

 ### Methods

 The aim of this power analysis is to estimate power to detect drivers of schistosomiasis infection. The association between the outcome (infection) and each driver is estimated and tested in a multivariable GLM. Power is defined as the proportion of drivers that are significantly associated with the outcome, averaged across 1000 simulated data analyses per scenario.

 The following assumptions are made:
- 10 drivers are associated with the outcome, of which 3 are binary and 7 continuous. The prevalences of the binary drivers are: 0.5, 0.2, 0.2, representing  malaria, soil-transmitted helminths, hybrid/resistance presence .
 - The drivers are correlated with each other, with a common correlation coefficient of 0.25.
 - In order to control inflation of the number of false positive results due to multiple testing of 10 drivers, the significance threshold of 0.05 is Bonferroni-adjusted to 0.005, i.e. a driver is significant if P < 0.005.

 We explore the effect on power of varying the following study design choices/assumptions:
 - Sample size: it is assumed that approximately 1600 infected individuals will be recruited (mean realised number of positives = 1460), while the number of negatives will be varied: 200, 375, 550, 725, 900, 1075, 1250, 1425, 1600.
 - The strength of association between each driver and failure to clear, defined as an odds ratio for binary drivers and as an odds ratio per standard deviation unit for continuous drivers: 1.25, 1.5, 1.75, 2. Full details are provided in the script [PowerAnalysis.R](https://github.com/pcdjohnson/SchistoDrivers/blob/main/PowerAnalysis.R). Results are output as CSV to the [results](https://github.com/pcdjohnson/SchistoDrivers/tree/main/results) directory and plotted to [schisto_power2bi.png](https://github.com/pcdjohnson/SchistoDrivers/blob/main/schisto_power2bi.png).

 ### Results
 ![Power2biCurve](schisto_power2bi.png) 


## Sample size calculation for Aim 2b-ii: identifying individual- and community-level drivers of re-infection following clearance

 ### Methods

 The aim of this power analysis is to estimate power to detect individual- and community-level drivers of schistosomiasis re-infection following clearance, and the expected margin of error around driver odds ratio estimates. The association between the outcome (infection) and each driver is estimated and tested in a multivariable GLMM. For this analysis, power is defined as the proportion of drivers that are significantly associated with the outcome, averaged across 500 simulated data analyses per scenario. Power and margin of error (half the width of a 95% confidence interval) in odds ratio estimation are presented across a range of intra-class correlation coefficient (ICC) values, where ICC = 0% represents drivers that vary between individuals within communities but not between communities, and ICC = 100% represents drivers that have the same value for all community members and differ between communities.

 The following assumptions are made:
- 10 drivers are associated with the outcome, of which 1 are binary and 9 continuous. The prevalences of the binary drivers are: 0.2, representing  hybrid/resistance presence .
 - The drivers are correlated with each other, with a common correlation coefficient of 0.25.
 - Log odds of re-infection varies among communities with a variance of 2.73.
 - In order to control inflation of the number of false positive results due to multiple testing of 10 drivers, the significance threshold of 0.05 was Bonferroni-adjusted to 0.005.

 We explore the effect on power of varying the following study design choices/assumptions:
 - Total sample size: 1600.
 - Community sample size (number of communities sampled): 22, 24, 26.
 - ICC: 0%, 25%, 50%, 75%, 100%.
 - Prevalence of re-infection: 0.5.
 - The strength of association between each driver and re-infection, defined as an odds ratio per standard deviation unit for continuous community-level drivers: 1.25, 1.5, 1.75, 2.

 Full details are provided in the script [PowerAnalysis.R](https://github.com/pcdjohnson/SchistoDrivers/blob/main/PowerAnalysis.R). Results are output as CSV to the [results](https://github.com/pcdjohnson/SchistoDrivers/tree/main/results) directory and plotted to [schisto_power2bii.png](https://github.com/pcdjohnson/SchistoDrivers/blob/main/schisto_power2bii.png) and [schisto_moe2bii.png](https://github.com/pcdjohnson/SchistoDrivers/blob/main/schisto_moe2bii.png).

 ### Results
 ![Power2biiCurve](schisto_power2bii.png) 


 ![MoE2biiCurve](schisto_moe2bii.png) 


## Sample size calculation for Aim 2c: identifying community-level drivers of schistosomiasis infection

 ### Methods

 The aim of this power analysis is to estimate power to detect community-level drivers of schistosomiasis infection, and the expected margin of error around community-level driver odds ratio estimates. The association between the outcome (infection) and each driver is estimated and tested in a multivariable GLMM. For this analysis, power is defined as the proportion of drivers that are significantly associated with the outcome, averaged across 500 simulated data analyses per scenario.

 The following assumptions are made:
- 4 continuous drivers are associated with the outcome.
 - The drivers are correlated with each other, with a common correlation coefficient of 0.25.
 - Log odds of infection prevalence varies among communities with a variance of 2.73.
 - In order to control inflation of the number of false positive results due to multiple testing of 4 drivers, the significance threshold of 0.05 is Bonferroni-adjusted to 0.0125, i.e. a driver is significant if P < 0.0125.

 We explore the effect on power of varying the following study design choices/assumptions:
 - Total sample size: 600, 1200, 1800, 2400, 3600, 4800, 6000, 7200, 8400, 9600.
 - Community sample size (number communities sampled): 22, 24, 26.
 - Prevalence of infection: 0.1, 0.5.
 - The strength of association between each driver and failure to clear, defined as an odds ratio per standard deviation unit for continuous community-level drivers: 1.25, 1.5, 1.75, 2.

 Full details are provided in the script [PowerAnalysis.R](https://github.com/pcdjohnson/SchistoDrivers/blob/main/PowerAnalysis.R). Results are output as CSV to the [results](https://github.com/pcdjohnson/SchistoDrivers/tree/main/results) directory and plotted to [schisto_power2c.png](https://github.com/pcdjohnson/SchistoDrivers/blob/main/schisto_power2c.png) and [schisto_moe2c.png](https://github.com/pcdjohnson/SchistoDrivers/blob/main/schisto_moe2c.png).

 ### Results
 ![Power2cCurve](schisto_power2c.png) 


 ![MoE2cCurve](schisto_moe2c.png) 


