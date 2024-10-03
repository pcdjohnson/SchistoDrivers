### Summary ----
# Script to estimate power to (1) identify multiple drivers of      
# praziquantel treatment failure in individuals infected with
# schistosomiasis, and (2) detect a difference in praziquantel 
# absorption between different instructions for taking food to
# aid absorption

# Load packages
library(ggplot2)
library(parallel)
library(MASS)
library(glmmTMB)

# Clear objects
rm(list = ls())

#### Power analysis 1 ----

# Study design options

# Total sample size
n <- seq(500, 3000, 500)

# Model parameters

# failure to clear proportion
p <- seq(0.05, 0.3, 0.05) 

# odds ratio per binary predictor, or per SD for continuous predictors
or <- seq(1, 2, 0.25) 

# how correlated are predictors (because abs(r) > 0 reduces power,
# and it would be unrealistic to assume zero correlation)
r <- c(0.25) 

# Choose candidate drivers:

# How many drivers?
# Continuous drivers x 6 (generic, but could be PK, genetic SNP score, etc)
# 4 binary drivers
# Total number of drivers:
n.x <- 10

# For the binary drivers, what are their proportions?
bin.x.p <- c(0.2, 0.5, 0.2, 0.2) # HIV 20%, malaria 50%, STH 20%, hybrid/resistance presence 20%

# Further simulation and analysis options

# number of data sets to simulate
n.sim <- 10

# adjust the significance thresholds for multiple testing (Bonferroni)
nominal.alpha <- 0.05
alpha <- c(nominal.alpha/n.x)

# Make table of all parameter and design combinations
par.tab <- expand.grid(n = n, p = p, or = or, r = r, alpha = alpha, n.sim = n.sim)

# Start the clock
start.time <- Sys.time()

# Loop over all scenarios = rows of par.tab
sim.res <- 
  sapply(1:nrow(par.tab), function(i) {
    
    # Print progress
    print(paste0(round(100*(i-1)/nrow(par.tab)), "% complete"))
    
    # Simulate Xs (the drivers)
    
    # Correlation among Xs
    sigma2 <- diag(n.x)
    sigma2[lower.tri(sigma2)] <- par.tab$r[i]
    sigma2[upper.tri(sigma2)] <- par.tab$r[i]
    
    # Simulate data and analysis n.sim times, returning the mean
    # number of drivers identified at P < alpha
    n.sim.out <-
      mclapply(1:n.sim, function(k) {
        
        # Simulate Xs before dichotomising (so the binary Xs are derived from latent scales)
        X.raw <- mvrnorm(par.tab$n[i], mu = rep(0, n.x), Sigma = sigma2)
        
        # Dichotomise the binary Xs, but centre on zero by subtracting 0.5
        # and add a column of 1s for the intercept
        X <-
          cbind(1, 
                sapply(1:n.x, function(j) {
                  if(j <= length(bin.x.p)) {
                    X.raw[, j] <- as.integer(X.raw[, j] < qnorm(bin.x.p[j])) - 0.5
                  }
                  X.raw[, j]
                }))
        
        # vector of intercept (log odds) and log odds ratios
        b <- c(qlogis(par.tab$p[i]), log(rep(par.tab$or[i], n.x)))
        
        # Simulate failure to clear from linear predictor X %*% b
        response <- rbinom(par.tab$n[i], 1, plogis(X %*% b))
        
        # Put everything into a data frame
        dat <- data.frame(X[, -1], response = response)
        
        # fit GLM
        fit <- glm(response ~ X[, -1], family = binomial)
        res.tab <- coef(summary(fit))
        
        # How many drivers were detected (P < alpha)?
        sum(res.tab[-1, "Pr(>|z|)"] < par.tab$alpha[i])
        
      }, mc.cores = detectCores())
    
    # take mean number of drivers (Xs) detected across simulated data sets
    mean(unlist(n.sim.out))
    
  })

# Stop the clock
run.time <- Sys.time() - start.time
print(run.time)

# Estimate power as the proportion of the n.x drivers with p < alpha
par.tab$prop.drivers <- sim.res/n.x

# Export results to CSV file with time stamp in file name
file.name <- paste0("results/schisto_power", 
                    substr(gsub(":", "", (gsub(" ", "-", Sys.time()))), 1, 15), ".csv")
write.csv(par.tab, file.name, row.names = FALSE, quote = FALSE)

# Make plot of results
par.tab$alpha <- factor(paste("alpha =", round(par.tab$alpha, 4)))
par.tab$n <- factor(par.tab$n, rev(n))
par.tab$r <- factor(paste("Correlation among drivers (r) =", par.tab$r))
par.tab$p <- factor(paste("P(failure to clear) =", par.tab$p))

power.plot <-
  ggplot(data = par.tab, aes(x = or, y = prop.drivers, color = n, shape = n, group = n)) + 
  geom_hline(yintercept = 0.8, linewidth = 0.3, linetype = 2) +
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  facet_wrap(~ p, ncol = 2) +
  xlab("Odds ratio") +
  ylab("Power") +
  labs(caption = file.name) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
power.plot
plot.file.name <- "schisto_power1.png"
ggsave(plot.file.name, width = 6, height = 6)



# output methods to README.md
readme.file <- "README.md"

cat("# SchistoDrivers\n\n",
    "## Power analysis 1: identifying drivers of schistosomiasis praziquantel treatment failure\n\n",
    "### Methods\n\n",
    "The aim of power analysis 1 is to estimate power to detect drivers of praziquantel treatment",
    "failure in individuals infected with schistosomiasis. This is a simulation-based power analysis,",
    "where the study data is simulated and analysed multiple times under varying study design",
    "scenarios, and power is estimated as the proportion of simulated analyses that achieve the",
    "desired outcome (detecting a true driver of treatment failure). The association between the",
    "outcome (treatment failure) and each driver is estimated and tested in a multivariable GLM.",
    "For this analysis, power is defined as the proportion of drivers that are significantly associated",
    paste("with the outcome, averaged across", n.sim, "simulated data analyses per scenario.\n\n"),
    "The following assumptions are made:\n-",
    paste(n.x, "drivers are associated with the outcome, of which", length(bin.x.p), "are binary and", 
          n.x - length(bin.x.p), "continuous."),
    paste0("The prevalences of the binary drivers are: ", paste(sort(bin.x.p), collapse = ", "), ","),
    "representing drivers such as co-infection (soil-transmitted helminths, malaria) and",
    "hybrid/resistance presence).\n",
    paste0("- The drivers are correlated with each other, with a common correlation coefficient of ", r, "."),
    "We don’t know what the true correlation is among drivers, but moderate correlations are likely",
    "and neglecting them will give optimistic power estimates.\n",
    "- In order to control inflation of the number of false positive results due to multiple testing of",
    n.x, "drivers, the significance threshold of", nominal.alpha,
    paste0("is Bonferroni-adjusted to ", alpha, ", i.e. a driver is significant if P < ", alpha, ".\n\n"),
    "We explore the effect on power of varying the following study design choices/assumptions:\n",
    paste0("- Sample size (number of infected and treated individuals): ", paste(n, collapse = ", "), ".\n"),  
    paste0("- The proportion of these that fail to clear: ", paste(p, collapse = ", "), ".\n"),
    "- The strength of association between each driver and failure to clear,",
    "defined as an odds ratio for binary drivers and as an odds ratio per standard deviation",
    paste0("unit for continuous drivers: ", paste(or[or != 1], collapse = ", "), "."),
    "To give a sense of what these effect sizes mean:\n",    
    "  - If a binary driver has an odds ratio of 1.5, then if 5% of people fail to clear", 
    "in the absence of a driver, that proportion will be 7.3% among those who are exposed to the driver.",
    "If the prevalence of failure to clear is 25% without the driver, it will be 33% with the driver.\n",
    "  - If we raise the odds ratio from 1.5 to 2, then the prevalences of failure to clear in the",
    "exposed population will be 9.5% relative to 5% in the unexposed population,",
    "and 40% relative to 25% in the unexposed population.\n",
    "Full details are provided in the script PowerAnalysis.R.",
    "Results are output as CSV to the results",
    paste0("directory and plotted to ", plot.file.name, ".\n\n"),
    "### Results\n",
    paste0("![PowerCurve](", plot.file.name, ")"),
    file = readme.file)

#### Power analysis 2 ----

# Trial to compare absorption of praziquantel between people who have
#   - received food at home
#   - brought in food
#   - food has been provided on site

#  Assumptions:
#    No cluster effect (conditional independence between observations)
#    Equal numbers in each group
#    AUC is normally distributed within each group (it is)
#    We know the SD of the AUC within each group (it's 403 units)
#    We know what our target effect size is, i.e. the smallest effect that we want to detect
#    alpha = 0.05 and target power = 90%.  

delta.AUC <- 0.25
target.power <- 0.9
n <- 
  power.anova.test(groups = 3, power = target.power, 
                   between.var = var(c(0, delta.AUC, delta.AUC * 2)),
                   within.var = 1, sig.level = nominal.alpha)$n

cat("\n\n## Power analysis 2: Randomised control trial protocol for the effect of",
    "food prior to treatment on praziquantel absorption\n\n",
    "### Methods\n\n",
    "The aim of power analysis 2 is to estimate power to detect a difference in praziquantel absorption",
    "(measured as area under the curve [AUC] of praziquantel metabolites) between three groups:\n",
    "- people who have taken food at home prior to praziquantel treatment;\n",
    "- people who have brought in food to be taken prior to praziquantel treatment;\n",
    "- people for whom food has been provided on site prior to praziquantel treatment.\n\n",
    "The null hypothesis is that mean AUC is equal across the three groups.",
    "The alternative hypothesis that mean AUC differs between the three groups.",
    "The effect size assumed here is that the group AUC means differ by",
    delta.AUC, "standard deviations from the lowest mean to the intermediate mean,",
    "and by", delta.AUC, "standard deviations from the intermediate mean to the highest mean.",
    "Target power is", paste0(target.power * 100, "%"), "and the significance threshold is",
    paste0(nominal.alpha, "."),
    "The required sample size per group was calculated using the R function *power.anova.test*.",
    "\n\n",
    "### Results\n",
    n, "people would be required per group in order to achieve",
    paste0(target.power * 100, "%"), 
    "to detect a", delta.AUC, "standard deviation difference between each of the three groups.",
    "\n\n",
    file = readme.file, append = TRUE)



if(FALSE) {

  
    
  title.text <-
    paste0("Sample size required per group for ", target.power * 100, 
           "% power at alpha = 0.05\ngiven target delta.AUC in SD units")
  title(title.text)
  cbind(delta.AUC, N.per.group)
  
}
