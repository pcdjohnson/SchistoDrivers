###########################################################################
# Script to estimate power to identify multiple drivers of praziquantel   #
# treatment failure in individuals infected with schistosomiasis.         #
#                                                                         #
# This is a simulation-based power analysis, where the project data is    #
# simulated and analysed multiple times under different study design      #
# scenarios, and power is estimated as the proportion of simulated        #
# analyses that achieve the desired outcome (e.g. detecting a true        # 
# driver of treatment failure).                                           #
###########################################################################

# Load packages
library(ggplot2)
library(parallel)
library(MASS)
library(glmmTMB)

# Clear objects
rm(list = ls())

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
n.sim <- 1000

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
plot.file.name <- "schisto_power.png"
ggsave(plot.file.name, width = 6, height = 6)


### Appendix: Additional power analyses ### 

if(FALSE) {
  
  # Additional power analysis 
  # Trial to compare absorption of praziquantel between people who have received food beforehand and those who haven’t
  
  #  Assumptions:
  #    No cluster effect (conditional independence between observations)
  #    Equal numbers in each group (this can be relaxed but usually isn’t)
  #    AUC is normally distributed within each group
  #    We know the SD of the AUC (or some transformation – lab measures are often positively skewed) within each group
  #    We know what our target effect size is, i.e. the smallest effect that we want to detect = the largest effect that we’d be relaxed about not detecting.
  #    alpha = 0.05 and target power = 80%.  
  
  delta.AUC <- seq(0.05, 1, 0.05)
  N.per.group <- ceiling(sapply(delta.AUC, function(delta) power.t.test(delta = delta, sd = 1, sig.level = 0.05, power = 0.8)$n))
  plot(delta.AUC, N.per.group, log = "y", type = "b")
  title("Sample size required per group for 80% power at alpha = 0.05\ngiven target delta.AUC in SD units")
  cbind(delta.AUC, N.per.group)
  
  # Additional power analysis for comparison of eggs per gram between two groups:
  #   - food then praziquental
  #   - no food then praziquental - predicting 3x higher EPG
  unfed.effect <- 3
  # Based on existing data mean and SD of EPG in fed subjects following treatment
  scaling.factor <- 24 # because EPG is calculated from >1 gram of faeces
  lambda <- 0.2128099 * scaling.factor
  s <- 2.45669 * scaling.factor
  theta <- lambda^2 / (s^2 - lambda) # calculate dispersion parameter from SD and mean
  # check we get the right mean and SD
  mean(rnbinom(100000, size = theta, mu = lambda))
  sd(rnbinom(100000, size = theta, mu = lambda))
  
  # estimate power
  power.estimate <- 
    unlist(mclapply(N.per.group, function(N) {
      # N <- 17
      p.value <-
        replicate(n.sim, {
          dat <- expand.grid(id = 1:N, group = c("fed", "unfed"))
          dat$lambda <- lambda * ifelse(dat$group == "fed", 1, unfed.effect)
          dat$epg <- rnbinom(nrow(dat), size = theta, mu = dat$lambda)
          (tab <- table(dat$epg, dat$group))
          # if the data is unanalysable (e.g. all zeros) then p.val=1 (effectively)
          if(any(tapply(dat$epg, dat$group, function(x) length(unique(x))) == 1)) return(1) 
          
          #fit <- glm.nb(epg ~ group, data = dat)
          #coef(summary(fit))["groupunfed", "Pr(>|z|)"]
          # ...switch from glm.nb to glmmTMB as it converges better at low N
          fit <- glmmTMB(epg ~ group, data = dat, family = nbinom2)
          p.val <- coef(summary(fit))$cond["groupunfed", "Pr(>|z|)"]
          if(is.nan(p.val)) return(1) else return(p.val)
        })
      mean(p.value < 0.05)
    }, mc.cores = detectCores()))
  
  plot(N.per.group, power.estimate, type = "b", log = "x", ylim = 0:1)
  abline(h = 0.8, lty = 3)
  title(paste0("Power to detect an unfed:fed ratio in EPG of ", 
               unfed.effect, ":1\nassuming mean (SD) fed EPG = ", 
               round(lambda, 2), " (", round(s, 2), ")\nestimated from ", n.sim, " simulated data sets"))
  cbind(N.per.group, power.estimate)
  
}


# output methods to README.md
readme.file <- "README.md"

cat("# SchistoDrivers\n\n",
    "## Power analysis for identifying drivers of schistosomiasis praziquantel treatment failure\n\n",
    "### Summary\n\n",
    "The script PowerAnalysis.R estimates power across a range of model parameter assumptions and",
    "sample sizes, detailed in comments in the R script. Results are output as CSV to the results",
    paste0("directory and plotted to ", plot.file.name, ".\n\n"),
    "### Methods\n\n",
    "The aim of the power analysis is to estimate power to detect drivers of praziquantel treatment",
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
    "representing drivers such as co-infection (HIV, soil-transmitted helminths, malaria) and",
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
    "and 40% relative to 25% in the unexposed population.\n\n",
    "### Results\n",
    paste0("![PowerCurve](", plot.file.name, ")"),
    file = readme.file)

