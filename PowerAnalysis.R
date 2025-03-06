### Overview ----
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
library(GLMMmisc) # from https://github.com/pcdjohnson/GLMMmisc

# Clear objects
rm(list = ls())

# Global settings
readme.file <- "README.md" # methods and results output file
nominal.alpha <- 0.05 # significance threshold
n.sim <- 50 # number of data sets to simulate (divided by 2 for the GLMM analysis, because it's slow)

#### Sample size for Aim 2 (individual clearance) ----

# Study design options

# Total sample size
n <- seq(1000, 3000, 500)

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
# Continuous drivers (generic, but could be PK, genetic SNP score, age, immune factors, AUC, etc)
n.cont <- 7

# 4 binary drivers
# For the binary drivers, what are their proportions?
# HIV 20%, malaria 50%, STH 20%, hybrid/resistance presence 20%
bin.x.p <- 
  c(HIV = 0.2, 
    malaria = 0.5, 
    `soil-transmitted helminths` = 0.2, 
    `hybrid/resistance presence` = 0.2)[c(2, 3, 4)] 

# Total number of drivers:
n.x <- n.cont + length(bin.x.p)

# Further simulation and analysis options

# adjust the significance thresholds for multiple testing (Bonferroni)
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
      mclapply(1:par.tab$n.sim[i], function(k) {
        
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
file.name.power2 <- 
  paste0("results/schisto_power2_", 
         substr(gsub(":", "", (gsub(" ", "-", Sys.time()))), 1, 15), ".csv")
write.csv(par.tab, file.name.power2, row.names = FALSE, quote = FALSE)
rm(par.tab)
par.tab <- read.csv(file.name.power2) 

# Make factors for plotting
par.tab$alpha <- factor(paste("alpha =", round(par.tab$alpha, 4)))
par.tab$n <- factor(par.tab$n, rev(n))
par.tab$r <- factor(paste("Correlation among drivers (r) =", par.tab$r))
par.tab$p <- factor(paste("P(failure to clear) =", par.tab$p))

# Make plot of results
power.plot <-
  ggplot(data = par.tab, aes(x = or, y = prop.drivers, color = n, shape = n, group = n)) + 
  geom_hline(yintercept = 0.8, linewidth = 0.3, linetype = 2) +
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  facet_wrap(~ p, ncol = 2) +
  xlab("Odds ratio") +
  ylab("Power") +
  labs(caption = file.name.power2) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
power.plot
plot.file.name <- "schisto_power2.png"
ggsave(plot.file.name, width = 6, height = 6)

# output methods and results to README.md

cat("## Sample size for Aim 2: identifying drivers of schistosomiasis praziquantel treatment failure\n\n",
    "### Methods\n\n",
    "The aim of this power analysis is to estimate power to detect drivers of praziquantel treatment",
    "failure in individuals infected with schistosomiasis. This is a simulation-based power analysis,",
    "where the study data is simulated and analysed multiple times under varying study design",
    "scenarios, and power is estimated as the proportion of simulated analyses that achieve the",
    "desired outcome (detecting a true driver of treatment failure). The association between the",
    "outcome (treatment failure) and each driver is estimated and tested in a multivariable GLM.",
    "For this analysis, power is defined as the proportion of drivers that are significantly associated",
    paste("with the outcome, averaged across", unique(par.tab$n.sim), "simulated data analyses per scenario.\n\n"),
    "The following assumptions are made:\n-",
    paste(n.x, "drivers are associated with the outcome, of which", length(bin.x.p), "are binary and", 
          n.cont, "continuous."),
    paste0("The prevalences of the binary drivers are: ", paste(bin.x.p, collapse = ", "), ","),
    "representing ", paste(names(bin.x.p), collapse = ", "), ".\n",
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
    "Full details are provided in the script",
    "[PowerAnalysis.R](https://github.com/pcdjohnson/SchistoDrivers/blob/main/PowerAnalysis.R).",
    "Results are output as CSV to the [results](https://github.com/pcdjohnson/SchistoDrivers/tree/main/results)",
    paste0("directory and plotted to [", plot.file.name, 
           "](https://github.com/pcdjohnson/SchistoDrivers/blob/main/", plot.file.name, ").\n\n"),
    "### Results\n",
    paste0("![Power2Curve](", plot.file.name, ")"),
    "\n\n\n",
    file = readme.file, append = FALSE)

#### Sample size for Aim 3 (individual reinfection) ----

# Remove objects except those still required
keep.obj <- c("n.sim", "readme.file", "r", "nominal.alpha", "n", "or")
rm(list = ls()[!ls() %in% keep.obj])

# Study design options

# Total sample size - 70% clear so are at risk of reinfection
# round up to the nearest 100 so they can be divided among up to 100 communities
n3 <- round(1 + n * 0.7, -2)

# ...divided among n.communities
n.communities <- c(25, 50, 100)

# Model parameters

# Prevalence of reinfection
p <- 0.5

# normal variance in logit prevalence among communities
community.var <- 2.73 # estimated from RP's data

# odds ratio per binary predictor, or per SD for continuous predictors
# already defined:
or

# how correlated are predictors (because abs(r) > 0 reduces power,
# and it would be unrealistic to assume zero correlation)
# already defined:
r 

# Choose candidate drivers:

# How many drivers?
# individual continuous drivers x 4 (Age, Immunology, Immunology, Baseline infection intensity)
n.cont <- 4
# individual binary drivers x 1 (Hybrid/resistance – 20%)
bin.x.p <- c(`hybrid/resistance presence` = 0.2)

# Total number of individual drivers:
n.x <- n.cont + length(bin.x.p)

# Community level drivers x 4 (generic, but standing for Snail infection rates,
# MDA coverage, Distance to water sites, How many produce sensible numbers?)
# All community level drivers are continuous
n.x.c <- 4

# Further simulation and analysis options

# adjust the significance thresholds for multiple testing (Bonferroni)
alpha.i <- nominal.alpha/n.x   # individual  drivers
alpha.c <- nominal.alpha/n.x.c # community drivers

# Make table of all parameter and design combinations
par.tab <- 
  expand.grid(n = n3, n.communities = n.communities, p = p, 
              community.var = community.var, or = or, 
              r = r, alpha.i = alpha.i, alpha.c = alpha.c, 
              n.sim = round(n.sim/2))

# No of subjects per community
par.tab$n.per.community <- ceiling(par.tab$n / par.tab$n.communities)

# Start the clock
start.time <- Sys.time()

# Loop over all scenarios = rows of par.tab
sim.res <- 
  sapply(1:nrow(par.tab), function(i) {
    
    # Print progress
    print(paste0(round(100*(i-1)/nrow(par.tab)), "% complete"))
    
    # Simulate Xs (the drivers)
    
    # Correlation among individual level drivers
    sigma2 <- diag(n.x)
    sigma2[lower.tri(sigma2)] <- par.tab$r[i]
    sigma2[upper.tri(sigma2)] <- par.tab$r[i]
    
    # Correlation among community level drivers
    sigma2.c <- diag(n.x.c)
    sigma2.c[lower.tri(sigma2.c)] <- par.tab$r[i]
    sigma2.c[upper.tri(sigma2.c)] <- par.tab$r[i]
    
    # Simulate data and analysis n.sim times, returning the mean
    # number of drivers identified at P < alpha
    n.sim.out <-
      mclapply(1:par.tab$n.sim[i], function(k) {
        
        # Make data frame storing each individual ID and their community ID
        dat <- expand.grid(id = factor(1:par.tab$n.per.community[i]), 
                           community = factor(1:par.tab$n.communities[i]))
        dim(dat)
        
        # Simulate individual Xs before dichotomising (so the binary Xs are derived from latent scales)
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
        names(b) <- c("intercept", paste0("X", 1:(ncol(X) - 1)))
        
        # Simulate community level Xs 
        X.c <- mvrnorm(par.tab$n.communities[i], mu = rep(0, n.x.c), Sigma = sigma2.c)
        
        # vector of log odds ratios
        b.c <- log(rep(par.tab$or[i], n.x.c))
        names(b.c) <- paste0("X", n.x + (1:ncol(X.c)))
        b.all <- c(b, b.c)
        rownames(X.c) <- levels(dat$community)
        
        # Add fixed effects
        X.all <- cbind(X, X.c[as.character(dat$community), ])
        colnames(X.all) <- names(b.all)
        dat <- data.frame(dat, X.all)
        head(dat)

        # Add column of no of trials (binary data so n_trials = 1)
        dat$n <- 1
        
        # Simulate infected status
        simdat <- 
          sim.glmm(design.data = dat, 
                   fixed.eff = b.all,
                   rand.V = list(community = par.tab$community.var[i]),  
                   distribution = "binomial")
        
        # fit GLMM and get z-test p-values for Xs
        form <- paste("response ~ (1 | community) +",  paste(names(b.all[-1]), collapse = " + "))
        fit <- glmmTMB(formula(form), family = binomial, data = simdat)
        res.tab <- coef(summary(fit))$cond
        
        # estimate (geometric) mean margin of error of odds ratio estimates across drivers
        # individual drivers
        MoE.i <- exp(1.96 * mean(res.tab[names(b[-1]), "Std. Error"])) - 1
        # community drivers
        MoE.c <- exp(1.96 * mean(res.tab[names(b.c), "Std. Error"])) - 1
        
        # How many drivers were detected (P < alpha)?
        # individual drivers
        n.drivers.sig.i <- sum(res.tab[names(b[-1]), "Pr(>|z|)"] < par.tab$alpha.i[i])
        # community drivers
        n.drivers.sig.c <- sum(res.tab[names(b.c), "Pr(>|z|)"] < par.tab$alpha.c[i])
        
        # output results
        c(n.drivers.sig.i = n.drivers.sig.i, n.drivers.sig.c = n.drivers.sig.c, MoE.i = MoE.i, MoE.c = MoE.c)
        
      }, mc.cores = detectCores())
    
    # take mean number of drivers (Xs) detected across simulated data sets
    # and the geometric mean margin of error
    n.sim.out.tab <- do.call("rbind", n.sim.out)
    c(n.drivers.sig.i = mean(n.sim.out.tab[, "n.drivers.sig.i"]),
      n.drivers.sig.c = mean(n.sim.out.tab[, "n.drivers.sig.c"]),
      MoE.i = exp(mean(log(n.sim.out.tab[, "MoE.i"]))),
      MoE.c = exp(mean(log(n.sim.out.tab[, "MoE.c"]))))
  })

# Stop the clock
run.time <- Sys.time() - start.time
print(run.time)

# Estimate power as the proportion of the n.x drivers with p < alpha
# For the individual drivers
par.tab$prop.drivers.i <- round(sim.res["n.drivers.sig.i", ]/n.x, 5)
# For the community drivers
par.tab$prop.drivers.c <- round(sim.res["n.drivers.sig.c", ]/n.x.c, 5)

# Estimate mean margin of error of odds ratio estimates 
# For the individual drivers
par.tab$or.margin.of.error.i <- round(sim.res["MoE.i", ], 5)
# For the community drivers
par.tab$or.margin.of.error.c <- round(sim.res["MoE.c", ], 5)

# Export results to CSV file with time stamp in file name
file.name.power3 <- 
  paste0("results/schisto_power3_", 
         substr(gsub(":", "", (gsub(" ", "-", Sys.time()))), 1, 15), ".csv")
write.csv(par.tab, file.name.power3, row.names = FALSE, quote = FALSE)
rm(par.tab)
par.tab <- read.csv(file.name.power3)

# Make plots of results
par.tab$alpha.i <- factor(paste("alpha =", round(par.tab$alpha.i, 4)))
par.tab$alpha.c <- factor(paste("alpha =", round(par.tab$alpha.c, 4)))
par.tab$n <- factor(par.tab$n, n3)
par.tab$n.communities.fac <- 
  factor(paste("N communities =", par.tab$n.communities), paste("N communities =", n.communities))
par.tab$r <- factor(paste("Correlation among drivers (r) =", par.tab$r))
par.tab$p <- factor(paste("Prevalence =", par.tab$p))

# plot power
# For the individual drivers
power.plot.i <-
  ggplot(data = par.tab, aes(x = or, y = prop.drivers.i, color = n, shape = n, group = n)) + 
  geom_hline(yintercept = 0.8, linewidth = 0.3, linetype = 2) +
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  facet_wrap(~ p + n.communities.fac) +
  xlab("Odds ratio") +
  ylab("Power to detect individual drivers") +
  labs(caption = file.name.power3) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
power.plot.i
power.plot.file.name.i <- "schisto_power3.i.png"
ggsave(power.plot.file.name.i, width = 6, height = 6)

# plot power
# For the community drivers
power.plot.c <-
  ggplot(data = par.tab, aes(x = or, y = prop.drivers.c, color = n, shape = n, group = n)) + 
  geom_hline(yintercept = 0.8, linewidth = 0.3, linetype = 2) +
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  facet_wrap(~ p + n.communities.fac) +
  xlab("Odds ratio") +
  ylab("Power to detect community drivers") +
  labs(caption = file.name.power3) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
power.plot.c
power.plot.file.name.c <- "schisto_power3.c.png"
ggsave(power.plot.file.name.c, width = 6, height = 6)


# plot margin of error
par.tab$or <- factor(paste("Odds ratio =", par.tab$or), 
                     paste("Odds ratio =", or))
# For the individual drivers
moe.plot.i <-
  ggplot(data = par.tab, aes(x = n.communities, y = 100 * or.margin.of.error.i, color = n, shape = n, group = n)) + 
  geom_point() +
  geom_line() +
  ylim(1, NA) +
  facet_wrap(~ or + p, ncol = 2) +
  xlab("N communities") +
  ylab("Margin of error (%) for individual driver OR") +
  scale_x_continuous(breaks = c(0, n.communities), 
                     limits = c(min(n.communities) - 5, max(n.communities) + 5)) + 
  labs(caption = file.name.power3) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
moe.plot.i
moe.plot.file.name.i <- "schisto_moe3.i.png"
ggsave(moe.plot.file.name.i, width = 6, height = 9)

# For the community drivers
moe.plot.c <-
  ggplot(data = par.tab, aes(x = n.communities, y = 100 * or.margin.of.error.c, color = n, shape = n, group = n)) + 
  geom_point() +
  geom_line() +
  ylim(1, NA) +
  facet_wrap(~ or + p, ncol = 2) +
  xlab("N communities") +
  ylab("Margin of error (%) for community driver OR") +
  scale_x_continuous(breaks = c(0, n.communities), 
                     limits = c(min(n.communities) - 5, max(n.communities) + 5)) + 
  labs(caption = file.name.power3) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
moe.plot.c
moe.plot.file.name.c <- "schisto_moe3.c.png"
ggsave(moe.plot.file.name.c, width = 6, height = 9)



# output methods and results to README.md

cat("## Sample size calculation for Aim 3: identifying individual- and community-level drivers of re-infection following clearance\n\n",
    "### Methods\n\n",
    "The aim of this power analysis is to estimate power to detect individual- and community-level drivers of",
    "schistosomiasis re-infection following clearance, and the expected margin of error around",
    "community-level driver odds ratio estimates. This is a simulation-based power analysis,",
    "where the study data is simulated and analysed multiple times under varying study design",
    "scenarios, and power is estimated as the proportion of simulated analyses that achieve the",
    "desired outcome (detecting a true driver of infection). The association between the",
    "outcome (infection) and each driver is estimated and tested in a multivariable GLMM.",
    "For this analysis, power is defined as the proportion of drivers that are significantly associated",
    paste("with the outcome, averaged across", unique(par.tab$n.sim), "simulated data analyses per scenario."),
    "Power and margin of error (i.e. 95% confidence interval) in odds ratio estimation are presented",
    "separately for individual and community-level drivers.\n\n",
    "The following assumptions are made:\n-",
    paste(n.x, "drivers are associated with the outcome, of which", length(bin.x.p), "are binary and", 
          n.cont, "continuous."),
    paste0("The prevalences of the binary drivers are: ", paste(bin.x.p, collapse = ", "), ","),
    "representing ", paste(names(bin.x.p), collapse = ", "), ".\n",
    paste("- A further", n.x.c, "continuous community-level drivers are associated with the outcome.\n"),
    paste0("- The drivers are correlated with each other, with a common correlation coefficient of ", r, "."),
    "We don’t know what the true correlation is among drivers, but moderate correlations are likely",
    "and neglecting them will give optimistic power estimates.\n",
    paste0("- Log odds of re-infection varies among communities with a variance of ", 
           community.var, ".\n"),
    "- In order to control inflation of the number of false positive results due to multiple testing of",
    n.x, "individual-level drivers and", n.x.c, "community-level drivers, the significance threshold of", 
    nominal.alpha,
    paste0("was Bonferroni-adjusted to ", alpha.i, " and ", alpha.c, 
           " respectively.\n\n"),
    "We explore the effect on power of varying the following study design choices/assumptions:\n",
    paste0("- Total sample size: ", paste(n3, collapse = ", "), ".\n"),  
    paste0("- Community sample size (number communities sampled): ", paste(n.communities, collapse = ", "), ".\n"),  
    paste0("- Prevalence of re-infection: ", paste(p, collapse = ", "), ".\n"),
    "- The strength of association between each driver and re-infection,",
    "defined as an odds ratio per standard deviation",
    paste0("unit for continuous community-level drivers: ", paste(or[or != 1], collapse = ", "), ".\n\n"),
    "Full details are provided in the script",
    "[PowerAnalysis.R](https://github.com/pcdjohnson/SchistoDrivers/blob/main/PowerAnalysis.R).",
    "Results are output as CSV to the [results](https://github.com/pcdjohnson/SchistoDrivers/tree/main/results)",
    paste0("directory and plotted to [", 
           power.plot.file.name.i, 
           "](https://github.com/pcdjohnson/SchistoDrivers/blob/main/", power.plot.file.name.i, 
           "), [",
           power.plot.file.name.c, 
           "](https://github.com/pcdjohnson/SchistoDrivers/blob/main/", power.plot.file.name.c, 
           "), [",
           moe.plot.file.name.i, 
           "](https://github.com/pcdjohnson/SchistoDrivers/blob/main/", moe.plot.file.name.i, 
           ") and [",
           moe.plot.file.name.c, 
           "](https://github.com/pcdjohnson/SchistoDrivers/blob/main/", moe.plot.file.name.c,").\n\n"),
    "### Results\n",
    paste0("![Power3CurveInd](", power.plot.file.name.i, ")"),
    "\n\n\n",
    paste0("![Power3CurveCom](", power.plot.file.name.c, ")"),
    "\n\n\n",
    paste0("![MoE3CurveInd](", moe.plot.file.name.i, ")"),
    "\n\n\n",
    paste0("![MoE3CurveCom](", moe.plot.file.name.c, ")"),
    "\n\n\n",
    file = readme.file, append = TRUE)



#### Sample size for Aim 4 (community-level drivers) ----

# Remove objects except those still required
keep.obj <- c("n.sim", "readme.file", "r", "nominal.alpha", "n", "or")
rm(list = ls()[!ls() %in% keep.obj])

# Study design options

# Total sample size - already defined
n

# ...divided among n.communities
n.communities <- c(25, 50, 100)

# Model parameters

# Prevalence
p <- c(0.1, 0.5)

# normal variance in logit prevalence among communities
community.var <- 2.73 # estimated from RP's data

# odds ratio per binary predictor, or per SD for continuous predictors
# already defined:
or

# how correlated are predictors (because abs(r) > 0 reduces power,
# and it would be unrealistic to assume zero correlation)
# already defined:
r 

# Choose candidate drivers:

# How many drivers?
# Continuous drivers x 4 (generic, but could be snail abundance,  infectivity combined risk measure,
# MDA coverage, estimate of % systematic non treated

n.cont <- 4

# Total number of drivers:
n.x <- n.cont

# Further simulation and analysis options

# adjust the significance thresholds for multiple testing (Bonferroni)
alpha <- c(nominal.alpha/n.x)

# Make table of all parameter and design combinations
par.tab <- 
  expand.grid(n = n, n.communities = n.communities, p = p, 
              community.var = community.var, or = or, 
              r = r, alpha = alpha, n.sim = round(n.sim/2))

# No of subjects per community
par.tab$n.per.community <- round(par.tab$n / par.tab$n.communities)

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
        
        # Simulate Xs 
        X <- mvrnorm(par.tab$n.communities[i], mu = rep(0, n.x), Sigma = sigma2)
        
        # vector of intercept (log odds) and log odds ratios
        b <- c(qlogis(par.tab$p[i]), log(rep(par.tab$or[i], n.x)))
        names(b) <- c("intercept", paste0("X", 1:ncol(X)))
        
        # Make data frame storing each individual ID and their community ID
        dat <- expand.grid(id = factor(1:par.tab$n.per.community[i]), 
                           community = factor(1:par.tab$n.communities[i]))
        
        # Add fixed effects
        rownames(X) <- levels(dat$community)
        dat <- data.frame(dat, X[as.character(dat$community), ])
        
        # Add column of no of trials (binary data so n_trials = 1)
        dat$n <- 1
        
        # Simulate infected status
        simdat <- 
          sim.glmm(design.data = dat, 
                   fixed.eff = b,
                   rand.V = list(community = par.tab$community.var[i]),  
                   distribution = "binomial")
        
        # fit GLMM and get z-test p-values for Xs
        form <- paste("response ~ (1 | community) +",  paste(names(b[-1]), collapse = " + "))
        fit <- glmmTMB(formula(form), family = binomial, data = simdat)
        res.tab <- coef(summary(fit))$cond
        
        # estimate (geometric) mean margin of error of odds ratio estimates across all drivers
        # MoE = exp(1.96 * SE) - 1
        MoE <- exp(1.96 * mean(res.tab[names(b[-1]), "Std. Error"])) - 1
        
        # How many drivers were detected (P < alpha)?
        n.drivers.sig <- sum(res.tab[names(b[-1]), "Pr(>|z|)"] < par.tab$alpha[i])
        
        # output results
        c(n.drivers.sig = n.drivers.sig, MoE = MoE)
        
      }, mc.cores = detectCores())
    
    # take mean number of drivers (Xs) detected across simulated data sets
    # and the geometric mean margin of error
    n.sim.out.tab <- do.call("rbind", n.sim.out)
    c(n.drivers.sig = mean(n.sim.out.tab[, "n.drivers.sig"]),
      MoE = exp(mean(log(n.sim.out.tab[, "MoE"]))))
  })

# Stop the clock
run.time <- Sys.time() - start.time
print(run.time)

# Estimate power as the proportion of the n.x drivers with p < alpha
par.tab$prop.drivers <- round(sim.res["n.drivers.sig", ]/n.x, 5)

# Estimate mean margin of error of odds ratio estimates 
par.tab$or.margin.of.error <- round(sim.res["MoE", ], 5)

# Export results to CSV file with time stamp in file name
file.name.power4 <- 
  paste0("results/schisto_power4_", 
         substr(gsub(":", "", (gsub(" ", "-", Sys.time()))), 1, 15), ".csv")
write.csv(par.tab, file.name.power4, row.names = FALSE, quote = FALSE)
rm(par.tab)
par.tab <- read.csv(file.name.power4)

# Make plots of results
par.tab$alpha <- factor(paste("alpha =", round(par.tab$alpha, 4)))
par.tab$n <- factor(par.tab$n, n)
par.tab$n.communities.fac <- 
  factor(paste("N communities =", par.tab$n.communities), paste("N communities =", n.communities))
par.tab$r <- factor(paste("Correlation among drivers (r) =", par.tab$r))
par.tab$p <- factor(paste("Prevalence =", par.tab$p))

# plot power
power.plot <-
  ggplot(data = par.tab, aes(x = or, y = prop.drivers, color = n, shape = n, group = n)) + 
  geom_hline(yintercept = 0.8, linewidth = 0.3, linetype = 2) +
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  facet_wrap(~ p + n.communities.fac) +
  xlab("Odds ratio") +
  ylab("Power") +
  labs(caption = file.name.power4) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
power.plot
power.plot.file.name <- "schisto_power4.png"
ggsave(power.plot.file.name, width = 6, height = 6)


# plot margin of error
par.tab$or <- factor(paste("Odds ratio =", par.tab$or), 
                     paste("Odds ratio =", or))

moe.plot <-
  ggplot(data = par.tab, aes(x = n.communities, y = 100 * or.margin.of.error, color = n, shape = n, group = n)) + 
  geom_point() +
  geom_line() +
  ylim(1, NA) +
  facet_wrap(~ or + p, ncol = nlevels(par.tab$p)) +
  xlab("N communities") +
  ylab("Margin of error (%)") +
  scale_x_continuous(breaks = c(0, n.communities), 
                     limits = c(min(n.communities) - 5, max(n.communities) + 5)) + 
  labs(caption = file.name.power4) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
moe.plot
moe.plot.file.name <- "schisto_moe4.png"
ggsave(moe.plot.file.name, width = 6, height = 9)

# output methods and results to README.md

cat("## Sample size calculation for Aim 4: identifying community-level drivers of schistosomiasis infection\n\n",
    "### Methods\n\n",
    "The aim of this power analysis is to estimate power to detect community-level drivers of",
    "schistosomiasis infection, and the expected margin of error around",
    "community-level driver odds ratio estimates. This is a simulation-based power analysis,",
    "where the study data is simulated and analysed multiple times under varying study design",
    "scenarios, and power is estimated as the proportion of simulated analyses that achieve the",
    "desired outcome (detecting a true driver of infection). The association between the",
    "outcome (infection) and each driver is estimated and tested in a multivariable GLMM.",
    "For this analysis, power is defined as the proportion of drivers that are significantly associated",
    paste("with the outcome, averaged across", unique(par.tab$n.sim), "simulated data analyses per scenario.\n\n"),
    "The following assumptions are made:\n-",
    paste(n.x, "continuous drivers are associated with the outcome.\n"),
    paste0("- The drivers are correlated with each other, with a common correlation coefficient of ", r, "."),
    "We don’t know what the true correlation is among drivers, but moderate correlations are likely",
    "and neglecting them will give optimistic power estimates.\n",
    paste0("- Log odds of infection prevalence varies among communities with a variance of ", 
           community.var, ".\n"),
    "- In order to control inflation of the number of false positive results due to multiple testing of",
    n.x, "drivers, the significance threshold of", nominal.alpha,
    paste0("is Bonferroni-adjusted to ", alpha, ", i.e. a driver is significant if P < ", alpha, ".\n\n"),
    "We explore the effect on power of varying the following study design choices/assumptions:\n",
    paste0("- Total sample size: ", paste(n, collapse = ", "), ".\n"),  
    paste0("- Community sample size (number communities sampled): ", paste(n.communities, collapse = ", "), ".\n"),  
    paste0("- Prevalence of infection: ", paste(p, collapse = ", "), ".\n"),
    "- The strength of association between each driver and failure to clear,",
    "defined as an odds ratio per standard deviation",
    paste0("unit for continuous community-level drivers: ", paste(or[or != 1], collapse = ", "), ".\n\n"),
    "Full details are provided in the script",
    "[PowerAnalysis.R](https://github.com/pcdjohnson/SchistoDrivers/blob/main/PowerAnalysis.R).",
    "Results are output as CSV to the [results](https://github.com/pcdjohnson/SchistoDrivers/tree/main/results)",
    paste0("directory and plotted to [", power.plot.file.name, 
           "](https://github.com/pcdjohnson/SchistoDrivers/blob/main/", power.plot.file.name, ") and [",
           moe.plot.file.name, 
           "](https://github.com/pcdjohnson/SchistoDrivers/blob/main/", moe.plot.file.name,").\n\n"),
    "### Results\n",
    paste0("![Power4Curve](", power.plot.file.name, ")"),
    "\n\n\n",
    paste0("![MoE4Curve](", moe.plot.file.name, ")"),
    "\n\n\n",
    file = readme.file, append = TRUE)


