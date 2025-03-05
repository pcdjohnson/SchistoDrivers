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
n.sim <- 1000 # number of data sets to simulate (divided by 2 for the GLMM analysis, because it's slow)

#### Sample size calculation 1 ----

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
target.power.AUC <- 0.9
n.AUC <- 
  power.anova.test(groups = 3, power = target.power.AUC, 
                   between.var = var(c(0, delta.AUC, delta.AUC * 2)),
                   within.var = 1, sig.level = nominal.alpha)$n

# output methods and results to README.md

cat("# SchistoDrivers\n\n",
    "## Sample size calculation 1: randomised controlled trial for the effect of",
    "food prior to treatment on praziquantel absorption (Rapid Answer Project 1) \n\n",
    "### Methods\n\n",
    "The aim of this power analysis is to estimate power to detect a difference in praziquantel absorption",
    "(measured as area under the curve [AUC] of praziquantel metabolites) between three groups:\n",
    "- people who have taken food at home prior to praziquantel treatment;\n",
    "- people who have brought in food to be taken prior to praziquantel treatment;\n",
    "- people for whom food has been provided on site prior to praziquantel treatment.\n\n",
    "The null hypothesis is that mean AUC is equal across the three groups.",
    "The alternative hypothesis is that mean AUC differs between the three groups.",
    "The effect size assumed here is that the group AUC means differ by",
    delta.AUC, "standard deviations from the lowest mean to the intermediate mean,",
    "and by", delta.AUC, "standard deviations from the intermediate mean to the highest mean.",
    "Target power is", paste0(target.power.AUC * 100, "%,"), 
    "which is set higher that the standard minimal threshold of 80%.",
    "Lowering the risk of a type II error in this way is justified by the likely",
    "public health benefit of finding an optimal strategy, if it exists, for",
    "taking food prior to praziquantel treatment.",
    "The null hypothesis is rejected if P <", nominal.alpha,
    "from a one-way ANOVA. The required sample",
    "size per group was calculated using the R function *power.anova.test*.",
    "Full details are provided in the script",
    "[PowerAnalysis.R](https://github.com/pcdjohnson/SchistoDrivers/blob/main/PowerAnalysis.R).",
    "\n\n",
    "### Results\n",
    ceiling(n.AUC), "people would be required per group in order to achieve",
    paste0(target.power.AUC * 100, "%"), 
    "power to detect a", delta.AUC, 
    "standard deviation difference between each of the three groups.",
    "For context, we note that [Kovač et al. (2018)](https://doi.org/10.1128/aac.02253-17)",
    "detected significant differences between groups of children treated with different",
    "praziquantel doses in both metabolite AUC and schistosome clearance rate with",
    "per-group sample sizes ranging from 29 to 47. Therefore we expect that", ceiling(n.AUC),
    "per group will give high sensitivity to detect the effects of different strategies for eating prior to treatment.",
    "\n\n\n",
    file = readme.file, append = FALSE)

#### Sample size calculation 2 ----
# RAP2: Miracidial trial to enable time and money saving mass storage in the field
# (Rapid Answer Project 2), unlocking greater flexibility for picking specific miracidia
# for individual and community-level genetic follow ups. 

n.mira <- 28
p.haplo <- 0.05

cat("## Sample size calculation 2: Miracidial trial to enable time and money",
    "saving mass storage in the field (Rapid Answer Project 2) \n\n",
    "The aim of this power analysis is to estimate the power to detect hybrid or resistant haplotypes",
    "present at at least", paste0(p.haplo * 100, "%"), "frequency in a population of miracidia.",
    "Assuming complete admixture within each community, by sampling", 
    n.mira, "miracidia ", paste0("(", 2 * n.mira, " haplotypes)"),
    "we have a", paste0(round(100 * (1 - (1 - p.haplo)^(2 * n.mira)), 1), "%"), 
    "chance of observing hybrid or resistant haplotypes present at", paste0(p.haplo * 100, "%"), 
    "frequency in the population.",
    "See the script [PowerAnalysis.R](https://github.com/pcdjohnson/SchistoDrivers/blob/main/PowerAnalysis.R) for details.",
    "\n\n\n",
    file = readme.file, append = TRUE)


#### Sample size calculation 3 ----

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
file.name.power3 <- 
  paste0("results/schisto_power3_", 
         substr(gsub(":", "", (gsub(" ", "-", Sys.time()))), 1, 15), ".csv")
write.csv(par.tab, file.name.power3, row.names = FALSE, quote = FALSE)
rm(par.tab)
par.tab <- read.csv(file.name.power3) 

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
  labs(caption = file.name.power3) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
power.plot
plot.file.name <- "schisto_power3.png"
ggsave(plot.file.name, width = 6, height = 6)

# output methods and results to README.md

cat("## Sample size calculation 3: identifying drivers of schistosomiasis praziquantel treatment failure (main study)\n\n",
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
           "](https://github.com/pcdjohnson/SchistoDrivers/blob/main/schisto_power3.png).\n\n"),
    "### Results\n",
    paste0("![PowerCurve](", plot.file.name, ")"),
    "\n\n\n",
    file = readme.file, append = TRUE)

#### Sample size calculation 4 ----

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
par.tab$or <- factor(paste("Odds ratio =", par.tab$or), 
                     paste("Odds ratio =", or))
par.tab$r <- factor(paste("Correlation among drivers (r) =", par.tab$r))
par.tab$p <- factor(paste("Prevalence =", par.tab$p))

# plot power
power.plot <-
  ggplot(data = par.tab, aes(x = or, y = prop.drivers, color = n, shape = n, group = n)) + 
  geom_hline(yintercept = 0.8, linewidth = 0.3, linetype = 2) +
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  facet_wrap(~ p + n.communities) +
  xlab("Odds ratio") +
  ylab("Power") +
  labs(caption = file.name.power4) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
power.plot
power.plot.file.name <- "schisto_power4.png"
ggsave(power.plot.file.name, width = 6, height = 6)

# plot margin of error
#par.tab$n <- factor(par.tab$n, rev(levels(par.tab$n)))
moe.plot <-
  ggplot(data = par.tab, aes(x = n.communities, y = 100 * or.margin.of.error, color = n, shape = n, group = n)) + 
  geom_point() +
  geom_line() +
  ylim(1, NA) +
  facet_wrap(~ p + or, nrow = nlevels(par.tab$p)) +
  xlab("N communities") +
  ylab("Margin of error (%)") +
  scale_x_continuous(breaks = c(0, n.communities), 
                     limits = c(min(n.communities) - 5, max(n.communities) + 5)) + 
  labs(caption = file.name.power4) + # link results filename to plot
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
moe.plot
moe.plot.file.name <- "schisto_moe4.png"
ggsave(moe.plot.file.name, width = 6, height = 6)

# output methods and results to README.md

cat("## Sample size calculation 4: identifying community-level drivers of schistosomiasis infection\n\n",
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
    paste0("![PowerCurve](", power.plot.file.name, ")"),
    "\n\n\n",
    paste0("![PowerCurve](", moe.plot.file.name, ")"),
    "\n\n\n",
    file = readme.file, append = TRUE)


