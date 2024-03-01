# Load packages
library(ggplot2)
library(parallel)
library(MASS)

# Clear out objects
rm(list = ls())


# Study design options
n <- seq(500, 3000, 500)

# Model parameters
p <- 0.3 # failure to clear proportion
or <- seq(1, 2, 0.25)
r <- c(0, 0.25)

# Choose candidate drivers:

# How many drivers?
# Generic continuous drivers x 4 (PK, genetic SNP score [????], etc)
# 2 binary drivers
n.x <- 10

# For the binary drivers, what proportions?
bin.x.p <- c(0.2, 0.5, 0.2, 0.2) # HIV 20%, malaria 50%, STH 20%, hybrid/resistance presence 20%

# simulation and analysis options
n.sim <- 1000
alpha <- c(0.05, 0.05/n.x)



# Make table of all parameter and design combinations
par.tab <- expand.grid(n = n, p = p, or = or, r = r, alpha = alpha, n.sim = n.sim)

start.time <- Sys.time()
sim.res <- 
  sapply(1:nrow(par.tab), function(i) {
    
    # i <- nrow(par.tab)
    
    print(paste0(round(100*(i-1)/nrow(par.tab)), "% complete"))
    
    # Simulate Xs
    
    # Correlation among Xs
    sigma2 <- diag(n.x)
    sigma2[lower.tri(sigma2)] <- par.tab$r[i]
    sigma2[upper.tri(sigma2)] <- par.tab$r[i]
    
    # Simulate data and analysis n.sim times, returning the mean
    # number of drivers identified at P < alpha
    n.sim.out <-
      mclapply(1:n.sim, function(k) {
        # Simulate Xs before dichotomising
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
        
        
        # scale all Xs to have zero mean and SD=1
        #X[, -1] <- apply(X[, -1], 2, scale)
        #apply(X[, -1], 2, sd)
        
        # vector or intercept (log odds) and log odds ratios
        b <- c(qlogis(par.tab$p[i]), log(rep(par.tab$or[i], n.x)))
        
        # Simulate failure to clear
        response <- rbinom(par.tab$n[i], 1, plogis(X %*% b))
        
        # Put everything into a data frame
        dat <- data.frame(X[, -1], response = response)
        
        # fit GLM
        fit <- glm(response ~ X[, -1], family = binomial)
        res.tab <- coef(summary(fit))
        
        
        # How many Xs were "detected" (P < alpha)?
        hits <- sum(res.tab[-1, "Pr(>|z|)"] < par.tab$alpha[i])
        hits
      }, mc.cores = detectCores())
    mean(unlist(n.sim.out))
  })
run.time <- Sys.time() - start.time
print(run.time)

#' Estimate power as the proportion of the n.x drivers with p < alpha
par.tab$prop.drivers <- sim.res/n.x


#' Export results to CSV file with time stamp in file name
file.name <- paste0("results/schisto_power", 
                    substr(gsub(":", "", (gsub(" ", "-", Sys.time()))), 1, 15), ".csv")
write.csv(par.tab, file.name, row.names = FALSE, quote = FALSE)


#' Make plot of results
par.tab$alpha <- factor(paste("alpha =", round(par.tab$alpha, 4)))
par.tab$n <- factor(par.tab$n, rev(n))
par.tab$r <- factor(paste("Correlation among drivers (r) =", par.tab$r))

power.plot <-
  ggplot(data = par.tab, aes(x = or, y = prop.drivers, color = n, shape = n, group = n)) + 
  geom_hline(yintercept = 0.8, linewidth = 0.3, linetype = 2) +
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  facet_wrap(~ r + alpha, ncol = 2) +
  xlab("Odds ratio") +
  ylab("Power") +
  labs(caption = file.name) +
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
power.plot
ggsave("schisto_power.pdf", width = 6, height = 6)


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
  lambda <- 0.2128099
  s <- 2.45669
  theta <- lambda^2 / (s^2 - lambda) # calculate dispersion parameter from SD and mean
  # check we get the right SD
  sd(rnbinom(100000, size = theta, mu = lambda))
  
  # estimate power
  power.estimate <- 
    unlist(mclapply(N.per.group, function(N) {
      # N <- 100
      p.value <-
        replicate(n.sim, {
          dat <- expand.grid(id = 1:N, group = c("fed", "unfed"))
          dat$lambda <- lambda * ifelse(dat$group == "fed", 1, unfed.effect)
          dat$epg <- rnbinom(nrow(dat), size = theta, mu = dat$lambda)
          (tab <- table(dat$epg, dat$group))
          if(nrow(tab) <= 2) return(1) # if the data is unanalysable (e.g. all zeros) P=1 (effectively)
          coef(summary(glm.nb(epg ~ group, data = dat)))["groupunfed", "Pr(>|z|)"]
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
