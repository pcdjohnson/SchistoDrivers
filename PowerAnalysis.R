# Load packages
library(ggplot2)
library(parallel)
library(MASS)

# Clear out objects
rm(list = ls())


# Study design options
n <- c(400, 600, 800, 1000, 1200)

# Model parameters
p <- 0.3 # failure to clear proportion
or <- seq(1, 2, 0.25)
r <- c(0, 0.5)

# Choose candidate drivers:

# How many drivers?
# Generic continuous drivers x 4 (PK, genetic SNP score [????], etc)
# 2 binary drivers
n.x <- 6

# For the binary drivers, what proportions?
bin.x.p <- c(0.2, 0.5) # HIV 20%, malaria 50%

# simulation and analysis options
n.sim <- 5000
alpha <- c(0.05, 0.05/n.x)



# Make table of all parameter and design combinations
par.tab <- expand.grid(n = n, p = p, or = or, r = r, alpha = alpha, n.sim = n.sim)

start.time <- Sys.time()
sim.res <- 
  sapply(1:nrow(par.tab), function(i) {
    # i <- nrow(par.tab)
    
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

par.tab$prop.drivers <- sim.res/n.x



par.tab$alpha <- factor(paste("alpha =", round(par.tab$alpha, 4)))
par.tab$n <- factor(par.tab$n, rev(n))
par.tab$r <- factor(paste("Correlation among drivers (r) =", par.tab$r))

power.plot <-
  ggplot(data = par.tab, aes(x = or, y = prop.drivers, color = n, shape = n, group = n)) + 
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  facet_wrap(~ r + alpha, ncol = 2) +
  xlab("Odds ratio") +
  ylab("Power") +
  labs(caption = paste(paste0("no of simulations = ", n.sim), substr(Sys.time(), 1, 16), sep = "\n")) +
  theme(plot.caption = element_text(colour = "grey60", size = rel(0.75)),
        plot.caption.position = "plot")
power.plot
ggsave("schisto_power.pdf", width = 6, height = 6)

