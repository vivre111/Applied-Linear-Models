#--- categorical predictors ------------------------------------------------
# data
mercury <- read.csv("fishermen_mercury.csv")
head(mercury)

fplevels <- c("none", "muscle", "muscle_whole", "whole")
mercury$fishpart <- factor(mercury$fishpart, levels=c(0,1,2,3),labels = fplevels) # custom ordering


# Question: How to test H0: fishpart has no effect in presence of weight?

# full model: 
Mfull <- lm(MeHg ~ factor(fishpart) + weight, data = mercury) 

# reduced model: all levels of fishpart have same effect (intercept)
Mred <- lm(MeHg ~ weight, data = mercury)

# degrees of freedom in full and reduced models
dff <- Mfull$df
dfr <- Mred$df

# F-statistic, using LR definition
SSRes_full <- sum(residuals(Mfull)^2) # sum-of-squares for Mfull
SSRes_red <- sum(residuals(Mred)^2) # sum-of-squares for Mred
Fobs <- (SSRes_red-SSRes_full)/(dfr-dff) / (SSRes_full/dff)
c(Fobs, (SSRes_red-SSRes_full)/(dfr-dff) / sigma(Mfull)^2) # equivalent denominator

# F-statistic using subset of MLE of full model (gamma)
gam.ind <- 2:4 # indices corresponding to beta's for fishpart
gam.hat <- coef(Mfull)[gam.ind] # subset of MLE
gam.ve <- vcov(Mfull)[gam.ind,gam.ind] # gam.ve = sigma.hat^2 * V
Fobs2 <- t(gam.hat) %*% solve(gam.ve, gam.hat) / length(gam.hat)

# check two version are identical
c(Fobs, Fobs2)

# p-value: large values of Fobs are more evidence against H0
pf(Fobs, df1 = (dfr-dff), df2 = dff, lower.tail = FALSE)

# fully R version:
anova(Mred, Mfull)








#--- colinearity simulation -----------------------------------------------------------

# true model: y ~ b0 + b1*x1 + b2*x2 + eps

beta0 <- 1
beta1 <- 2
beta2 <- 0 ## x2 has no effect
sigma <- 1.5

# simulate data
n <- 100 # sample size

# let x1 and x2 be extremely correlated
# covariate matrix
rho <- 0.99
varX <- rbind(c(1, rho), c(rho, 1))
colnames(varX) <- c("x1", "x2")
rownames(varX) <- colnames(varX)
varX # extremely correlated x1 and x2


nreps <- 1000 # replicate the whole experiment nreps times
sim.out <- matrix(NA, nreps, 4) # store beta.hats and their p-values
colnames(sim.out) <- c("beta.hat1", "beta.SE1", "beta.hat2", "beta.SE2")
rownames(sim.out) <- paste0("rep", 1:nreps)

library(mvtnorm)
set.seed(6)
for(ii in 1:nreps) {
  # generate MVN covariates via cholesky decomp:
  # Z <- matrix(rnorm(n*2), n, 2) # n samples from 2 iid N(0,1)
  # X <- Z %*% chol(varX) # n samples from N(0, varX)
  # generate MVN covariates via mvtnorm package:
  X <- rmvnorm(n,c(0,0),varX) # covariates
  x1 <- X[,1]
  x2 <- X[,2]
  eps <- rnorm(n, sd = sigma) # error
  y <- beta0 + beta1*x1 + beta2*x2 + eps # response
  M1 <- lm(y ~ x1 ) # fit regression on x1 only
  M2 <- lm(y ~ x1 + x2) # fit regression on both x1 and x2
  sim.out[ii,] <- round(c(coef(M1)[2], sqrt(diag(vcov(M1))[2]), coef(M2)[2], sqrt(diag(vcov(M2))[2])),2)
}
round(apply(sim.out,2,mean),2)

round(apply(sim.out[,c(1,3)],2,sd),2)












#--- variance inflation factor PRACTICE --------------------------------------------


library(readr)
real_estate <- read_csv("real_estate.txt")

## fit model
g <- lm(SalePrice~SqrFeet+Bedrooms+Lotsize,data=real_estate)
summary(g)$coef[2,]

## regress Bedrooms on SqrFeet and get residuals
g_bed <- lm(Bedrooms~SqrFeet,data=real_estate)
real_estate$e_bed <- residuals(g_bed)
## regress Lotsize on SqrFeet and get residuals
g_lot <- lm(Lotsize~SqrFeet,data=real_estate)
real_estate$e_lot <- residuals(g_lot)

# fit model with idealized design matrix (replacing other covariates with residuals)
g_ideal <- lm(SalePrice~SqrFeet+e_bed+e_lot,data=real_estate)
summary(g_ideal)$coef[2,]

# compute ratio of variance estimates from original model fit and idealized model fit 
(summary(g)$coef[2,2]/summary(g_ideal)$coef[2,2])^2

# regress SqrFeet on other covariates
g_feet <- lm(SqrFeet~Bedrooms+Lotsize,data=real_estate)
## compute R-squared 
R2 <- cor(real_estate$SqrFeet,fitted(g_feet))^2 
R2 <- summary(g_feet)$r.squared ## other way
VIF_1 <- 1/(1-R2) ## compute VIF for SqrFeet


## compare fitted values from regular and idealized model fit
pred_g <- predict(g,interval="confidence")
pred_ideal <- predict(g_ideal,interval="confidence")
## same fitted values
plot(pred_g[,1]~pred_ideal[,1])
cor(pred_g[,1],pred_ideal[,1])

## compare CI widths 
widths_g <- pred_g[,3]-pred_g[,2]
widths_ideal <- pred_ideal[,3]-pred_ideal[,2]
## same fitted values
plot(widths_g~widths_ideal)
cor(widths_g,widths_ideal)


