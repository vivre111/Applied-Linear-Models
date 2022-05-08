#--- criteria illustration -----------------------------------------

satis <- read.csv("satisfaction.csv") # load the data
head(satis)
M <- lm(Satisfaction ~ ., data = satis) # regress on all variables

# throw in garbage predictors: iid N(0,1) not related to anything
n <- nrow(satis)
ng <- 41 # max number of predictors
# set.seed(20)
satis2 <- matrix(rnorm(n*ng), n, ng)
colnames(satis2) <- paste0("Grbg_", 1:ng)
satis2 <- cbind(satis, satis2)
#satis2[1:6,1:10]
## calculate R^2, R^2_adj, etc. upon adding in garbage predictors one at a time
R2 <- rep(NA, ng+1)
R2a <- rep(NA, ng+1)
aic <- rep(NA, ng+1)
bic <- rep(NA, ng+1)
mse <- rep(NA, ng+1)
for(ii in 0:ng) {
  M <- lm(Satisfaction ~ ., data = satis2[1:(4+ii)])
  R2[ii+1] <- summary(M)$r.squared
  R2a[ii+1] <- summary(M)$adj.r.squared
  aic[ii+1] <- AIC(M)
  bic[ii+1] <- BIC(M)
  mse[ii+1] <- sigma(M)^2 
}

# plot
plot(0:ng, mse, xlab = "# of Garbage Predictors", ylab = "",
     main = "MSE",pch = 16, type = "b")

plot(0:ng, R2, xlab = "# of Garbage Predictors", ylab = "",
     main = "Coefficient of Determination",
     pch = 16, type = "b", ylim = range(R2, R2a))
points(0:ng, R2a, pch = 16, type = "b", col = "red")
legend(x = "topleft", legend = expression(R^2, R[adj]^2),
       fill = c("black", "red"))

plot(0:ng, aic, xlab = "# of Garbage Predictors", ylab = "",
     main = "Criteria",pch = 16, type = "b", ylim = range(aic, bic))
points(0:ng, bic, pch = 16, type = "b", col = "blue")
legend(x = "topleft", legend = expression(AIC, BIC),
       fill = c("black", "blue"))






##---same thing in larger simulated sample---------------------------

library(mvtnorm)
nn <- 200 # n obs
n.cov <- 20 # 50 covs
n.cov.plot <- 20
X <- rmvnorm(nn,rep(0,n.cov),diag(rep(1,n.cov)))
## generate y as function of first 4 covs
y <- X[,1]+X[,2]+X[,3]+X[,4]+rnorm(nn,0,1)
dat <- data.frame(y,X)
# calculate R^2 and R^2_adj upon adding in garbage predictors one at a time
R2 <-  rep(NA,n.cov.plot-1)
R2a <- rep(NA, n.cov.plot-1)
aic <- rep(NA, n.cov.plot-1)
bic <- rep(NA, n.cov.plot-1)
mse <- rep(NA, n.cov.plot-1)
for(ii in 2:n.cov.plot) {
  M <- lm(y ~ ., data = dat[,1:ii])
  R2[ ii-1] <- summary(M)$r.squared
  R2a[ii-1] <- summary(M)$adj.r.squared
  aic[ii-1] <- AIC(M)
  bic[ii-1] <- BIC(M)
  mse[ii-1] <- sigma(M)  
}



plot(2:n.cov.plot, mse, xlab = "# of Predictors", ylab = "",
     main = "MSE",
     pch = 16, type = "b")

plot(2:n.cov.plot, R2, xlab = "# of Predictors", ylab = "",
     main = "Coefficient of Determination",
     pch = 16, type = "b", ylim = range(R2, R2a))
points(2:n.cov.plot, R2a, pch = 16, type = "b", col = "red")
legend(x = "bottomright", legend = expression(R^2, R[adj]^2),
       fill = c("black", "red"))

plot(2:n.cov.plot, aic, xlab = "# of Predictors", ylab = "",
     main = "Criteria",pch = 16, type = "b", ylim = range(aic, bic))
points(2:n.cov.plot, bic, pch = 16, type = "b", col = "blue")
legend(x = "bottomright", legend = expression(AIC, BIC),
       fill = c("black", "blue"))





#--- model selection: overfitting ------------------------------

## seed <- sample(1000, 1)
par(mfrow = c(1,1), mar = c(4,4,1,1))
cex <- 1
lwd <- 2
clrs <- c("blue", "red")
seed <- 704
set.seed(seed)
n <- 20
sig <- 2
x <- runif(n, -3, 3)
y <- exp(x) + sig * rnorm(n)
M1 <- lm(y ~ x)
form2 <- paste0("I(x^", 2:12, ")", collapse = " + ")
## M2 <- lm(y ~ poly(x, degree = 12, raw = TRUE) - 1)
## form2 <- paste0("sin(", 1:8, "*x)", collapse = " + ")
form2 <- as.formula(paste0("y ~ ", form2, " - 1"))
M2 <- lm(form2)
xlim <- c(-3.1, 3.1)
ylim <- range(y)
ylim <- (ylim - mean(ylim) * 1.5) + mean(ylim)
xseq <- seq(xlim[1], xlim[2], len = 1000)
ypred <- predict(M2, newdata = data.frame(x = xseq))
plot(x, y, xlim = xlim, ylim = ylim, pch = 16,
     xlab = "", ylab = "", cex.axis = cex, cex = cex)
title(xlab = "x", ylab = "y", cex.lab = cex, line = 2.5)
abline(reg = M1, col = "blue", lwd = lwd)
lines(xseq, ypred, col = "red", lwd = lwd)
lgd <- expression(M[1]:E*group("[",y*" | "*x,"]")==beta[0]+beta[1]*x,
                  M[2]:E*group("[",y*" | "*x,"]")==beta[2]*x^2+cdots+beta[12]*x^{12})
legend("bottomright", legend = lgd,
       fill = clrs, title = "Model")





# ------ automatic selection ------------------


# forward, backward, and stepwise selection

# high school and beyond dataset
hsb <- read.csv("HSB2.csv")
hsbm <- hsb[,c(13, 1:6, 10, 7:9)] # math response only

summary(hsbm)

# bounds for model selection
M0 <- lm(math ~ 1, data = hsbm) # minimal model

# maximal model: all main effect and interactions, and a couple non-linearities
Mmax <- lm(math ~ (.)^2 +
             I(locus^2) + I(concept^2) + I(mot^2), data = hsbm)
beta.max <- coef(Mmax)
length(beta.max) # number of coefficients
names(beta.max)[is.na(beta.max)] # coefficients that couldn't be estimated

table(hsbm[c("lang", "minor")]) # number of observations in each category


# full model:
# 1. remove all interactions with career, but leave as main effect
# 2. remove the interaction between lang and minor
# 3. add nonlinear effects for continuous variables locus, concept, mot
Mfull <- lm(math ~ (.-career)^2 + career - lang:minor +    
              I(locus^2) + I(concept^2) + I(mot^2),
            data= hsbm)
length(coef(Mfull))
anyNA(coef(Mfull))

df.penalty <- 2 # this is the k penalty
## 2 corresponds to AIC
## log(n) corresponds to BIC

help(system.time)

# forward
system.time({
  Mfwd <- step(object = M0, # base model
               scope = list(lower = M0, upper = Mfull), # smallest and largest model
               direction = "forward",
               trace = 1, # trace prints out information
               k = df.penalty
  )
})

Mfwd
Mfwd$coefficients

Mfwd$model

str(Mfwd)

# backward
system.time({
  Mback <- step(object = Mfull, # base model
                scope = list(lower = M0, upper = Mfull),
                direction = "backward", trace = 0, k = df.penalty)
})

# stepwise (both directions)
Mstart <- lm(math ~ ., data = hsbm) # starting point model: main effects only
system.time({
  Mstep <- step(object = Mstart,
                scope = list(lower = M0, upper = Mfull),
                direction = "both", trace = 1, k = df.penalty)
})

# compare the three different models
beta.fwd <- coef(Mfwd)
beta.back <- coef(Mback)
beta.step <- coef(Mstep)
c(fwd = length(beta.fwd), back = length(beta.back),
  step = length(beta.step)) # number of coefficients in each
# check if models are nested
names(beta.fwd)[!names(beta.fwd) %in% names(beta.back)]
names(beta.back)[!names(beta.back) %in% names(beta.fwd)]

# Mfwd and Mback different, but also each has some covariates the other does not.



## AIC
n <- nrow(hsbm)
ll_fwd <- -n/2 * (1 + log(sum(resid(Mfwd)^2)/n) + log(2*pi))
aic_fwd <- -2*ll_fwd + 2*(n - Mfwd$df + 1) # total number of parameters includes sigma
aic_fwd - AIC(Mfwd)
aic_step <- AIC(Mstep)
aic_back <- AIC(Mback)

aic_all <- round(c(aic_fwd,aic_step,aic_back),2)
names(aic_all) <- c("FWD","Step","Back")
aic_all

## BIC
bic_fwd <- -2*ll_fwd + log(n)*(n - Mfwd$df + 1) # total number of parameters includes sigma
bic_fwd - BIC(Mfwd)
bic_step <- BIC(Mstep)
bic_back <- BIC(Mback)

bic_all <- round(c(bic_fwd,bic_step,bic_back),2)
names(bic_all) <- c("FWD","Step","Back")
bic_all


