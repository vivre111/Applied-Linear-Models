### LECTURE 9

#--- categorical predictors ------------------------------------------------
# data
mercury <- read.csv("fishermen_mercury.csv")
head(mercury)

# predict MeHg as a function of weight and fishpart
M <- lm(MeHg ~ weight + fishpart, data = mercury)
summary(M)

# since fishpart is coded as:
# None = 0, Muscle = 1, Muscle/Whole = 2, Whole = 3
# model assumes increase in expected MeHg from None to Muscle is
# exactly equal to increase in expected MeHg from Muscle to Muscle/Whole.
# Imagine if we encoded: None = 0, M = 10, MW = 11, W = 20?

# the "design matrix" X for the problem
head(model.matrix(M))

#--- fishpart as categorical predictor -------------------------------------

# convert fishpart to non-numeric categorical variable
fishpart2 <- mercury$fishpart
table(fishpart2)

# "levels" of variable fishpart
fplevels <- c("none", "muscle", "muscle_whole", "whole")
# using for-loop
for(ii in 0:3) {
  fishpart2[fishpart2 == ii] <- fplevels[ii+1]
}
# without for-loop
all(fplevels[mercury$fishpart+1] == fishpart2)


# in R, categorical variables are coded as "factor" variables
tmp <- factor(fishpart2)
head(tmp)
levels(tmp) # alphabetical group ordering

fishpart2 <- factor(fishpart2, levels = fplevels) # custom ordering
levels(fishpart2)

mercury$fishpart <- fishpart2 # replace in dataset

# run regression in R:
Mb <- lm(MeHg ~ weight + fishpart, data = mercury) # intercept model (I)
Mc <- lm(MeHg ~ weight + fishpart - 1,         # no-intercept model (0)
         data = mercury)          # note the "-1" to remove intercept



# compare design matrices at selected values
ind <- c(102, 35, 2, 1, 14, 106)

mercury[ind, c("MeHg", "weight", "fishpart")] # original dataset
model.matrix(Mb)[ind,] # no intercept model
model.matrix(Mc)[ind,] # with intercept model

# interpretation of parameters:
coef(Mb)
# fishpartmuscle: estimate of E[MeHg | ...]
pbinom(0,100,2/100)

coef(Mc)
# fishpartmuscle: estimate of ...

# note that model predictions are identical:
range(predict(Mb) - predict(Mc))




# confidence interval for E[MeHg | x = xstar]
predict(Mb,
        newdata = data.frame(weight = c(70, 60),
                             fishpart = c("whole", "none")),
        interval = "confidence", level = .95)


# prediction interval for MeHg.star | x = xstar
predict(Mb, newdata = data.frame(weight = 70, fishpart = "whole"),
        interval = "prediction", level = .95)





#--- Hypothesis testing ----------------------------------------------------

# Question: How to test H0: fishpart not associated with mercury in presence of weight?
Mfull <- lm(MeHg ~ fishpart + weight, data = mercury) # full model
coef(Mfull)
# degrees of freedom 
dff <- Mfull$df

# F-statistic using subset of MLE of full model (betastar)
betastar.ind <- 2:4 # indices corresponding to beta's for fishpart
betastar.hat <- coef(Mfull)[betastar.ind] # subset of MLE
# vcov(Mfull) = sigma.hat^2 * (X'X)^{-1}
help(vcov)
betastar.ve <- vcov(Mfull)[betastar.ind,betastar.ind] # betastar.ve = sigma.hat^2 * V
# long hand calculation
X <- model.matrix(Mfull) # design matrix for full model

solve(t(X)%*%X)*  sum(Mfull$residuals^2) / Mfull$df


V <- solve(crossprod(X))[betastar.ind,betastar.ind] # crossprod(X) = t(X) %*% X
range(betastar.ve - sigma(Mfull)^2 * V)
# Fobs = (betastar.hat' V^{-1} betastar.hat / q) / sigma.hat^2
# note: betastar.ve already includes the denominator
Fobs <- t(betastar.hat) %*% solve(betastar.ve, betastar.hat) / length(betastar.hat)

# p-value: large values of Fobs are more evidence against H0
pf(Fobs, df1 = length(betastar.hat), df2 = dff, lower.tail = FALSE)











### LECTURE 10
Mb <- lm(MeHg ~ weight + fishpart, data = mercury) ## common slopes, different intercepts model

wgrid <- seq(min(mercury$weight),max(mercury$weight),length.out = 100)

Mbpred_n <- predict(Mb,newdata = data.frame(weight=wgrid,fishpart="none"))
Mbpred_m <- predict(Mb,newdata = data.frame(weight=wgrid,fishpart="muscle"))
Mbpred_mw <- predict(Mb,newdata = data.frame(weight=wgrid,fishpart="muscle_whole"))
Mbpred_w <- predict(Mb,newdata = data.frame(weight=wgrid,fishpart="whole"))

plot(mercury$MeHg~mercury$weight,
     ylab="Mean MeHg", ## axis labels
     xlab="Weight",
     pch=16) # dots

lines(Mbpred_n~wgrid,col="black")
lines(Mbpred_m~wgrid,col="blue")
lines(Mbpred_mw~wgrid,col="red")
lines(Mbpred_w~wgrid,col="gray")


##
Md <- lm(MeHg ~ weight*fishpart, data = mercury) ## different slopes, different intercepts model

Mdpred_n <- predict(Md,newdata = data.frame(weight=wgrid,fishpart="none"))
Mdpred_m <- predict(Md,newdata = data.frame(weight=wgrid,fishpart="muscle"))
Mdpred_mw <- predict(Md,newdata = data.frame(weight=wgrid,fishpart="muscle_whole"))
Mdpred_w <- predict(Md,newdata = data.frame(weight=wgrid,fishpart="whole"))

plot(mercury$MeHg~mercury$weight,
     ylab="Mean MeHg", ## axis labels
     xlab="Weight",
     pch=16) # dots

lines(Mdpred_n~wgrid,col="black")
lines(Mdpred_m~wgrid,col="blue")
lines(Mdpred_mw~wgrid,col="red")
lines(Mdpred_w~wgrid,col="gray")





#--- interactions  --------------------------------------------------

# load data
real <- read.csv("real_estate.txt")
real$SalePrice <- real$SalePrice/1000 # price in 1000$ of dollars
head(real)

# plot Sale Price (SalePrice) vs Floor Area (SqrFeet)
clrs <- c("blue", "red")
mrkr <- c(16,17)
cex <- 1
plot(y = real$SalePrice, x = real$SqrFeet,
     pch = mrkr[2-real$Air], col = clrs[2-real$Air], cex = cex,
     xlab = "Floor Area (square feet)", ylab = "Sale Price (in 1000$)")
legend(x = "topleft",
       legend = c("With AC", "Without AC"),
       col = c(clrs), lwd = c(2, 2), seg.len = 1)

# Fit Sale Price vs Floor Area, with common slope but different
# intercepts for AC
M1 <- lm(SalePrice ~ SqrFeet + Air, data = real)

# add lines to plot
beta1 <- coef(M1)
abline(a = beta1["(Intercept)"],
       b = beta1["SqrFeet"], col = clrs[2], lwd = 2) # no AC
abline(a = beta1["(Intercept)"] + beta1["Air"],
       b = beta1["SqrFeet"], col = clrs[1], lwd = 2) # with AC
legend(x = "topleft",
       legend = c("With AC", "Without AC"),
       col = c(clrs), lwd = c(2, 2), seg.len = 1)


# Fit Sale Price vs Floor Area, with common intercept but different
# slopes for AC
M2 <- lm(SalePrice ~ SqrFeet + SqrFeet:Air, data = real)

# Fit Sale Price vs Floor Area, with different intercept and slope for AC
M3 <- lm(SalePrice ~ Air + SqrFeet + SqrFeet:Air, data = real)
M3 <- lm(SalePrice ~ Air*SqrFeet, data = real) # shorthand

# add lines to plot
clrs <- c("blue", "red")
mrkr <- c(16,17)
plot(y = real$SalePrice, x = real$SqrFeet,
     pch = mrkr[2-real$Air], col = clrs[2-real$Air],
     xlab = "Floor Area (square feet)", ylab = "Sale Price (in 1000$)")
# no interactions
beta1 <- coef(M1)
abline(a = beta1["(Intercept)"],
       b = beta1["SqrFeet"], col = clrs[2], lwd = 2) # no AC
abline(a = beta1["(Intercept)"] + beta1["Air"],
       b = beta1["SqrFeet"], col = clrs[1], lwd = 2) # with AC
# with interactions
beta2 <- coef(M3)
abline(a = beta2["(Intercept)"],
       b = beta2["SqrFeet"], lty = 3, col = clrs[2], lwd = 2) # no AC
abline(a = beta2["(Intercept)"] + beta2["Air"],
       b = beta2["SqrFeet"] + beta2["Air:SqrFeet"],
       lty = 3, col = clrs[1], lwd = 2) # with AC
legend(x = "topleft",
       legend = c("With AC", "Without AC", "No Interaction", "With Interaction"),
       col = c(clrs, "black", "black"), pch = c(mrkr, NA, NA),
       lty = c(NA, NA, 1, 3), lwd = c(NA, NA, 2, 2))







#--- non-linear terms --------------------------------------------------
tortoise <- read.csv("tortoise.csv")  # load data 
head(tortoise)

# fit a linear model in CarapDiam
M1 <- lm(NumEggs ~ CarapDiam, data = tortoise) 
summary(M1)

# fit a quadratic model in CarapDiam
M2 <- lm(NumEggs ~ CarapDiam + I(CarapDiam^2), data = tortoise)
summary(M2)


# plots
cgrid <- seq(min(tortoise$CarapDiam),max(tortoise$CarapDiam),length.out = 100)
pred1 <- predict(M1,newdata = data.frame(CarapDiam=cgrid))
pred2 <- predict(M2,newdata = data.frame(CarapDiam=cgrid))

plot(tortoise$NumEggs~tortoise$CarapDiam,
     ylab="Number ofEggs", ## axis labels
     xlab="Diameter",
     pch=16) # dots

lines(pred1~cgrid,col="black")
lines(pred2~cgrid,col="blue")


