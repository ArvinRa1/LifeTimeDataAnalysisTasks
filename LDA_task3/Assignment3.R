#setwd("C:/Users/johan/Google Drive/UNI/UPC/Lifetime Data Analysis/Labs")


##### IMPORT ######
#install.packages("rms")
library(survival)
library(KMsurv)
library(rms)


##### EXERCISE 1 #####
data(tongue)
tongue$type=factor(tongue$type, labels=c("aneuploid","diploid"))
summary(tongue)
attach(tongue)

### 1a
ston <- with(tongue, Surv(time, delta))
svf.type=npsurv(ston~type, tongue)
par(font = 2, font.axis = 2, font.lab = 4, las = 1)
survplot(svf.type, ylab = expression(bold(hat(S)(t))), col = 2:3, lwd = 3, xlab="Time to death [Weeks]")
title("Survival functions according to Tumor DNA profile")

### 1b
ston.type<- with(tongue, Surv(time, delta) ~ type)
survdiff(ston.type)  
#H0: no diff
#H1: yes diff
#p=0.09 -> not reject

### 1c
svf.loglo=survreg(ston ~ type, tongue, dist = "loglo")
summary(svf.loglo)
#H0: type is not impacting the surv time param=0 for type
#H1: type is incluencing surv time
#p=0.051 -> not reject

### 1d

### 1e
lnopred <- predict(lnomod3, type = "linear")
residsLN <- (log(larynx$time) - lnopred) / lnomod3$scale
residsLN

loglo.pred <- predict(svf.loglo, type = "linear")
summary(loglo.pred)
resids <- (log(tongue$time) - loglo.pred) /svf.loglo$scale
resids
par(font = 2, font.lab = 4, font.axis = 2, las = 1, oma = c(0, 0, 1, 0),
    mar = c(5, 5, 4, 2))
plot(survfit(Surv(resids, tongue$delta) ~ 1), xlab = "Years", lwd = 3,
     ylab = expression(bold(hat(S)(t))), yaxs = "i")
title("Residuals of the Loglogistic regression model")
survgumb <- function(x) {
  return(exp(-exp(x)))
}
curve(survgumb(x), from = min(resids), to = max(resids), col = 2, lwd = 3,
      add = TRUE)


#### EXERCISE 2 ####
### 2a
survlog=rlnorm(300,2,1)
### 2b
lam=1/20
censexp=rexp(300,lam)

### 2c
Y=pmin(survlog,censexp)
del=as.numeric(survlog<=censexp)
#proportion of censoring
table(del)
1-(sum(del)/length(del))
head(cbind(survlog, censexp, Y, del))

### 2d
#stype=2 - cum.hazard
#ctype=1 - Nelson-Aalen Formula
svf <- survfit(Surv(Y, del) ~ 1, stype = 2, ctype = 1)
svf
summary(svf)
## Uncensored survival times
times <- summary(svf)$time

## Nelson-Aalen estimate of the cumulative hazard function
chaz <- -log(summary(svf)$surv)

## The probability plots
par(mfrow = c(1, 2), las = 1, font.lab = 4, font.axis = 2, pch = 16)
plot(log(chaz) ~ log(times), xlab = "Log(Time)",
     ylab = "Log(Cumulative hazard)", col = 4)
abline(lm(log(chaz) ~ log(times)), lwd = 3)
title("Check for Weibull distribution")

plot(log(exp(chaz) - 1) ~ log(times), xlab = "Log(Time)",
     ylab = "Log(Exp(Cumulative hazard) - 1)", cex = 1.3, col = 4)
abline(lm(log(exp(chaz) - 1) ~ log(times)), lwd = 3)
title("Check for Log-logistic distribution")
#exclude Weibull as valid model, choose log-logistic

#### EXERCISE 4 ####

data(hodg)
### 4a
hodg$gtype=factor(hodg$gtype, labels=c("allogenic", "autologous"))
hodg$dtype=factor(hodg$dtype, labels=c("Non Hodgkin lymphoma", "Hodgkins disease"))
summary(hodg)
### 4b
?hodg
hodg$months <- hodg$time / 30
shodg <- with(hodg, Surv(months, delta))
summary(shodg)
plot(shodg)
# (b) Draw the survival functions corresponding to the four combinations of graft and disease types
# measuring time until relapse or death in years. Comment on the graph.
par(mfrow = c(1, 2), font = 2, font.lab = 4, font.axis = 2, las = 1, mar = c(3, 3, 4, 2))

plot(survfit(shodg ~ gtype, hodg), col = 1:2, xlab = "Months", lwd = 2, lty = 2:3,
     ylab = expression(bolditalic(hat(S)(t))), bty = "n")
legend("bottomleft", levels(hodg$gtype), col = 1:2, lwd = 2, bty = "n", lty = 2:3,
       title = "Graft type")
title("Survival functions with different types of Graft")

plot(survfit(shodg ~ dtype, hodg), col = 3:4, xlab = "Months", lwd = 2, lty = 1:2,
     ylab = expression(bolditalic(hat(S)(t))), bty = "n")
legend("bottomleft", levels(hodg$dtype), col = 3:4, lwd = 2, bty = "n", lty = 1:2,
       title = "disease type")
title("Survival functions with different types of disease")

# (c) Fit the proportional hazards model that includes graft type, disease type, the interaction of both,
# and the Karnofsky index. Interpret the model fit.
## The fit of a Cox model
## ======================
(coxh <- coxph(shodg ~ gtype + dtype, hodg))


## Including the interaction between both variables and the Karnofsky index
## -------------------------------------------------------------------------
coxh <- update(coxh,  ~ . + gtype:dtype + score)
summary(coxh)
coxh$var

# (d) Estimate the hazard ratios associated to the graft type (comparing autologous to allogenic 
# transplantations) and interpret both values.
library(Epi)
ci.lin(coxh)
round(ci.lin(coxh, Exp = TRUE), 3)
(ctmat <- matrix(c(2, 0, 0, 0, 1, 0, 0, 0), byrow = TRUE, nr = 2))
round(ci.lin(coxh, ctr.mat = ctmat, Exp = TRUE), 3)

# A somewhat nicer presentation
HRmat <- round(ci.lin(coxh, ctr.mat = ctmat, Exp = TRUE), 3)[, c(1, 5:7)]
rownames(HRmat) <- c("HR| allogenic", "HR| autologous")
colnames(HRmat) <- c("logHR", "HR", "Lower 95%", "Upper 95%")
HRmat

# (e) Check the proportional hazards assumption.
## Checking proportional hazard assumption
## -----------------------------------------
residuals(coxh, "schoenfeld")
## (i) Use of function cox.zph
chpa =  cox.zph(coxh)
## (ii) Use of function plot.cox.zph
#windows(width = 12, height = 7)
par(mfrow = c(1, 1), font = 2, font.lab = 4, font.axis = 2, las = 1,
    cex.lab = 1.3, cex.axis = 1.2)
plot(chpa[1], lwd = 2)
abline(h=0, col= 2)
plot(chpa[2], lwd = 2)
abline(h=0, col= 2)
plot(chpa[3], lwd = 2)
abline(h=0, col= 2)
plot(chpa[4], lwd = 2)
abline(h=0, col= 2)
# (f ) Concerning the estimation of the four model parameters, are there any in
# influential observations? 
dfbet <- residuals(coxh, type = "dfbeta")
dim(dfbet)
par(mfrow = c(2, 2), font = 2, font.lab = 4, font.axis = 2, las = 1,
    cex.lab = 1.3, cex.axis = 1.2)
for (i in 1:5) {
  plot(dfbet[, i], pch = 16, ylab = "")
  title(names(coef(coxh))[i])
  axis(1, at = seq(5, 45, 5))
}

