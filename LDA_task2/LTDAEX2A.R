#install.packages("KMsurv")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("bshazard")
#install.packages("km.ci")
library(KMsurv)
library(ggplot2)
library(dplyr)
library(survival)
library(bshazard)
library(km.ci)
data(kidney)
kidney1 <- kidney[kidney$type == 1, ]
kidney2 <- kidney[kidney$type == 2, ]
#Draw the (estimated) survival functions in both study groups
#and interpret the resulting graph.
# Kaplan Meier Survival Curve
km1 <- with(kidney1, Surv(time, delta))
km2 <- with(kidney2, Surv(time, delta))
head(km1)
head(km2)
km_fit1 <- survfit(Surv(time, delta) ~ 1, data=kidney1)
km_fit2 <- survfit(Surv(time, delta) ~ 1, data=kidney2)
summary(km_fit1)
summary(km_fit2)
plot(km_fit1, xlab="Months", main = 'Kaplan Meyer Survival function Group One') 
#base graphics is always ready
plot(km_fit2, xlab="Months", main = 'Kaplan Meyer Survival function Group two')
#Estimate the median time until renal infection in both groups
median_times <- survfit(Surv(time, delta) ~ type, kidney)
#draw smoothed estimates of the hazard functions in both study groups
smoothed_estimate1 <- bshazard(Surv(time, delta) ~ 1, data=kidney1)
plot(smoothed_estimate1, xlab="Months", main = 'Smooth hazard functions treatment 1')
smoothed_estimate2 <- bshazard(Surv(time, delta) ~ 1, data=kidney2)
plot(smoothed_estimate2, xlab="Months", main = 'Smooth hazard functions treatment 2')
#Which are the lower and upper limits
#of the linear EP confidence bands with level 95%
#in both study groups after 6, 12, 18, 24, and 30 months.
svforall <- survfit(Surv(time, delta) ~ type, kidney)
# Linear EP confidence bands
summary(km.ci(svforall[1], method='linear'), times=c(6,12,18,24,30))
summary(km.ci(svforall[2], method='linear'), times=c(6,12,18,24,30))
#Draw the number of people at risk (of dying) as a function of age.
Elderly <- read.csv("/Volumes/GoogleDrive/My Drive/Assignment 2 Exercises and data sets-20211117/Elderly.txt", sep="")
#Lab3
## Number of individuals at risk over time (Slide 90/101)
## ======================================================
attach(Elderly)
Times <- seq(min(Elderly$entry), max(Elderly$exit))
Nrisk <- numeric(length(Times))
for (i in 0:length(Times)) {
  Nrisk[i] <- with(Elderly, sum(entry <= Times[i] & exit >= Times[i]))
}
plot(Times+65, Nrisk, type = "S", lwd = 3, xlab = "Age [Years]",
     ylab = "Individuals at risk", bty = "n", col = 2,
     xlim = c(65, 110), ylim = c(0, 75), yaxs = "i", xaxs = "i")

axis(1, at = seq(65, 110, 5))
axis(2, at = seq(0, 75, 10))
title("Number of people at risk")
#Draw the conditional survival functions for men aged 70 and 85 years,
#respectively. Which are the estimated probabilities
#of surviving 90 years in both cases?
## Conditional survival functions (Slide 93/101)
## ---------------------------------------------
# 1st step: Subsets of those people exactly in 70 and 85 years, respectively.
# (at study exit)
sch <- survfit(Surv(entry, exit, Death) ~ 1)
summary(sch)
#Which are the corresponding estimated probabilities of surviving 90 years,
#if left truncation was ignored? Comment on the results.
ignored <- survfit(Surv(exit, Death) ~ 1)
summary(ignored)
# 2nd step: Graph
plot(sch, col=2 ,lty = 2, xlab ="Years",
     ylab = expression(bolditalic(hat(S)(t))))
lines(ignored, col = 3, lty = 1)
legend("bottomleft",col = 2:3, lty = c(2, 1), lwd = 2, bty = "n",
       legend = paste0(c(": Considering ", ": Ignoring ")))
#Draw the survival functions of time until turnover for both women and men.
#What do you observe?
# $>$ Hmisc::Label(turnover)
# label(stag) $<-$ 'Time in company until turnover or end of study [months]'
# label(event) $<$ 'Turnover indicator'
# label(gender) $<-$ 'Gender'
# label(headgend) $<-$ 'Gender of the supervisor'
load("Assign2Exer4.RData")
head(turnover)
attach(turnover)
male <- turnover[gender=="Male",]
female <- turnover[gender=="Female",]
unturnmale <- survfit(Surv(stag, event) ~ 1, data = male)
unturnfemale <- survfit(Surv(stag, event) ~ 1, data = female)
summary(unturnmale)
summary(unturnfemale)
unturn <- survfit(Surv(stag, event) ~ gender)
plot(unturn, col = 3:2, lwd = 3, xlab = "Time",
     ylab = expression(bolditalic(hat(S)(t))))
#lines(unturnfemale[1], col = 2)
legend("topright", c("Women", "Men"), col = 2:3, lwd = 3, bty = "n")
#Draw the survival functions separately for both genders of the supervisor.
#What do you observe?
unturnhead <- survfit(Surv(stag, event) ~ headgend)
plot(unturnhead, col = 3:2, lwd = 3, xlab = "Time",
     ylab = expression(bolditalic(hat(S)(t))))
#lines(unturnfemale[1], col = 2)
legend("topright", c("Women supervisor", "Men supervisor"), col = 2:3, lwd = 3, bty = "n")
#Test the hypothesis that time to turnover does not depend
#on the employee’s gender using any
#test of the Fleming Harrington family of tests. 
#What do you conclude?




