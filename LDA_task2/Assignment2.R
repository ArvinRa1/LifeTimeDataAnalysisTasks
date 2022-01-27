#EXERCISE 2

install.packages("KMsurv")
library(KMsurv)
data(kidney)
kidney

#EXERCISE 3

setwd("C:/Users/johan/Google Drive/UNI/UPC/Lifetime Data Analysis/Labs")
Elderly=read.table("Elderly.txt", skip=8)
summary(Elderly)

names(Elderly)[1] <- "entry"
names(Elderly)[2] <- "exit"
names(Elderly)[3] <- "death"


Elderly["entry"]=65+Elderly["entry"]
Elderly["exit"]=65+Elderly["exit"]
attach(Elderly)


#3b
Ages <- seq(min(Elderly$entry), max(Elderly$exit), 1)
Nrisk <- numeric(length(Ages))
for (i in 1:length(Ages)) {
  Nrisk[i] <- with(Elderly, sum(Elderly$entry <= Ages[i] & Elderly$exit >= Ages[i]))
}

par(font = 2, font.axis = 2, font.lab = 4, las = 1)
plot(Ages, Nrisk, type = "S", lwd = 3, xlab = "Age [Years]",
     ylab = "Individuals at risk", bty = "n", col = 2,xlim = c(65, 105), ylim = c(0, 50),
     yaxs = "i", xaxs = "i")
axis(1, at = seq(65, 110, 5))
axis(2, at = seq(0, 50, 5))
title("Elderly data: number of persons at risk")
#3c

sub70 <- subset(Elderly, entry > 70)
sub85 <- subset(Elderly, exit > 85)
sel70 <- survfit(Surv(entry,exit, death)~1, sub70)
sel85 <- survfit(Surv(entry, exit, death)~1, sub85)

# Plotting Conditional Survival Function
par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sel70, lwd = 2, col = 3, lty = 1, xlab = "Age [Years]",
     xlim = c(70, 102), ylab = expression(bolditalic(hat(S)(t))))
lines(sel85, lwd = 2, col = 2, lty = 1)
title("Conditional survival functions of people older than 70 and 85 years")
legend("bottomleft", c(" > 85", " > 70"),col = 2:3, lty = 1, lwd = 3, bty = "n")

## prob to survive 90 years: 
## summary(sel70): 0.1103  
## summary(sel85): 0.444   

#3d
sel70i <- survfit(Surv(exit, death)~1, sub70)
sel85i <- survfit(Surv(exit, death)~1, sub85)


par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(sel70i, lwd = 2, col = 3, lty = 1, xlab = "Years",
     xlim = c(70, 100), ylab = expression(bolditalic(hat(S)(t))))
lines(sel85i, lwd = 2, col = 2, lty = 1)
axis(1, at = seq(70, 100, 5))
title("Conditional survival functions of people older than 68 years\n Ignoring Truncation")
legend("bottomleft", c(" > 85", " > 70"),col = 2:3, lty = 1, lwd = 3, bty = "n")

## prob to survive 90 years: 
## summary(sel70i): 0.301    
## summary(sel85i): 0.625   

#EXERCISE 4
load("Assign2Exer4.RData")
attach(turnover)
summary(turnover)

#turnover$gender <- factor(turnover$gender, labels = c("Allogeneic", "Autologous"))


#4a
stu=survfit(Surv(stag,event)~gender, turnover)
summary(stu)
par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(stu, lwd = 2, col = 2:3, lty = 1, xlab = "Time in Months",
     xlim = c(0, 180), ylab = expression(bolditalic(hat(S)(t))))
title("Survival function until turnover")
legend("bottomleft", c("Men", "Women"),col = 2:3, lty = 1, lwd = 3, bty = "n")


#4b
stuh=survfit(Surv(stag,event)~headgend, turnover)
summary(stuh)
par(las = 1, font = 2, font.axis = 2, font.lab = 4, xaxs = "i", yaxs = "i",
    mar = c(5, 5, 4, 2), bty = "l", cex.lab = 1.5, cex.axis = 1.25)
plot(stuh, lwd = 2, col = 2:3, lty = 1, xlab = "Time in Months",
     xlim = c(0, 180), ylab = expression(bolditalic(hat(S)(t))))
title("Survival function until turnover\n with respect to the supervisor's gender")
legend("bottomleft", c("Men", "Women"),col = 2:3, lty = 1, lwd = 3, bty = "n")

#4c
#install.packages("FHtest")
library(FHtest)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install("Icens")
stu <- with(kidney, Surv(stag, event) ~ gender)
survfit(stu)

survdiff(stu)
survdiff(stu, rho = 1)


FHtestrcc(stu)
FHtestrcc(stu, rho = 1)
FHtestrcc(stu, lambda = 1)
FHtestrcc(stu, rho = 1, lambda = 1)
FHtestrcc(stu, rho = 0.5, lambda = 2)


#4d


survdiff(Surv(stag, event) ~gender+strata(gender), turnover)
