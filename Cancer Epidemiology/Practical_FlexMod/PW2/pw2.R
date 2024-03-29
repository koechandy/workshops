#*------------------------------------------------------*
#|                             |
#|                                                      |
#| source("C:/Users/Aurel/Dropbox/WorkStat/Teaching/Nairobi-2019_ISCB_ShortCourse/Material/Practical_FlexMod/PW2/pw2.R", echo=T)
#|                                                      |
#-------------------------------------------------------*
# 
# 		 ===========================
#   		|	Pratical work n�2   |
# 		 ===========================

setwd("C:/Users/Aurel/Dropbox/WorkStat/Teaching/Nairobi-2019_ISCB_ShortCourse/Material/Practical_FlexMod/")


sink("PW2/pw2.lis")

date()
system("hostname", intern=T)
library("survival")
library("mexhaz")


sessionInfo()

#------------------------------------------------------------------------------------------------------------------------
## 2.1 ==> Get the population hazard and plot hazard vs age for year=2000 and females&Males in log scale

load("data/popex1.RData")
summary(popex1)
dim(popex1)
head(popex1)
plot(popex1$age[popex1$sex=="Females"&popex1$year==2000], popex1$rate[popex1$sex=="Females"&popex1$year==2000], 
     xlab="Age", ylab="Population Hazard", log="y", type="l", col="red")
lines(popex1$age[popex1$sex=="Males"&popex1$year==2000], popex1$rate[popex1$sex=="Males"&popex1$year==2000], col="green")
legend("topleft",col=c("red","green"),lty=1,legend=c("Females","Males"))

#------------------------------------------------------------------------------------------------------------------------
## 2.2 ==> Add the population hazard at death in datacolon ===> save as datacolon.v2

load("data/datacolon.RData")


datacolon$fu=as.numeric((datacolon$finmdy-datacolon$diagmdy)/365.25)
datacolon$agedeath=trunc(datacolon$agediag+datacolon$fu)


datacolon$yeardeath=as.numeric(format(datacolon$finmdy, format = "%Y"))
datacolon$agedeath=ifelse(datacolon$agedeath>100, 100, datacolon$agedeath)


datacolon.v2=merge(datacolon, popex1[popex1$sex=="Males", c("age","year","rate")], 
           by.x=c("agedeath","yeardeath"), by.y=c("age","year")
          )

save(datacolon.v2,file = "data/datacolon.v2.RData")


#------------------------------------------------------------------------------------------------------------------------
## 2.3 ==> fit a excess hazard model where the baseline is piecewise constant and plot the hazard
#	   with intervas: "[0,0.25]"   "(0.25,0.5]" "(0.5,1]"    "(1,2]"      "(2,3]"      "(3,4]"      "(4,5]"
# 	   save the values of the hazard in /PW2/hazard23.Rdata

stage="St4"
temp=datacolon.v2[datacolon.v2$stage==stage,]
interval=c(0.25, 0.5, 1:5)
model23=mexhaz(Surv(time=fu,event=dead)~ 1, data=temp, base="pw.cst", knots=interval, expected="rate")

# Note the WARNING message: "Some observations had a follow-up time of length 0. ..."

# summary, coefficients and variance-covariance matrix
summary(model23)
model23$coef
round(model23$vcov,4)


# With mexhaz
# pred.model23=predict(model23, time.pts=c(0.0001, seq(0.1,5, by=0.0001)))
# note that the above "by" parameters must be low
# plot(pred.model23, which="hazard", conf.int=TRUE)


begin.interval=c(0, 0.25, 0.5, 1:4)
end.interval=c(0.25, 0.5, 1:5)
hazard23=exp(model23$coef[-match("(5,7.96]", names(model23$coef))])


plot(0,0, type="n", xlim=c(0,5), ylim=c(0,2),xlab="Time since diagnosis", ylab="hazard ")
segments(begin.interval, hazard23, end.interval, hazard23, lwd=3,col=3)
title( paste("\nExcess hazard, stage=",stage ), cex=0.6)

save(hazard23,file = "PW2/hazard23.RData"   )

#------------------------------------------------------------------------------------------------------------------------
## 2.4 ==> with model defined in 2.3, calculate the net survival at 1 year and 1.37 years; check the results with predict.mexhaz

# cumulative hazard at 1 year then NS
ch1=0.25*exp(model23$coef["[0,0.25]"])+ 0.25*exp(model23$coef["(0.25,0.5]"]) + 0.5*exp(model23$coef["(0.5,1]"])
exp(-ch1)


# cumulative hazard at 1.37 year then NS
ch137=0.25*exp(model23$coef["[0,0.25]"])+ 0.25*exp(model23$coef["(0.25,0.5]"]) + 0.5*exp(model23$coef["(0.5,1]"])+ 0.37*exp(model23$coef["(1,2]"])
exp(-ch137)

predict(model23, time.pts=c(1, 1.37) )$results


#------------------------------------------------------------------------------------------------------------------------------
## 2.5 ==> calculate the CI of these net survivals (at 1 year and 1.37 years), using a Wald-type CI for the cumulative hazard
# check the results with predict.mexhaz

# var(exp(beta))# [exp(beta)]^2 * var(beta)

var.ch1=0.25^2*exp(model23$coef["[0,0.25]"])^2*  model23$vcov["[0,0.25]","[0,0.25]"] +
        0.25^2*exp(model23$coef["(0.25,0.5]"])^2*model23$vcov["(0.25,0.5]","(0.25,0.5]"] +
        0.50^2*exp(model23$coef["(0.5,1]"])^2*model23$vcov["(0.5,1]","(0.5,1]"]

binf.ch1=ch1-qnorm(0.975)*sqrt(var.ch1)
bsup.ch1=ch1+qnorm(0.975)*sqrt(var.ch1)
c(exp(-bsup.ch1), exp(-binf.ch1))


var.ch137=0.25^2*exp(model23$coef["[0,0.25]"])^2*  model23$vcov["[0,0.25]","[0,0.25]"] +
        0.25^2*exp(model23$coef["(0.25,0.5]"])^2*model23$vcov["(0.25,0.5]","(0.25,0.5]"] +
        0.50^2*exp(model23$coef["(0.5,1]"])^2*model23$vcov["(0.5,1]","(0.5,1]"]+
        0.37^2*exp(model23$coef["(1,2]"])^2*model23$vcov["(1,2]","(1,2]"]

binf.ch137=ch137-qnorm(0.975)*sqrt(var.ch137)
bsup.ch137=ch137+qnorm(0.975)*sqrt(var.ch137)
c(exp(-bsup.ch137), exp(-binf.ch137))


predict(model23, time.pts=c(1, 1.37), delta.type.s="log" )$results

#------------------------------------------------------------------------------------------------------------------------------
##  2.6: 
# Using predict.mexhaz, compute CI for NS at the same times assuming normality on log-log scale (i.e. : on the log-cumulative scale)

predict(model23, time.pts=c(1, 1.37), delta.type.s="log-log" )$results


date()


sink()





























