#*------------------------------------------------------*
#|                             |
#|                                                      |
#| source("C:/Users/Aurel/Dropbox/WorkStat/Teaching/Nairobi-2019_ISCB_ShortCourse/Material/Practical_FlexMod/PW4/pw4.R", echo=T)
#|                                                      |
#-------------------------------------------------------*
# 
# 		 ===========================
#   		|	Pratical work n�4   |
# 		 ===========================


setwd("C:/Users/Aurel/Dropbox/WorkStat/Teaching/Nairobi-2019_ISCB_ShortCourse/Material/Practical_FlexMod/")


sink("PW4/pw4.lis")

date()
system("hostname", intern=T)
library(survival)
library(mexhaz)
library(splines)


sessionInfo()

#----------------------------------------------------------------------------------------------------------------------------------------
## 4.1 ==> Using datacolon.v2, stage 4, age [30;90]  fit an excess hazard model with
#	-  a cubic B-spline (one knot at 1 year) as baseline
#	-  a PH effect for covariable age in 5 groups: 
#		"[30;45[", "[45;55[", "[55;65[","[65;75[", [75;90]" (reference="[65;75[")


load("data/datacolon.v2.RData")

stage="St4"

rm(temp)
temp=datacolon.v2[datacolon.v2$agediag>=30&datacolon.v2$agediag<=90&datacolon.v2$stage==stage,]


temp$agegrp <- rep("[30;45[")
temp$agegrp <- ifelse(temp$agediag>=45 & temp$agediag<55, "[45;55[", temp$agegrp)
temp$agegrp <- ifelse(temp$agediag>=55 & temp$agediag<65, "[55;65[", temp$agegrp)
temp$agegrp <- ifelse(temp$agediag>=65 & temp$agediag<75, "[65;75[", temp$agegrp)
temp$agegrp <- ifelse(temp$agediag>=75 & temp$agediag<=90, "[75;90]", temp$agegrp)

# check the code:
tapply(temp$agediag, temp$agegrp, range) 

temp$Iagegrp3045 <- ifelse(temp$agegrp=="[30;45[",1,0) 
temp$Iagegrp4555 <- ifelse(temp$agegrp=="[45;55[",1,0) 
temp$Iagegrp5565 <- ifelse(temp$agegrp=="[55;65[",1,0) 
temp$Iagegrp6575 <- ifelse(temp$agegrp=="[65;75[",1,0) 
temp$Iagegrp7590 <- ifelse(temp$agegrp=="[75;90]",1,0) 


# [65;75[ is the reference
#
model41=mexhaz(Surv(time=fu,event=dead)~Iagegrp3045+Iagegrp4555+Iagegrp5565+Iagegrp7590 , 
               data=temp, base="exp.bs", knots=c(1), bound=c(0,8), expected="rate")


#----------------------------------------------------------------------------------------------------------------------------------------
## 4.2 ==> Compute the hazard at time 0:5 for the 5 age groups; check the results for age group [30;45[ and [65;75[ with predict
# for help type "? predict.mexhaz"

library(splines)
formula41=as.formula("~bs(fu, k=1, B=c(0,8)) + Iagegrp3045+Iagegrp4555+Iagegrp5565+Iagegrp7590" )
beta41=as.numeric(model41$coef)

my.data41=data.frame(fu=0:5)
my.data41=merge(my.data41, data.frame(agegrp=c("[30;45[", "[45;55[", "[55;65[", "[65;75[", "[75;90]")   ))
my.data41$Iagegrp3045 <- ifelse(my.data41$agegrp=="[30;45[",1,0) 
my.data41$Iagegrp4555 <- ifelse(my.data41$agegrp=="[45;55[",1,0) 
my.data41$Iagegrp5565 <- ifelse(my.data41$agegrp=="[55;65[",1,0) 
my.data41$Iagegrp6575 <- ifelse(my.data41$agegrp=="[65;75[",1,0) 
my.data41$Iagegrp7590 <- ifelse(my.data41$agegrp=="[75;90]",1,0) 


my.mat41=model.matrix(formula41, data=my.data41)
my.data41$hazard=exp(my.mat41%*%beta41)
head(my.data41)
my.data41[my.data41$agegrp=="[30;45[",]
my.data41[my.data41$agegrp=="[65;75[",]


predict(model41, time.pts=c(0.0000001, 1:5), 
               data.val = data.frame(Iagegrp3045=1,Iagegrp4555=0,Iagegrp5565=0,Iagegrp7585=0,Iagegrp7590=0) 
                )

predict(model41, time.pts=c(0.0000001, 1:5), 
               data.val = data.frame(Iagegrp3045=0,Iagegrp4555=0,Iagegrp5565=0,Iagegrp7585=0,Iagegrp7590=0) 
                )



#----------------------------------------------------------------------------------------------------------------------------------------
## 4.3 ==> Compute the HR for each age groups (ref=[65;75[) ; given the model, show that HR do not depend on fu (PH assumption)
#		check that HR correspond to exp(model41$coef) for the corresponding age group

# get the hazards of the reference group
#
ref=my.data41[my.data41$agegrp=="[65;75[",c("fu","hazard")]
ref$hz.ref=ref$hazard
my.data41=merge(my.data41, ref[,c("fu","hz.ref")])
my.data41$HR=my.data41$hazard/my.data41$hz.ref
my.data41=my.data41[order(my.data41$agegrp, my.data41$fu), c("fu","agegrp","hazard","hz.ref","HR")]

my.data41

my.data41[my.data41$agegrp=="[30;45[",]
exp(model41$coef)


#----------------------------------------------------------------------------------------------------------------------------------------
## 4.4 ==> plot the HR vs agediag

my.data41$begin.interval=as.numeric(substring(my.data41$agegrp,2,3))
my.data41$end.interval=  as.numeric(substring(my.data41$agegrp,5,6))

plot(0,0, type="n", xlim=c(30,90), ylim=c(0.5,2),xlab="Age at diagnosis", ylab="HR (ref=[65;75[) ")
segments(my.data41$begin.interval, my.data41$HR, my.data41$end.interval, my.data41$HR, lwd=3,col=3)
title("Excess Hazard Ratio for each age group")


#----------------------------------------------------------------------------------------------------------------------------------------
## 4.5 ==> model the age effect with a cubic spline with a knot at 70 years in truncated powers basis, centered and reduced agediag ie (agediag-70)/10 
#      ==> show that the log-likehood and predicted are the same with bs basis (using bs(agediag,k=70))
# note that 70 # median(temp$agediag[temp$dead==1])


temp$agecr=(temp$agediag-70)/10
temp$agecr2=temp$agecr^2
temp$agecr3=temp$agecr^3
posk=0
temp$agecr3plus=(temp$agecr-posk)^3*(temp$agecr>posk)

model45=mexhaz(Surv(time=fu,event=dead)~ agecr+ agecr2 +agecr3+ agecr3plus, 
               data=temp, base="exp.bs", knots=c(1), bound=c(0,8), expected="rate")
summary(model45)


model45.bs=mexhaz(Surv(time=fu,event=dead)~ bs(agediag, k=70, B=c(20,120)), data=temp, base="exp.bs", knots=c(1), bound=c(0,8), expected="rate")
summary(model45.bs)

c(model45$loglik, model45.bs$loglik)

# predicted for an age like the age of idenpat=8 for example
predict(model45, time.pts=1:5, data.val =temp[temp$idenpat==8,] )
predict(model45.bs, time.pts=1:5, data.val =temp[temp$idenpat==8,] )

#----------------------------------------------------------------------------------------------------------------------------------------
## 4.6 ==> add to the current graph the HR from model45


beta45.HR=as.numeric(model45$coef[c("agecr","agecr2","agecr3","agecr3plus")])


my.data45=data.frame(age=30:100)
my.data45$agecr=(my.data45$age-70)/10
my.data45$agecr2=my.data45$agecr^2
my.data45$agecr3=my.data45$agecr^3
my.data45$agecr3plus=(my.data45$agecr-posk)^3*(my.data45$agecr>posk)

my.mat45=as.matrix(my.data45[,c("agecr","agecr2","agecr3","agecr3plus")])
my.data45$HR=exp(my.mat45%*%beta45.HR)
lines(my.data45$age, my.data45$HR, col="blue")


#----------------------------------------------------------------------------------------------------------------------------------------
## 4.7 ==> add a non-proportional effect of age to model45 and 
#           predicted the hazard(t,age) for each combinaison of t=seq(0,5,by=0.1), age=30:90


model47=mexhaz(Surv(time=fu,event=dead)~ agecr+ agecr2 +agecr3+ agecr3plus + nph(agecr), 
               data=temp, base="exp.bs", knots=c(1), bound=c(0,8), expected="rate" )
summary(model47)


formula47=as.formula("~bs(fu, k=1, B=c(0,8)) + agecr+ agecr2 +agecr3+ agecr3plus + agecr:(bs(fu, k=1, B=c(0,8)))  " )
beta47=as.numeric(model47$coef)

my.data47=data.frame(fu=seq(0,5,by=0.1))
my.data47=merge(my.data47, data.frame(age=30:90)   )
my.data47$agecr=(my.data47$age-70)/10
my.data47$agecr2=my.data47$agecr^2
my.data47$agecr3=my.data47$agecr^3
my.data47$agecr3plus=(my.data47$agecr-posk)^3*(my.data47$agecr>posk)


# my.mat47 is the design matrix
#
my.mat47=model.matrix(formula47, data=my.data47)
my.data47$hazard=exp(my.mat47%*%beta47)
head(my.data47)

# Check that the "by-hand" design matrix is correct 
names(model47$coef)
dimnames(my.mat47)[[2]]




#----------------------------------------------------------------------------------------------------------------------------------------
## 4.8 ==> from model47 make figures:
#	- hazard vs t for age=30,60,90 ==> fig1
#	- hazard vs age for t=0.2, 0.5 1, 5 ==> fig2
#	- a 3D plot hazard(t,a) vs (t,a) ==> fig3

par(mfrow=c(1,2))

# Fig1
#
plot(0,0, type="n", xlim=c(0,5), ylim=c(0,2.5),xlab="time since diagnosis", ylab="excess mortality hazard ")
lines( my.data47$fu[my.data47$age==30], my.data47$hazard[my.data47$age==30], col=1)
lines( my.data47$fu[my.data47$age==60], my.data47$hazard[my.data47$age==60], col=2)
lines( my.data47$fu[my.data47$age==90], my.data47$hazard[my.data47$age==90], col=3)
axis(1, label=F, at=c(0.2, 0.5, 1, 5), tck=1,lty=8, lwd=0.1)
axis(2, label=F, at=seq(0,2.5, by=0.25), tck=1,lty=8, lwd=0.1)
legend(2,2, c("age=30","age=60","age=90"), text.col=1:3)
title("hazard vs time for given age")



# Fig2
#
plot(0,0, type="n", xlim=c(30,90), ylim=c(0,2.5),xlab="age at diagnosis", ylab="excess mortality hazard ")
lines( my.data47$age[my.data47$fu==0.2], my.data47$hazard[my.data47$fu==0.2], col=1)
lines( my.data47$age[my.data47$fu==0.5], my.data47$hazard[my.data47$fu==0.5], col=2)
lines( my.data47$age[my.data47$fu==1], my.data47$hazard[my.data47$fu==1], col=3)
lines( my.data47$age[my.data47$fu==5], my.data47$hazard[my.data47$fu==5], col=4)
axis(1, label=F, at=c(30,60,90), tck=1,lty=8, lwd=0.1)
axis(2, label=F, at=seq(0,2.5, by=0.25), tck=1,lty=8, lwd=0.1)
legend(40,2.5, c("t=0.2","t=0.5","t=1","t=5"), text.col=1:4)
title("hazard vs age for given time")


# Fig3
#
par(mfrow=c(1,1))
z=matrix(my.data47$hazard, ncol=length(30:90))
persp(seq(0,5,by=0.1), 30:90, z, xlab="time since diagnosis", ylab="age at diagnosis", zlab="excess mortality hazard",theta=20)
title("Hazard(t,age) vs (t,age)")


#----------------------------------------------------------------------------------------------------------------------------------------
## 4.9 ==> from model47 make figures of the HR of age at time 0.2, 0.5, 1, 5

# The HR are calculted by omiting the columns and parameters corresponding to the baseline

my.data47$HR=exp(my.mat47[,-(1:5)]%*%beta47[-(1:5)])
head(my.data47)

plot(0,0, type="n", xlim=c(30,90), ylim=c(0.5,2),xlab="Age at diagnosis", ylab="HR (ref=70) ")
lines( my.data47$age[my.data47$fu==0.2], my.data47$HR[my.data47$fu==0.2], col=1)
lines( my.data47$age[my.data47$fu==0.5], my.data47$HR[my.data47$fu==0.5], col=2)
lines( my.data47$age[my.data47$fu==1], my.data47$HR[my.data47$fu==1], col=3)
lines( my.data47$age[my.data47$fu==5], my.data47$HR[my.data47$fu==5], col=4)
axis(1, label=F, at=70, tck=1,lty=8, lwd=0.1)
axis(2, label=F, at=1, tck=1,lty=8, lwd=0.1)
legend(30,2, c("t=0.2","t=0.5","t=1","t=5"), text.col=1:4)
title("Hazard Ratio vs age for given time  \n reference= 70 years old")


#----------------------------------------------------------------------------------------------------------------------------------------
## 4.10 ==> from model47, calculate "by hand" (ie like exercise 3.4) the NS at 2.1 years for individuals with the same age as idenpat n�8
#	    check the result with predict	

library(statmod)
GL <- gauss.quad(n=20,kind="legendre")
Rescale <- function(gl,a,b){gl$nodes <- gl$nodes*(b-a)/2+(a+b)/2 ; gl$weights <- gl$weights*(b-a)/2; return(gl) }
gg <- Rescale(GL,0,2.1)
gg

my.data410=data.frame(fu=gg$nodes)
my.data410=merge(my.data410, temp[temp$idenpat==8,c("agecr","agecr2","agecr3","agecr3plus")] )

my.mat410=model.matrix(formula47, data=my.data410 )


# Compute of cumulative hazard ch as a simple weighted sum of hazard then ns=exp(-ch)
exp(-sum(gg$weights*exp(my.mat410%*%beta47)))


# check with mexhaz
predict(model47, time.pts=c(2.1), data.val=temp[temp$idenpat==8,] )




#----------------------------------------------------------------------------------------------------------------------------------------
## 4.11 ==> calculate the populational NS by averaging the individual NS

mean(predict(model47, time.pts=c(2.1), data.val=temp)$results$surv)


date()


sink()






