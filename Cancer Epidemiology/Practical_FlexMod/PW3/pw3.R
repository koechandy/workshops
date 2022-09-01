#*------------------------------------------------------*
#|                             |
#|                                                      |
#| source("C:/Users/Aurel/Dropbox/WorkStat/Teaching/Nairobi-2019_ISCB_ShortCourse/Material/Practical_FlexMod/PW3/pw3.R", echo=T)
#|                                                      |
#-------------------------------------------------------*
# 
# 		 ===========================
#   		|	Pratical Work n°3   |
# 		 ===========================

setwd("C:/Users/Aurel/Dropbox/WorkStat/Teaching/Nairobi-2019_ISCB_ShortCourse/Material/Practical_FlexMod/")


sink("PW3/pw3.lis")

date()
system("hostname", intern=T)
library(survival)
library(mexhaz)
library(splines)

sessionInfo()

#----------------------------------------------------------------------------------------------------------------------------------------
## 3.1 ==> Using datacolon.v2, fit an excess hazard model with a cubic B-spline (one knot at 1 year) as baseline and 
#	   make hazard predictions for time t=seq(0,5,by=0.01) with predict.mexhaz

load("data/datacolon.v2.RData")

stage="St4"
temp=datacolon.v2[datacolon.v2$stage==stage,]

model31=mexhaz(Surv(time=fu,event=dead)~ 1, data=temp, base="exp.bs", knots=c(1), bound=c(0,8), expected="rate")

pred.model31=predict(model31, time.pts=c(0.0001, seq(0.01,5, by=0.01)))
head(pred.model31$results)


#----------------------------------------------------------------------------------------------------------------------------------------
## 3.2 ==> check the fit done in 3.1 by re-drawing the plot of 2.3 and adding the predicted hazard from model31

load("PW2/hazard23.RData"   )

begin.interval=c(0, 0.25, 0.5, 1:4)
end.interval=c(0.25, 0.5, 1:5)

plot(0,0, type="n", xlim=c(0,5), ylim=c(0,2),xlab="Time since diagnosis", ylab="Excess hazard ")
segments(begin.interval, hazard23, end.interval, hazard23, lwd=3,col=3)
lines(pred.model31$results$time.pts, pred.model31$results$hazard, col="red")


#----------------------------------------------------------------------------------------------------------------------------------------
## 3.3 ==> calculate the predicted hazard from model31 "by hand" ie helps are:
#		help 1) write the formulae of model31 in R terms (ie write the model as if you use bs() of library splines,  )
#		help 2) use model.matrix and matrix product with coefficients of model31
#	   add this predicted hazard on the current graph (with lty=8)

library(splines)

# Note the default values for intercept and degree of the bs function !
formula31=as.formula("~bs(fu, k=1, B=c(0,8))" )
beta31=as.numeric(model31$coef)

my.data31=data.frame(fu=seq(0,5,by=0.01))

# my.mat31 is the design matrix
#
my.mat31=model.matrix(formula31, data=my.data31)
my.data31$hazard=exp(my.mat31%*%beta31)
head(my.data31)

lines(my.data31$fu, my.data31$hazard, col="blue", lty=8, lwd=2)


#----------------------------------------------------------------------------------------------------------------------------------------
## 3.4 ==> compute NS at 1.5 years from model31 using Gauss-Legendre approximation with 100 nodes (try also with 20 nodes) and compare the results with predict.mexhaz;
#	    
#	help 1: GL <- gauss.quad(n=100,kind="legendre")
#	help 2: use fn Rescale <- function(gl,a,b){gl$nodes <- gl$nodes*(b-a)/2+(a+b)/2 ; gl$weights <- gl$weights*(b-a)/2; return(gl) }
#		gg <- Rescale(GL,0,1.5)		


library(statmod)
GL <- gauss.quad(n=100,kind="legendre")

## ‘Rescale’ to modify the nodes and weight of Gauss-Legendre approximation according to the integration interval:
Rescale <- function(gl,a,b){gl$nodes <- gl$nodes*(b-a)/2+(a+b)/2 ; gl$weights <- gl$weights*(b-a)/2; return(gl) }
gg <- Rescale(GL,0,1.5)
gg

# Compute the design matrix evaluated  at the 100 integration points (100 nodes)
#
rm(my.mat31)
my.mat31=model.matrix(formula31,data.frame(fu=gg$nodes))


# Compute of cumulative hazard ch as a simple weighted sum of hazard
# then NS=exp(-ch)
ch=sum(gg$weights*exp(my.mat31%*%beta31))
exp(-ch)


# check with mexhaz
predict(model31, time.pts=c(1.5))

# Approximation RECTANGLE ==> NOT TO DO IN REAL LIFE !!
rm(my.mat31)
my.mat31=model.matrix(formula31,data.frame(fu=seq(0,1.5, by=0.0005)))
hazard=exp(my.mat31%*%beta31)
ch2=sum(0.0005*hazard)
exp(-ch2)

# with 20 nodes
GL <- gauss.quad(n=20,kind="legendre")
gg <- Rescale(GL,0,1.5)
rm(my.mat31)
my.mat31=model.matrix(formula31,data.frame(fu=gg$nodes))
exp(-sum(gg$weights*exp(my.mat31%*%beta31)))



#----------------------------------------------------------------------------------------------------------------------------------------
## 3.5 ==> plot the hazard and NS in black color using plot.predMexhaz; add CI in red

par(mfrow=c(1,2))
plot(pred.model31, which="hazard", conf.int=F)
lines(pred.model31$results$time.pts, pred.model31$results$hazard.inf, col="red")
lines(pred.model31$results$time.pts, pred.model31$results$hazard.sup, col="red")
#
plot(pred.model31, which="surv", conf.int=F)
lines(pred.model31$results$time.pts, pred.model31$results$surv.inf, col="red")
lines(pred.model31$results$time.pts, pred.model31$results$surv.sup, col="red")
#

#---------------------------------------------------------------------------------------------------------------
## 3.6 ==> check that the value of the likelihood do not depend on population mortality hazard of alive patients
#	   	

# initial log-likelihood
model31$loglik

temp$rate0=ifelse(temp$dead==0, 0, temp$rate)
temp$rate10=ifelse(temp$dead==0, 10*temp$rate, temp$rate)
temp$rate500=ifelse(temp$dead==0, 500*temp$rate, temp$rate)

ll0=mexhaz(Surv(time=fu,event=dead)~ 1, data=temp, base="exp.bs", knots=c(1), bound=c(0,8), expected="rate0", verbose=0)$loglik
ll10=mexhaz(Surv(time=fu,event=dead)~ 1, data=temp, base="exp.bs", knots=c(1), bound=c(0,8), expected="rate10", verbose=0)$loglik
ll500=mexhaz(Surv(time=fu,event=dead)~ 1, data=temp, base="exp.bs", knots=c(1), bound=c(0,8), expected="rate500")$loglik


c(model31$loglik, ll0, ll10, ll500)



date()

sink()






