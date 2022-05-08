library(relsurv)
library(survival)
library(statmod)
library(mexhaz)
data(colrec)
source("frpop.r")

#put two graphs on the same figure:
par(mfrow=c(1,2))
#(i) age range:
hist(colrec$age/365.241,main="Age distribution",xlab="Age")
range(colrec$age/365.241)
#sex distribution
table(colrec$sex)
#(ii) year of diagnosis range
range(colrec$diag)
hist(colrec$time/365.241, main="Follow-up time distribution",xlab="Time")
#(iii) follow-up times range
range(colrec$time/365.241)
#number of deaths
sum(colrec$stat)

#frpop <- transrate.hmd("mltper_1x1.txt","fltper_1x1.txt")
#attributes(frpop)
attributes(frpop)$dimid
attributes(frpop)$dim
frpop["80", "2003",]
exp(-frpop["80", "2003", "male"]*365)

####2.0 Overall and expected survival of the patients
overall_surv <- survfit(Surv(time, stat)~1, data = colrec)
plot(overall_surv, xlab = "Time (years)", ylab = "Survival", xscale = 365.241)
#(ii)
summary(overall_surv, times=c(5,10)*365.241)
#(iii)
exp.surv <- survexp( ~ 1, rmap=list(sex=sex, year=diag,age=age),times=(0:22)*365.241,data=colrec, ratetable=frpop)
lines(exp.surv, col=2)
#(iv)
summary(exp.surv,times=c(5,10)*365.241)

#(i) Censor all individuals after 5 years of follow-up
colrec$time5 <- pmin(colrec$time,5*365.241) #limit to 5 years
colrec$stat5 <- ifelse(colrec$time<5*365.241,colrec$stat,0) #censor those with longer follow-> #(ii)
cru <- cmp.rel(Surv(time5, stat5) ~ 1, rmap=list(age = age, sex = sex, year = diag), data = colrec, ratetable = frpop)
summary(cru,times=5*365.241)

#(iii) Estimate 5-year net survival
net_surv <- rs.surv(Surv(time5, stat5) ~ 1, rmap=list(age = age, sex = sex,year = diag),data = colrec, ratetable = frpop, method = "pohar-perme")
summary(net_surv, times = 5*365.241)

#put two graphs on the same figure:
par(mfrow=c(1,2))
plot(cru, xlab = "Time (years)", ylab = "Crude mortality", xscale = 365.241)
plot(net_surv, xlab = "Time (years)", ylab = "Net survival", xscale = 365.241)

###4.0 Net survival with respect to covariates
colrec$age_over65 <- ifelse(colrec$age<=65*365.241,0,1)
net_surv_sex <- rs.surv(Surv(time5, stat5) ~ sex, rmap=list(age = age, sex = sex, year = diag),data = colrec, ratetable = frpop, method = "pohar-perme")
summary(net_surv_sex,times=5*365.241)

####Differences 
rs.diff(Surv(time5, stat5) ~ sex, rmap=list(age = age, sex = sex,year = diag),data = colrec, ratetable = frpop)
net_surv_age <- rs.surv(Surv(time5, stat5) ~ age_over65, rmap=list(age = age,sex = sex, year = diag),data = colrec, ratetable = frpop, method = "pohar-perme")
summary(net_surv_age,times=5*365.241)

rs.diff(Surv(time5, stat5) ~ age_over65, rmap=list(age = age, sex = sex, year = diag),data = colrec, ratetable = frpop)
#The prognosis of this tumor is much worse for the old compared to the young maybe because of biology of the elderly, treatment and maybe late diagnosis of the patinets or cancers
#put two graphs on the same figure:
par(mfrow=c(1,2))
plot(net_surv_sex, xlab = "Time (years)", ylab = "Net survival",xscale = 365.241,col=c(4,2))
legend("topright",fill=c(4,2),legend=c("men","women"))
plot(net_surv_age, xlab = "Time (years)", ylab = "Net survival",xscale = 365.241,col=c(1,"grey"))
legend("topright",fill=c(1,"grey"),legend=c("young","old"))
###Using net survival for chldren????It depends:-childhood mortality of different cancers due to their diff natures in severity

#(i)
breaks <- c(0, seq(from = 45, to = 90, by = 5), Inf)
colrec$agegr <- cut(colrec$age / 365.241, breaks)
nessie(Surv(time, stat) ~ agegr+sex, data = colrec, ratetable = frpop,times = c(0,2,5,10,15), rmap = list(age = age, sex = sex, year = diag))


#(ii)
colrec$time15 <- pmin(colrec$time,15*365.241) #limit to 15 years
#censor those with longer follow-up times
colrec$stat15 <- ifelse(colrec$time<15*365.241,colrec$stat,0)
net_surv <- rs.surv(Surv(time15, stat15)~ 1 , rmap=list(age = age, sex = sex,year = diag), data = colrec[colrec$age <= 80*365.241,],ratetable = frpop,method = "pohar-perme")
plot(net_surv, xlab = "Time (years)", ylab = "Net survival", xscale = 365.241)