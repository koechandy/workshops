# Load necessary package
library(mexhaz)
mypath <-"F:/International Cancer Institute-ICI/Cancer Epidemiology/ForStudent/Practical_FlexMixMod/"
mydat <- read.delim(paste0(mypath, "fakeLOCP.dat"))
# select only men
mydat <- subset(mydat, sex==1)
summary(mydat)
# create survival time and corresponding vital status
mydat$timesurv10y <- pmin(mydat$timesurv, 10)
mydat$status10y <- ifelse(mydat$timesurv10y==mydat$timesurv,mydat$status,0)
 # Creation of useful variables

mydat$Iagecat1545 <- ifelse(mydat$agediag>=15 & mydat$agediag<45,1,0)
mydat$Iagecat4555 <- ifelse(mydat$agediag>=45 & mydat$agediag<55,1,0)
mydat$Iagecat5565 <- ifelse(mydat$agediag>=55 & mydat$agediag<65,1,0)
mydat$Iagecat6575 <- ifelse(mydat$agediag>=65 & mydat$agediag<75,1,0)
mydat$Iagecat75pp <- ifelse(mydat$agediag>=75 ,1,0)
# Alternative coding:
# mydat$agecat <- cut(mydat$agediag, breaks=c(15, 45, 55, 65, 75, 150), right = FALSE)
# table(mydat$agecat)

# Creation of the variable for a spline of age (deg 3, knot 70) in a truncated power basis
# not reduced
mydat$agediagc=mydat$agediag-70
mydat$agediagc2 <- mydat$agediagc^2
mydat$agediagc3 <- mydat$agediagc^3
mydat$agediagc3plus <- (mydat$agediagc-0)^3*(mydat$agediagc>0)
 # redcued
mydat$agediagcr=(mydat$agediag-70)/10
mydat$agediagcr2 <- mydat$agediagcr^2
mydat$agediagcr3 <- mydat$agediagcr^3
mydat$agediagcr3plus <- (mydat$agediagcr-0)^3*(mydat$agediagcr>0)

# Number of events by deprivation quintiles
with(mydat, table(quintile, status10y))

# Distribution of EDI and agediag in the population
 summary(mydat$EDI) 
 summary(mydat$agediag)
 
 
 FPM1 <- mexhaz(Surv(timesurv10y, status10y)~ Iagecat1545 + Iagecat4555 + Iagecat5565 + Iagecat75pp,
                expected="expectedrate", degree=3, knots = c(1,5), base="exp.bs", data=mydat, verbose = 0)
summary(FPM1)
 
mytime <- seq(0.01,10,0.01)
predFPM1 <- predict(FPM1, time.pts = mytime,
                       data.val = data.frame(Iagecat1545=0, Iagecat4555=0, Iagecat5565=0, Iagecat75pp=0))
par(mfrow=c(1,2), oma = c(0, 0, 2, 0), mgp=c(1.5,0.4,0), mar=c(3, 3, 2.5, 1))
plot(predFPM1, which="hazard", xlim=c(0,10), ylim=c(0,0.5),conf.int = T, main="Baseline excess \n mortality hazard", lwd=2,panel.first=abline(h=seq(0,0.5,by=0.1), tck=1, lty=8, col="grey"))
plot(predFPM1, which="surv", xlim=c(0,10), ylim=c(0,1),conf.int = T, main="Net survival \n (Reference group)", lwd=2,panel.first=abline(h=seq(0,1,by=0.1), tck=1, lty=8, col="grey"))
mtext("Model 1 (Age-category)", outer = TRUE, cex = 1.5)
 
 