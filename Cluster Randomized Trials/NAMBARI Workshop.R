# R code for regression review

library ("arm")
library("foreign")
library("readstata13")
library("geepack")

# you will have to change this to point to the directory where you save the dataset
setwd("C:/Users/ANDREW K/Desktop/Advanced Biostatistical Methods/Cluster Randomized Trials")


## This section reads in the birthweight data and does some pre-processing
bwt = read.dta13("birthwt_long.dta")

## partial listing of birthweight data
head(bwt, n=10)

## restrict to first births (order=1)
bwt.first = bwt[bwt$order==1,]

# restrict to first births (order=1)
head(bwt.first, n=10)

# scatterplot
plot(bwt.first$age, bwt.first$weight)

# regression of weight on age at first birth (baseline)
M0 = lm( weight ~ 1, data=bwt.first )
display(M0)

# regression of weight on age at first birth (baseline)
M1 = lm( weight ~ 1 + age, data=bwt.first )
display(M1)

plot(bwt.first$age, bwt.first$weight, pch="*", xlab = "Age",ylab = "Weight")
abline(M1, col="green",lty=4,lwd=3)

##########################################################################
#Radon Analysis
## Analysis of Radon Data
## NAMBARI Short Course 2019
## Adapted from Gelman and Hill, 2006

library ("arm")

# set your own directory here
setwd("C:/Users/ANDREW K/Desktop/Advanced Biostatistical Methods/Cluster Randomized Trials")


### This section reads in and re-organizes radon data
srrs2 <- read.table ("srrs2.dat", header=T, sep=",")
idnum <- srrs2$idnum
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon
x <- floor

# get county index variable
county.name <- as.vector(srrs2$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}

# look at data
data = cbind(idnum, county, radon, log.radon, floor)

head(data, n=20)

# no predictors
ybarbar = mean(y)

sample.size <- as.vector (table (county))
sample.size.jittered <- sample.size*exp (runif (J, -.1, .1))
cty.mns = tapply(y,county,mean)
cty.vars = tapply(y,county,var)
cty.sds = mean(sqrt(cty.vars[!is.na(cty.vars)]))/sqrt(sample.size)
cty.sds.sep = sqrt(tapply(y,county,var)/sample.size)


par(mfrow=c(1,1))
plot (sample.size.jittered, cty.mns, cex.lab=.9, cex.axis=1,
      xlab="sample size in county j",
      ylab="avg. log radon in county j",
      pch=20, log="x", cex=.3, mgp=c(1.5,.5,0),
      ylim=c(0,3.2), yaxt="n", xaxt="n")
axis (1, c(1,3,10,30,100), cex.axis=.9, mgp=c(1.5,.5,0))
axis (2, seq(0,3), cex.axis=.9, mgp=c(1.5,.5,0))
for (j in 1:J){
  lines (rep(sample.size.jittered[j],2),
         cty.mns[j] + c(-1,1)*cty.sds[j], lwd=.5)
  #         cty.mns[j] + c(-1,1)*mean(cty.sds[!is.na(cty.sds)]), lwd=.5)
}
abline(h=ybarbar)
title("No pooling",cex.main=.9, line=1)

#abline(h=ybarbar)
points(sample.size.jittered[36],cty.mns[36],cex=4)

### pooled estimate
M0 = lm( log.radon ~ 1 )
display(M0)

### county-specific intercepts (mean for each county)
M1 = lm( log.radon ~ floor+factor(county)-1 )
display(M1)

### fit linear mixed model with random intercept
M2 = lmer(log.radon ~ floor+ (1 | county))
display(M2)

fixef(M2)
ranef(M2)

### generate predictions for each county
beta.j = as.matrix( ranef(M2)$county, ncol=1 )  # random effects
mu.j    = as.matrix( coef(M2)$county, ncol=1 )
se.mu.j = as.matrix( se.ranef(M2)$county, ncol=1 )
lower.ci.mu.j = mu.j - 2*se.mu.j
upper.ci.mu.j = mu.j + 2*se.mu.j

plot (sample.size.jittered, mu.j, cex.lab=.9, cex.axis=1,
      xlab="sample size in county j",
      ylab="avg. log radon in county j",
      pch=20, log="x", cex=.3, mgp=c(1.5,.5,0),
      ylim=c(0,3.2), yaxt="n", xaxt="n")
axis (1, c(1,3,10,30,100), cex.axis=.9, mgp=c(1.5,.5,0))
axis (2, seq(0,3), cex.axis=.9, mgp=c(1.5,.5,0))
for (j in 1:J){
  lines (rep(sample.size.jittered[j],2),
         mu.j[j] + c(-1,1)*se.mu.j[j], lwd=.5)
  #         cty.mns[j] + c(-1,1)*mean(cty.sds[!is.na(cty.sds)]), lwd=.5)
}
abline(h=ybarbar)
title("Pooling with multilevel model",cex.main=.9, line=1)
points(sample.size.jittered[36],mu.j[36],cex=4)



###############################################################################
#Radon 2 analysis
## Read & clean the data
# get radon data
# Data are at http://www.stat.columbia.edu/~gelman/arm/examples/radon
library ("arm")

setwd("/Users/jhogan/GoogleDrive/TEACHING/PHP2517-MultilevelModeling/Analyses/")


### This section reads in and re-organizes radon data
srrs2 <- read.table ("srrs2.dat", header=T, sep=",")
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon
x <- floor

# get county index variable
county.name <- as.vector(srrs2$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}




## this part looks at individual-level predictor (floor)

## Complete pooling regression
lm.pooled <- lm (y ~ x)
display (lm.pooled)

## No pooling regression
lm.unpooled <- lm (y ~ x + factor(county) -1)
display (lm.unpooled)

## Comparing-complete pooling & no-pooling (Figure 12.2)
x.jitter <- x + runif(n,-.05,.05)
display8 <- c (36, 1, 35, 21, 14, 71, 61, 70)  # counties to be displayed
y.range <- range (y[!is.na(match(county,display8))])

par (mfrow=c(2,4), mar=c(4,4,3,1), oma=c(1,1,2,1))
for (j in display8){
  plot (x.jitter[county==j], y[county==j], xlim=c(-.05,1.05), ylim=y.range,
        xlab="floor", ylab="log radon level", cex.lab=1.2, cex.axis=1.1,
        pch=20, mgp=c(2,.7,0), xaxt="n", yaxt="n", cex.main=1,
        main=uniq[j])
  axis (1, c(0,1), mgp=c(2,.7,0), cex.axis=1.1)
  axis (2, seq(-1,3,2), mgp=c(2,.7,0), cex.axis=1.1)
  curve (coef(lm.pooled)[1] + coef(lm.pooled)[2]*x, lwd=.5, lty=2, add=TRUE)
  curve (coef(lm.unpooled)[j+1] + coef(lm.unpooled)[1]*x, lwd=.5, col="blue", add=TRUE)
}

## No-pooling ests vs. sample size (plot on the left on figure 12.3)
sample.size <- as.vector (table (county))
sample.size.jittered <- sample.size*exp (runif (J, -.1, .1))


par(mfrow=c(1,1))
par (mar=c(5,5,4,2)+.1)
plot (sample.size.jittered, coef(lm.unpooled)[-1], cex.lab=1.2, cex.axis=1.2,
      xlab="sample size in county j", ylab=expression (paste
                                                       ("est. intercept, ", alpha[j], "   (no pooling)")),
      pch=20, log="x", ylim=c(.15,3.5), yaxt="n", xaxt="n")
axis (1, c(1,3,10,30,100), cex.axis=1.1)
axis (2, seq(0,3), cex.axis=1.1)
for (j in 1:J){
  lines (rep(sample.size.jittered[j],2),
         coef(lm.unpooled)[j+1] + c(-1,1)*se.coef(lm.unpooled)[j+1], lwd=.5)
}



## Varying-intercept model w/ no predictors
M0 <- lmer (y ~ 1 + (1 | county))
display (M0)

## Including x as a predictor
M1 <- lmer (y ~ x + (1 | county))
display (M1)

# estimated regression coefficicents
coef (M1)

# fixed and random effects
fixef (M1)
ranef (M1)

# uncertainties in the estimated coefficients
se.fixef (M1)
se.ranef (M1)

# 95% CI for the slope
fixef(M1)["x"] + c(-2,2)*se.fixef(M1)["x"]
#or
fixef(M1)[2] + c(-2,2)*se.fixef(M1)[2]

# 95% CI for the intercept in county 26
coef(M1)$county[26,1] + c(-2,2)*se.ranef(M1)$county[26]

# 95% CI for the error in the intercept in county 26
as.matrix(ranef(M1)$county)[26] + c(-2,2)*se.ranef(M1)$county[26]

# to plot Figure 12.4
a.hat.M1 <- coef(M1)$county[,1]                # 1st column is the intercept
b.hat.M1 <- coef(M1)$county[,2]                # 2nd element is the slope

par (mfrow=c(2,4))
for (j in display8){
  plot (x.jitter[county==j], y[county==j], xlim=c(-.05,1.05), ylim=y.range,
        xlab="floor", ylab="log radon level", main=uniq[j],cex.lab=1.2,
        cex.axis=1.1, pch=20, mgp=c(2,.7,0), xaxt="n", yaxt="n", cex.main=1.1)
  axis (1, c(0,1), mgp=c(2,.7,0), cex.axis=1)
  axis (2, c(-1,1,3), mgp=c(2,.7,0), cex.axis=1)
  curve (coef(lm.pooled)[1] + coef(lm.pooled)[2]*x, lty=2, col="black", add=TRUE)
  curve (coef(lm.unpooled)[j+1] + coef(lm.unpooled)[1]*x, col="blue", add=TRUE)
  curve (a.hat.M1[j] + b.hat.M1[j]*x, lwd=1, col="red", add=TRUE)
}  

## Multilevel model ests vs. sample size (plot on the right on figure 12.3)
a.se.M1 <- se.coef(M1)$county

par(mfrow=c(1,1))
par (mar=c(5,5,4,2)+.1)
plot (sample.size.jittered, t(a.hat.M1), cex.lab=1.2, cex.axis=1.1,
      xlab="sample size in county j", ylab=expression (paste
                                                       ("est. intercept, ", alpha[j], "   (multilevel model)")),
      pch=20, log="x", ylim=c(.15,3.5), yaxt="n", xaxt="n")
axis (1, c(1,3,10,30,100), cex.axis=1.1)
axis (2, seq(0,3), cex.axis=1.1)
for (j in 1:J){
  lines (rep(sample.size.jittered[j],2),
         as.vector(a.hat.M1[j]) + c(-1,1)*a.se.M1[j], lwd=.5, col="gray10")
}
abline (coef(lm.pooled)[1], 0, lwd=.5)













# this part looks at group level predictor (uranium level)

## Get the county-level predictor
srrs2.fips <- srrs2$stfips*1000 + srrs2$cntyfips
cty <- read.table ("cty.dat", header=T, sep=",")
usa.fips <- 1000*cty[,"stfips"] + cty[,"ctfips"]
usa.rows <- match (unique(srrs2.fips[mn]), usa.fips)
uranium <- cty[usa.rows,"Uppm"]
u <- log (uranium)


## Varying-intercept model w/ group-level predictors
u.full <- u[county]
M2 <- lmer (y ~ x + u.full + (1 | county))
display (M2)

coef (M2)
fixef (M2)
ranef (M2)

## Plots on Figure 12.5
M1 <- lmer (y ~ x + (1 | county))
a.hat.M1 <- fixef(M1)[1] + ranef(M1)$county                
b.hat.M1 <- fixef(M1)[2]

a.hat.M2 <- fixef(M2)[1] + fixef(M2)[3]*u + ranef(M2)$county
b.hat.M2 <- fixef(M2)[2]

par (mfrow=c(2,4), mar=c(4,4,3,1), oma=c(1,1,2,1))
for (j in display8){
  plot (x.jitter[county==j], y[county==j], xlim=c(-.05,1.05), ylim=y.range,
        xlab="floor", ylab="log radon level", cex.lab=1.2, cex.axis=1.1,
        pch=20, mgp=c(2,.7,0), xaxt="n", yaxt="n", cex.main=1.1, main=uniq[j])
  axis (1, c(0,1), mgp=c(2,.7,0), cex.axis=1.1)
  axis (2, seq(-1,3,2), mgp=c(2,.7,0), cex.axis=1.1)
  curve (a.hat.M1[j,] + b.hat.M1*x, lwd=.5, col="black", add=TRUE)
  curve (a.hat.M2[j,] + b.hat.M2*x, lwd=1, col="red", add=TRUE)
}

# Plot of ests & se's vs. county uranium (Figure 12.6)
a.se.M2 <- se.coef(M2)$county

par(mfrow=c(1,1))
par (mar=c(5,5,4,2)+.1)
plot (u, t(a.hat.M2), cex.lab=1.2, cex.axis=1.1,
      xlab="county-level uranium measure", ylab="est. regression intercept", pch=20,
      ylim=c(0.5,2.0), yaxt="n", xaxt="n", mgp=c(3.5,1.2,0))
axis (1, seq(-1,1,.5), cex.axis=1.1, mgp=c(3.5,1.2,0))
axis (2, cex.axis=1.1, mgp=c(3.5,1.2,0))
curve (fixef(M2)["(Intercept)"] + fixef(M2)["u.full"]*x, lwd=1, col="black", add=TRUE)
for (j in 1:J){
  lines (rep(u[j],2), a.hat.M2[j,] + c(-1,1)*a.se.M2[j,], lwd=.5, col="gray10")
}


## Predictions


x.tilde <- 1 
sigma.y.hat <- sigma.hat(M2)$sigma$data
coef.hat <- as.matrix(coef(M2)$county)[26,]
y.tilde <- rnorm (1, coef.hat %*% c(1, x.tilde, u[26]), sigma.y.hat)

n.sims <- 1000 
coef.hat <- as.matrix(coef(M2)$county)[26,]
y.tilde <- rnorm (1000, coef.hat %*% c(1, x.tilde, u[26]), sigma.y.hat)
quantile (y.tilde, c(.25,.5,.75))


u.tilde <- mean (u)

# extract model components from level 2 model
g.0.hat <- fixef(M2)["(Intercept)"]
g.1.hat <- fixef(M2)["u.full"]
tau.hat <- sigma.hat(M2)$sigma$county

# coefficients from level 1 model
b.hat <- coef.hat[2]

# Simulate possible intercepts for the new county, then possible values
# of radon level for the new house in the county

a.tilde <- rnorm (n.sims, g.0.hat + g.1.hat*u.tilde, tau.hat)
y.tilde <- rnorm (n.sims, a.tilde + b.hat*x.tilde, sigma.y.hat)

