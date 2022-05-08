## Analysis of Radon Data
## NAMBARI Short Course 2019
## Adapted from Gelman and Hill, 2006

library ("arm")

# set your own directory here
setwd("/Users/jhogan/GoogleDrive/TEACHING/PHP2517-MultilevelModeling/Analyses/")


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
M0 = lm( y ~ 1 )
display(M0)

### county-specific intercepts (mean for each county)
M1 = lm( y ~ factor(county) - 1 )
display(M1)

### fit linear mixed model with random intercept
M2 = lmer(y ~  (1 | county))
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







