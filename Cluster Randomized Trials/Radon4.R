# This analysis goes with Regression for Multilevel Data

library ("arm")

setwd("C:/Users/ANDREW K/Desktop/Advanced Biostatistical Methods/Cluster Randomized Trials")


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


# this part looks at group level predictor (uranium level)

## Get the county-level predictor
srrs2.fips <- srrs2$stfips*1000 + srrs2$cntyfips
cty <- read.table ("cty.dat", header=T, sep=",")
usa.fips <- 1000*cty[,"stfips"] + cty[,"ctfips"]
usa.rows <- match (unique(srrs2.fips[mn]), usa.fips)
uranium <- cty[usa.rows,"Uppm"]
u <- log (uranium)
u.full <- u[county]

## Varying intercept & slopes w/ no group level predictors
M3 <- lmer (y ~ x + (1 + x | county))
display (M3)

coef (M3)
fixef (M3)
ranef (M3)


# plots on Figure 13.1

x.jitter <- x + runif(n,-.05,.05)
display8 <- c (36, 1, 35, 21, 14, 71, 61, 70)  # counties to be displayed
y.range <- range (y[!is.na(match(county,display8))])

a.hat.M3 <- fixef(M3)[1] + ranef(M3)$county[,1] 
b.hat.M3 <- fixef(M3)[2] + ranef(M3)$county[,2]

b.hat.unpooled.varying <- array (NA, c(J,2))
for (j in 1:J){
  lm.unpooled.varying <- lm (y ~ x, subset=(county==j))
  b.hat.unpooled.varying[j,] <- coef(lm.unpooled.varying)
}

lm.pooled <- lm (y ~ x)

par (mfrow=c(2,4), mar=c(4,4,3,1), oma=c(1,1,2,1))
for (j in display8){
  plot (x.jitter[county==j], y[county==j], xlim=c(-.05,1.05), ylim=y.range,
        xlab="floor", ylab="log radon level", cex.lab=1.2, cex.axis=1.1,
        pch=20, mgp=c(2,.7,0), xaxt="n", yaxt="n", cex.main=1.1, main=uniq[j])
  axis (1, c(0,1), mgp=c(2,.7,0), cex.axis=1.1)
  axis (2, seq(-1,3,2), mgp=c(2,.7,0), cex.axis=1.1)
  curve (coef(lm.pooled)[1] + coef(lm.pooled)[2]*x, lwd=1, lty=2, col="black", add=TRUE)
  curve (b.hat.unpooled.varying[j,1] + b.hat.unpooled.varying[j,2]*x, lwd=1, col="blue", add=TRUE)
  curve (a.hat.M3[j] + b.hat.M3[j]*x, lwd=1, col="red", add=TRUE)
}


## Including group level predictors
M4 <- lmer (y ~ x + u.full + x:u.full + (1 + x | county))
display (M4)

coef (M4)
fixef (M4)
ranef (M4)

a.hat.M4 <- fixef(M4)[1] + fixef(M4)[3]*u + ranef(M4)$county[,1]
b.hat.M4 <- fixef(M4)[2] + fixef(M4)[4]*u + ranef(M4)$county[,2]
a.se.M4 <- se.ranef(M4)$county[,1]
b.se.M4 <- se.ranef(M4)$county[,2]

par(mfrow=c(1,1))
# plot on Figure 13.2(a)
lower <- a.hat.M4 - a.se.M4
upper <- a.hat.M4 + a.se.M4
par (mar=c(5,5,4,2)+.1)
plot (u, a.hat.M4, cex.lab=1.2, cex.axis=1.1, ylim=range(lower,upper), 
      xlab="county-level uranium measure", ylab="regression intercept", 
      pch=20, yaxt="n")
axis (2, c(0,1,1.5,2))
curve (fixef(M4)[1] + fixef(M4)[3]*x, lwd=1, col="black", add=TRUE)
segments (u, lower, u, upper, lwd=.5, col="gray10")
mtext ("Intercepts", line=1)

# plot on Figure 13.2(b)
lower <- b.hat.M4 - b.se.M4
upper <- b.hat.M4 + b.se.M4
par (mar=c(5,5,4,2)+.1)
plot (u, b.hat.M4, cex.lab=1.2, cex.axis=1.1, ylim=range(lower,upper),
      xlab="county-level uranium measure", ylab="regression slope", pch=20)
curve (fixef(M4)[2] + fixef(M4)[4]*x, lwd=1, col="black", add=TRUE)
segments (u, lower, u, upper, lwd=.5, col="gray10")
mtext ("Slopes", line=1)


