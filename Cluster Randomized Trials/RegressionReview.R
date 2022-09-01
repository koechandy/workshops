#Regression Review
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

plot(bwt.first$age, bwt.first$weight)
abline(M1)


