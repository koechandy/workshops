######################
# sample size for cRCT section
# note that all these results are PER GROUP

# calculate the sample size per group if there were no clustering
ss.simple = power.t.test(n=NULL, delta=20, sd=30, power=0.8, sig.level=0.05)
ss.simple

# Assume cluster size of 10 and ICC of 0.05
# Calculate the design effect
J = 10
ICC = 0.05
deff = 1 + (J-1)*ICC
deff

# calculate the inflated sample size
ss.star = ceiling(ss.simple$n)*deff
ss.star

# calculate the minimum number of clusters
library(clusterPower)
clusters.simple = ceiling(ss.star/J)
clusters.simple

# calculate minimum number of clusters using clusterPower routine
cluster.R1 = crtpwr.2mean(alpha=0.05, 
                          power=0.8, 
                          m=NA, 
                          n=10, 
                          cv=0, 
                          d=20, 
                          icc=0.05, 
                          varw=900)
cluster.R1

# calculate total needed sample size
nstar.R1 = cluster.R1*J
nstar.R1
