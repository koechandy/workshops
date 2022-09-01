library(genetics)
data1=c(rep("A/A",6),rep("A/B",25),rep("B/B",19))
g1=genotype(data1)
test=HWE.exact(g1)
test

data2=c(rep("A/A",5),rep("A/B",29),rep("B/B",14))
g1=genotype(data2)
test=HWE.exact(g2)
test

data3=c(rep("A/A",10),rep("A/B",218),rep("B/B",929))
g3=genotype(data3)
test=HWE.exact(g3)
test
