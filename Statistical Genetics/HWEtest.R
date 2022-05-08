#install packages the first time only
#install.packages("genetics")

library(genetics)

#SNP from 17snps2.csv where there is disparity
data1=c(rep("A/A",249),rep("A/B",496),rep("B/B",396))
g1=genotype(data1)
#give full output
HWE.test(g1)
#give only test output
test=HWE.exact(g1)

#SNP from 17snps2.csv where there is no disparity
data2=c(rep("A/A",26),rep("A/B",220),rep("B/B",862))
g2=genotype(data2)
HWE.exact(g2)

#setwd("C:/Users/mmlan/OneDrive/Documents/CaseControls")
#data=read.csv("17snps2.csv")
library(readxl)
mydata <- read_excel("C:/Users/USER/Desktop/Statistical genetics workshop/mydata.xlsx")
View(mydata)
HWE_p=matrix(mydata=NA,nrow=nrow(mydata),ncol=1)
for (i in (1:nrow(mydata))){
  data3=c(rep("A/A",mydata[i,5]),rep("A/B",mydata[i,6]),rep("B/B",mydata[i,7]))
  g3=genotype(data3)
  test=HWE.exact(g3)
  HWE_p[i,]=unlist(test["p.value"])
}

mydata=data.frame(mydata,HWE_p)

write.table(data,file="HWEtest17snps.csv",row.names=FALSE,na="",col.names=TRUE,sep=",")

