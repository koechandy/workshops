library(readr)
SGdata <- read_csv("C:/Users/USER/Desktop/Statistical genetics workshop/17snps2.csv")
View(SGdata)
SGdata2 <- SGdata[,-1]
SGdataT <- as.data.frame(t(SGdata2))
snp1 <- t(matrix(SGdataT$V1, nrow = 3, ncol = 2, dimnames = list(c("AA","AB","BB"), c("Case", "Control"))))
print(snp1)
chisq.test(snp1) # gives the Pearson (genotype) chi-square test statistic
library(genetics)
SNP1 <- genotype(snp1)
HWE.chisq(SNP1) # gives the HWE test statistic
#HWE.exact(SNP1) --can only be computed for 2 markers with 2 alleles
#This gives the list of the 17 2x3 tables and their chi-square statistics
s <- vector("list")  
ChisqT <- double(length = 17)
for (i in colnames(SGdataT)) {
  s[[i]] <- t(matrix(eval(parse(text = paste0("SGdataT", "$", i))), nrow = 3, ncol = 2, dimnames = list(c("AA","AB","BB"), c("Case", "Control"))))
}

for (i in seq(1, length(SGdataT))) {
  ChisqT[[i]] <- chisq.test(s[[i]])$p.value
}

#library(genetics)
#sgt<-genotype(SGdataT$)
#for (i in seq(1, length(SGdataT))) {
 # ChisqT[[i]] <- HWE.exact(s[[i]])$p.value
#}
# Calculates the HWE p-values for the 17 SNPs based on the data for "controls" only
HWE_p=matrix(data=NA,nrow=nrow(SGdata),ncol=1)
for (i in (1:nrow(SGdata))){
  data3=c(rep("A/A",SGdata[i,5]),rep("A/B",SGdata[i,6]),rep("B/B",SGdata[i,7]))
  g3=genotype(data3)
  test=HWE.exact(g3)
  HWE_p[i,]=unlist(test["p.value"])
}

SGdata=data.frame(SGdata,HWE_p)

#write.table(data,file="HWEtest17snps.csv",row.names=FALSE,na="",col.names=TRUE,sep=",")


#*********************************
colnames(HWE_p) <- c("pval_HWE")
combined_SNdata=data.frame(SGdata,HWE_p)

s2 <- NULL
for (i in colnames(SGdataT)) {
  snp <- t(matrix(SGdataT[,i], nrow = 3, ncol = 2, dimnames = list(c("AA","AB","BB"), c("Case", "Control"))))#General Chi-Square
  snp
  snpd=matrix(c(snp[,1]+snp[,2],snp[,3]), nrow = 2, ncol = 2, dimnames = list(c("Case", "Control"),c("AA or AB","BB")))#Dominant Chi-Square
  snpd
  snpr=matrix(c(snp[,1],snp[,2]+snp[,3]), nrow = 2, ncol = 2, dimnames = list(c("Case", "Control"),c("AA","BB or AB")))#Recessive Chi-Square
  snpr
  snpAllele=matrix(c(2*(snp[,1])+snp[,2],(2*(snp[,3])+snp[,2])), nrow = 2, ncol = 2, dimnames = list(c("Case", "Control"),c("A","B")))#Allelic Chi-Square
  snpAllele
  
  s1 <- cbind(chisq.test(snp)$p.value,chisq.test(snpd)$p.value,chisq.test(snpr)$p.value,chisq.test(snpAllele)$p.value)
  
  s2= rbind(s2,s1)
  
}

colnames(s2) <- c("pval_general", "pval_dominant", "pval_recessive", "pval_allelic")
s1
s2

combined_SNdata=data.frame(combined_SNdata,s2)


