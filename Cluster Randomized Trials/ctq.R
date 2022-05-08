library ("arm")
library("foreign")
library("readstata13")
library("geepack")
setwd("C:/Users/ANDREW K/Desktop/Advanced Biostatistical Methods/Cluster Randomized Trials")

ctq = read.dta13("ctq_mlm.dta")

ctq$include = rep(1, nrow(ctq))
for (i in unique(ctq$id))
{ 
  nmiss                   = sum(is.na( ctq$Y[(ctq$id == i)] ) )
  ctq$include[ctq$id==i]  = rep( ifelse(nmiss<9, 1, 0), 9 )
}

ctq = ctq[(ctq$include==1),]


#include = 1 - is.na(ctq$Y)
#ctq = ctq[(include==1),]
#ctq = ctq[(include==1 & ctq$target_wk==1),]
#ctq$week = ctq$week - 12

# intercept only model
M0 = glmer( Y ~ 1 + (1 | id) , family=binomial, data=ctq)
display(M0)
alpha.hat = coef(M0)$id
prob.hat  = exp(alpha.hat) / (1 + exp(alpha.hat))

cbind(alpha.hat, prob.hat)


M1 = glmer( Y ~ totfager + (1 | id), family=binomial(link=logit), data=ctq)
display(M1)

totfager = (0:10)
unique.id = unique(ctq$id)

## logit scale
plot(totfager, ((totfager/10)* 15 - 8), pch="", ylab="log[ p / (1-p) ]" )
for (j in 1:length(unique.id))
{
  beta.hat = coef(M1)$id[j,]
  X        = cbind( rep(1,11), totfager)
  linpred  = X %*% t(beta.hat)
  lines(totfager, linpred, col="grey")
  
  avepred  = X %*% fixef(M1)
  lines(totfager, avepred, lwd=3)
}

## probability scale

plot( totfager, totfager/10, pch="", ylab="Pr(Yj = 1)" )
phat = matrix(0, nrow=length(unique.id), ncol=11)

for (j in 1:length(unique.id))
{
  beta.hat = coef(M1)$id[j,]
  X        = cbind( rep(1,11), totfager)
  linpred  = X %*% t(beta.hat)
  linpred  = matrix(linpred, nrow=1)
  phat[j,]     = exp(linpred) / (1 + exp(linpred))
  lines(totfager, phat[j,], col="grey")
}

lines(totfager, apply(phat, 2, mean), lwd=3, col="red")
cond.logit = apply( (X %*% t(coef(M1)$id)), 1, mean) 
cond.prob  = exp(cond.logit) / (1 + exp(cond.logit))
lines(totfager, cond.prob, lwd=4)


ctq$week_centered = ctq$week - 12

M2 = glmer( Y ~ target_wk + target_wk:Z + (1  | id ), family=binomial(link=logit), data=ctq)
display(M2)

logit.hat = cbind(unique(ctq$id), matrix(0,nrow=length(unique(ctq$id)), ncol=10) )
prob.hat = cbind(unique(ctq$id), matrix(0,nrow=length(unique(ctq$id)), ncol=10) )
                                    
j=0

par(mfrow=c(2,3))
for (i in unique(ctq$id))
{ j = j+1
  select = (i == ctq$id)
  
  X.i = cbind( rep(1,9), ctq$target_wk[select], ctq$Z_target[select])
  beta.hat.i = t( coef(M2)$id[(i == as.numeric(rownames(coef(M2)$id))), ] )
  
  logit.hat[j,2]    = ctq$Z[select][1]
  prob.hat[j,2]     = logit.hat[j,2]
  
  logit.hat[j,3:11] = t( X.i %*% beta.hat.i )
  prob.hat[j,3:11]  = exp(logit.hat[j,3:11]) / (1+exp(logit.hat[j,3:11]))
  
  
  plot( ctq$week[select], ctq$Y[select], xlim=c(4,12), ylim=c(0,1) )
  color = ifelse(ctq$Z[select][1]==1, "red", "blue")
  lines(ctq$week[select==1], prob.hat[j,3:11], col=color, lwd=2)
}

par(mfrow=c(1,1))

prob.0 = apply(prob.hat[prob.hat[,2]==0,], 2, mean)
prob.1 = apply(prob.hat[prob.hat[,2]==1,], 2, mean)

plot( 4:12, prob.0[3:11], col="blue", lwd=2, ylim=c(0,.5), pch="")
lines( 4:12, prob.0[3:11], col="blue", lwd=2)
lines( 4:12, prob.1[3:11], col="red", lwd=2)

logit.0 = apply(logit.hat[logit.hat[,2]==0,3:11], 2, mean)
logit.1 = apply(logit.hat[logit.hat[,2]==1,3:11], 2, mean)

ss.prob.0 = exp(logit.0) / (1 + exp(logit.0))
ss.prob.1 = exp(logit.1) / (1 + exp(logit.1))

plot( 4:12, prob.0[3:11], col="blue", lwd=2, ylim=c(0,.5), pch="")
lines( 4:12,ss.prob.0, col="blue", lwd=2, lty=2)
lines( 4:12, ss.prob.1, col="red", lwd=2, lty=2)



# repeat these analyses using GEE
G0.exch = geeglm(Y ~ totfager, family=binomial("logit"), data=ctq, 
            id=id, corstr="exchangeable", waves=week)
summary(G0.exch)

G0.indep = geeglm(Y ~ totfager, family=binomial("logit"), data=ctq, 
                    id=id, corstr="independence", waves=week)
summary(G0.indep)

G0.unst = geeglm(Y ~ totfager, family=binomial("logit"), data=ctq, 
                  id=id, corstr="unstructured", waves=week)
summary(G0.unst)


G1 = geeglm(Y ~ target_wk + target_wk:Z, family=binomial("logit"), data=ctq, 
                  id=id, corstr="exchangeable", waves=week)
summary(G1)

G2 = geeglm(Y ~ target_wk + target_wk:Z, family=binomial("logit"), data=ctq, 
                id=id, corstr="unstructured", waves=week)
summary(G2)


