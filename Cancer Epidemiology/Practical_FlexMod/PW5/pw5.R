#*------------------------------------------------------*
#|                             |
#|                                                      |
#| source("C:/Users/Aurel/Dropbox/WorkStat/Teaching/Nairobi-2019_ISCB_ShortCourse/Material/Practical_FlexMod/PW5/pw5.R", echo=T)
#|                                                      |
#-------------------------------------------------------*
# 
# 		 ===========================
#   		|	Pratical work n°5   |
# 		 ===========================

setwd("C:/Users/Aurel/Dropbox/WorkStat/Teaching/Nairobi-2019_ISCB_ShortCourse/Material/Practical_FlexMod/")


sink("PW5/pw5.lis")

date()
rm(list=ls())
library(mexhaz)

sessionInfo()


#-------------------------------------------------------------------------------
## 5.0 ==> load data datacolon.v2

load("data/datacolon.v2.RData")

data5=datacolon.v2
end.time=max(data5$fu) # last day of follow-up
  
#		        ==========================================
#   		   |	     model 1: log(h.e) = beta0      |
# 		        =========================================

#-------------------------------------------------------------------------------
## 5.1 ==> write log-likelihood function for model 1 

  log.lik=function(beta0, data){
    
	t0=rep(0,length(data$fu)) # initial time is zero
	t1=data$fu # time since diagnostic
	delta=data$dead # 1 if the patient died, 0 if he was censored
	h.p=data$rate # population hazard
	
 
	
    # excess hazard is
    h.e = exp(beta0)
    
    # Cumulative excess hazard is
    cumul.h.e=(t1-t0)*h.e
 
    # Vector of indiviual log-likelihoods
    ll.indiv=-cumul.h.e + delta*log(h.e+h.p)
    
    # Log-likelihood of the model is the sum of the indiviual ones
    ll=sum(ll.indiv)
    
    # returns the value of the log-likelihood for a given vector beta
    return(ll)
  }
  
  
#-------------------------------------------------------------------------------
## 5.2 ==> Maximize the log-likelihood function for model 1 then get estimation of 
## beta and the covariance matrix (and so the standard-error)

  # initialize parameter estimate
  beta.ini=log(sum(data5$dead)/sum(data5$fu))
 
	# sum(delta) gives the number of deaths and sum(t) the number
	# of person-days
  
  
  # proceed to maximization
  result=optim(beta.ini,log.lik,data=data5, hessian=TRUE, control=list(fnscale=-1),
                method="Brent",lower=-100,upper=0)

				
# N.B : control=list(fnscale=-1) transforms the default 
# minimization problem into a maximization one.
# We could have also minimized minus the log-likelihood
  
  # get estimation of beta
  beta.estim=result$par
  beta.estim
 
  # get hessian matrix of the log-likelihood
  Hess=result$hessian
  Hess
  
  # get covariance matrix for beta.estim (and standard error)
  Covariance=solve(-Hess)
  Covariance
  sqrt(Covariance)


  # get the log-likelihood at beta.estim
  log.lik(beta.estim, data=data5)
 	
  result$value
 	

#-------------------------------------------------------------------------------
## 5.3 ==> Check results with mexhaz package for model 1 and make plot of the 
## hazard and the net survival
  
  model.check=mexhaz(formula=Surv(time=fu,event=dead)~1,data5,
                      expected="rate",base="pw.cst",print.level=0,verbose=0)
  

  summary(model.check)
  
  # log-likelihood
  model.check$loglik

  # predictions
  surv.time=seq(0.001,end.time,length=50)
  

  res=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agediag=50)) # the value for age doesn't
  # matter since the excess hazard doesn't depend on age in model 1 
  
   

  ##### Plot the results
  
  par(mfrow=c(1,2)) # we want two graphs sides by side
  # Excess hazard
  # plot(res$results$time.pts,res$results$hazard,type="l",
       # ylab="Excess hazard",xlab="Time since diagnosis (years)",
	   # main="Model 1 : Excess hazard")
  plot(res$results$time.pts,res$results$hazard,type="l",ylim=c(0,0.15),
	ylab="Excess hazard",xlab="Time since diagnosis (years)",main="Model 1 : Excess hazard")
   
  # Net survival
  plot(res$results$time.pts,res$results$surv,type="l", ylim=c(0,1),
       ylab="Net survival",xlab="Time since diagnosis (years)",
	   main="Model 1 : Net survival")
    
   
#  =======================================================================
# |	      model 2: log(h.e) = beta0 + beta1*time + beta2*age         |
#  =======================================================================

   
  # we can substract its mean to agediag in order to center the data
  data5$agec=data5$agediag-mean(data5$agediag)
 

#-------------------------------------------------------------------------------
## 5.4 ==> write log-likelihood function for model 2
  
# remove the last function
rm(log.lik)

  log.lik=function(beta, data){
    
		# data<-data5
		# beta0<- -0.1
	
	t0=rep(0,length(data$fu)) # initial time is zero
	t1=data$fu # time since diagnostic
	delta=data$dead # 1 if the patient died, 0 if he was censored
	h.p=data$rate # population hazard
	end.time=max(t1) # last day of follow-up
	age=data5$agec
	
    # parameter vector beta has three components (beta0, slope for 
    # time and slope for age)
    beta0=beta[1]

    # Technical point for the line below: our reference time is "standardardized"
	# so our parameter estimated are comparable with the ones of mexhaz which uses 
	# for the baseline bs(time, degree=1, Boundary.knots=c(0, end.time) )  
	beta1=beta[2]/end.time 
		
    beta2=beta[3]
    
    # excess hazard is
    h.e = exp(beta0+beta1*t1+beta2*age)
    
    # Cumulative excess hazard is
    if (beta1==0){
      cumul.h.e=(t1-t0)*exp(beta0+beta2*age)
    }else{

      cumul.h.e=exp(beta0+beta2*age)*(exp(beta1*t1)-exp(beta1*t0))/beta1
      
    }
	# note: t0=0 but the above writing is more general and allows delayed entry (ie t0>0)

    # Vector of indiviual log-likelihoods
    ll.indiv=-cumul.h.e + delta*log(h.e+h.p)
    
    # Log-likelihood of the model is the sum of the indiviual ones
    ll=sum(ll.indiv)
    
    # returns the value of the log-likelihood for a given vector beta
    return(ll)
  }

#-------------------------------------------------------------------------------
## 5.5 ==> Maximize the log-likelihood function for model 2 then get estimation of beta, the covariance matrix, and  the standard-error of beta

 
  # initialize parameter estimates
  beta.ini=c(log(sum(data5$dead)/sum(data5$fu)),0,0)
  
  # proceed to maximization
  result=optim(beta.ini,log.lik,data=data5, hessian=TRUE, control=list(fnscale=-1))

  # get estimation of beta
  beta.estim=result$par
  beta.estim
  
  # get hessian matrix of the log-likelihood
  Hess=result$hessian
  Hess
  
  # get covariance matrix for beta.estim and standard-errors
  Covariance=solve(-Hess)
  Covariance
  

  # standard errors of the parameters
  sqrt(diag(Covariance))
 

  # get the log-likelihood at beta.estim
  log.lik(beta.estim, data=data5)
  
  result$value
  

  # predictions of the hazard at 0.5 year (agec=20)
  # see technical point 
  my.t=0.5/end.time
  exp(beta.estim[1]+my.t*beta.estim[2]+20*beta.estim[3])
  


#-------------------------------------------------------------------------------
## 5.6 ==> Check results with mexhaz package for model 2 and make plot of the hazard and the net survival for 3 age: 50, 60 and 70 years old

  model.check=mexhaz(formula=Surv(time=fu,event=dead)~agec,data5,
                      expected="rate",base="exp.bs",degree=1,print.level=0,verbose=0)
					  
  
  summary(model.check)
  
  # log-likelihood
  model.check$loglik
  
 
  # check predictions predictions of the hazard at 0.5 year (agec=20) with mexhaz
  predict(model.check, time.pts=c(0.5), data.val=data.frame(agec=20))

  # predictions for plots
  
  surv.time=seq(0.001,end.time,length=50)
  
  res50=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agec=50-mean(data5$agediag)) ) # age=50
  res60=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agec=60-mean(data5$agediag)) ) # age=60
  res70=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agec=70-mean(data5$agediag)) ) # age=70
  
  
  
  ##### Plot the results
  par(mfrow=c(1,2))
  # Excess hazard
  plot(res50$results$time.pts,res50$results$hazard,type="l",
       ylab="Excess hazard",xlab="Time since diagnosis (years)",ylim=range(c(res50$results$hazard,
                                                     res60$results$hazard,
                                                     res70$results$hazard)),
	   main="Model 2 : Excess hazard")
  lines(res60$results$time.pts,res60$results$hazard,col="blue")
  lines(res70$results$time.pts,res70$results$hazard,col="red")
  legend(2.5, 0.3, c("age = 50", "age = 60", "age = 70"), col = c("black", "blue", "red"),
       text.col = c("black", "blue", "red"),lty=c(1,1,1))
  # Net survival
  plot(res50$results$time.pts,res50$results$surv,type="l", ylim=c(0,1),
       ylab="Net survival",xlab="Time since diagnosis (years)",
	   main="Model 2 : Net survival")
  lines(res60$results$time.pts,res60$results$surv,col="blue")
  lines(res70$results$time.pts,res70$results$surv,col="red")
  legend(2, 0.5, c("age = 50", "age = 60", "age = 70"), col = c("black", "blue", "red"),
       text.col = c("black", "blue", "red"),lty=c(1,1,1))
  #!!!!!!!!!!!!!!!!
  




#  =======================================================================
# |    model 3: log(h.e) = beta0 + Bspline(time,degree=3) + beta*age   |
#  =======================================================================
# contrary to model 1 and 2, the cumulative hazard of the model3 does not have a 
# closed form and thus requires numerical approximation. 
#-------------------------------------------------------------------------------
## 5.7 ==> Define model 3 with a formula
 

  library(splines) # for using splines
  library(statmod) # for numerical integration (Gauss-Legendre)
  # write the model as a formula object

	
  # Technical point: the baseline of model 3 is simply a classical polynomial; however, 
  # in order to be comparable to mexhaz, last day of follow-up is used to parametrize 
  # the baseline as a BS splines.
  formula1=as.formula(~bs(t1,Boundary.knots=c(0,end.time))+agec)

#-------------------------------------------------------------------------------
## 5.8 ==> Write the design matrices that will be used for Gauss-Legendre integration
  
  # table that gives the nodes and weights for the 10 points Gauss-Legendre 
  # quadrature 
  n.legendre=10
  
  leg1=gauss.quad(n=n.legendre,kind="legendre")
  
	#-------------------------------------------------------------------
	# We calculate the n.legendre design matrices that we will need to 
	# get the cumulative hazard
	#-------------------------------------------------------------------
	X_func<-function(t1,data,formula){
      
      data.t<-data
      data.t$t1<-t1
      model.matrix(formula,data.t)
      
    }
	
	X_GL1<-list()
	
	t0=rep(0,length(data5$fu)) # initial time is zero
	t1=data5$fu # time since diagnostic
	
	for(i in 1:n.legendre){
    
	X_GL1[[i]]<-X_func((t1-t0)/2*leg1$nodes[i]+(t0+t1)/2,data=data5,formula=formula1)

	
	}
	
	
  
#-------------------------------------------------------------------------------
## 5.9 ==> write log-likelihood function for model 3

# remove the last function
rm(log.lik)
	

log.lik=function(beta, data, formul, leg, X_GL){
    
	t0=rep(0,length(data$fu)) # initial time is zero
	t1=data$fu # time since diagnostic
	delta=data$dead # 1 if the patient died, 0 if he was censored
	h.p=data$rate # population hazard
	
	
	# construct the model matrix
	X<-model.matrix(formul,data)
  
	# excess hazard is
    h.e = as.vector(exp(X%*%beta))
    
    # Cumulative excess hazard is
    
    cumul.h.e<-0
    
	for(i in 1:length(leg$nodes)){
    
    cumul.h.e<-cumul.h.e+as.vector(exp((X_GL[[i]]%*%beta)))*leg$weights[i]
    
	}
  #!!!!!!!!!!!!!!!!
	# cumul.h.e<-(t1-t0)/2*cumul.h.e
	cumul.h.e<-cumul.h.e*(t1-t0)/2

    # Vector of indiviual log-likelihoods
    ll.indiv=-cumul.h.e + delta*log(h.e+h.p)
    
    # Log-likelihood of the model is the sum of the indiviual ones
    ll=sum(ll.indiv)
    
    # returns the value of the log-likelihood for a given vector beta
    return(ll)
  }
  
  
#-------------------------------------------------------------------------------
## 5.10 ==> Maximize the log-likelihood function for model 3 then get estimation 
# of beta and the covariance matrix (and so the standard-error)

  
  # initialize parameter estimates
  beta.ini=c(log(sum(data5$dead)/sum(data5$fu)),0,0,0,0)
  
  # proceed to maximization
  result=optim(beta.ini,log.lik,data=data5,formul=formula1,leg=leg1,X_GL=X_GL1, hessian=TRUE, control=list(fnscale=-1))
  
  # get estimation of beta
  beta.estim=result$par
  beta.estim
  
  # get hessian matrix of the log-likelihood
  Hess=result$hessian
  Hess
  
  # get covariance matrix and standard errors for beta.estim
  Covariance=solve(-Hess)
  Covariance
  sqrt(diag(Covariance))

  # get the log-likelihood at beta.estim
  log.lik(beta.estim, data=data5, formul=formula1,leg=leg1,X_GL=X_GL1)
  
  result$value
  
#-------------------------------------------------------------------------------
## 5.11 ==> Check results with mexhaz package for model 3 and make plot of the 
## hazard and the net survival for 3 ages: 50, 60 and 70 years old

  
  
  model.check=mexhaz(formula=Surv(time=fu,event=dead)~agec,data5,
                      expected="rate",base="exp.bs",print.level=0,
                      bound=c(0,end.time),verbose=0)
  
  summary(model.check)
  
  # log-likelihood
  model.check$loglik
  
  
  # predictions
  surv.time=seq(0.001,end.time,length=50)
  
  res50=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agec=50-mean(data5$agediag)) ) # age=50
  res60=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agec=60-mean(data5$agediag)) ) # age=60
  res70=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agec=70-mean(data5$agediag)) ) # age=70
  
  ##### Plot the results
  par(mfrow=c(1,2))
   # Excess hazard
  plot(res50$results$time.pts,res50$results$hazard,type="l",
       ylab="Excess hazard",xlab="Time since diagnosis (years)",ylim=range(c(res50$results$hazard,
                                                     res60$results$hazard,
                                                     res70$results$hazard)),
	   main="Model 3 : Excess hazard")
  lines(res60$results$time.pts,res60$results$hazard,col="blue")
  lines(res70$results$time.pts,res70$results$hazard,col="red")
  legend(2.5, 0.3, c("age = 50", "age = 60", "age = 70"), col = c("black", "blue", "red"),
       text.col = c("black", "blue", "red"),lty=c(1,1,1))
  # Net survival
  plot(res50$results$time.pts,res50$results$surv,type="l", ylim=c(0,1),
       ylab="Net survival",xlab="Time since diagnosis (years)",
	   main="Model 3 : Net survival")
  lines(res60$results$time.pts,res60$results$surv,col="blue")
  lines(res70$results$time.pts,res70$results$surv,col="red")
  legend(2, 0.5, c("age = 50", "age = 60", "age = 70"), col = c("black", "blue", "red"),
       text.col = c("black", "blue", "red"),lty=c(1,1,1))
  
 
 

#  =======================================================================
# |	  model 4: log(h.e) = beta0 + Bspline(time,degree=3) +               |
# |     beta*age +  Bspline(time,degree=3)*age                          |
#  =======================================================================

#-------------------------------------------------------------------------------
## 5.12 ==> Define model 4 with a formula
  
  formula2=as.formula(~bs(t1,Boundary.knots=c(0,end.time)) + agec + 
  agec: bs(t1,Boundary.knots=c(0,end.time))) 
  
 
 
#-------------------------------------------------------------------------------
## 5.13 ==> Calculate the new 10 design matrices for Gauss Legendre quadrature

	X_GL2<-list()
	
	t0=rep(0,length(data5$fu)) # initial time is zero
	t1=data5$fu # time since diagnostic
	
	for(i in 1:n.legendre){
    
	X_GL2[[i]]<-X_func((t1-t0)/2*leg1$nodes[i]+(t0+t1)/2,data=data5,formula=formula2)

	
	}
	
  
#-------------------------------------------------------------------------------
## 5.14 ==> Maximize the log-likelihood function for model 4
  
  
  # initialize parameter estimates
  beta.ini=c(log(sum(data5$dead)/sum(data5$fu)),0,0,0,0,0,0,0)
  
  # proceed to maximization (we use BFGS method because Nelder-Mead is quite unstable here)
  result=optim(beta.ini,log.lik,data=data5,formul=formula2,leg=leg1,X_GL=X_GL2, hessian=TRUE,method="BFGS",control=list(fnscale=-1))
  
  # get estimation of beta
  beta.estim=result$par
  beta.estim
  
  # get hessian matrix of the log-likelihood
  Hess=result$hessian
  Hess
  
  # get covariance matrix for beta.estim
  Covariance=solve(-Hess)
  Covariance
  sqrt(diag(Covariance))
  
  
  # get the log-likelihood at beta.estim
  log.lik(beta.estim, data=data5, formul=formula2,leg=leg1,X_GL=X_GL2)

#-------------------------------------------------------------------------------
## 5.15 ==> Check results with mexhaz package for model 4 and plot them
  
  model.check=mexhaz(formula=Surv(time=fu,event=dead)~agec+nph(agec),data5,
                      expected="rate",base="exp.bs",print.level=0,
                      bound=c(0,end.time),verbose=0)
  
  
  summary(model.check)
  
  # log-likelihood
  model.check$loglik
  
  # predictions
  surv.time=seq(0.001,end.time,length=50)
  
  res50=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agec=50-mean(data5$agediag)) ) # age=50
  res60=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agec=60-mean(data5$agediag)) ) # age=60
  res70=predict(model.check, time.pts=surv.time, 
                 data.val=data.frame(agec=70-mean(data5$agediag)) ) # age=70
  
  ##### Plot the results
  par(mfrow=c(1,2))
  # Excess hazard
  plot(res50$results$time.pts,res50$results$hazard,type="l",
       ylab="Excess hazard",xlab="Time since diagnosis (years)",ylim=range(c(res50$results$hazard,
                                                     res60$results$hazard,
                                                     res70$results$hazard)),
	   main="Model 4 : Excess hazard")
  lines(res60$results$time.pts,res60$results$hazard,col="blue")
  lines(res70$results$time.pts,res70$results$hazard,col="red")
  legend(2, 0.30, c("age = 50", "age = 60", "age = 70"), col = c("black", "blue", "red"),
       text.col = c("black", "blue", "red"),lty=c(1,1,1))
  # Net survival
  plot(res50$results$time.pts,res50$results$surv,type="l", ylim=c(0,1),
       ylab="Net survival",xlab="Time since diagnosis (years)",
	   main="Model 4 : Net survival")
  lines(res60$results$time.pts,res60$results$surv,col="blue")
  lines(res70$results$time.pts,res70$results$surv,col="red")
  legend(2, 0.5, c("age = 50", "age = 60", "age = 70"), col = c("black", "blue", "red"),
       text.col = c("black", "blue", "red"),lty=c(1,1,1))
  
  
# 		               ===========================
#   		          |	   End of Practical n°5   |
# 		               ===========================

date()


sink()



####################################################################################

