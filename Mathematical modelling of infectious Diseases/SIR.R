#example1(No interaction between Susceptibles and infectives)
parameters<-c(lambda=0.2,v=36.5)
parameters
state<-c(x=4999,y=1,z=0)
state
SIR<-function(t,state,parameters)
{
  with(as.list(c(state,parameters)),
       {
         dx<--lambda*x
         dy<-lambda*x-v*y
         dz<-v*y
         list(c(dx,dy,dz))
         
       })
}
times<-seq(0,40,by=0.01)
times
require(deSolve)
out<-as.data.frame(ode(y=state,times =times,func = SIR,parms = parameters ))
plot(times,out$x,type="l",lwd=3)
lines(times,out$y,col="red",lwd=3)
lines(times,out$z,col="blue",lwd=3)
legend(30,3500,c("S","I","R"),lty ="l",lwd=3,col = c("black","red","blue"))
#example2(interaction between suseptibles and infectives)
parameters<-c(beta=0.0085,v=36.5)
parameters
state<-c(x=4999,y=1,z=0)
state
SIR<-function(t,state,parameters)
{
  with(as.list(c(state,parameters)),
       {
         dx<--beta*x*y
         dy<-beta*y*x-v*y
         dz<-v*y
         list(c(dx,dy,dz))
         
       })
}
times<-seq(0,40,by=0.01)
times
require(deSolve)
out<-as.data.frame(ode(y=state,times =times,func = SIR,parms = parameters ))
plot(times,out$x,type="l",lwd=3)
lines(times,out$y,col="red",lwd=3)
lines(times,out$z,col="blue",lwd=3)
legend(30,3500,c("S","I","R"),lty ="l",lwd=2,col = c("black","red","blue"))
#example3(Assuming demographics)
parameters<-c(mu=0.017,beta=0.0003,sig=0.03,v=0.1)
parameters
state<-c(x=999,y=1,z=0)
state
SIR<-function(t,state,parameters)
{
  with(as.list(c(state,parameters)),
       {
         dx<-1000*mu-beta*x*y-mu*x
         dy<-beta*y*x-v*y-mu*y-sig*y
         dz<-v*y-mu*z
         list(c(dx,dy,dz))
         
       }) 
}
times<-seq(0,50,by=0.01)
times
require(deSolve)
out<-as.data.frame(ode(y=state,times =times,func = SIR,parms = parameters ))
plot(times,out$x,type="l",lwd=3 ,ylab="Freq",ylim=c(0,1000),xlab="Times in years")
lines(times,out$y,col="red",lwd=3)
lines(times,out$z,col="blue",lwd=3)

legend(30,3500,c("S","I","R"),lty ="l",lwd=2,col = c("black","red","blue"))

