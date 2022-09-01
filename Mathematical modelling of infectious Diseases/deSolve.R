#ODE solution in R
#deSolve package
##################################################

#The Lorenz function
#######################

#Parameters
parameters <- c(a = -8/3, b = -10, c = 28)

#Initial conditions
state <- c(X = 1, Y = 1, Z = 1)

#Model

Lorenz<-function(t, state, parameters) {
   with(as.list(c(state, parameters)),{
   # rate of change
   dX <- a*X + Y*Z
   dY <- b * (Y-Z)
   dZ <- -X*Y + c*Y - Z
  
   # return the rate of change
   list(c(dX, dY, dZ))
    }) # end with(as.list...
}

#Time specification
######################

#We run the model for 100 days, and give output at 0.01 daily intervals. R's function seq() creates the time sequence.
times <-seq(0,100,by=0.01)


#Model Intergration
########################

#install.packages('deSolve')
require(deSolve)
out <- as.data.frame(ode(y=state,times=times,func=Lorenz,parms=parameters))
head(out)


#Plotting Results

par(mfrow=c(2,2), oma=c(0,0,3,0))
plot (times,out$X ,type="l",main="X", xlab="time", ylab="-")
plot (times,out$Y ,type="l",main="Y", xlab="time", ylab="-")
plot (out$X,out$Y, pch=".")
plot (out$X,out$Z, pch=".")
mtext(outer=TRUE,side=3,"Lorenz model",cex=1.5)

?plot   # to view the plotting parameter (pch, col, bg, cex and lwd) options.

#Different ODE solvers. Time elapsed to solution covergence.
########################################

print(system.time(out1 <- rk4 (state, times, Lorenz, parameters)))
print(system.time(out2 <- lsode (state, times, Lorenz, parameters)))
print(system.time(out <- lsoda (state, times, Lorenz, parameters)))
print(system.time(out <- lsodes(state, times, Lorenz, parameters)))
print(system.time(out <- daspk (state, times, Lorenz, parameters)))
print(system.time(out <- vode (state, times, Lorenz, parameters)))












