library(deSolve)
library(MASS)
library(phaseR)

## Vector Field for SIR model
SIR.vector.field <- function(t, vars=c(S,I,R), parms=c(R_0, epsilon)) {
  with(as.list(c(parms, vars)), {
    dS <- epsilon - R_0*S*I - epsilon*S # dS/dtao
    dI <- R_0*S*I - (1-epsilon)*I - epsilon*I # dI/dtao
    dR <- (1-epsilon)*I - epsilon*R #dR/dtao
    vec.fld <- c(dS=dS, dI=dI, dR=dR)
    return(list(vec.fld))
  })
}

## Plot solutions of the SIR model
tmax <- 71 # end time for numerical integration of the ODE

draw.soln <- function(ic=c(S=1,I=0), tmax=1,
                      times=seq(0,tmax,by=tmax/500),
                      func, parms,...){
  soln <- ode(ic, times, func, parms)
  lines(soln[,"S"], soln[,"I"], col=Rep_nums[i], lwd=3,... )
}#plots prevalence vs. susceptibles phase diagram

## Initial conditions:
I0 <- 0.001
S0 <- 1 - I0
R0 <- 0

par(mfrow=c(2,2))

## draw box for plot:
plot(0,0,xlim=c(0,1),ylim=c(0,1),
     type="n",xlab="Susceptibles (S)",ylab="Prevelance (I)",las=1)

## Draw phase diagrams for different R_0 values:
Rep_nums <- c(0.9)
epsilon <- 8/9
for (i in 1:length(Rep_nums)) {
  draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
            func=SIR.vector.field,
            parms=c(R_0=Rep_nums[i], epsilon),
            lty=i #Different line style for each solution
  )
}
##Legend for R_0:
text(0.95,0.9,expression(paste(italic("Ro" == 0.9))))
##Legend for gamma:
text(0.95,0.6, (expression(epsilon == frac(8,9))))

## Draw solutions for several values of parameter R_0:
## draw box for plot:
plot(0,0,xlim=c(0,1),ylim=c(0,1),
     type="n",xlab="Susceptibles (S)",ylab="Prevelance (I)",las=1)

Rep_nums <- c(1.4)
epsilon <- 8/9
for (i in 1:length(Rep_nums)) {
  draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
            func=SIR.vector.field,
            parms=c(R_0=Rep_nums[i], epsilon),
            lty=i #Different line style for each solution
  )
}
##Legend for R_0:
text(0.95,0.9,expression(paste(italic("Ro" == 1.4))))
##Legend for gamma:
text(0.95,0.6, (expression(epsilon == frac(8,9))))

## Draw solutions for several values of parameter R_0:
## draw box for plot:
plot(0,0,xlim=c(0,1),ylim=c(0,1),
     type="n",xlab="Susceptibles (S)",ylab="Prevelance (I)",las=1)

Rep_nums <- c(2.5)
epsilon <- 8/9
for (i in 1:length(Rep_nums)) {
  draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
            func=SIR.vector.field,
            parms=c(R_0=Rep_nums[i], epsilon),
            lty=i #Different line style for each solution
  )
}
##Legend for R_0:
text(0.95,0.9,expression(paste(italic("Ro" == 2.5))))
##Legend for gamma:
text(0.95,0.6, (expression(epsilon == frac(8,9))))

## Draw solutions for several values of parameter R_0:
## draw box for plot:
plot(0,0,xlim=c(0,1),ylim=c(0,1),
     type="n",xlab="Susceptibles (S)",ylab="Prevelance (I)",las=1)
Rep_nums <- c(5)
epsilon <- 8/9
for (i in 1:length(Rep_nums)) {
  draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
            func=SIR.vector.field,
            parms=c(R_0=Rep_nums[i], epsilon),
            lty=i #Different line style for each solution
  )
}
##Legend for R_0:
text(0.95,0.9,expression(paste(italic("Ro" == 5))))
##Legend for gamma:
text(0.95,0.6, (expression(epsilon == frac(8,9))))