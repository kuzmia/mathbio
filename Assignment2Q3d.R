library(deSolve)

## Vector Field for SIR model
SIR.vector.field <- function(t, vars=c(S,I,R), parms=c(R_0,gamma)) {
  with(as.list(c(parms, vars)), {
    dS <- -gamma*R_0*S*I # dS/dt
    dI <- gamma*R_0*S*I - gamma*I # dI/dt
    dR <- gamma*I #dR/dt
    vec.fld <- c(dS=dS, dI=dI, dR=dR)
    return(list(vec.fld)) # ode() requires a list
  })
}

## Plot solutions of the SIR model
tmax <- 70 # end time for numerical integration of the ODE
## draw box for plot:
plot(0,0,xlim=c(0,tmax),ylim=c(-10,2),
     type="n",xlab="Time (t)",ylab="Log (I)",las=1)

draw.soln <- function(ic=c(S=1,I=0), tmax=1,
                      times=seq(0,tmax,by=tmax/500),
                      func, parms, legend,... ) {
  soln <- ode(ic, times, func, parms)
  lines(times+14, log(soln[,"I"]), col=Rep_nums[i], lwd=3,... )
}#translate times to point where pim data begins exponential growth

##Initial conditions:
I0 <- 0.001
S0 <- 1 - I0
R0 <- 1 - I0 - S0

##Draw solutions for several values of parameter R_0:
Rep_nums <- c(2.1)
gamma <- 1/4.1

for (i in 1:length(Rep_nums)) {
  draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
            func=SIR.vector.field,
            parms=c(R_0=Rep_nums[i],gamma),
            lty=i #Different line style for each solution
  )
}
##Legend for Ro
legend("topright",legend=Rep_nums, title = expression(paste(italic("Ro ="))))

##Legend for gamma
legend("topleft", legend=(1/gamma), title = expression(paste(gamma," =")))

## Plots log of mortality data and translates down 8 units
points(df$date,log(df$pim)-8)