library(deSolve)
library(MASS)

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
tmax <- 70 # end time for numerical integration of the ODE

## draw box for plot:
plot(0,0,xlim=c(0,tmax),ylim=c(0,0.4),
     type="n",xlab="Tao",ylab="Prevelance (I)",las=1, main="Undamped Oscillation with e-folding of Amplitude")

draw.soln <- function(ic=c(S=1,I=0), tmax=1,
                      times=seq(0,tmax,by=tmax/500),
                      period_times=seq(0,4*pi/(sqrt(4*(Rep_nums*epsilon-epsilon)-(Rep_nums*epsilon)^2))),
                      func, parms, ... ) {
  soln <- ode(ic, times, func, parms)
  lines(times, soln[,"I"], col=Rep_nums[i], lwd=3,... ) #plots prevalence
  lines(times, exp((-Rep_nums*epsilon/2)*(times+13.5))+0.048, col=4, lwd=2,...) #translated plot fo exponential decay
  abline(v=7.1,col=1,lty=2) #plots left dashed vertical line
  abline(v=4*pi/(sqrt(4*(Rep_nums*epsilon-epsilon)-(Rep_nums*epsilon)^2))+7.1, col=1, lty=2) #plots right dashed vertical line
}

## Initial conditions:
I0 <- 0.001
S0 <- 1 - I0
R0 <- 0

## Draw solutions for several values of parameter R_0:
Rep_nums <- c(2)
epsilon <- 0.1
for (i in 1:length(Rep_nums)) {
  draw.soln(ic=c(S=S0,I=I0,R=R0), tmax=tmax,
            func=SIR.vector.field,
            parms=c(R_0=Rep_nums[i], epsilon),
            lty=i #Different line style for each solution
            
  )
}


##Legend for R_0:
legend("topright",legend=Rep_nums,col=Rep_nums,lty=1:6, title = expression(paste(italic("Ro"))))

##Legend for gamma:
legend("right",legend=(epsilon), title = (expression(paste(epsilon," ="))))
