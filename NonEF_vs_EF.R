library(deSolve)

## Vector Field for SIR model
SIR.vector.field <- function(t, vars=c(S1,S2,I1,I2,R1,R2), parms=c(beta,gamma,p1,p2,c11,c12,c21,c22,ef1,ef2)) {
  with(as.list(c(parms, vars)), {
    dS1 <- -beta*S1*(c11*I1*ef1/p1+c12*I2*ef2/p2)
    dS2 <- -beta*S2*(c21*I1*ef1/p1+c22*I2*ef2/p2)
    dI1 <- beta*S1*(c11*I1*ef1/p1+c12*I2*ef2/p2) - gamma*I1
    dI2 <- beta*S2*(c21*I1*ef1/p1+c22*I2*ef2/p2) - gamma*I2
    dR1 <- gamma*I1
    dR2 <- gamma*I2
    vec.fld <- c(dS1=dS1, dS2=dS2, dI1=dI1, dI2=dI2, dR1=dR1, dR2=dR2)
    return(list(vec.fld))
  })
}

## Plot solutions of the SIR model
tmax <- 250 # end time for numerical integration of the ODE
## draw box for plot:
plot(0,0,xlim=c(0,tmax),ylim=c(0,0.03),
     type="n",xlab="Time (days)",ylab="Prevalence",las=1)

draw.soln <- function(ic=c(S1,S2,I1,I2,R1,R2), tmax=1,
                      times=seq(0,tmax,by=tmax/500),
                      func, parms, legend, ...) {
  soln <- ode(ic, times, func, parms)
  
  lines(times, soln[,"I1"]+soln[,"I2"], col = c("cyan","green")[i], lwd=2,... )
  
}#translate times to point where pim data begins exponential growth

##Initial conditions:

p1 <- 0.216
p2 <- 0.784

I10 <- 0.0001
I20 <- 0.0001
S10 <- p1 - I10
S20 <- p2 - I20
R10 <- 0
R20 <- 0 

c11 <- 10.6
c12 <- 4.6
c21 <- 1.3
c22 <- 10.4

ef1 <- c(1, 1.03618)
ef2 <- c(1, 1.0195)
##Draw solutions for several values of parameter R_0:
beta <- 0.0291
gamma <- 0.3

for(i in 1:2){  
  draw.soln(ic=c(S1=S10,S2=S20,I1=I10,I2=I20,R1=R10,R2=R20), tmax=tmax,
            func=SIR.vector.field,
            parms=c(beta,gamma,p1,p2,c11,c12,c21,c22,ef1=ef1[i],ef2=ef2[i])
  )
}
legend("topright", c("Non-EF Total", "EF Total"), col = c("cyan", "green"), lwd = c(2,2))