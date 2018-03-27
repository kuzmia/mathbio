#install.packages("deSolve")
#install.packages("dplyr")
#install.packages("tikzDevice")
library(deSolve)
library(dplyr)
library(tikzDevice)

###################################################################
####################HOMOGENEOUSLY-MIXED MODEL######################
###################################################################

out1 <- function(b,y=0.3){
  parameters1 <- function(b, y = 0.3){
    X<- c(b,y)
    return(X)
  }
  state1 <- c(S = 0.998, I = 0.002)
  Eq1 <- function(t, state1, parameters1){
    with(as.list(c(state1,parameters1)),{
      dS <- -b*S*I
      dI <- b*S*I - y*I
      list(c(dS,dI))
    })
  }
  times1 <- seq(0, 150, by = 0.01)
  out1 <- ode(y=state1, times = times1, 
              func = Eq1, parms = parameters1)
  return(out1)
}

####################BETA=0.5####################
b = 0.54

Min <- min(out1(b)[,2])

M <- max(out1(b)[,3]) ##maximum percentage of infected individuals (low)
a <- out1(b)[which(out1(b)[,3]> M - 0.000001 & out1(b)[,3] <  M + 0.000001),]

####################PLOTTING FOR BETA=0.5####################
par(mfrow=c(1,2))

plot(out1(b)[,1],out1(b)[,2], type = "l", ylim = c(0,1), xlab="Time (Days)", ylab = "Susceptible", col = 1, lwd = 2)

plot(out1(b)[,1],out1(b)[,3], type = "l", ylim = c(0, M+0.03), xlab = "Time (Days)", ylab = "Infected", lwd = 2, col = 1)
segments(x0=a[1,1], y0=0, x1=a[1,1], y1=a[1,3], col = 1)


###################################################################
####################HETEROGENEOUSLY-MIXED MODEL####################
###################################################################

out2 <- function(b1,b2,y=0.3){
  parameters2 <- function(b1, b2, y = 0.3){
    X<- c(b1,b2,y)
    return(X)
  }
  state2 <- c(S1 = 0.213, S2 = 0.785, I1 = 0.001, I2 = 0.001)
  Eq2 <- function(t, state2, parameters2){
    with(as.list(c(state2,parameters2)),{
      dS1 <- -b1*S1*(I1+I2)
      dS2 <- -b2*S2*(I1+I2)
      dI1 <- b1*S1*(I1+I2) - y*I1
      dI2 <- b2*S2*(I1+I2) - y*I2
      list(c(dS1, dS2, dI1, dI2))
    })
  }
  times2 <- seq(0, 150, by = 0.01)
  out2 <- ode(y=state2, times = times2, 
              func = Eq2, parms = parameters2(b1=b1,b2=b2))
  return(out2)
}

####################BETA1=0.6 & BETA2=0.4####################
b1 = 0.54*1.2
b2 = 0.54*0.8

MAXI <- max(out2(b1,b2)[,4] + out2(b1,b2)[,5])
MINS <- min(out2(b1,b2)[,2] + out2(b1,b2)[,3])

mins1 <-  min(out2(b1,b2)[,2])
mins2 <- min(out2(b1,b2)[,3])

M1 <- max(out2(b1,b2)[,4]) ##maximum percentage of infected individuals (high contact)
a1 <- out2(b1,b2)[which(out2(b1,b2)[,4]>M1-0.0000001 & out2(b1,b2)[,4] < M1 + 0.0000001),]

M2 <- max(out2(b1,b2)[,5]) ##maximum percentage of infected individuals (low contact)
a2 <- out2(b1,b2)[which(out2(b1,b2)[,5]>M2 -0.0000001 & out2(b1,b2)[,5] < M2 + 0.0000001),]

#################PLOTTING FOR BETA1=0.6 & BETA2=0.4#################

par(mfrow=c(1,2),oma = c(0,0,2,0)) #change b1 and b2 to include enhancement factor for second row of plots

b1 = 0.54*1.2
b2 = 0.54*0.8
plot(out2(b1,b2)[,1],out2(b1,b2)[,3], type = "l", ylim = c(0,1), xlab="Time (Days)", ylab = "Susceptible", col = "mediumorchid", lwd = 2)
lines(out2(b1,b2)[,1],out2(b1,b2)[,2], type="l", col = "mediumblue", lwd = 2)
b1 = 0.54*1.2*1.03618
b2 = 0.54*0.8*1.0195
lines(out2(b1,b2)[,1],out2(b1,b2)[,3], type = "l", ylim = c(0,1), xlab="Time (Days)", ylab = "Susceptible", col = "black", lwd = 2)
lines(out2(b1,b2)[,1],out2(b1,b2)[,2], type="l", col = "red", lwd = 2)
legend("topright", c("Low Contact", "Low Contact EF", "High Contact", "High Contact EF"), col = c("mediumorchid", "black", "blue", "red"), lwd = c(2,2,2,2))

b1 = 0.54*1.2
b2 = 0.54*0.8
plot(out2(b1,b2)[,1],out2(b1,b2)[,5], type = "l", ylim = c(0,max(M1,M2)+0.03), xlab = "Time (Days)", ylab = "Infected", col = "mediumorchid", lwd=2)
lines(out2(b1,b2)[,1],out2(b1,b2)[,4], type="l", col = "mediumblue", lwd = 2)
b1 = 0.54*1.2*1.03618
b2 = 0.54*0.8*1.0195
lines(out2(b1,b2)[,1],out2(b1,b2)[,5], type = "l", col = "black", lwd=2)
lines(out2(b1,b2)[,1],out2(b1,b2)[,4], type="l", col = "red", lwd = 2)
legend("topright", c("Low Contact", "Low Contact EF", "High Contact", "High Contact EF"), col = c("mediumorchid", "black", "blue", "red"), lwd = c(2,2,2,2))

#################PLOTTING FOR BOTH MODELS#################
plot(out2(b1,b2)[,1],out2(b1,b2)[,3] + out2(b1,b2)[,2], type = "l", ylim = c(0,1), xlab="Time (Days)", ylab = "Susceptible", col = "cyan3", lwd = 2)
legend("topright", c("Mixed", "Single"), col = c("cyan3", 1), lwd = c(2,2))

plot(out2(b1,b2)[,1],out2(b1,b2)[,5] + out2(b1,b2)[,4], type = "l", ylim = c(0,max(max(M1,M2)+0.03,M+0.03)), xlab = "Time (Days)", ylab = "Infected", col = "cyan3", lwd=2)
legend("topright", c("Mixed", "Single"), col = c("cyan3", 1), lwd = c(2,2))


###################################################################
###########################FINAL SIZE#############################
###################################################################
rm(list = ls(all = TRUE)) 
SIR <- function(time,state,parameters) {
  with(as.list(c(state,parameters)),{
    dS1 <- -b1*S1*(I1+I2)
    dS2 <- -b2*S2*(I1+I2)
    dI1 <- b1*S1*(I1+I2) - gamma*I1
    dI2 <- b2*S2*(I1+I2) - gamma*I2
    dR <- gamma*(I1+I2)
    return(list(c(dS1,dS2,dI1,dI2,dR)))
  })
}
init <-  c(S1 = 0.213, S2 = 0.785, I1 = 0.001, I2 = 0.001, R=0) #initial conditions for odes
time <- seq(0,150,by=0.01) #time period

betavec1 <- seq(0.1,3,by=0.1)

res1 <- vector(length(betavec1),mode="list")

for (k in seq_along(betavec1)){ #range of values for beta
  res1[[k]] <- ode(y=init,times=time,func=SIR,
                   parms=c(b1=betavec1[k], gamma=0.3, b2=0.54))
}

names(res1) <- betavec1  ## to get beta value incorporated in results
dd1 <- dplyr::bind_rows(lapply(res1,as.data.frame),.id="beta1")
dd1$b1 <- as.numeric(dd1$beta1)

newdata1 <-dd1[ which(dd1$time==150), ]

##final size:R_0= log(1-z)/-z
newdata1$R_0 <- log(1-newdata1$R)/-newdata1$R


##changing beta_2
betavec2 <- seq(0.1,3,by=0.1)

res2 <- vector(length(betavec2),mode="list")

for (k in seq_along(betavec2)){ #range of values for beta
  res2[[k]] <- ode(y=init,times=time,func=SIR,
                   parms=c(b2=betavec2[k], gamma=0.3, b1=0.54))
}

names(res2) <- betavec2  ## to get beta value incorporated in results
dd2 <- dplyr::bind_rows(lapply(res2,as.data.frame),.id="beta2")
dd2$b2 <- as.numeric(dd2$beta2)

newdata2 <-dd2[ which(dd2$time==150), ]

##final size:R_0= log(1-z)/-z
newdata2$R_0 <- log(1-newdata2$R)/-newdata2$R

plot(newdata2$beta2, newdata2$R_0, type = "l", col = "blue", xlab = "\beta", lwd = 2)
lines(newdata1$beta1, newdata1$R_0, type = "l", col = "red", lwd = 2)
legend("bottomright", c("\beta_1 = 0.54","\beta_2 = 0.54"), col = c("red", "blue"),lwd = c(2,2))
