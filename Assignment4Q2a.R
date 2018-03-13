SI.gillespie <- function(beta, N, I0, tmax){
  t0 <- 0
  times =(t0:tmax)
  x <-c(S=N-I0, I=I0, t=t0)

  res <- matrix(nrow=length(t0:tmax),ncol=length(x),
              dimnames = list(times,names(x))) #matrix to store values
  
  for (i in 1:(tmax+1)){
    res[i,] <- x
    rate <- with(as.list(x), beta*S*I) ## calculate current rate
    if(rate<=0) break #rate != 0; t_next would return NaN
    t_next <- rexp(1,rate)  #time to next event
    t0 <- t0 + t_next # update time
    x <- x+c(-1*rate*t_next, 1*rate*t_next, t0) #updates x <- c(S, I, t)
    if(x[1]<0) #if S is for some reason negative, change it to its previous positive
      x <- res[i,]
  }
  cbind(res[,3],res[,2]) #returns cbind(t, I)
}

## N=32
plot(0,0,xlim=c(0,10),ylim=c(0,32),
     type="n",xlab="Time (t)",ylab="Prevalence (I)",las=1)
for(i in 1:30){
  G.SI <- SI.gillespie(beta=1, N=32, I0=1, tmax=80)
  lines(G.SI, col=i)
}
N <- 32
beta <- 1
I0 <- 1
It <- I0*exp(N*beta*G.SI[,1])/(1+(I0/N)*(exp(N*beta*G.SI[,1])-1))
lines(G.SI[,1],It,lwd=3)

## N=100
plot(0,0,xlim=c(0,10),ylim=c(0,100),
     type="n",xlab="Time (t)",ylab="Prevalence (I)",las=1)
for(i in 1:30){
  G.SI <- SI.gillespie(beta=1, N=100, I0=1, tmax=300)
  lines(G.SI, col=i)
}
N <- 100
G.SI <- SI.gillespie(beta=1, N=100, I0=1, tmax=80)
It <- I0*exp(N*beta*G.SI[,1])/(1+(I0/N)*(exp(N*beta*G.SI[,1])-1))
lines(G.SI[,1],It,lwd=3)


## N=1000
plot(0,0,xlim=c(0,10),ylim=c(0,1000),
     type="n",xlab="Time (t)",ylab="Prevalence (I)",las=1)
for(i in 1:30){
  G.SI <- SI.gillespie(beta=1, N=1000, I0=1, tmax=1000)
  lines(G.SI, col=i)
}
N <- 1000
It <- I0*exp(N*beta*G.SI[,1])/(1+(I0/N)*(exp(N*beta*G.SI[,1])-1))
lines(G.SI[,1],It,lwd=3)

## N=10,000
plot(0,0,xlim=c(0,10),ylim=c(0,10000),
     type="n",xlab="Time (t)",ylab="Prevalence (I)",las=1)
for(i in 1:30){
  G.SI <- SI.gillespie(beta=1, N=10000, I0=1, tmax=10000)
  lines(G.SI, col=i)
}
N <- 10000
It <- I0*exp(N*beta*G.SI[,1])/(1+(I0/N)*(exp(N*beta*G.SI[,1])-1))
lines(G.SI[,1],It,lwd=3)
