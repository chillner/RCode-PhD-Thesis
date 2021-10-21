########################################
#Exit probabilities Magnusson & Turnbull
########################################
library(multcomp)
library(pracma)
library(xtable)
a <- 0.025
l1 <- qnorm(sqrt((1-a)/2))
f05 <- c(.5,.5)
f025 <- c(.25,.75)
simpsons.rule <- function(f, a, b){
  if (is.function(f) == FALSE) {
    stop('f must be a function with one parameter (variable)')
  }
  
  h <- (b - a) / 2
  x0 <- a
  x1 <- a + h
  x2 <- b
  
  s <- (h / 3) * (f(x0) + 4 * f(x1) + f(x2))
  
  return(s)
}

################################################
#Exit probabilities
################################################
#Stage 1
Zi0 <- function(l=l1, f=f025, theta=c(0,0), Imax=10){
  Delta1 <- c(f[1]*Imax/2, f[2]*Imax/2)
  pnorm(l-theta[1]*sqrt(Delta1[1]))*pnorm(l-theta[2]*sqrt(Delta1[2]))
}
Psi11 <- function(l=l1,u, f=f025, theta=c(0,0), Imax=10){
  S <- 1
  fS <- sum(f[S])
  Delta <- matrix(c(f[1]*Imax/2, f[2]*Imax/2, 
                    f[1]/fS*Imax/2, f[2]/fS*Imax/2), nr=2, byrow = TRUE)
  pnorm(l-theta[2]*sqrt(Delta[1,2]))*pnorm(theta[1]*sqrt(Delta[1,1])-u)
}
Psi12 <- function(l=l1,u, f=f025, theta=c(0,0), Imax=10){
  S <- 2
  fS <- sum(f[S])
  Delta <- matrix(c(f[1]*Imax/2, f[2]*Imax/2, 
                    f[1]/fS*Imax/2, f[2]/fS*Imax/2), nr=2, byrow = TRUE)
  pnorm(l-theta[1]*sqrt(Delta[1,1]))*pnorm(theta[2]*sqrt(Delta[1,2])-u)
}
Psi10 <- function(l=l1,u, f=f025, theta=c(0,0), Imax=10){
  S <- 1:2
  fS <- sum(f[S])
  Delta <- matrix(c(f[1]*Imax/2, f[2]*Imax/2, 
                    f[1]/fS*Imax/2, f[2]/fS*Imax/2), nr=2, byrow = TRUE)
  lv <- l*sqrt(c(Delta[1,1],Delta[1,2]))
  u1S <- u*sqrt(sum(Delta[1,S]))
  integ <- function(x11){
    1/sqrt(Delta[1,1])*dnorm((x11-theta[1]*Delta[1,1])/sqrt(Delta[1,1]))*pnorm((theta[2]*Delta[1,2]-(u1S-x11))/sqrt(Delta[1,2]))
  }
  rint <- pnorm(theta[1]*sqrt(Delta[1,1])-l)*pnorm((theta[1]*Delta[1,1]-(u1S-lv[2]))/sqrt(Delta[1,1]))
  rint+integrate(integ, lower = lv[1], upper = u1S-lv[2], abs.tol = 0L)$val
  #rint+ simpsons.rule(integ, lv[1], u1S-lv[2])
}
#find u1
u1.025 <- uniroot(function(x){Psi11(u=x, f=f025)+Psi12(u=x, f=f025)+Psi10(u=x, f=f025)-a/2}, interval = c(1,3))$root #2.278222 #2.538586
u1.05 <- uniroot(function(x){Psi11(u=x, f=f05)+Psi12(u=x, f=f05)+Psi10(u=x, f=f05)-a/2}, interval = c(1,3))$root #2.297639 #2.389095
#Stage 2
Psi21 <- function(l,u, f=f025, theta=c(0,0), Imax=10){
  S <- 1
  thetaS <- sum(f*theta)
  fS <- sum(f[S])
  Delta <- matrix(c(f[1]*Imax/2, f[2]*Imax/2, 
                    f[1]/fS*Imax/2, f[2]/fS*Imax/2), nr=2, byrow = TRUE)
  u1S <- u[1]*sqrt(sum(Delta[1,S]))
  integ <- function(y11){
    1/sqrt(Delta[1,1])*dnorm((y11-theta[1]*Delta[1,1])/sqrt(Delta[1,1]))*pnorm((y11+theta[1]*Delta[2,1]-u[2]*sqrt(Delta[1,1]+Delta[2,1]))/sqrt(Delta[2,1]))
  }
  integral <- integrate(integ, lower=l[1]*sqrt(Delta[1,1]),upper=u1S)$val
  #l[1]*sqrt(Delta[1,1])
  pnorm(l[1]-theta[2]*sqrt(Delta[1,2]))*integral
}
Zi21 <- function(l,u, f=f025, theta=c(0,0), Imax=10){
  S <- 1
  thetaS <- sum(f*theta)
  fS <- sum(f[S])
  Delta <- matrix(c(f[1]*Imax/2, f[2]*Imax/2, 
                    f[1]/fS*Imax/2, f[2]/fS*Imax/2), nr=2, byrow = TRUE)
  u1S <- u[1]*sqrt(sum(Delta[1,S]))
  integ <- function(y11){
    1/sqrt(Delta[1,1])*dnorm((y11-theta[1]*Delta[1,1])/sqrt(Delta[1,1]))*pnorm((u[2]*sqrt(Delta[1,1]+Delta[2,1])-y11-theta[1]*Delta[2,1])/sqrt(Delta[2,1]))
  }
  integral <- integrate(integ, lower=l[1]*sqrt(Delta[1,1]),upper=u1S)$val
  pnorm(l[1]-theta[2]*sqrt(Delta[1,2]))*integral
}
Psi22 <- function(l,u, f=f025, theta=c(0,0), Imax=10){
  S <- 2
  thetaS <- sum(f*theta)
  fS <- sum(f[S])
  Delta <- matrix(c(f[1]*Imax/2, f[2]*Imax/2, 
                    f[1]/fS*Imax/2, f[2]/fS*Imax/2), nr=2, byrow = TRUE)
  u1S <- u[1]*sqrt(sum(Delta[1,S]))
  integ <- function(y12){
    1/sqrt(Delta[1,2])*dnorm((y12-theta[2]*Delta[1,2])/sqrt(Delta[1,2]))*pnorm((y12+theta[2]*Delta[2,2]-u[2]*sqrt(Delta[1,2]+Delta[2,2]))/sqrt(Delta[2,2]))
  }
  integral <-integrate(integ, lower=l[1]*sqrt(Delta[1,2]),upper=u1S)$val
  #l[1]*sqrt(Delta[1,2])
  pnorm(l[1]-theta[1]*sqrt(Delta[1,1]))*integral
}
Zi22 <- function(l,u, f=f025, theta=c(0,0), Imax=10){
  S <- 2
  thetaS <- sum(f*theta)
  fS <- sum(f[S])
  Delta <- matrix(c(f[1]*Imax/2, f[2]*Imax/2, 
                    f[1]/fS*Imax/2, f[2]/fS*Imax/2), nr=2, byrow = TRUE)
  u1S <- u[1]*sqrt(sum(Delta[1,S]))
  integ <- function(y12){
    1/sqrt(Delta[1,2])*dnorm((y12-theta[2]*Delta[1,2])/sqrt(Delta[1,2]))*pnorm((u[2]*sqrt(Delta[1,2]+Delta[2,2])-y12-theta[2]*Delta[2,2])/sqrt(Delta[2,2]))
  }
  integral <-integrate(integ, lower=l[1]*sqrt(Delta[1,2]),upper=u1S)$val
  pnorm(l[1]-theta[1]*sqrt(Delta[1,1]))*integral
}
Psi20 <- function(l,u, f=f025, theta=c(0,0), Imax=10){
  S <- c(1,2)
  thetaS <- sum(f*theta)
  fS <- sum(f[S])
  Delta <- matrix(c(f[1]*Imax/2, f[2]*Imax/2, 
                    f[1]/fS*Imax/2, f[2]/fS*Imax/2), nr=2, byrow = TRUE)
  u1S <- u[1]*sqrt(sum(Delta[1,S]))
  #Density f1(y10|theta):
  f1 <- function(y10){
    integf1 <- function(x1.1){
      1/sqrt(Delta[1,1]*Delta[1,2])*dnorm((x1.1-theta[1]*Delta[1,1])/sqrt(Delta[1,1]))*dnorm((y10-x1.1-theta[2]*Delta[1,2])/sqrt(Delta[1,2]))
  }
    value <- integrate(integf1, lower=l*sqrt(Delta[1,1]), upper=y10-l*sqrt(Delta[1,2]), subdivisions = 1e4)$val
    value
  }
  #Main integral:
  integr <- function(y10){
    Vectorize(f1)(y10)*pnorm((y10+thetaS*(Delta[2,1]+Delta[2,2])-u[2]*sqrt(sum(Delta)))/sqrt(Delta[2,1]+Delta[2,2]))
  }
  integrate(integr, lower=l*(sqrt(Delta[1,1])+sqrt(Delta[1,2])), upper=u1S)$val
  #f <- function(x,y){
  #  1/sqrt(Delta[1,1]*Delta[1,2])*dnorm((y-theta[1]*Delta[1,1])/sqrt(Delta[1,1]))*dnorm((x-y-theta[2]*Delta[1,2])/sqrt(Delta[1,2]))*pnorm((x+thetaS*(Delta[2,1]+Delta[2,2])-u[2]*sqrt(sum(Delta)))/sqrt(Delta[2,1]+Delta[2,2]))
  #}
  #xmin <- l[1]*(sqrt(Delta[1,1])+sqrt(Delta[1,2])); xmax <- u1S
  #ymin <- l[1]*sqrt(Delta[1,1]); ymax <- function(x){x-l[1]*sqrt(Delta[1,2])}
  #integral2(f,xmin, xmax, ymin, ymax)$Q
}
Zi20 <- function(l,u, f=f025, theta=c(0,0), Imax=10){
  S <- c(1,2)
  thetaS <- sum(f*theta)
  fS <- sum(f[S])
  Delta <- matrix(c(f[1]*Imax/2, f[2]*Imax/2, 
                    f[1]/fS*Imax/2, f[2]/fS*Imax/2), nr=2, byrow = TRUE)
  u1S <- u[1]*sqrt(sum(Delta[1,S]))
  #Density f1(y10|theta):
  f1 <- function(y10){
    integf1 <- function(x11){
      1/sqrt(Delta[1,1]*Delta[1,2])*dnorm((x11-theta[1]*Delta[1,1])/sqrt(Delta[1,1]))*dnorm((y10-x11-theta[2]*Delta[1,2])/sqrt(Delta[1,2]))
    }
    integrate(integf1, lower=l*sqrt(Delta[1,1]), upper=y10-l*sqrt(Delta[1,2]))$val
  }
  #Main integral:
  integr <- function(z){
    Vectorize(f1)(z)*pnorm(-(z+thetaS*(Delta[2,1]+Delta[2,2])-u[2]*sqrt(sum(Delta)))/sqrt(Delta[2,1]+Delta[2,2]))
  }
  integrate(integr, lower=l*(sqrt(Delta[1,1])+sqrt(Delta[1,2])), upper=u1S)$val
  #f <- function(x,y){
  #  1/sqrt(Delta[1,1]*Delta[1,2])*dnorm((y-theta[1]*Delta[1,1])/sqrt(Delta[1,1]))*dnorm((x-y-theta[2]*Delta[1,2])/sqrt(Delta[1,2]))*pnorm((u[2]*sqrt(sum(Delta))-x-thetaS*(Delta[2,1]+Delta[2,2]))/sqrt(Delta[2,1]+Delta[2,2]))
  #}
  #xmin <- l[1]*sqrt(Delta[1,1]+Delta[1,2]); xmax <- u1S
  #ymin <- l[1]*sqrt(Delta[1,1]); ymax <- function(x){x-l[1]*sqrt(Delta[1,2])}
  #integral2(f,xmin, xmax, ymin, ymax)$Q
}
#find u2:
u2.025 <- uniroot(function(x){
  Psi21(l=l1, u=c(u1.025,x), f=f025)+Psi22(l=l1, u=c(u1.025,x), f=f025)+Psi20(l=l1, u=c(u1.025,x), f=f025)-a/2},interval = c(1,3))$root
l2.025 <- uniroot(function(x){
  Zi21(l=l1, u=c(u1.025,x), f=f025)+Zi22(l=l1, u=c(u1.025,x), f=f025)+Zi20(l=l1, u=c(u1.025,x), f=f025)-(1-a)/2},interval = c(1,3))$root
u2.05 <- uniroot(function(x){
  Psi21(l=l1, u=c(u1.05,x), f=f05)+Psi22(l=l1, u=c(u1.05,x), f=f05)+Psi20(l=l1, u=c(u1.05,x), f=f05)-a/2},interval = c(1,3))$root
#################
#Test simulation
#################
sim <- function(Nsim = 10^7, l, u, Imax=10, f=c(1/4,3/4), theta = c(0,0)){
  Delta1 <- Imax/2*f
  Delta2 <- Imax/2
  thetaS <- sum(theta*f)
  X1 <- mvrnorm(Nsim, mu=c(theta[1]*Delta1[1],theta[2]*Delta1[2]), Sigma=matrix(c(Delta1[1],0,0,Delta1[2]), nr=2, byrow=T))
  Psi11 <- mean(X1[,1]>l*sqrt(Delta1[1]) & X1[,2]<=l*sqrt(Delta1[2]) & X1[,1]>=u[1]*sqrt(Delta1[1]))
  Psi12 <- mean(X1[,2]>l*sqrt(Delta1[2]) & X1[,1]<=l*sqrt(Delta1[1]) & X1[,2]>=u[1]*sqrt(Delta1[2]))
  Psi10 <- mean(X1[,1]>l*sqrt(Delta1[1]) & X1[,2]>l*sqrt(Delta1[2]) & X1[,1]+X1[,2]>u[1]*sqrt(sum(Delta1)))
  
  #X2 <- mvrnorm(Nsim, mu=c(theta[1]*Delta2,theta[2]*Delta2), Sigma = matrix(c(Delta2,0,0,Delta2), nr=2, byrow=T))
  #Y20 <- rnorm(Nsim, mean = (f[1]*theta[1]+f[2]*theta[2])*Imax, sd = sqrt(Imax))
  #Psi21 <- mean(X1[,1]>l*sqrt(Delta1[1]) & X1[,2]<=l*sqrt(Delta1[2]) & X1[,1]<u[1]*sqrt(Delta1[1]) &
  #              X2[,1]+X1[,1]>=u[2]*sqrt(Imax/2+Imax/2*f[1]))
  #Psi22 <- mean(X1[,2]>l*sqrt(Delta1[2]) & X1[,1]<=l*sqrt(Delta1[1]) & X1[,2]<u[1]*sqrt(Delta1[2]) &
  #              X2[,2]+X1[,2]>=u[2]*sqrt(Imax/2+Imax/2*f[2]))
  #Psi20 <- mean(X1[,1]>l*sqrt(Delta1[1]) & X1[,2]>l*sqrt(Delta1[2]) & X1[,1]+X1[,2]<=u[1]*sqrt(sum(Delta1)) & 
  #              X1[,1]+X1[,2]>l*(sqrt(Delta1[1])+sqrt(Delta1[2])) & Y20>= u[2]*sqrt(Imax))
  return(c(Psi11,Psi12,Psi10))
}
#extended design probability:
Psi10sim <- function(Nsim = 10^7, l=l1, u, Imax = 10.31, f=f025, theta=c(0,0)){
  Delta1 <- Imax/2*f
  X1 <- mvrnorm(Nsim, mu=c(theta[1]*Delta1[1],theta[2]*Delta1[2]), 
                Sigma=matrix(c(Delta1[1],0,0,Delta1[2]), nr=2, byrow=T))
  mean(X1[,1]>l*sqrt(Delta1[1]) & X1[,2]>l*sqrt(Delta1[2]) & X1[,1]+X1[,2]>u[1]*sqrt(sum(Delta1)))
}
Psi10j <- function(Nsim = 10^7, l=l1, u, Imax, f=f025, theta = c(0,0)){
  Delta <- Imax/2*f
  X11 <- rnorm(Nsim, mean = theta[1]*Delta[1], sd = sqrt(Delta[1]))
  X12 <- rnorm(Nsim, mean = theta[2]*Delta[2], sd = sqrt(Delta[2]))
  Y10 <- X11+X12
  
  A11 <- mean(X11>l*sqrt(Delta[1]) & X12>l*sqrt(Delta[2]) & (Y10>=u[1]*sqrt(sum(Delta)) | X11>=u[1]*sqrt(Delta[1])))
  A12 <- mean(X11>l*sqrt(Delta[1]) & X12>l*sqrt(Delta[2]) & (Y10>=u[1]*sqrt(sum(Delta)) | X12>=u[1]*sqrt(Delta[2])))
  return(list(A11=A11, A12=A12))
}
Psi20j <- function(Nsim = 10^7, l=l1, u, Imax, f=f025, theta = c(0,0)){
  Delta <- Imax/2*f
  X11 <- rnorm(Nsim, mean = theta[1]*Delta[1], sd = sqrt(Delta[1]))
  X12 <- rnorm(Nsim, mean = theta[2]*Delta[2], sd = sqrt(Delta[2]))
  X21 <- rnorm(Nsim, mean = theta[1]*Delta[1], sd = sqrt(Delta[1]))
  X22 <- rnorm(Nsim, mean = theta[2]*Delta[2], sd = sqrt(Delta[2]))
  Y10 <- X11+X12
  Y20 <- X11+X12+X21+X22
  Y21 <- X21+X11
  Y22 <- X22+X21
  
  A21 <- mean(X11>l*sqrt(Delta[1]) & X12>l*sqrt(Delta[2]) & Y10<u[1]*sqrt(sum(Delta)) & X11<u[1]*sqrt(Delta[1]) & (Y20>=u[2]*sqrt(Imax) | Y21>=u[2]*sqrt(Imax*f[1])))
  A22 <- mean(X11>l*sqrt(Delta[1]) & X12>l*sqrt(Delta[2]) & Y10<u[1]*sqrt(sum(Delta)) & X12<u[1]*sqrt(Delta[2]) & (Y20>=u[2]*sqrt(Imax) | Y22>=u[2]*sqrt(Imax*f[2])))
  return(list(A21=A21, A22=A22))
}
PowersimR <- function(Nsim = 10^7, l=l1, u, Imax, f=f025, theta=c(0,0)){
  Delta <- Imax/2*f
  X11 <- rnorm(Nsim, mean = theta[1]*Delta[1], sd = sqrt(Delta[1]))
  X12 <- rnorm(Nsim, mean = theta[2]*Delta[2], sd = sqrt(Delta[2]))
  X21 <- rnorm(Nsim, mean = theta[1]*Delta[1], sd = sqrt(Delta[1]))
  X22 <- rnorm(Nsim, mean = theta[2]*Delta[2], sd = sqrt(Delta[2]))
  Y10 <- X11+X12
  Y20 <- X11+X12+X21+X22
  Y21 <- X21+X11
  Y22 <- X22+X21
  
  R1 <- mean(X11>l*sqrt(Delta[1]) & X12>l*sqrt(Delta[2]) & (Y10>=u[1]*sqrt(sum(Delta)) | X11>=u[1]*sqrt(Delta[1]) | X12>=u[1]*sqrt(Delta[2])))
  R2 <- mean(X11>l*sqrt(Delta[1]) & X12>l*sqrt(Delta[2]) & Y10<u[1]*sqrt(sum(Delta)) & X11<u[1]*sqrt(Delta[1]) & X12<u[1]*sqrt(Delta[2]) & 
             (Y20>=u[2]*sqrt(Imax) | Y21>=u[2]*sqrt(Imax*f[1]) | Y22>=u[2]*sqrt(Imax*f[2])))
  return(list(R1=R1, R2=R2, R=R1+R2))
}
PowersimA <- function(Nsim = 10^7, l=l1, u, Imax, f=f025, theta=c(0,0)){
  Delta <- Imax/2*f
  X11 <- rnorm(Nsim, mean = theta[1]*Delta[1], sd = sqrt(Delta[1]))
  X12 <- rnorm(Nsim, mean = theta[2]*Delta[2], sd = sqrt(Delta[2]))
  X21 <- rnorm(Nsim, mean = theta[1]*Delta[1], sd = sqrt(Delta[1]))
  X22 <- rnorm(Nsim, mean = theta[2]*Delta[2], sd = sqrt(Delta[2]))
  Y10 <- X11+X12
  Y20 <- X11+X12+X21+X22
  Y21 <- X21+X11
  Y22 <- X22+X21
  
  A <- mean(X11>l*sqrt(Delta[1]) & X12>l*sqrt(Delta[2]) & (Y20<u[2]*sqrt(Imax) | Y21<u[2]*sqrt(Imax*f[1]) | Y22<u[2]*sqrt(Imax*f[2])))
  return(list(A=A))
}

#########################################
#PWER
#########################################
PWERacc_1 <- function(l, f=f025){
  f[1]*(pnorm(-l)*pnorm(l)+pnorm(l)^2)+f[2]*(pnorm(-l)*pnorm(l)+pnorm(l)^2)
}
PWERrej_1 <- function(l, u, f=f025, theta=c(0,0), Imax = 10){
  f[1]*(Psi11(l=l, u=u, f=f, theta=theta, Imax=Imax)+Psi10(l=l, u=u, f=f, theta=theta, Imax=Imax))+
  f[2]*(Psi12(l=l, u=u, f=f, theta=theta, Imax=Imax)+Psi10(l=l, u=u, f=f, theta=theta, Imax=Imax))
}
PWERrej_2 <- function(l, u, f=f025, theta=c(0,0), Imax = 10){
  f[1]*(Psi21(l=l, u=u, f=f, theta=theta, Imax=Imax)+Psi20(l=l, u=u, f=f, theta=theta, Imax=Imax))+
  f[2]*(Psi22(l=l, u=u, f=f, theta=theta, Imax=Imax)+Psi20(l=l, u=u, f=f, theta=theta, Imax=Imax))
}
PWERacc_2 <- function(l, u, f=f025, theta=c(0,0), Imax = 10){
  f[1]*(Zi21(l=l, u=u, f=f, theta=theta, Imax=Imax)+Zi20(l=l, u=u, f=f, theta=theta, Imax=Imax))+
    f[2]*(Zi22(l=l, u=u, f=f, theta=theta, Imax=Imax)+Zi20(l=l, u=u, f=f, theta=theta, Imax=Imax))
}
#PWER extended design:
PWERext1 <- function(l, u, f=f025, theta=c(0,0), Imax = 10, seed = 42){
  set.seed(seed)
  A <- Psi10j(l=l, u=u, Imax=Imax, theta=theta, f=f)
  f[1]*(Psi11(l=l, u=u, f=f, theta=theta, Imax=Imax)+A$A11)+
  f[2]*(Psi12(l=l, u=u, f=f, theta=theta, Imax=Imax)+A$A12)
}
PWERext2 <- function(l, u, f=f025, theta=c(0,0), Imax = 10, seed = 42){
  set.seed(seed)
  A <- Psi20j(l=l, u=u, f=f, theta=theta, Imax=Imax)
  f[1]*(Psi21(l=l, u=u, f=f, theta=theta, Imax=Imax)+A$A21)+
  f[2]*(Psi22(l=l, u=u, f=f, theta=theta, Imax=Imax)+A$A22)
  
}

###critical values:
##stage 1:
u1PWER.025 <- uniroot(function(x){PWERrej_1(l=l1, u=x, f=f025)-a/2}, interval=c(1,3))$root #2.125422 #2.401283
u1PWER.05 <- uniroot(function(x){PWERrej_1(l=l1, u=x, f=f05)-a/2}, interval=c(1,3))$root #2.155697 #
#extended:
u1PWERext.025 <- uniroot(function(x){PWERext1(l=l1, u=x, f=f025)-a/2}, interval=c(1,3))$root #2.136224
##stage 2:
u2PWER.025 <- uniroot(function(x){PWERrej_2(l=l1, u=c(u1PWER.025,x), f=f025)-a/2}, interval=c(1,3))$root #1.854863 #2.205478
u2PWER.05 <- uniroot(function(x){PWERrej_2(l=l1, u=c(u1PWER.05,x), f=f05)-a/2}, interval=c(1,3))$root # 1.870743 
#extend:
u2PWERext.025 <- uniroot(function(x){PWERext2(l=l1, u=c(u1PWERext.025,x), f=f025)-a/2}, interval=c(1,3))$root #1.871693
##Power and Information:
powerfct <- function(theta_star = c(1,1), l, u, f=f025, Imax, powercrit = "GSDS", extended = FALSE){
  if(powercrit == "GSDS"){
    if(!extended){
      rejprobH0 <- Psi10sim(l=l, u=u[1], f=f, Imax=Imax, theta = theta_star)+
                   Psi20(l=l, u=u, f=f, Imax=Imax, theta = theta_star)
      rejprobH1 <- Psi11(l=l, u=u[1], f=f, Imax=Imax, theta = theta_star)+
                   Psi21(l=l, u=u, f=f, Imax=Imax, theta = theta_star)
      rejprobH2 <- Psi12(l=l, u=u[1], f=f, Imax=Imax, theta = theta_star)+
                   Psi22(l=l, u=u, f=f, Imax=Imax, theta = theta_star)
      rejprobTotal <- rejprobH0+rejprobH1+rejprobH2
  
      return(list(H0=rejprobH0, H1=rejprobH1, H2=rejprobH2, Total=rejprobTotal))
    }
    else{
      rejprob <- PowersimR(l=l, u=u, Imax=Imax, f=f, theta=theta_star)$R
      return(list(Total = rejprob))
    }
  }
}
Imaxfct <- function(beta = .1, theta_star = c(1,1), l=l1, u=c(u1.025,u2.025), f=f025, powercrit = "GSDS"){
  if(powercrit == "PWER"){
    uniroot(function(x){powerfct(theta_star=theta_star, l=l, u=u, f=f, Imax=x, powercrit = powercrit)-(1-beta)},
            interval = c(1,20))$root
  }
  else{
    uniroot(function(x){powerfct(theta_star=theta_star, l=l, u=u, f=f, Imax=x, powercrit = powercrit)$Total-(1-beta)},
            interval = c(1,20))$root
  }
}
E_IT <- function(l, u, theta, Imax, f=f025, extended = FALSE){
  if(!extended){
    1/2*Imax*(Psi10sim(l=l, u=u[1], Imax = Imax, f=f, theta=theta)+
              Psi11(l=l, u=u[1], f=f, theta=theta, Imax=Imax)+
              Psi12(l=l, u=u[1], f=f, theta=theta, Imax=Imax)+ 
              pnorm(l-theta[1]*sqrt(Imax/2*f[1]))*pnorm(l-theta[2]*sqrt(Imax/2*f[2])))+
    Imax*(Psi20(l=l, u=u, f=f, theta=theta, Imax=Imax)+ Zi20(l=l, u=u, f=f, theta=theta, Imax=Imax)+
          Psi21(l=l, u=u, f=f, theta=theta, Imax=Imax)+ Zi21(l=l, u=u, f=f, theta=theta, Imax=Imax)+
          Psi22(l=l, u=u, f=f, theta=theta, Imax=Imax)+ Zi22(l=l, u=u, f=f, theta=theta, Imax=Imax))
  }
  else{
    Rej <- PowersimR(l=l, u=u, Imax=Imax, f=f, theta=theta_star)
    Acc <- PowersimA(l=l, u=u, Imax=Imax, f=f, theta=theta_star)$A
    1/2*Imax*(Psi11(l=l, u=u[1], f=f, theta=theta, Imax=Imax)+
              Psi12(l=l, u=u[1], f=f, theta=theta, Imax=Imax)+ 
              pnorm(l-theta[1]*sqrt(Imax/2*f[1]))*pnorm(l-theta[2]*sqrt(Imax/2*f[2]))+
              Rej$R1)+
    Imax*(Psi21(l=l, u=u, f=f, theta=theta, Imax=Imax)+ Zi21(l=l, u=u, f=f, theta=theta, Imax=Imax)+
            Psi22(l=l, u=u, f=f, theta=theta, Imax=Imax)+ Zi22(l=l, u=u, f=f, theta=theta, Imax=Imax)+
            Acc+Rej$R2)
  }
}



Psi11(l=l, u=u[1], f=f, Imax=10.31, theta = theta_star)+
Psi21(l=l, u=u, f=f, Imax=10.31, theta = theta_star)



##############################
#Tables:
##############################
tablefct <- function(l=l1, uPWER, uGSDS, powercrit = "GSDS", f = f025){
  Imax_PWER1 <- Imaxfct(l=l, u=uPWER, powercrit = powercrit, f=f)
  Imax_GSDS1 <- Imaxfct(l=l, u=uGSDS, powercrit = powercrit, f=f)
  theta_configs <- matrix(c(0,0,1,1,1,0,2,0), nr=4, byrow=TRUE)
  table1 <- as.data.frame(matrix(0, nr=4, nc=13))
  colnames(table1) <- c("theta0","theta1", "theta2", "H0PWER", "H1PWER", "H2PWER", "TotalPWER", "H0GSDS", "H1GSDS", "H2GSDS","TotalGSDS", "E(I_PWER)", "E(I_GSDS)")
  table1$theta0 <- f[1]*theta_configs[,1]+f[2]*theta_configs[,2]
  table1$theta1 <- theta_configs[,1]; table1$theta2 <- theta_configs[,2]
  for(i in 1:4){
    if(i > 1){
      table1[i,4:7] <- as.numeric(powerfct(f=f, theta_star = theta_configs[i,], l=l, u=uPWER, Imax = Imax_PWER1, powercrit = powercrit)) 
      table1[i,8:11] <- as.numeric(powerfct(f=f, theta_star = theta_configs[i,], l=l, u=uGSDS, Imax = Imax_PWER1, powercrit = powercrit)) 
    }
    else{
      table1[i,4] <- Psi10sim(l=l, u=uPWER[1], f=f, theta=c(0,0), Imax=Imax_PWER1)+
                     Psi20(l=l, u=uPWER, f=f, theta=c(0,0), Imax=Imax_PWER1)
      table1[i,5] <- f[1]*Psi11(l=l, u=uPWER[1], f=f, theta=c(0,0), Imax=Imax_PWER1)+
                     f[1]*Psi21(l=l, u=uPWER, f=f, theta=c(0,0), Imax=Imax_PWER1)
      table1[i,6] <- f[2]*Psi12(l=l, u=uPWER[1], f=f, theta=c(0,0), Imax=Imax_PWER1)+
                     f[2]*Psi22(l=l, u=uPWER, f=f, theta=c(0,0), Imax=Imax_PWER1)
      table1[i,7] <- PWERrej_1(l=l, u=uPWER[1], f=f)+PWERrej_2(l=l, u=uPWER, f=f)
      
      table1[i,8] <- Psi10sim(l=l, u=uGSDS[1], f=f, theta=c(0,0), Imax=Imax_PWER1)+
                     Psi20(l=l, u=uGSDS, f=f, theta=c(0,0), Imax=Imax_PWER1)
      table1[i,9] <- f[1]*Psi11(l=l, u=uGSDS[1], f=f, theta=c(0,0), Imax=Imax_PWER1)+
                     f[1]*Psi21(l=l, u=uGSDS, f=f, theta=c(0,0), Imax=Imax_PWER1)
      table1[i,10] <- f[2]*Psi12(l=l, u=uGSDS[1], f=f, theta=c(0,0), Imax=Imax_PWER1)+
                      f[2]*Psi22(l=l, u=uGSDS, f=f, theta=c(0,0), Imax=Imax_PWER1)
      table1[i,11] <- PWERrej_1(l=l, u=uGSDS[1], f=f)+PWERrej_2(l=l, u=uGSDS, f=f)
    }
  }
  for(i in 1:4){
    table1[i,12] <- E_IT(theta = theta_configs[i,], l=l, u=uPWER, Imax = Imax_PWER1, f=f)
    table1[i,13] <- E_IT(theta = theta_configs[i,], l=l, u=uGSDS, Imax = Imax_PWER1, f=f)
  }
  return(list(table = table1, Imax_reldiff = abs((Imax_GSDS1-Imax_PWER1)/Imax_GSDS1)))
}
#Power criterion: P(reject an HS) = 1-beta under (1,1):
tab025 <- tablefct(l=l1, uPWER = c(u1PWER.025, u2PWER.025), uGSDS = c(u1.025,u2.025), powercrit = "GSDS", f=f025)
tab05 <- tablefct(l=l1, uPWER = c(u1PWER.05, u2PWER.05), uGSDS = c(u1.05,u2.05), powercrit = "GSDS", f=f05)
xtable(tab025$table, digits = 3)
xtable(tab05$table, digits = 3)

#####stage 1 table#####
tab025s1 <- matrix(0, nr = 3, nc = 8)
colnames(tab025s1) <- c("H_0", "H_1", "H_2", "sum", "H_0", "H_1", "H_2", "sum")
#PWER
tab025s1[1,1] <- Psi10sim(l=l1, u=u1PWER.025, theta = c(1,1), Imax = 8.99)
tab025s1[1,2] <- Psi11(l=l1, u=u1PWER.025, theta = c(1,1), Imax = 8.99)
tab025s1[1,3] <- Psi12(l=l1, u=u1PWER.025, theta = c(1,1), Imax = 8.99)
tab025s1[2,1] <- Psi10sim(l=l1, u=u1PWER.025, theta = c(1,0), Imax = 8.99)
tab025s1[2,2] <- Psi11(l=l1, u=u1PWER.025, theta = c(1,0), Imax = 8.99)
tab025s1[2,3] <- Psi12(l=l1, u=u1PWER.025, theta = c(1,0), Imax = 8.99)
tab025s1[3,1] <- Psi10sim(l=l1, u=u1PWER.025, theta = c(2,0), Imax = 8.99)
tab025s1[3,2] <- Psi11(l=l1, u=u1PWER.025, theta = c(2,0), Imax = 8.99)
tab025s1[3,3] <- Psi12(l=l1, u=u1PWER.025, theta = c(2,0), Imax = 8.99)
#FWER
tab025s1[1,5] <- Psi10sim(l=l1, u=u1.025, theta = c(1,1), Imax = 8.99)
tab025s1[1,6] <- Psi11(l=l1, u=u1.025, theta = c(1,1), Imax = 8.99)
tab025s1[1,7] <- Psi12(l=l1, u=u1.025, theta = c(1,1), Imax = 8.99)
tab025s1[2,5] <- Psi10sim(l=l1, u=u1.025, theta = c(1,0), Imax = 8.99)
tab025s1[2,6] <- Psi11(l=l1, u=u1.025, theta = c(1,0), Imax = 8.99)
tab025s1[2,7] <- Psi12(l=l1, u=u1.025, theta = c(1,0), Imax = 8.99)
tab025s1[3,5] <- Psi10sim(l=l1, u=u1.025, theta = c(2,0), Imax = 8.99)
tab025s1[3,6] <- Psi11(l=l1, u=u1.025, theta = c(2,0), Imax = 8.99)
tab025s1[3,7] <- Psi12(l=l1, u=u1.025, theta = c(2,0), Imax = 8.99)

xtable(tab025s1, digits = 3)

Psi10sim(u=u[1], theta = c(1,1), Imax = 8.99)+Psi20(l = l1, u=u, theta = c(1,1), Imax = 8.99)
Psi10sim(u=u[1], theta = c(1,0), Imax = 8.99)+Psi20(l = l1, u=u, theta = c(1,0), Imax = 8.99)
Psi10sim(u=u[1], theta = c(2,0), Imax = 8.99)+Psi20(l = l1, u=u, theta = c(2,0), Imax = 8.99)
