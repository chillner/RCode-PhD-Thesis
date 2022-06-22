############################################
#Simulation estimation of pi (Section 5.3)
############################################
library(multcomp)
#PWER-function:
pwerfct <- function(piv, corr, crit){
  Sigma <- matrix(c(1,corr,corr,1),nr=2)
  (1-piv[3])*(1-pnorm(crit)) + piv[3]*(1-pmvnorm(upper=rep(crit,2), corr = Sigma)[1])
}
#Correlation-function:
corrpi<- function(piv, corrfct){
  if(corrfct == 1){
    piv[3]/sqrt((piv[1]+piv[3])*(piv[2]+piv[3]))  #one treatment case
  }
  else if(corrfct == 2){
    pi1 <- piv[1]+piv[3]; pi2 <- piv[2]+piv[3]
    V1 <- sqrt((2*piv[1]/pi1^2)+(3*piv[3]/pi1^2))
    V2 <- sqrt((2*piv[2]/pi2^2)+(3*piv[3]/pi2^2))
    3*piv[3]/(2*pi1*pi2*V1*V2) #multiple treatments case
  }
}
#critical values for PWER control:
critpwer <- function(alpha, piv, corrfct = 1){
  if(piv[1]+piv[2] == 1){
    #disjoint populations
    return(qnorm(1-alpha))
  }
  else if(piv[3] == 1){
    #one population p12 = 1
    return(ifelse(corrfct == 1, qnorm(1-alpha), qmvnorm(1-alpha, tail="lower.tail", corr = matrix(c(1,0.5,0.5,1), nr=2))$quantile))
  }
  else{
    corr <- corrpi(piv=piv, corrfct=corrfct)
    Sigma <- matrix(c(1, corr, corr, 1), nr = 2, byrow = T)
    res <- function(crit){(1-piv[3])*(1-pnorm(crit)) + piv[3]*(1-pmvnorm(upper=rep(crit,2), corr = Sigma)[1])-alpha}
    return(uniroot(res, interval = c(1,5))$root)
  }
}
#Simulation:
sim_piest <- function(Nsim = 10^4, samplesize = 100, stepsize = 0.05, seed = 423, corrfct = 1){
  #sequence of true prevalences:
  pi1.true <- seq(0,1,stepsize)
  #matrix of PWER "estimates":
  z <- matrix(0, nr=length(pi1.true), nc=length(pi1.true))
  set.seed(seed)
  #Simulation:
  for(i in 1:length(pi1.true)){
    #For each entry in pi1.true create values of pi2.true such that pi1.true + pi2.true <= 1:
    jstar <- which(round(pi1.true[i]+1-sort(pi1.true, decreasing = TRUE),2) == 1)
    for(j in 1:jstar){
      #Exploit symmetry of the PWER w.r.t. pi1 and pi2:
      if(z[j,i] != 0){
        z[i,j] <- z[j,i]
      }
      else{
        #simulate observed pis und true pi-vector :
        pitrue <- c(pi1.true[i], pi1.true[j], 1-(pi1.true[i]+pi1.true[j]))
        pi.est.matrix <- rmultinom(n=Nsim, size=samplesize, prob=pitrue)
        pwer.vector <- numeric(Nsim)
        #For each generated pi compute the critical value under this estimated pi and see how much changes when plugging this critical value
        #in the PWER with the true pis:
        for(k in 1:(dim(pi.est.matrix)[2])){
          pi.est <- pi.est.matrix[,k]
          pi.est <- pi.est/sum(pi.est)
          #critical value from estimated PWER:
          crit.1 <- critpwer(alpha=0.025, piv=pi.est, corrfct = corrfct)
          if((pitrue[1]+pitrue[2] == 1)|pitrue[3] == 1 | (pi.est[1]+pi.est[2] == 1)|pi.est[3] == 1){
            pwer.vector[k] <- alpha
          }
          else{
            #"true" PWER using estimated correlations (conditional on sample sizes)
            pwer.vector[k] <- pwerfct(piv=pitrue, corr=corrpi(piv=pi.est, corrfct=corrfct), crit=crit.1)
          }
        }
        z[i,j] <- mean(pwer.vector)
      }
    }
  }
  return(z)
}

#RUN SIMULATIONS:
#equal treatment case (corrfct == 1):
z100 <- sim_piest(samplesize = 100)
z50 <- sim_piest(samplesize = 50)
#unequal treatment case (corrfct == 2):
y100 <- sim_piest(samplesize = 100, corrfct = 2)
y50 <- sim_piest(samplesize = 50, corrfct = 2)

#function to find minimum without 0 values in matrix
min0 <- function(x){y <- as.numeric(x); min(y[y>0])}

##############################
#Contour plots:
##############################
#Equal treatment case:
#n = 100:
filled.contour(x=pi1.true, y=pi1.true, z = z100, zlim = c(min0(z100), max(z100)), 
               color.palette = function(n) gray.colors(n), 
               main = "n=100", xlab = expression(pi["{1}"]), ylab = expression(pi["{2}"]),
               plot.axes = contour(x=pi1.true, y=pi1.true, z = z100, levels = seq(min0(z100), max(z100), 0.00001), add=T))

#n = 50:
filled.contour(x=pi1.true, y=pi1.true, z = z50, zlim = c(min0(z50), max(z50)), 
               color.palette = function(n) terrain.colors(n), 
               main = "n=50", xlab = expression(pi["{1}"]), ylab = expression(pi["{2}"]),
               plot.axes = contour(x=pi1.true, y=pi1.true, z = z50, levels = seq(min0(z50), max(z50), 0.00002), add=T))
#Unequal treatment case:
#n = 100:
filled.contour(x=pi1.true, y=pi1.true, z = y100, zlim = c(min0(y100), max(y100)), 
               color.palette = function(n) terrain.colors(n), 
               main = "n=100", xlab = expression(pi["{1}"]), ylab = expression(pi["{2}"]),
               plot.axes = contour(x=pi1.true, y=pi1.true, z = y100, levels = seq(min0(y100), max(y100), 0.00001), add=T))
#n = 50:
filled.contour(x=pi1.true, y=pi1.true, z = y50, zlim = c(min0(y50), max(y50)), 
               color.palette = function(n) terrain.colors(n), 
               main = "n=50", xlab = expression(pi["{1}"]), ylab = expression(pi["{2}"]),
               plot.axes = contour(x=pi1.true, y=pi1.true, z = y50, levels = seq(min0(y50), max(y50), 0.00002), add=T))

# values of 0.00001 and 0.00002 can be changed to adjust the look of the plot

#######################################################################################################
#computationally more efficient way to conduct the simulation that also calculates standard deviations
#######################################################################################################
#expand.grid unique 
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}

#Simulation:
sim_piest <- function(Nsim = 10^4, samplesize = 100, stepsize = 0.05, seed = 423, corrfct = 1){
  #all combinations of pi1, pi2 true:
  pi1.true <- seq(0,1,stepsize)
  pi.true <- expand.grid.unique(pi1.true, pi1.true, include.equals = T)
  pi.true <- pi.true[pi.true[,1]+pi.true[,2]<=1,]
  pi.true <- cbind(pi.true, 1-rowSums(pi.true))
  set.seed(seed)
  #function conducting the simulation:
  f <- function(p){
    pi.est.matrix <- rmultinom(n=Nsim, size=samplesize, prob=p)/samplesize
    #pwer.vector <- numeric(Nsim)
    g <- function(k){
      pi.est <- pi.est.matrix[,k]
      #critical value estimated PWER:
      crit.1 <- critpwer(alpha=0.025, piv=pi.est, corrfct = corrfct)
      if((p[1]+p[2] == 1)|p[3] == 1 | (pi.est[1]+pi.est[2] == 1)|pi.est[3] == 1){
        return(0.025)
      }
      else{
        #true PWER conditional on sample sizes
        return(pwerfct(piv=p, corr=corrpi(piv=pi.est, corrfct=corrfct), crit=crit.1))
      }
    }
    pwer.vector <- sapply(1:ncol(pi.est.matrix), g)
    return(c(mean(pwer.vector), sd(pwer.vector)/sqrt(Nsim)))
  }
  cbind(pi.true, t(apply(pi.true, 1, f)))
}

