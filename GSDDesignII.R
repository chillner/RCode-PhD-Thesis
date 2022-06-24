#################################################
#Design I - H_1 = {h_1, h_2}, H_2 = {H_{1,2}}
#################################################
library(multcomp)
corrmat <- function(prev, w){
  piv <- prev[1:2]+prev[3]
  rho12 <- prev[3]/sqrt(piv[1]*piv[2])
  rho1_12 <- sqrt(prev[3]/piv[1])
  rho2_12 <- sqrt(prev[3]/piv[2])
  sigma <- matrix(c(1, rho12, w*rho1_12,
                    0, 1, w*rho2_12,
                    0, 0, 1), nr = 3, byrow = T)
  sigma <- sigma + t(sigma) - diag(3)
}
#PWER-function:
#crit matrix:
#  P1  P2
#1 c11 c12
#2 c21 c22
pwer <- function(crit, prev, nu = c(0,0,0), sigma){
  names(prev) <- c("P1", "P2", "P12")
  if(length(crit)==2){
    crit <- c(crit[1],crit[1],crit[2])
  }
  #PWER:
  P <- numeric(3)
  P[1] <- 1-pnorm(crit[1], mean = nu[1], sd = 1)
  P[2] <- 1-pnorm(crit[2], mean = nu[2], sd = 1)
  P[3] <- 1-pmvnorm(upper= crit, mean = nu, sigma = sigma, algorithm = Miwa())[1]
  PWER <- sum(prev*P)
  #PWER^(k):
  PWER1 <- prev[1]*P[1]+ prev[2]*P[2]+
           prev[3]*(1-pmvnorm(upper = crit[1:2], mean = nu[1:2],
                       sigma = sigma[1:2,1:2], algorithm = Miwa())[1])
  
  PWER2 <- prev[3]*(pmvnorm(lower=c(-Inf,-Inf,crit[3]),upper=c(crit[1:2],Inf),
                     mean=nu, sigma=sigma, algorithm = Miwa())[1])
  #results list:
  res <- list(PWER = PWER, P1 = P[1], P2 = P[2], P12 = P[3], PWER_stage1 = PWER1, PWER_stage2 = PWER2)
  return(res)
}

#Wang & Tsiatis critical values
critWT <- function(p, prev, sigma, alpha = 0.025, tau){
  f <- function(cr){
    crit <- c(cr,cr,cr*(1/tau)^(p-0.5))
    pwer(crit=crit,prev=prev, nu = rep(0,3), sigma=sigma)$PWER-alpha
  }
  cr_const <- uniroot(f, interval = c(qnorm(1-alpha),5))$root
  res <- c(cr_const,cr_const,cr_const*(1/tau)^(p-0.5))
  names(res) <- c("c_P1", "c_P2", "c_P12")
  return(res)
}
#error spending critical values
critES <- function(p, prev, sigma, alpha=0.025, tau, aspend){
    if(aspend == "Kim"){
      a <- function(l, par = p){alpha*l^par}
    }
    else if(aspend == "Hwang"){
      a <- function(l, par = p){
        alpha*ifelse(par != 0, (1-exp(-par*l))/(1-exp(-par)), l)
      }
    }
    #first stage critical value
    f1 <- function(cr){
      crit = c(cr,cr,0)
      pwer(crit =  crit, prev = prev, sigma = sigma, nu = rep(0,3))$PWER_stage1 - a(l=tau[1])
    }
    cr1 <- uniroot(f1, interval = c(qnorm(1-alpha), 5))$root
    #second stage critical value
    f2 <- function(cr){
      pwer(crit =  c(cr1,cr1,cr), prev = prev, sigma = sigma, nu = rep(0,3))$PWER_stage2 - 
          (alpha-a(l=tau[1]))
    }
    cr2 <- uniroot(f2, interval = c(1, 5))$root
    res <- c(cr1,cr1,cr2)
    names(res) <-  c("c_P1", "c_P2", "c_P12")
    return(res)
}

#POWER
pwp <- function(p, prev, delta, N, gamma = 1, spendingfct = NULL, tau, alpha = 0.025){
  n1 <- prev*N
  n2 <- n1[3]*gamma
  nu <- delta*sqrt(c(n1[1:2]+n1[3], n1[3]+n2))
  w <- sqrt(n1[3]/(n1[3]+n2))
  sigma <- corrmat(prev = prev, w=w)
  #critical values:
  if(is.null(spendingfct)){
    crit <- critWT(p=p, prev=prev, sigma=sigma, alpha=alpha, tau = tau)
  }
  else{
    crit <- critES(p=p, prev=prev, sigma=sigma, alpha=alpha, tau=tau, aspend = spendingfct)
  }
  ##power value
  #denominator in PWP
  den <- prev[1]*ifelse(delta[1]>0,1,0) + prev[2]*ifelse(delta[2]>0,1,0) + 
       prev[3]*ifelse(delta[1]>0 | delta[2]>0 | delta[3]>0,1,0)
  #numerator in PWP
  upper12 <- c(ifelse(delta[1]>0,1,Inf)*crit[1], ifelse(delta[2]>0,1,Inf)*crit[2], 
  (ifelse(delta[1]>0 | delta[2]>0, 1,Inf)+ifelse(delta[1]<=0 & delta[2]<=0 & delta[3]>0,1,Inf))*crit[3])
  
  num <- prev[1]*(1-pnorm(crit[1], mean = nu[1]))*ifelse(delta[1]>0,1,0)+
    prev[2]*(1-pnorm(crit[2], mean = nu[2]))*ifelse(delta[2]>0,1,0)+
    prev[3]*(1-pmvnorm(lower=rep(-Inf,3), upper=upper12, mean=nu, sigma=sigma)[1])
  return(ifelse(den>0,num/den,0))
}

pow1 <- function(p, prev, delta, N, gamma = 1, spendingfct = NULL, tau, alpha = 0.025){
  n1 <- prev*N
  n2 <- n1[3]*gamma
  nu <- delta*sqrt(c(n1[1:2]+n1[3], n1[3]+n2))
  w <- sqrt(n1[3]/(n1[3]+n2))
  sigma <- corrmat(prev = prev, w=w)
  if(is.null(spendingfct)){
    crit <- critWT(p=p, prev=prev, sigma=sigma, alpha=alpha, tau = tau)
  }
  else{
    crit <- critES(p=p, prev=prev, sigma=sigma, alpha=alpha, tau=tau, aspend = spendingfct)
  }
  pos <- which(nu > 0)
  1-pmvnorm(upper=crit[pos], mean = nu[pos], sigma = sigma[pos,pos], algorithm = Miwa())[1]
}

power <- function(p, prev, delta, N, gamma = 1, spendingfct = NULL, tau, alpha = 0.025, powertype){
  if(powertype == "pwp"){
    pwp(p=p, prev=prev, delta=delta, N=N, gamma=gamma, 
        spendingfct=spendingfct, tau=tau, alpha=alpha)
  }
  else if(powertype == "pow1"){
    pow1(p=p, prev=prev, delta=delta, N=N, gamma=gamma, 
         spendingfct=spendingfct, tau=tau, alpha=alpha)
  }
}

#Function for finding N:
Nroot <- function(p, prev, delta, gamma=1, alpha = 0.025, 
                  beta = 0.2, spendingfct = NULL, tau, powertype,
                  search.interval = c(1,2000)){
  fzero <- function(x){
    power(N = x, p=p, prev=prev, delta = delta, gamma=gamma, powertype = powertype,
          spendingfct = spendingfct, tau=tau, alpha = alpha)-(1-beta)
  }
  uniroot(fzero, interval = search.interval)$root
  
}

ASNroot <-function(p, prev, delta, gamma=1, rho, alpha = 0.025, 
                   beta = 0.2, spendingfct = NULL, tau, powertype,
                   search.interval = c(1,2000)){
  #find N such that power = 1-beta:
  fNzero <- function(x){
    power(N=x, p=p, prev=prev, delta=delta, gamma=gamma, alpha=alpha, 
          spendingfct = spendingfct, tau=tau, powertype = powertype) - (1-beta)
  }
  N0 <- uniroot(fNzero, interval = search.interval)$root
  #critical values:
  n1 <- prev*N0
  n2 <- n1[3]*gamma
  nu <- delta*sqrt(c(n1[1:2]+n1[3], n1[3]+n2))
  w <- sqrt(n1[3]/(n1[3]+n2))
  sigma <- corrmat(prev=prev, w=w)
  if(is.null(spendingfct)){
    crit <- critWT(p=p, prev=prev, sigma=sigma, alpha=alpha, tau = tau)
  }
  else{
    crit <- critES(p=p, prev=prev, sigma=sigma, alpha=alpha, tau=tau, aspend = spendingfct)
  }
  #ASN
  N0 + n2*(1-pmvnorm(upper=crit[1:2], mean = nu[1:2], sigma = sigma[1:2,1:2])[1])
}

opt_par <- function(prev, delta, tau, gamma = 1, spendingfct = NULL, alpha = 0.025,
                    beta = 0.2, powertype, search.inteval.N, search.interval.p,
                    opt.criterion){
  if(opt.criterion == "Nmin"){
    fzero <- function(x){
      Nroot(p=x, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, tau=tau,
            spendingfct = spendingfct, powertype = powertype, search.interval = search.inteval.N)
    }
  }
  if(opt.criterion == "ASNmin"){
    fzero <- function(x){
      ASNroot(p=x, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, tau=tau,
              spendingfct = spendingfct, powertype = powertype, search.interval = search.inteval.N)
    }
  }
  optimize(fzero, interval = search.interval.p)$min
}



##########
#tests
##########
prev = c(.4,.4,.2)
w = 1/sqrt(2)
sigma = corrmat(prev=prev, w=w)
crit <- critWT(p=0.5, prev=prev, sigma=sigma, tau = 1/2)
crit2 <- critES(p=0.5, prev=prev, sigma=sigma, tau = 1/2, aspend = "Hwang")



###################
#Table Section 5
###################
tabD2 <- as.data.frame(matrix(0, nr = 6, nc = 9))
colnames(tabD2) <- c("pi{1}","pi{2}", "pi{1,2}", "c0", "c05", "N_PWP0", "N_PWP05", "N_Pow1_0", "N_Pow1_05")
pr <- rbind(c(0.3,0.3,0.4),
            c(0.35,0.35,0.3),
            c(0.4, 0.4, 0.2),
            c(0.4,0.2,0.4),
            c(0.4,0.3,0.3),
            c(0.6,0.2,0.2))
tabD2[,1:3] <- pr 
w <- 1/sqrt(2)
delta <- c(0.3,0.3,0.3)
for(i in 1:6){
  pri <- pr[i,]
  corr <- corrmat(prev = pri, w = w)
  tabD2[i, 4] <- critWT(p=0, prev=pri, sigma = corr, tau =0.5)[1]
  tabD2[i, 5] <- critWT(p=.5, prev=pri, sigma = corr, tau =0.5)[1]
  tabD2[i, 6] <- Nroot(p=0, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau =0.5, powertype = "pwp")
  tabD2[i, 7] <- Nroot(p=.5, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau =0.5, powertype = "pwp")
  tabD2[i, 8] <- Nroot(p=0, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau =0.5, powertype = "pow1")
  tabD2[i, 9] <- Nroot(p=.5, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau =0.5, powertype = "pow1")
}
tabD2[,4:5] <- round(tabD2[,4:5], 3)
tabD2[,6:9] <- ceiling(tabD2[,6:9])
tabD2

xtable(x=tabD2, caption = "A caption", label = "tab: tabD1", digits = 3)


