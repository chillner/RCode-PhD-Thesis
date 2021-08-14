#######################################################################
#Design I - H_1 = {h_1, h_2}, H_2 = {H_1, H_2, H_{1}, H_{2}, H_{1,2}}
#######################################################################
library(multcomp)
#PWER-function:
#crit vector:
#(c_1^, c_1^1, c_1^2, c_2^2, c_{1}^2, c_{2}^2, c_{1,2}^2)
#prev: (pi_J)_{J\in I}
#sigma: correlation matrix
#nu: non-centrality parameter (7-dim.)
#w: weights for inverse normal combination
#tau.option: "overall": classic treatment of information rates
#            "individual": population-wise information rates
#n: matri x of theform
#        P{1} P{2} P{1,2}
#stage 1 n{1}^(1) n{2}^(1) n{1,2}^(1)
#stage 2 n{1}^(2) n{2}^(2) n{1,2}^(2)
#aspend: type of spending function: "Hwang" or "Kim"
#spendingfct: should the spendingfunction approach be used (logical)?
corrmat <- function(prev, w){
  rho12 <- prev[3]/sqrt(piv[1]*piv[2])
  rho1_1 <- sqrt(prev[1]/piv[1])
  rho2_2 <- sqrt(prev[2]/piv[2])
  rho1_12 <- sqrt(prev[3]/piv[1])
  rho2_12 <- sqrt(prev[3]/piv[2])
  corr12 <- rho12*(w[1]*w[2]+sqrt((1-w[1]^2)*(1-w[2]^2)))
  corr1_1 <- rho1_1*(w[1]*w[3]+sqrt((1-w[1]^2)*(1-w[3]^2)))
  corr2_2 <- rho2_2*(w[2]*w[4]+sqrt((1-w[2]^2)*(1-w[4]^2)))
  corr1_12 <- rho1_12*(w[1]*w[5]+sqrt((1-w[1]^2)*(1-w[5]^2)))
  corr2_12 <- rho2_12*(w[2]*w[5]+sqrt((1-w[2]^2)*(1-w[5]^2)))
  
  sigma <- matrix(c(1, rho12, w[1], w[2]*rho12, w[3]*rho1_1, 0, w[5]*rho1_12,
                    0, 1, w[1]*rho12, w[2], 0, w[4]*rho2_2, w[5]*rho2_12,
                    0, 0, 1, corr12,corr1_1, 0, corr1_12,
                    0, 0, 0, 1, 0, corr2_2, corr2_12,
                    0, 0, 0, 0, 1, 0, 0,
                    0, 0, 0, 0, 0, 1, 0, 
                    0, 0, 0, 0, 0, 0, 1), nr = 7, byrow = T)
  sigma <- sigma + t(sigma) - diag(7)
  return(sigma)
}

pwer <- function(crit, prev, nu = rep(0,7), sigma){
  names(prev) <- c("P1", "P2", "P12")
  #PWER:
  P <- numeric(3)
  P[1] <- 1-pmvnorm(upper = crit[c(1,3,5)], mean = nu[c(1,3,5)], 
                    sigma = sigma[c(1,3,5),c(1,3,5)], algorithm = Miwa())[1]
  P[2] <- 1-pmvnorm(upper = crit[c(2,4,6)], mean = nu[c(2,4,6)], 
                    sigma = sigma[c(2,4,6),c(2,4,6)], algorithm = Miwa())[1]
  P[3] <- 1-pmvnorm(upper = crit[c(1:4,7)], mean = nu[c(1:4,7)], 
                    sigma = sigma[c(1:4,7), c(1:4,7)], algorithm = Miwa())[1]
  PWER <- sum(prev*P)
  #PWER^(k):
  PWER1 <- prev[1]*(1-pnorm(crit[1], mean = nu[1],  sd = 1))+
    prev[2]*(1-pnorm(crit[2], mean = nu[2],  sd = 1))+
    prev[3]*(1-pmvnorm(upper = crit[1:2], mean = nu[1:2],
                       sigma = sigma[1:2,1:2], algorithm = Miwa())[1])
  
  PWER2 <- prev[1]*(pmvnorm(lower=c(-Inf,crit[3]),upper=c(crit[1],Inf),mean=nu[c(1,3)], 
                           sigma=sigma[c(1,3),c(1,3)], algorithm = Miwa())[1]+
                    pmvnorm(lower=c(-Inf,-Inf,crit[5]),upper=c(crit[c(1,3)],Inf),
                            mean=nu[c(1,3,5)],sigma=sigma[c(1,3,5),c(1,3,5)]))+
    
           prev[2]*(pmvnorm(lower=c(-Inf,crit[4]),upper=c(crit[2],Inf),mean=nu[c(2,4)], 
                           sigma=sigma[c(2,4),c(2,4)], algorithm = Miwa())[1]+ 
                    pmvnorm(lower=c(-Inf,-Inf,crit[6]),upper=c(crit[c(2,4)],Inf),
                            mean=nu[c(2,4,6)],sigma=sigma[c(2,4,6),c(2,4,6)]))+
    
           prev[3]*(pmvnorm(lower=c(-Inf,-Inf,crit[3]),upper=c(crit[1:2],Inf),
                            mean=nu[1:3], sigma=sigma[1:3,1:3], algorithm = Miwa())[1]+
                    pmvnorm(lower=c(rep(-Inf,3),crit[4]),upper=c(crit[1:3],Inf),mean=nu[1:4],
                            sigma=sigma[1:4,1:4])+
                    pmvnorm(lower=c(rep(-Inf,4),crit[7]),upper=c(crit[1:4],Inf),mean=nu[c(1:4,7)],
                            sigma=sigma[c(1:4,7),c(1:4,7)]))
  #results list:
  res <- list(PWER = PWER, P1 = P[1], P2 = P[2], P12 = P[3], PWER_stage1 = PWER1, PWER_stage2 = PWER2)
  return(res)
}
#WT critical values
critWT <- function(p, prev, sigma, alpha = 0.025, n, tau.option){
  if(tau.option == "overall"){
    tau <- sum(n[1,])/sum(n)
    tau <- tau*rep(1,5)
  }
  else if(tau.option == "population-wise"){
    tau1 <- (n[1,1]+n[1,3])/sum(n[,c(1,3)])
    tau2 <- (n[1,2]+n[1,3])/sum(n[,c(2,3)])
    tau_indpop <- n[1,]/c(sum(n[,1]),sum(n[,2]),sum(n[,3]))
    tau <- c(tau1, tau2, tau_indpop)
  }
  f <- function(cr){
    crit <- cr*(1/c(1,1,tau1,tau2,tau_indpop))^(p-0.5)
    pwer(crit=crit,prev=prev, sigma=sigma)$PWER-alpha
  }
  cr_const <- uniroot(f, interval = c(qnorm(1-alpha),5))$root
  res <- cr_const*(1/c(1,1,tau))^(p-0.5)
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
    crit = c(cr,cr,rep(0,5))
    pwer(crit =  crit, prev = prev, sigma = sigma)$PWER_stage1 - a(l=tau[1])
  }
  cr1 <- uniroot(f1, interval = c(qnorm(1-alpha), 5))$root
  #second stage critical value
  f2 <- function(cr){
    pwer(crit =  c(cr,cr,rep(cr1,5)), prev = prev, sigma = sigma)$PWER_stage2 - 
      (alpha-a(l=tau[1]))
  }
  cr2 <- uniroot(f2, interval = c(qnorm(1-alpha), 5))$root
  return(c(rep(cr1,2),rep(cr2,5)))
}
#Powers
#power-functions:
pwp <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025){
  if(nrow(n) == 2 & ncol(n) == 3){
    piv <- prev[1:2]+prev[3]
    nu <- delta*sqrt(n[1,1:2]+n[1,3], n[1,1:2]+n[2,1:2]+sum(n[,3]),
                     n[1,]+n[2,])
    w <- c((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])), #w1,w2
           n[1,]/(n[1,]+n[2,]))
    sigma <- corrmat(prev=prev, w=w)
    #critical values:
    if(is.null(spendingfct)){
      crit <- critWT(p=p, prev=prev, sigma=sigma, n=n,
                     alpha=alpha, tau.option = tau.option)
    }
    else{
      tau <- sum(n[1,])/sum(n)
      crit <- critES(p=p, prev=prev, sigma=sigma, alpha=alpha, 
                     tau = tau, aspend = spendingfct)
    }
    ##power value
    #denominator
    den <- prev[1]*ifelse(delta[1]>0 | delta[3]>0)+
      prev[2]*ifelse(delta[2]>0 | delta[4]>0)+
      prev[3]*ifelse(delta[1]>0 | delta[2]>0 | delta[5]>0)
    #numerator
    upper1 <- upper2 <- upper12 <- rep(Inf,7)
    upper1[c(1,3,5)] <- c(crit[c(1,3)]*ifelse(delta[1]>0,1,Inf), crit[5]*ifelse(delta[3]>0,1,Inf))
    upper2[c(2,4,6)] <- c(crit[c(2,4)]*ifelse(delta[2]>0,1,Inf), crit[6]*ifelse(delta[4]>0,1,Inf))
    upper12[c(1:4,7)] <- c(crit[c(1,3)]*ifelse(delta[1]>0,1,Inf), crit[c(2,4)]*ifelse(delta[2]>0,1,Inf),
                           crit[7]*ifelse(delta[1]>0|delta[2]>0|delta[5]>0))[c(1,3,2,4,5)]
    num <- prev[1]*(1-pmvnorm(upper=upper1, mean = nu, sigma=sigma))+
      prev[2]*(1-pmvnorm(upper=upper2, mean = nu, sigma=sigma))+
      prev[3]*(1-pmvnorm(upper=upper12, mean = nu, sigma=sigma))
    
    return(ifelse(den>0,num/den,0))
  }
}
pow1 <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025){
  piv <- prev[1:2]+prev[3]
  nu <- delta*sqrt(n[1,1:2]+n[1,3], n[1,1:2]+n[2,1:2]+sum(n[,3]),
                   n[1,]+n[2,])
  w <- c((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])), #w1,w2
         n[1,]/(n[1,]+n[2,]))
  sigma <- corrmat(prev=prev, w=w)
  #critical values:
  if(is.null(spendingfct)){
    crit <- critWT(p=p, prev=prev, sigma=sigma, n=n,
                   alpha=alpha, tau.option = tau.option)
  }
  else{
    tau <- sum(n[1,])/sum(n)
    crit <- critES(p=p, prev=prev, sigma=sigma, alpha=alpha, 
                   tau = tau, aspend = spendingfct)
  }
  pos <- which(nu > 0)
  1-pmvnorm(upper=crit[pos], mean = nu[pos], sigma = sigma[pos,pos], algorithm = Miwa())[1]
}

power <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025, powertype){
  if(powertype == "pwp"){
    pwp(p=p, prev=prev, delta=delta, n=n, spendingfct=spendingfct, tau.option=tau.option, alpha=alpha)
  }
  else if(powertype == "pow1"){
    pow1(p=p, prev=prev, delta=delta, N=N, spendingfct=spendingfct, tau.option=tau.option, alpha=alpha)
  }
}

#N root
Nroot <- function(p, prev, delta, gamma=c(1,1,1), alpha = 0.025, 
                  beta = 0.2, spendingfct = NULL, tau, powertype,
                  search.interval = c(1,2000)){
  fzero <- function(N){
    n <- rbind(N*prev, gamma*N*prev)
    power(n=n, p=p, prev=prev, delta = delta, gamma=gamma, 
          spendingfct = spendingfct, tau=tau, alpha = alpha)-(1-beta)
  }
  uniroot(fzero, interval = search.interval)$root
}

ASNroot <-function(p, prev, delta, gamma=c(1,1,1), alpha = 0.025, 
                   beta = 0.2, spendingfct = NULL, tau.option, powertype,
                   search.interval = c(1,2000)){
  #find N such that power = 1-beta:
  N0 <- Nroot(p=p, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, 
              spendingfct=spendingfct, tau.option=tau.option, powertype=powertype,
              search.interval=search.interval)
  piv <- prev[1:2]+prev[3]
  n <- rbind(N0*prev, N0*prev*gamma)
  nu <- delta*sqrt(n[1,1:2]+n[1,3], n[1,1:2]+n[2,1:2]+sum(n[,3]),
                   n[1,]+n[2,])
  w <- c((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])), #w1,w2
         n[1,]/(n[1,]+n[2,]))
  sigma <- corrmat(prev=prev, w=w)
  #critical values
  if(is.null(spendingfct)){
    crit <- critWT(p=p, prev=prev, sigma=sigma, alpha=alpha, tau.option=tau.option, n=n)
  }
  else{
    crit <- critES(p=p, prev=prev, sigma=sigma, alpha=alpha, tau=tau, aspend = spendingfct)
  }
  #ASN
  ASN <- N0 +
         (n[2,1]+n[2,3])*pmvnorm(lower=c(-Inf,crit[2]), upper=c(crit[1],Inf), 
                                 mean=nu[1:2], sigma=sigma[1:2,1:2])[1]+
         (n[2,2]+n[2,3])*pmvnorm(lower=c(crit[1],-Inf), upper=c(Inf,crit[2]), 
                            mean=nu[1:2], sigma=sigma[1:2,1:2])[1]+
         sum(n[2,])*pmvnorm(upper=crit[1:2], mean=nu[1:2], sigma=sigma[1:2,1:2])[1]
  return(ASN)
}

opt_par <- function(prev, delta, tau.option, gamma = c(1,1,1), spendingfct = NULL, alpha = 0.025,
                    beta = 0.2, powertype, search.inteval.N, search.interval.p, opt.criterion){
  if(opt.criterion == "Nmin"){
    fzero <- function(x){
      Nroot(p=x, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, tau.option=tau.option,
            spendingfct=spendingfct, powertype=powertype, search.interval=search.inteval.N)
    }
  }
  if(opt.criterion == "ASNmin"){
    fzero <- function(x){
      ASNroot(p=x, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, tau.option=tau.option,
              spendingfct=spendingfct, powertype=powertype, search.interval=search.inteval.N)
    }
  }
  optimize(fzero, interval = search.interval.p)$min
}