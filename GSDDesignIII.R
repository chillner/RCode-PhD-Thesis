#######################################################################
#Design III - H_1 = {h_1, h_2}, H_2 = {H_1, H_2, H_{1}, H_{2}, H_{1,2}}
#######################################################################
library(multcomp)
#PWER-function:
#crit list
#  P1  P2
#(c_1^1, c_2^1, c_1^2, c_2^2, c_{1}^2, c_{2}^2, c_{1,2}^2)
corrmat <- function(prev, w){
  piv <- prev[1:2]+prev[3]
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
    tau <- c(1,1,tau*rep(1,5))
  }
  else if(tau.option == "population-wise"){
    tau1 <- (n[1,1]+n[1,3])/sum(n[,c(1,3)])
    tau2 <- (n[1,2]+n[1,3])/sum(n[,c(2,3)])
    tau_indpop <- n[1,]/c(sum(n[,1]),sum(n[,2]),sum(n[,3]))
    tau <- c(1,1, tau1, tau2, tau_indpop)
  }
  f <- function(cr){
    crit <- cr*(1/tau)^(p-0.5)
    pwer(crit=crit,prev=prev, sigma=sigma)$PWER-alpha
  }
  cr_const <- uniroot(f, interval = c(qnorm(1-alpha),5))$root
  res <- cr_const*(1/tau)^(p-0.5)
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
    nu <- c(delta[1],delta[2],delta[1],delta[2],delta[3:5])*sqrt(c(n[1,1:2]+n[1,3], n[1,1:2]+n[2,1:2]+sum(n[,3]), n[1,]+n[2,]))
    w <- sqrt(c((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])),
           n[1,]/(n[1,]+n[2,])))
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
    den <- prev[1]*(delta[1]>0 | delta[3]>0)+
      prev[2]*(delta[2]>0 | delta[4]>0)+
      prev[3]*(delta[1]>0 | delta[2]>0 | delta[5]>0)
    #numerator
    upper1 <- upper2 <- upper12 <- rep(Inf,7)
    upper1[c(1,3,5)] <- c(crit[c(1,3)]*ifelse(delta[1]>0,1,Inf), crit[5]*ifelse(delta[3]>0,1,Inf))
    upper2[c(2,4,6)] <- c(crit[c(2,4)]*ifelse(delta[2]>0,1,Inf), crit[6]*ifelse(delta[4]>0,1,Inf))
    upper12[c(1:4,7)] <- c(crit[c(1,3)]*ifelse(delta[1]>0,1,Inf), crit[c(2,4)]*ifelse(delta[2]>0,1,Inf),
                           crit[7]*ifelse(delta[1]>0|delta[2]>0|delta[5]>0,1,Inf))[c(1,3,2,4,5)]
    num <- prev[1]*(1-pmvnorm(upper=upper1, mean = nu, sigma=sigma))+
      prev[2]*(1-pmvnorm(upper=upper2, mean = nu, sigma=sigma))+
      prev[3]*(1-pmvnorm(upper=upper12, mean = nu, sigma=sigma))
    
    return(ifelse(den>0,num/den,0))
  }
}

pow1 <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025){
  piv <- prev[1:2]+prev[3]
  nu <- c(delta[1],delta[2],delta[1],delta[2],delta[3:5])*sqrt(c(n[1,1:2]+n[1,3], n[1,1:2]+n[2,1:2]+sum(n[,3]), n[1,]+n[2,]))
  w <- sqrt(c((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])), n[1,]/(n[1,]+n[2,])))
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
  1-pmvnorm(upper=crit[pos], mean = nu[pos], sigma = sigma[pos,pos])[1]
}

power <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025, powertype){
  if(powertype == "pwp"){
    pwp(p=p, prev=prev, delta=delta, n=n, spendingfct=spendingfct, tau.option=tau.option, alpha=alpha)
  }
  else if(powertype == "pow1"){
    pow1(p=p, prev=prev, delta=delta, n=n, spendingfct=spendingfct, tau.option=tau.option, alpha=alpha)
  }
}

#N root
Nroot <- function(p, prev, delta, gamma=c(1,1,1), alpha = 0.025, 
                  beta = 0.2, spendingfct = NULL, tau.option, powertype,
                  search.interval = c(1,2000)){
  fzero <- function(N){
    n <- rbind(N*prev, gamma*N*prev)
    power(n=n, p=p, prev=prev, delta = delta, powertype = powertype,
          spendingfct = spendingfct, tau.option=tau.option, alpha = alpha)-(1-beta)
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


######################
#Table Section 5
######################
tabD3 <- as.data.frame(matrix(0, nr = 6, nc = 9))
colnames(tabD3) <- c("pi{1}","pi{2}", "pi{1,2}", "c0", "c05", "N_PWP0", "N_PWP05", "N_Pow1_0", "N_Pow1_05")
pr <- rbind(c(0.3,0.3,0.4),
            c(0.35,0.35,0.3),
            c(0.4, 0.4, 0.2),
            c(0.4,0.2,0.4),
            c(0.4,0.3,0.3),
            c(0.6,0.2,0.2))
tabD3[,1:3] <- pr 
w <- 1/sqrt(c(2,2,2,2,2))
delta <- c(0.3,0.3,0.3,0.3,0.3)
for(i in 1:6){
  pri <- pr[i,]
  n <- rbind(pri,pri)
  corr <- corrmat(prev = pri, w = w)
  tabD3[i, 4] <- critWT(p=0, prev=pri, sigma = corr, n = n, tau.option = "overall")[1]
  tabD3[i, 5] <- critWT(p=.5, prev=pri, sigma = corr, n = n, tau.option = "overall")[1]
  tabD3[i, 6] <- Nroot(p=0, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau.option ="overall", powertype = "pwp")
  tabD3[i, 7] <- Nroot(p=.5, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau.option ="overall", powertype = "pwp")
  tabD3[i, 8] <- Nroot(p=0, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau.option ="overall", powertype = "pow1")
  tabD3[i, 9] <- Nroot(p=.5, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau.option ="overall", powertype = "pow1")
}
tabD3[,4:5] <- round(tabD3[,4:5], 3)
tabD3[,6:9] <- ceiling(tabD3[,6:9])
tabD3

xtable(x=tabD3, caption = "A caption", label = "tab: tabD3", digits = 3)
