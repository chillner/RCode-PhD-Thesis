#####################################
#Design I - H_1 = H_2 = {h_1, h_2}
#####################################
library(multcomp)
library(xtable)
corrmat <- function(prev, w){
  piv <- prev[1:2]+prev[3]
  rho <- prev[3]/sqrt(piv[1]*piv[2])
  sigma <- matrix(c(1, w[1], rho, w[2]*rho,
                    0, 1, w[1]*rho, rho,
                    0, 0, 1, w[2],
                    0, 0, 0, 1), nr = 4, byrow = T)
  sigma <- sigma + t(sigma) - diag(4)
  return(sigma)
}
#PWER-function:
#crit matrix:
#  P1  P2
#1 c11 c12
#2 c21 c22
pwer <- function(crit, prev, nu = c(0,0,0,0), sigma){
  names(prev) <- c("P1", "P2", "P12")
  if(!is.matrix(crit)){
    crit <- matrix(c(crit,crit), byrow = FALSE, nr=2, nc=2)
  }
  #PWER:
  P <- numeric(3)
  P[1] <- 1-pmvnorm(upper= crit[,1], mean = nu[1:2], sigma = sigma[1:2,1:2], algorithm = Miwa())[1]
  P[2] <- 1-pmvnorm(upper= crit[,2], mean = nu[3:4], sigma = sigma[3:4,3:4], algorithm = Miwa())[1]
  P[3] <- 1-pmvnorm(upper= c(crit[,1], crit[,2]), mean = nu, sigma = sigma, algorithm = Miwa())[1]
  PWER <- sum(prev*P)
  #PWER^(k):
  PWER1 <- prev[1]*(1-pnorm(crit[1,1], mean = nu[1],  sd = 1))+
           prev[2]*(1-pnorm(crit[1,2], mean = nu[3],  sd = 1))+
           prev[3]*(1-pmvnorm(upper = crit[1,], mean = nu[c(1,3)],
                              sigma = sigma[c(1,3),c(1,3)], algorithm = Miwa())[1])
  
  PWER2 <- prev[1]*pmvnorm(lower=c(-Inf, crit[2,1]), upper=c(crit[1,1],Inf), mean = nu[1:2], 
                           sigma=sigma[1:2,1:2], algorithm = Miwa())[1]+
           prev[2]*pmvnorm(lower=c(-Inf, crit[2,2]), upper=c(crit[1,2],Inf), mean = nu[3:4],
                           sigma=sigma[3:4,3:4], algorithm = Miwa())[1]+
           prev[3]*(pmvnorm(lower=c(-Inf,crit[2,1],-Inf),upper=c(crit[1,1],Inf,crit[1,2]),
                            mean=nu[1:3], sigma=sigma[1:3,1:3], algorithm = Miwa())[1]+
                    pmvnorm(lower=c(-Inf,-Inf,-Inf,crit[2,2]), upper=c(crit[1,1],crit[2,1],crit[1,2],Inf),
                            mean=nu, sigma=sigma, algorithm = Miwa())[1])
  #results list:
  res <- list(PWER = PWER, P1 = P[1], P2 = P[2], P12 = P[3], PWER_stage1 = PWER1, PWER_stage2 = PWER2)
  return(res)
}

#crit Wang Tsiatis:
#tau matrix:
#  P1  P2
#1 n11 n12
#2 n21 n22
critWT <- function(p, prev, sigma, alpha = 0.025, n, tau.option){
  if(!is.matrix(n)){
    n <- rbind(n,n)
  }
  if(tau.option == "overall"){
    tau <- sum(n[1,])/sum(n)
    tau <- matrix(c(tau,tau,1,1), nr = 2, byrow = T)
  }
  else if(tau.option == "population-wise"){
    tau <- (n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)]))
    tau <- rbind(tau,c(1,1))
  }
  f <- function(cr){
    crit <- matrix(0, nr=2, nc=2)
    crit <- cr*sweep(tau,2,tau[1,],'/')^(p-0.5)
    pwer(crit=crit,prev=prev, nu = rep(0,4), sigma=sigma)$PWER-alpha
  }
  cr_const <- uniroot(f, interval = c(qnorm(1-alpha),5))$root
  res <- cr_const*sweep(tau,2,tau[1,],'/')^(p-0.5)
  colnames(res) <- c("c_P1", "c_P2")
  rownames(res) <- c("stage1", "stage2")
  return(res)
}

#crit error spending:
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
    pwer(crit =  cr, prev = prev, sigma = sigma, nu = rep(0,4))$PWER_stage1 - a(l=tau[1])
  }
  cr1 <- uniroot(f1, interval = c(qnorm(1-alpha), 5))$root
  #second stage critical value
  f2 <- function(cr){
    pwer(crit =  matrix(c(cr1,cr1,cr,cr), byrow=T, nr=2), 
         prev = prev, sigma = sigma, nu = rep(0,4))$PWER_stage2 - (alpha-a(l=tau[1]))
  }
  cr2 <- uniroot(f2, interval = c(qnorm(1-alpha), 5))$root
  res <- rbind(c(cr1, cr1),c(cr2, cr2))
  colnames(res) <- c("c_P1", "c_P2")
  rownames(res) <- c("stage1", "stage2")
  return(res)
}

#power-functions:
pwp <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025){
  if(nrow(n) == 2 & ncol(n) == 3){
    nu <- c(delta*sqrt(n[1,1:2]+n[1,3]), delta*sqrt(n[2,1:2]+n[2,3]))[c(1,3,2,4)]
    w <- sqrt((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])))
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
    #power value
    if(nu[1]<= 0 & nu[3] > 0){
      pwer(crit = crit, prev = prev, nu=nu, sigma=sigma)$P2
    }
    else if(nu[3]<= 0 & nu[1] > 0){
      pwer(crit = crit, prev = prev, nu=nu, sigma=sigma)$P1
    }
    else if(nu[1]>0 & nu[3]>0){
      pwer(crit = crit, prev = prev, nu=nu, sigma=sigma)$PWER
    }
    else{
      0
    }
  }
}

pow1 <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025){
  nu <- c(delta*sqrt(n[1,1:2]+n[1,3]), delta*sqrt(n[2,1:2]+n[2,3]))[c(1,3,2,4)]
  w <- sqrt((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])))
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
  crit.vec <- as.vector(crit)
  1-pmvnorm(upper=crit.vec[pos], mean = nu[pos], sigma = sigma[pos,pos], algorithm = Miwa())[1]
}

power <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025, powertype){
  if(powertype == "pwp"){
    pwp(p=p, prev=prev, delta=delta, n=n, spendingfct=spendingfct, tau.option=tau.option, alpha=alpha)
  }
  else if(powertype == "pow1"){
    pow1(p=p, prev=prev, delta=delta, n=n, spendingfct=spendingfct, tau.option=tau.option, alpha=alpha)
  }
}

#function for finding N:
Nroot <- function(p, prev, delta, gamma=c(1,1,1), alpha = 0.025, 
              beta = 0.2, spendingfct = NULL, tau.option, powertype,
              search.interval = c(1,2000)){
  fzero <- function(N){
    n <- rbind(N*prev, gamma*N*prev)
    power(p=p, prev=prev, delta = delta, n=n, powertype = powertype,
          spendingfct=spendingfct, tau.option=tau.option, alpha=alpha)-(1-beta)
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
  #critical values:
  n <- rbind(N0*prev, N0*prev*gamma)
  nu <- c(delta*sqrt(n[1,1:2]+n[1,3]), delta*sqrt(n[1,1:2]+n[2,1:2]+sum(n[,3])))[c(1,3,2,4)]
  w <- sqrt((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])))
  sigma <- corrmat(prev=prev, w=w)
  if(is.null(spendingfct)){
    crit <- critWT(p=p, prev=prev, sigma=sigma, alpha=alpha, tau.option=tau.option, n=n)
  }
  else{
    tau <- sum(n[1,])/sum(n)
    crit <- critES(p=p, prev=prev, sigma=sigma, alpha=alpha, tau=tau, aspend = spendingfct)
  }
  #ASN
  N0 + (n[2,1]+n[2,3])*pmvnorm(lower=c(-Inf,crit[1,2]), upper=c(crit[1,1],Inf), 
                           mean = nu[c(1,3)], sigma = sigma[c(1,3),c(1,3)])[1] +
  (n[2,2]+n[2,3])*pmvnorm(lower=c(crit[1,1], -Inf), upper=c(Inf, crit[1,2]), 
                                  mean = nu[c(1,3)], sigma = sigma[c(1,3),c(1,3)])[1] +
  sum(n[2,])*pmvnorm(lower=c(-Inf,-Inf), upper=c(crit[1,1],crit[1,2]), 
                       mean = nu[c(1,3)], sigma = sigma[c(1,3),c(1,3)])[1]
}

#find optimal critical values (W&T or errspend parameter):
opt_par <- function(prev, delta, tau.option, gamma = c(1,1,1), spendingfct = NULL, alpha = 0.025,
                    beta = 0.2, powertype, search.inteval.N = c(1,2000), search.interval.p, opt.criterion,
                    grid.method = FALSE, grid.steps = 0.001){
  if(!grid.method){
    if(opt.criterion == "Nmin"){
      fzero <- function(x){
        Nroot(p=x, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, tau.option=tau.option,
              spendingfct = spendingfct, powertype = powertype, search.interval = search.inteval.N)
      }
    }
    if(opt.criterion == "ASNmin"){
      fzero <- function(x){
        ASNroot(p=x, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, tau.option=tau.option,
              spendingfct = spendingfct, powertype = powertype, search.interval = search.inteval.N)
      }
    }
    p_opt <- optimize(fzero, interval = search.interval.p)$min
  }
  else{
    s <- search.interval.p
    lspace <- seq(s[1],s[2], grid.steps)
    Nstored <- numeric(length(lspace))
    if(opt.criterion == "Nmin"){
      for(i in 1:length(lspace)){
        #skip_to_next <- FALSE
        Nstored[i] <- tryCatch(Nroot(p=lspace[i], prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, 
                            tau.option=tau.option, spendingfct = spendingfct, powertype = powertype,
                            search.interval = search.interval.N), error = function(e){Inf})
      }
    }
    else if(opt.criterion == "ASNmin"){
      for(i in 1:length(lspace)){
        Nstored[i] <- tryCatch(ASNroot(p=lspace[i], prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, 
                              tau.option=tau.option, spendingfct = spendingfct, powertype = powertype,
                              search.interval = search.interval.N), error = function(e){Inf})
      }
    }
    p_opt <- lspace[which.min(Nstored)]
  }
  return(p_opt)
}



################
#tests:
##################
n <- rbind(c(40,40,20),c(50,50,25))
w <- sqrt((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])))
prev <- c(.4,.4,.2)
sigma <- corrmat(prev = prev, w=w)

#pwer test:
pwer(crit = c(2,2.1), prev=prev, sigma=sigma)

#critWT test:
p = 0
critWT(p=p, prev=prev, sigma=sigma, n=n, tau.option="population-wise")

#critES test:
critES(p=0.5, prev=prev, sigma=sigma, tau = 0.5, aspend = "Kim")

#power test:
pwp(p = 0.5, prev = prev, delta = c(0.3,0.3), n=n, spendingfct = NULL, tau.option="overall", alpha=0.025)
pow1(p = 0.5, prev = prev, delta = c(0.3,0.3), n=n, spendingfct = "Kim", tau.option="overall", alpha=0.025)

#Nroot test:
Nroot(p = 0.5, prev = prev, gamma = c(1,1,1), delta = c(0.3,0.3),
      spendingfct = "Kim", tau.option="overall", powertype = "pwp")
Nroot(p = 0.5, prev = prev, gamma = c(1,1,1), delta = c(0.2,0.2),
      spendingfct = NULL, tau.option="overall", powertype = "pow1")

#ASNroot test:
ASNroot(p = 0.5, prev = prev, gamma = c(1,1,1), delta = c(0.3,0.3),
        spendingfct = NULL, tau.option="overall", powertype = "pwp")

#opt_par test:
opt_par(prev=prev, delta=c(0.2,0.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
        search.interval.p = c(-.8,0.1), opt.criterion = "Nmin")



######Example values:
prev1 = c(.4,.4,.2); prev2 = c(.6,.3,.1)
block.N.ASN <- data.frame(N = numeric(3), ASN = numeric(3))
optWT <- list(PWP.prev1 = block.N.ASN, PWP.prev2 = block.N.ASN, 
              Pow1.prev1 = block.N.ASN, Pow1.prev2 = block.N.ASN)

####################
#####critWT
####################
####powertype PWP:
###prev = c(.4,.4,.2)
#####################
##gamma=c(1,1,1)
#####################
#N
optWT$PWP.prev1[1,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pwp", 
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT$PWP.prev1[2,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pwp", 
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT$PWP.prev1[3,1] <- optWT$PWP.prev1[2,1]
#ASN
optWT$PWP.prev1[1,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pwp",
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT$PWP.prev1[2,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pwp",
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT$PWP.prev1[3,2] <- optWT$PWP.prev1[2,2] 
###prev = c(.6,.3,.1)
#N
optWT$PWP.prev2[1,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pwp", 
        search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT$PWP.prev2[2,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pwp",
        search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT$PWP.prev2[3,1] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pwp",
        search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
#ASN
optWT$PWP.prev2[1,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pwp",
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT$PWP.prev2[2,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pwp",
                                search.inteval.N = c(1,2000), search.interval.p =c(-1,1), opt.criterion = "ASNmin")
optWT$PWP.prev2[3,2] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pwp",
                                search.inteval.N = c(1,2000), search.interval.p =c(-1,1), opt.criterion = "ASNmin")

####powertype Pow1:
###prev = c(.4,.4,.2)
#N
optWT$Pow1.prev1[1,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pow1", 
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT$Pow1.prev1[2,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pow1",
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT$Pow1.prev1[3,1] <- optWT$Pow1.prev1[2,1]
#ASN
optWT$Pow1.prev1[1,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pow1",
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT$Pow1.prev1[2,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pow1",
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT$Pow1.prev1[3,2] <- optWT$Pow1.prev1[2,2]
###prev = c(.6,.3,.1)
#N
optWT$Pow1.prev2[1,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pow1",
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT$Pow1.prev2[2,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pow1", 
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT$Pow1.prev2[3,1] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pow1",
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
#ASN
optWT$Pow1.prev2[1,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pow1",
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT$Pow1.prev2[2,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pow1",
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT$Pow1.prev2[3,2] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pow1",
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")


#########################
#gamma=c(1.25,1.25,1.25)
#########################
optWT125 <- list(PWP.prev1 = block.N.ASN, PWP.prev2 = block.N.ASN, 
                 Pow1.prev1 = block.N.ASN, Pow1.prev2 = block.N.ASN)
gamma125 <- c(1.25,1.25,1.25)
#N
optWT125$PWP.prev1[1,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pwp", gamma=gamma125,
                                   search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT125$PWP.prev1[2,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pwp", gamma=gamma125,
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT125$PWP.prev1[3,1] <- optWT125$PWP.prev1[2,1]
#ASN
optWT125$PWP.prev1[1,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pwp",gamma=gamma125,
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT125$PWP.prev1[2,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pwp",gamma=gamma125,
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT125$PWP.prev1[3,2] <- optWT125$PWP.prev1[2,2] 
###prev = c(.6,.3,.1)
#N
optWT125$PWP.prev2[1,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pwp", gamma=gamma125,
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT125$PWP.prev2[2,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pwp",gamma=gamma125,
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT125$PWP.prev2[3,1] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pwp",gamma=gamma125,
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
#ASN
optWT125$PWP.prev2[1,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pwp",gamma=gamma125,
                                search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT125$PWP.prev2[2,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pwp",gamma=gamma125,
                                search.inteval.N = c(1,2000), search.interval.p =c(-1,1), opt.criterion = "ASNmin")
optWT125$PWP.prev2[3,2] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pwp",gamma=gamma125,
                                search.inteval.N = c(1,2000), search.interval.p =c(-1,1), opt.criterion = "ASNmin")

####powertype Pow1:
###prev = c(.4,.4,.2)
#N
optWT125$Pow1.prev1[1,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pow1", gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT125$Pow1.prev1[2,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pow1",gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT125$Pow1.prev1[3,1] <- optWT125$Pow1.prev1[2,1]
#ASN
optWT125$Pow1.prev1[1,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pow1",gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT125$Pow1.prev1[2,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pow1",gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT125$Pow1.prev1[3,2] <- optWT125$Pow1.prev1[2,2]
###prev = c(.6,.3,.1)
#N
optWT125$Pow1.prev2[1,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pow1",gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT125$Pow1.prev2[2,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pow1", gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
optWT125$Pow1.prev2[3,1] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pow1",gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "Nmin")
#ASN
optWT125$Pow1.prev2[1,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pow1",gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT125$Pow1.prev2[2,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pow1",gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")
optWT125$Pow1.prev2[3,2] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pow1",gamma=gamma125,
                                 search.inteval.N = c(1,2000), search.interval.p = c(-1,1), opt.criterion = "ASNmin")



######################
#####critES
######################
######################
#gamma = c(1,1,1)
######################
#####critES
####powertype PWP:
###prev = c(.4,.4,.2)
optES <- list(PWP.prev1 = block.N.ASN, PWP.prev2 = block.N.ASN, 
              Pow1.prev1 = block.N.ASN, Pow1.prev2 = block.N.ASN)
#N
optES$PWP.prev1[1,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                             search.interval.p = c(-1,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES$PWP.prev1[2,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                search.interval.p = c(0,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES$PWP.prev1[3,1] <- optES$PWP.prev1[2,1]  
#ASN
optES$PWP.prev1[1,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES$PWP.prev1[2,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES$PWP.prev1[3,2] <- optES$PWP.prev1[2,2] 

###prev = c(.6,.3,.1)
#N
optES$PWP.prev2[1,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES$PWP.prev2[2,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES$PWP.prev2[3,1] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "Nmin", spendingfct = "Hwang")
#ASN
optES$PWP.prev2[1,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES$PWP.prev2[2,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES$PWP.prev2[3,2] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")

####powertype Pow1:
###prev = c(.4,.4,.2)
#N
optES$Pow1.prev1[1,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 search.interval.p = c(-1,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES$Pow1.prev1[2,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 search.interval.p = c(-1,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES$Pow1.prev1[3,1] <- optES$Pow1.prev1[2,1]  
#ASN
optES$Pow1.prev1[1,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES$Pow1.prev1[2,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES$Pow1.prev1[3,2] <- optES$Pow1.prev1[2,2]

###prev = c(.6,.3,.1)
#N
optES$Pow1.prev2[1,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES$Pow1.prev2[2,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES$Pow1.prev2[3,1] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                search.interval.p = c(-1,1), opt.criterion = "Nmin", spendingfct = "Hwang")
#ASN
optES$Pow1.prev2[1,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES$Pow1.prev2[2,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES$Pow1.prev2[3,2] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 search.interval.p = c(-1,1), opt.criterion = "ASNmin", spendingfct = "Hwang")

#########################
#gamma=c(1.25,1.25,1.25)
#########################
optES125 <- list(PWP.prev1 = block.N.ASN, PWP.prev2 = block.N.ASN, 
              Pow1.prev1 = block.N.ASN, Pow1.prev2 = block.N.ASN)
#N
optES125$PWP.prev1[1,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES125$PWP.prev1[2,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES125$PWP.prev1[3,1] <- optES125$PWP.prev1[2,1]  
#ASN
optES125$PWP.prev1[1,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES125$PWP.prev1[2,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES125$PWP.prev1[3,2] <- optES125$PWP.prev1[2,2] 

###prev = c(.6,.3,.1)
#N
optES125$PWP.prev2[1,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES125$PWP.prev2[2,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES125$PWP.prev2[3,1] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
#ASN
optES125$PWP.prev2[1,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES125$PWP.prev2[2,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES125$PWP.prev2[3,2] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
                                gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")

####powertype Pow1:
###prev = c(.4,.4,.2)
#N
optES125$Pow1.prev1[1,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES125$Pow1.prev1[2,1] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES125$Pow1.prev1[3,1] <- optES125$Pow1.prev1[2,1]  
#ASN
optES125$Pow1.prev1[1,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES125$Pow1.prev1[2,2] <- opt_par(prev=c(.4,.4,.2), delta=c(.2,0), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES125$Pow1.prev1[3,2] <- optES125$Pow1.prev1[2,2]

###prev = c(.6,.3,.1)
#N
optES125$Pow1.prev2[1,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES125$Pow1.prev2[2,1] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
optES125$Pow1.prev2[3,1] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "Nmin", spendingfct = "Hwang")
#ASN
optES125$Pow1.prev2[1,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES125$Pow1.prev2[2,2] <- opt_par(prev=c(.6,.3,.1), delta=c(.2,0), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")
optES125$Pow1.prev2[3,2] <- opt_par(prev=c(.6,.3,.1), delta=c(0,.2), tau.option="overall", powertype = "pow1", search.inteval.N = c(1,2000),
                                 gamma=gamma125, search.interval.p = c(-2,1), opt.criterion = "ASNmin", spendingfct = "Hwang")

###################
#Latex-Table
###################
optWTtex<-round(rbind(cbind(optWT$PWP.prev1, optWT$PWP.prev2), cbind(optWT$Pow1.prev1, optWT$Pow1.prev2)),2)
optWT125tex<-round(rbind(cbind(optWT125$PWP.prev1, optWT125$PWP.prev2), cbind(optWT125$Pow1.prev1, optWT125$Pow1.prev2)),2)
optEStex<-round(rbind(cbind(optES$PWP.prev1, optES$PWP.prev2), cbind(optES$Pow1.prev1, optES$Pow1.prev2)),2)
optES125tex<-round(rbind(cbind(optES125$PWP.prev1, optES125$PWP.prev2), cbind(optES125$Pow1.prev1, optES125$Pow1.prev2)),2)
optcrittable <- rbind(cbind(optWTtex, optWT125tex), cbind(optEStex, optEStex))
xtable(optcrittable)
