#####################################
#Design I - H_1 = H_2 = {h_1, h_2}
#####################################
library(multcomp)
library(xtable)
corrmat <- function(prev, w){
  piv <- prev[1:2]+prev[3]
  rho <- prev[3]/sqrt(piv[1]*piv[2])
  sigma <- matrix(c(1, w[1], rho, w[2]*rho,
                    0, 1, w[1]*rho, rho*(w[1]*w[2]+sqrt((1-w[1]^2)*(1-w[2]^2))),
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
pwer <- function(crit, prev, nu = c(0,0,0,0), sigma, uF = -Inf){
  names(prev) <- c("P1", "P2", "P12")
  if(!is.matrix(crit)){
    crit <- matrix(c(crit,crit), byrow = FALSE, nr=2, nc=2)
  }
  #PWER:
  P <- numeric(3)
  P[1] <- 1-pmvnorm(upper= crit[,1], mean = nu[1:2], sigma = sigma[1:2,1:2], algorithm = Miwa())[1]
  P[2] <- 1-pmvnorm(upper= crit[,2], mean = nu[3:4], sigma = sigma[3:4,3:4], algorithm = Miwa())[1]
  P[3] <- 1-pmvnorm(upper= c(crit[,1], crit[,2]), mean = nu, sigma = sigma, algorithm = Miwa())[1]
  #PWER^(k):
  P11 <- 1-pnorm(crit[1,1], mean = nu[1],  sd = 1)
  P12 <- 1-pnorm(crit[1,2], mean = nu[3],  sd = 1)
  P13 <- 1-pmvnorm(upper = crit[1,], mean = nu[c(1,3)],
                   sigma = sigma[c(1,3),c(1,3)], algorithm = Miwa())[1]
  P21 <-pmvnorm(lower=c(uF, crit[2,1]), upper=c(crit[1,1],Inf), mean = nu[1:2], 
                sigma=sigma[1:2,1:2], algorithm = Miwa())[1]
  P22 <- pmvnorm(lower=c(uF, crit[2,2]), upper=c(crit[1,2],Inf), mean = nu[3:4],
                 sigma=sigma[3:4,3:4], algorithm = Miwa())[1]
  P23 <- pmvnorm(lower=c(uF,crit[2,1],uF),upper=c(crit[1,1],Inf,crit[1,2]),
                 mean=nu[1:3], sigma=sigma[1:3,1:3], algorithm = Miwa())[1]+
         pmvnorm(lower=c(uF,uF,uF,crit[2,2]), upper=c(crit[1,1],crit[2,1],crit[1,2],Inf),
                 mean=nu, sigma=sigma, algorithm = Miwa())[1]
  PWER1 <- prev[1]*P11+prev[2]*P12+prev[3]*P13
  PWER2 <- prev[1]*P21+prev[2]*P22+prev[3]*P23
  
  if(is.null(uF)){
    PWER <- sum(prev*P)
  }
  else{
    PWER <- PWER1 + PWER2
  }
  #results list:
  res <- list(PWER = PWER, P1 = P11+P21, P2 = P12+P22, P12 = P13+P23, PWER_stage1 = PWER1, PWER_stage2 = PWER2)
  return(res)
}

#crit Wang Tsiatis:
#tau matrix:
#  P{1}  P{2} P{1,2}
#1 n{1}^(1)  n{2}^(1) n{1,2}^(1)
#2 n{1}^(2)  n{2}^(2) n{1,2}^(2)
critWT <- function(p, prev, sigma, alpha = 0.025, n, tau.option, uF = -Inf){
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
    pwer(crit=crit,prev=prev, nu = rep(0,4), sigma=sigma, uF = uF)$PWER-alpha
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
pwp <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025, uF = -Inf){
  if(nrow(n) == 2 & ncol(n) == 3){
    nu <- c(delta*sqrt(n[1,1:2]+n[1,3]), delta*sqrt(c(sum(n[,c(1,3)]),sum(n[,c(2,3)]))))[c(1,3,2,4)]
    w <- sqrt((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])))
    sigma <- corrmat(prev=prev, w=w)
    #critical values:
    if(is.null(spendingfct)){
      crit <- critWT(p=p, prev=prev, sigma=sigma, n=n, uF = uF,
                     alpha=alpha, tau.option = tau.option)
    }
    else{
      tau <- sum(n[1,])/sum(n)
      crit <- critES(p=p, prev=prev, sigma=sigma, alpha=alpha, 
                     tau = tau, aspend = spendingfct)
    }
    #power value
    if(nu[1]<= 0 & nu[3] > 0){
      pwer(crit = crit, prev = prev, nu=nu, sigma=sigma, uF=uF)$P2
    }
    else if(nu[3]<= 0 & nu[1] > 0){
      pwer(crit = crit, prev = prev, nu=nu, sigma=sigma, uF=uF)$P1
    }
    else if(nu[1]>0 & nu[3]>0){
      pwer(crit = crit, prev = prev, nu=nu, sigma=sigma, uF=uF)$PWER
    }
    else{
      0
    }
  }
}

pow1 <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025){
  nu <- c(delta*sqrt(n[1,1:2]+n[1,3]), delta*sqrt(c(sum(n[,c(1,3)]),sum(n[,c(2,3)]))))[c(1,3,2,4)]
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

power <- function(p, prev, delta, n, spendingfct = NULL, tau.option, alpha = 0.025, powertype, uF = -Inf){
  if(powertype == "pwp"){
    pwp(p=p, prev=prev, delta=delta, n=n, spendingfct=spendingfct, tau.option=tau.option, alpha=alpha, uF = uF)
  }
  else if(powertype == "pow1"){
    pow1(p=p, prev=prev, delta=delta, n=n, spendingfct=spendingfct, tau.option=tau.option, alpha=alpha)
  }
}

#function for finding N:
Nroot <- function(p, prev, delta, gamma=c(1,1,1), alpha = 0.025, uF = -Inf,
              beta = 0.2, spendingfct = NULL, tau.option, powertype,
              search.interval = c(1,2000)){
  fzero <- function(N){
    n <- rbind(N*prev, gamma*N*prev)
    power(p=p, prev=prev, delta = delta, n=n, powertype = powertype, uF = uF,
          spendingfct=spendingfct, tau.option=tau.option, alpha=alpha)-(1-beta)
  }
  uniroot(fzero, interval = search.interval)$root
  
}

ASNroot <-function(p, prev, delta, gamma=c(1,1,1), alpha = 0.025, uF = -Inf,
                   beta = 0.2, spendingfct = NULL, tau.option, powertype,
                   search.interval = c(1,2000)){
  #find N such that power = 1-beta:
  N0 <- Nroot(p=p, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, uF = uF,
              spendingfct=spendingfct, tau.option=tau.option, powertype=powertype,
              search.interval=search.interval)
  #critical values:
  n <- rbind(N0*prev, N0*prev*gamma)
  nu <- c(delta*sqrt(n[1,1:2]+n[1,3]), delta*sqrt(n[1,1:2]+n[2,1:2]+sum(n[,3])))[c(1,3,2,4)]
  w <- sqrt((n[1,1:2]+n[1,3])/c(sum(n[,c(1,3)]),sum(n[,c(2,3)])))
  sigma <- corrmat(prev=prev, w=w)
  if(is.null(spendingfct)){
    crit <- critWT(p=p, prev=prev, sigma=sigma, alpha=alpha, tau.option=tau.option, n=n, uF = uF)
  }
  else{
    tau <- sum(n[1,])/sum(n)
    crit <- critES(p=p, prev=prev, sigma=sigma, alpha=alpha, tau=tau, aspend = spendingfct)
  }
  #ASN
  N0 + (n[2,1]+n[2,3])*pmvnorm(lower=c(uF,crit[1,2]), upper=c(crit[1,1],Inf), 
                           mean = nu[c(1,3)], sigma = sigma[c(1,3),c(1,3)])[1] + 
    (n[2,2]+n[2,3])*pmvnorm(lower=c(crit[1,1], uF), upper=c(Inf, crit[1,2]), 
                            mean = nu[c(1,3)], sigma = sigma[c(1,3),c(1,3)])[1] + 
    sum(n[2,])*pmvnorm(lower=c(uF,uF), upper=c(crit[1,1],crit[1,2]), mean = nu[c(1,3)], sigma = sigma[c(1,3),c(1,3)])[1]
}

#find optimal critical values (W&T or errspend parameter):
opt_par <- function(prev, delta, tau.option, gamma = c(1,1,1), spendingfct = NULL, alpha = 0.025, uF = -Inf, 
                    beta = 0.2, powertype, search.inteval.N = c(1,2000), search.interval.p, opt.criterion,
                    grid.method = FALSE, grid.steps = 0.001){
  if(!grid.method){
    if(opt.criterion == "Nmin"){
      fzero <- function(x){
        Nroot(p=x, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, tau.option=tau.option, uF = uF, 
              spendingfct = spendingfct, powertype = powertype, search.interval = search.inteval.N)
      }
    }
    if(opt.criterion == "ASNmin"){
      fzero <- function(x){
        ASNroot(p=x, prev=prev, delta=delta, gamma=gamma, alpha=alpha, beta=beta, tau.option=tau.option, uF = uF, 
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
      spendingfct = NULL, tau.option="overall", powertype = "pwp")
Nroot(p = 0.5, prev = prev, gamma = c(1,1,1), delta = c(0.2,0.2),
      spendingfct = NULL, tau.option="overall", powertype = "pow1")

#ASNroot test:
ASNroot(p = 0.5, prev = prev, gamma = c(1,1,1), delta = c(0.3,0.3),
        spendingfct = NULL, tau.option="overall", powertype = "pwp")

#opt_par test:
opt_par(prev=prev, delta=c(0.2,0.2), tau.option="overall", powertype = "pwp", search.inteval.N = c(1,2000),
        search.interval.p = c(-.8,0.1), opt.criterion = "Nmin")

###################
#Table Section 5
###################
tabD1 <- as.data.frame(matrix(0, nr = 6, nc = 9))
colnames(tabD1) <- c("pi{1}","pi{2}", "pi{1,2}", "c0", "c05", "N_PWP0", "N_PWP05", "N_Pow1_0", "N_Pow1_05")
pr <- rbind(c(0.3,0.3,0.4),
            c(0.35,0.35,0.3),
            c(0.4, 0.4, 0.2),
            c(0.4,0.2,0.4),
            c(0.4,0.3,0.3),
            c(0.6,0.2,0.2))
tabD1[,1:3] <- pr 
w <- 1/sqrt(c(2,2))
delta <- c(0.3,0.3)
for(i in 1:6){
  pri <- pr[i,]
  n <- rbind(pri,pri)
  corr <- corrmat(prev = pri, w = w)
  tabD1[i, 4] <- critWT(p=0, prev=pri, sigma = corr, n = n, tau.option = "overall")[1,1]
  tabD1[i, 5] <- critWT(p=.5, prev=pri, sigma = corr, n = n, tau.option = "overall")[1,1]
  tabD1[i, 6] <- Nroot(p=0, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau.option ="overall", powertype = "pwp")
  tabD1[i, 7] <- Nroot(p=.5, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau.option ="overall", powertype = "pwp")
  tabD1[i, 8] <- Nroot(p=0, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau.option ="overall", powertype = "pow1")
  tabD1[i, 9] <- Nroot(p=.5, prev=pri, delta=delta, beta = 0.1, spendingfct = NULL, tau.option ="overall", powertype = "pow1")
}
tabD1[,4:5] <- round(tabD1[,4:5], 3)
tabD1[,6:9] <- ceiling(tabD1[,6:9])
tabD1

xtable(x=tabD1, caption = "A caption", label = "tab: tabD1", digits = 3)






