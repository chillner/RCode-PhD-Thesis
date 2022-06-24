##############################
#PWER nested pop
##############################
corrmatnest <- function(prev){
  m <- length(prev)
  rho <- matrix(0, nr = m , nc = m)
  for(j in 1:m){
    for(i in 1:j){
      rho[i,j] <- sqrt(i*prev[j]/(j*prev[i]))
    }
  }
  rho <- t(rho) + rho - diag(m)
  rho
}
#prev: prevalences pi[i] of the disjoint populations P[i]
#crit: vector of critical values (same length as prev and nu)
#nu: non-centrality parameter/ expectation of the multi-variate normal Z = (Z_1,...,Z_m)
#sigma: m x m Correlation matrix of Z
pwernest <- function(crit, prev, nu, corrmatr){
  m <- length(prev)
  #PWER
  Pstage <- numeric(m)
  for(j in 1:m){
    lowerj <- c(rep(-Inf,j-1), crit[j]) 
    if(j == 1){
      upperj <- Inf
    } else { upperj <- c(crit[1:(j-1)], Inf)}
    Pstage[j] <- prev[j]*pmvnorm(lower = lowerj, upper = upperj, mean = nu[1:j], sigma = corrmatr[1:j,1:j])[1]
  }
  res <- list(PWER = sum(Pstage), PWERj = Pstage/prev, PWERstages = Pstage)
  return(res)
}
#critical values
critWTnest <- function(p, prev, corrmatr, alpha = 0.025){
  m <- length(prev)
  tau <- numeric(m)
  for(j in 1:m){
    tau[j] <- sum(prev[1:j])/sum(prev)
  }
  f <- function(cr){
    crit <- cr*(tau/tau[1])^(p-0.5)
    pwernest(crit=crit,prev=prev, nu = rep(0,m), corrmatr = corrmatr)$PWER-alpha
  }
  const <- uniroot(f, interval = c(0,5))$root
  return(const*(tau/tau[1])^(p-0.5))
}
critESnest <- function(p, prev, corrmatr, alpha=0.025, aspend){
  if(aspend == "Kim"){
    a <- function(l, par = p){alpha*l^par}
  }
  else if(aspend == "Hwang"){
    a <- function(l, par = p){
      alpha*ifelse(par != 0, (1-exp(-par*l))/(1-exp(-par)), l)
    }
  }
  m <- length(prev)
  tau <- numeric(m)
  for(j in 1:m){
    tau[j] <- sum(prev[1:j])/sum(prev)
  }
  crits <- numeric(m)
  for(j in 1:m){
    f <- function(cr){
      if(j == 1){
        pwernest(crit = c(cr,rep(0, m-1)), prev=prev, corrmatr = corrmatr)$Pstage[1]- a(tau[1])
      }else{
        pwernest(crit = c(crits[1:(j-1)],cr,rep(0, m-j)), prev=prev, corrmatr = corrmatr)$Pstage[1]- (a(tau[j])-a(tau[j-1]))
      }
    }
    crits[j] <- uniroot(f, interval = c(0,5))$root
  }
  return(crits)
}
#power
pwpnest <- function(p, prev, delta, N, spendingfct = NULL, alpha = 0.025){
  m <- length(prev)
  nu <- delta*sqrt((1:m)*prev*N)
  corrmatr <- corrmatnest(prev)
  if(is.null(spendingfct)){
    crits <- critWTnest(p=p, prev=prev, corrmatr=corrmatr, alpha = alpha)
  }
  else{
    crits <- critESnest(p=p, prev=prev, corrmatr=corrmatr, alpha=alpha, aspend = spendingfct)
  }
  #under assumption that all delta_j > 0:
  pwernest(crit=crits, prev=prev, nu=nu, corrmatr = corrmatr)$PWER
}
pow1nest <- function(p, prev, delta, N, spendingfct = NULL, alpha = 0.025){
  m <- length(prev)
  nu <- delta*sqrt((1:m)*prev*N)
  corrmatr <- corrmatnest(prev)
  if(is.null(spendingfct)){
    crits <- critWTnest(p=p, prev=prev, corrmatr=corrmatr, alpha = alpha)
  }
  else{
    crits <- critESnest(p=p, prev=prev, corrmatr=corrmatr, alpha=alpha, aspend = spendingfct)
  }
  #under assumption that all delta_j > 0:
  1-pmvnorm(upper = crits, mean = nu, sigma = corrmatr)
}
powernest <- function(p, prev, delta, N, spendingfct = NULL, alpha = 0.025, aspend = NULL, powertype = "pwp"){
  if(powertype == "pwp"){
    pwpnest(p=p, prev=prev, delta=delta, N=N, spendingfct=spendingfct, alpha=alpha)
  }
  else if(powertype == "pow1"){
    pow1nest(p=p, prev=prev, delta=delta, N=N, spendingfct=spendingfct, alpha=alpha)
  }
}
#Nroot:
Nrootnest <- function(p, prev, delta, alpha = 0.025, beta = 0.1, spendingfct = NULL, powertype = "pwp",
                      search.interval = c(1,2000)){
  fzero <- function(N){
    powernest(p=p, prev=prev, delta = delta, N=N, powertype = powertype,
          spendingfct=spendingfct, alpha=alpha)-(1-beta)
  }
  uniroot(fzero, interval = search.interval)$root
}
ASNnest <- function(p, prev, delta, N, alpha, spendingfct = NULL){
  m <- length(prev)
  nu <- delta*sqrt((1:m)*prev*N)
  corrmatr <- corrmatnest(prev)
  if(is.null(spendingfct)){
    crits <- critWTnest(p=p, prev=prev, corrmatr=corrmatr, alpha = alpha)
  }
  else{
    crits <- critESnest(p=p, prev=prev, corrmatr=corrmatr, alpha=alpha, aspend = spendingfct)
  }
  Pj <- numeric(m-1)
  for(j in 2:m){
    Pj<-pmvnorm(upper = crits[1:(j-1)], mean = nu[1:(j-1)], sigma = corrmatr[1:(j-1),1:(j-1)])[1]
  }
  N+N*sum(prev[2:m]*Pj)
}

########################
#Example
########################
#m = 3 nested populations
alpha = 0.025; beta = 0.1
p <- c(0, 0.5)
prev = c(1, 0.6, 0.2)
corrmatr <- corrmatnest(prev)
delta = c(0.2, 0.25, 0.3)
#critical values:
cP <- critWTnest(p=p[2], prev=prev, corrmatr = corrmatr, alpha = alpha) #2.168274 2.168274 2.168274
cOBF <- critWTnest(p=p[1], prev=prev, corrmatr = corrmatr, alpha = alpha) #2.494216 1.971851 1.859079
NP <- Nrootnest(p = cP, prev=prev, delta=delta) #253
NOBF <- Nrootnest(p = cOBF, prev=prev, delta=delta) #245
ASNP <- ASNnest(p=p[2], prev = prev, delta=delta, N = NP, alpha = alpha) #255.1742
ASNOBF <- ASNnest(p=p[1], prev = prev, delta=delta, N = NOBF, alpha = alpha) #246.8021
