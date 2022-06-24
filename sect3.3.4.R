#################################################
#Simulation in Section 3.3.4
#################################################
library(multcomp)
library(gtools)
library(rje)
library(mvtnorm)
##########################################
#Contrast matrix all subsets
##########################################
contrastMatrix <- function(n, frac){
  m <- (permutations(2,n,rep=T)-1)[-1,]
  m1 <- lapply(seq_along(1:nrow(m)), function(i){
    l <- m[i,]*frac
    l
  })
  m2 <- do.call(rbind ,m1)
  C <- cbind(-m2,m2)
  rownames(C) <- apply(m,1, function(l){paste0('g',paste0(which(l>0), collapse=''))})
  C
}
##########################################
#PWER and FWER functions
##########################################
critpwerfct <- function(K, df, piv, R, rnames, alpha = 0.025){
  pSetK <- powerSet(1:K)[-1] #set of all subsets without empty set
  #probabilities for each addend in the PWER-expression:
  pj <- function(crit, j){
    S <- pSetK[unlist(lapply(pSetK, function(x){j %in% x}))] #list of subsets containing j
    names(S) <- unlist(lapply(S, function(w){paste0('g', paste0(w, collapse = ''))}))
    posj <- match(names(S), rnames)
    upper <- rep(crit, length(posj))
    1-pmvt(upper=upper, df=df, corr=R[posj,posj], abseps = 0.001)[1]
  }
  #pwer:
  pwer_zero <- function(crit){
    probs <- unlist(lapply(1:K, function(z){pj(crit=crit, j=z)}))
    return(sum(piv*probs)-alpha)
  }
  uniroot(pwer_zero, interval = c(qnorm(1-alpha),qnorm(1-alpha/2^K)))$root
}
critfwerfct <- function(df, R, alpha = 0.025){
  K = log2(dim(R)[1]+1)
  fwer_zero <- function(crit){
    1-pmvt(upper=rep(crit, dim(R)[1]), corr=R, df = df, abseps = 0.001)[1]-alpha
  }
  uniroot(fwer_zero, interval = c(qnorm(1-alpha),qnorm(1-alpha/2^K)))$root
}
##########################################
#Simulation
##########################################
simul_sun <- function(Nsim = 10^4, cp, cf, alpha = 0.025, theta_A = 0.15, P, K, N = 1056, tau, seed = 4328){
  #Determination of n and pi:
  n <- rep(N/K, K)
  piv <- n/N
  #Defining true effect sizes theta_g, g=1,...,K:
  theta <- numeric(K)
  f <- function(g){(1-tau)+2*tau*(g-1)/(K-K0-1)}
  #K0 (number of nulls):
  K0 <- ifelse((P*K %% 1) == 0, P*K, sample(x = c(floor(P*K), ceiling(P*K)), size = 1, prob = c(ceiling(P*K)-P*K, P*K-floor(P*K))))
  piGplus <- sum(piv[1:(K-K0)])
  thetaGplus <- theta_A/piGplus
  if(K0>=1){theta[(K-K0+1):K] <- 0}
  if(K-K0 > 1){
    thetastar <- thetaGplus/sum(piv[1:(K-K0)]/piGplus*f(1:(K-K0)))
    theta[1:(K-K0)] <- thetastar*f(1:(K-K0))
  }
  else if(K-K0 == 1){
    theta[1:(K-K0)] <- thetaGplus
  }
  #contrast matrix:
  C <- contrastMatrix(n=K, frac=piv)
  Diagonal <- diag(1/rowSums(C[,(K+1):(2*K)])); rownames(Diagonal) <- rownames(C)
  C <- Diagonal%*%C
  #thetaS:
  thetaS <- as.numeric(C[,(K+1):(2*K)]%*%theta); names(thetaS) <- rownames(C)
  #design matrix:
  v <- as.factor(rep(1:(2*K), each = N/(2*K)))
  X <- model.matrix(~v+0)
  ##simulation of t test statistics:
  nu <- N-2*K #df of multivar. t-distr.
  D <- C%*%solve(t(X)%*%X)%*%t(C)
  D <- diag(1/sqrt(diag(D)))
  R <- D%*%C%*%solve(t(X)%*%X)%*%t(C)%*%D #Correlation matrix of the multivar. t-distr.
  ncp <- thetaS*diag(D) #non-centrality-parameter
  set.seed(seed)
  t <- rmvt(n=Nsim, sigma = R, df = nu, delta = ncp); colnames(t) <- rownames(C)
  #critical value for PWER-control:
  crit.pwer <- cp
  #critical value for FWER-control:
  crit.fwer <- cf
  
  ##Measured quantitites:
  cnt0 <- c("g0", colnames(t))
  g <- function(x){as.numeric(strsplit(sub(".", "", x), "")[[1]])}  #creates vector out of names ("g12" turns into c(1,2))
  Sstar_fwer <- sapply(cnt0[apply(t, 1, function(x){ifelse(max(x)>=crit.fwer, which.max(x)+1, 1)})], 
                       FUN = g, USE.NAMES = FALSE) #list of chosen populations S* under FWER-control
  Sstar_pwer <- sapply(cnt0[apply(t, 1, function(x){ifelse(max(x)>=crit.pwer, which.max(x)+1, 1)})], 
                       FUN = g, USE.NAMES = FALSE) #list of chosen populations S* under FWER-control
  piSstar_fwer <- unlist(lapply(Sstar_fwer, FUN = function(x){sum(piv[x])}))
  piSstar_pwer <- unlist(lapply(Sstar_pwer, FUN = function(x){sum(piv[x])}))
  thetaSstar_fwer <- unlist(lapply(Sstar_fwer, FUN = function(x){ifelse(any(x==0), 0, thetaS[names(thetaS) == paste0("g",paste0(x, collapse=""))])}))
  thetaSstar_pwer <- unlist(lapply(Sstar_pwer, FUN = function(x){ifelse(any(x==0), 0, thetaS[names(thetaS) == paste0("g",paste0(x, collapse=""))])}))
  #success rate
  success_fwer <-  unlist(lapply(Sstar_fwer, FUN = function(x){ifelse(any(x==0), FALSE, TRUE)}))
  success_pwer <-  unlist(lapply(Sstar_pwer, FUN = function(x){ifelse(any(x==0), FALSE, TRUE)}))
  #power (probability of rejecting at least one null where thetaS > 0): 
  pow_fwer <- apply(t, 1, function(x){max(x[thetaS > 0]) >= crit.fwer})
  pow_pwer <- apply(t, 1, function(x){max(x[thetaS > 0]) >= crit.pwer})
  #impact:
  #impact_fwer <- piSstar_fwer*thetaSstar_fwer
  impact_fwer <- unlist(lapply(Sstar_fwer, FUN = function(x){ifelse(any(x==0), 0, sum(piv[x]*theta[x]))}))
  impact_pwer <- piSstar_pwer*thetaSstar_pwer
  #%correct with respect to S*:
  percentage_correct_fwer <- unlist(lapply(Sstar_fwer, FUN = function(x){ifelse(any(x==0), 0, sum(piv[intersect(x, which(theta>0))])/sum(piv[x]))}))
  percentage_correct_pwer <- unlist(lapply(Sstar_pwer, FUN = function(x){ifelse(any(x==0), 0, sum(piv[intersect(x, which(theta>0))])/sum(piv[x]))}))
  #%correct with respect to full pop.:
  percentage_correct_fwer_fp <- percentage_correct_fwer*piSstar_fwer
  percentage_correct_pwer_fp <- percentage_correct_pwer*piSstar_pwer
  #%false with respect to S*:
  percentage_false_fwer <- unlist(lapply(Sstar_fwer, FUN = function(x){ifelse(any(x==0), 0, sum(piv[intersect(x, which(theta<=0))])/sum(piv[x]))}))
  percentage_false_pwer <- unlist(lapply(Sstar_pwer, FUN = function(x){ifelse(any(x==0), 0, sum(piv[intersect(x, which(theta<=0))])/sum(piv[x]))}))
  #%false with respect to full pop.:
  percentage_false_fwer_fp <- percentage_false_fwer*piSstar_fwer
  percentage_false_pwer_fp <- percentage_false_pwer*piSstar_pwer
  
  ##Results:
  reslist <- list()
  results_params <- c(K, tau, P, alpha)
  names(results_params) <- c("Number of pop.", "tau", "fraction true nulls", "alpha")
  results_sim <- matrix(0, nr = 7, nc = 2)
  colnames(results_sim) <- c("PWER", "FWER")
  rownames(results_sim) <- c("Success rate", "Power", "% correct", "% correct (full)", "% false", "% false (full)", "I(S*)")
  results_sim[1,1] <- mean(success_pwer)
  results_sim[1,2] <- mean(success_fwer)
  results_sim[2,1] <- mean(pow_pwer)
  results_sim[2,2] <- mean(pow_fwer)
  results_sim[3,1] <- mean(percentage_correct_pwer)
  results_sim[3,2] <- mean(percentage_correct_fwer)
  results_sim[4,1] <- mean(percentage_correct_pwer_fp)
  results_sim[4,2] <- mean(percentage_correct_fwer_fp)
  results_sim[5,1] <- mean(percentage_false_pwer)
  results_sim[5,2] <- mean(percentage_false_fwer)
  results_sim[6,1] <- mean(percentage_false_pwer_fp)
  results_sim[6,2] <- mean(percentage_false_fwer_fp)
  results_sim[7,1] <- mean(impact_pwer)
  results_sim[7,2] <- mean(impact_fwer)
  
  reslist$parameters <- results_params
  reslist$results <- results_sim
  return(reslist)
}
