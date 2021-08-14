############################
#Design objects:
############################
library(multcomp)
library(gtools)

#help functions:
getPops <- function(object){
  return(union(colnames(object$n),colnames(object$test)))
}
intersectUV <- function(u,v){
  res <- intersect(unlist(strsplit(u, split = "v")), unlist(strsplit(v, split = "v")))
  paste(res, collapse = "v")
}
#does pop have no overlap with any population in popvec?
noints <- function(pop, popvec){
  tmplog <- sapply(popvec, function(x){intersectUV(u=x, v=pop)})
  return(all(tmplog==""))
}
#rft: relevant for theta
rft <- function(x,y){
  s <- stringi::stri_detect_fixed(str=y, pattern=x)
  return(x[s])
}
vecgrepl <- function(a,b){
  res <- logical(length(b))
  for(i in 1:length(b)){
    res[i] <- grepl(a,b[i], fixed = TRUE)
  }
  return(res)
}
breakU <- function(u){
  unlist(strsplit(u, split = "v"))
}
gets <- function(u, cond=FALSE){
  if(!cond){
    as.numeric(substring(u, first = nchar(u)-1, last = nchar(u)-1))
  }
  else{
    as.numeric(substring(u, first = nchar(u)-3, last = nchar(u)-3))
  }
}
gets2 <- function(u){
  if(is.cond(u)){
    as.numeric(substring(u, first = nchar(u)-3, last = nchar(u)-3))
  }
  else{
    as.numeric(substring(u, first = nchar(u)-1, last = nchar(u)-1))
  }
}
getcond <- function(u){
  if(is.cond(u)){
    as.numeric(substring(u, first = nchar(u)-1, last = nchar(u)-1))
  }
}
rmcond <- function(u){
  if(is.cond(u)){
    paste0(substring(text=u, first=1, last=nchar(u)-3),")")
  }
}
#function for pi_U, U union of partitions:
piU <- function(object, u){
  sum(object$prev[unlist(strsplit(u, split = "v"))])
}
rms <- function(u, cond = FALSE){
  if(!cond){
    substring(u, first = 1, last = nchar(u)-3)
  }
  else{
    substring(u, first = 1, last = nchar(u)-5)
  }
}
rms2 <- function(u){
  if(!is.cond(u)){
    substring(u, first = 1, last = nchar(u)-3)
  }
  else{
    substring(u, first = 1, last = nchar(u)-5)
  }
}
nU <- function(object, u){
  l <- object$n[,unlist(strsplit(rms2(u), split = "v"))]
  if(is.vector(l)){
    return(l)
  }
  else{
    return(rowSums(l))
  }
}
nU2 <- function(object, u){
  l <- object$n[,unlist(strsplit(u, split = "v"))]
  if(is.vector(l)){
    return(l)
  }
  else{
    return(rowSums(l))
  }
}
wUs <- function(object, u, s, from = 1, remove.s = TRUE){
  if(!remove.s){
    w <- sqrt(nU2(object=object, u=u)[from:s]/sum(nU2(object=object, u=u)[from:s]))
  }
  else{
    w <- sqrt(nU(object=object, u=u)[from:s]/sum(nU(object=object, u=u)[from:s]))
  }
  w/sqrt(sum(w^2))
}
rhoUV <- function(object, u,v){
  piU(object=object, u=intersectUV(u,v))/sqrt(piU(object=object,u=u)*piU(object=object, u=v))
}
is.UinVector <- function(u, vec){
  which(rms2(u) == names(vec))
}
testasvec <- function(testmatrix, from){
  stgs <- from:nrow(testmatrix)
  tnames <- colnames(testmatrix)
  res <- c(testmatrix[stgs,])
  names(res) <- c(adds(u=tnames, s=stgs))
  return(res)
}
adds <- function(u,s){
  f <- function(u){sapply(s, function(x){paste0(u, "(", x, ")")})}
  sapply(u, f)
}
#function getweights(design, from, to):
getweights <- function(design, from = 1){
  cnames <- names(testasvec(testmatrix = design$test, from = from))
  f <- function(x){
    wUs(object = design, u = x, s = gets(x), from = from)
  }
  L <- lapply(cnames, f)
  names(L) <- cnames
  return(L)
}
changetoconds <- function(u, s){
  left <- substring(u, first=1, last = nchar(u)-1)
  paste0(left, "|",s, ")")
}
extendn <- function(object){
  notinn <- !is.element(colnames(object$test), colnames(object$n))
  testv <- colnames(object$test[,notinn])
  L <- setNames(vector("list", length(testv)), testv)
  for(i in 1:length(testv)){
    L[[i]] <- nU2(object=object, u=testv[i])
  }
  cbind(object$n, do.call(cbind, L))
}
is.cond <- function(x){
  unlist(lapply(strsplit(x=x, split=""), function(y){"|" %in% y}))
}
nstages <- function(d){
  nrow(d$test)
}
#constructor function for a group sequential design consisting of m populations and K stages:
design <- function(m, I, K){
  if(!is.matrix(I) | dim(I)[1]!=m | dim(I)[2]!=m){
    return("ERROR: I is not an (m x m) - matrix.")
  }
  else if(K < 1 | !is.integer(K)){
    return("ERROR: K >= 1 must be an integer.")
  }
  else if(m < 1 | !is.integer(m)){
    return("ERROR: m >= 1 must be an integer.")
  }
  else if(length(piv)<m | lenght(piv)>2^m-1){
    return("ERROR: piv must be a vector with m <= length(piv) <= 2^m - 1.")
  }
  if(is.null(colnames(I))& !is.null(rownames(I))){
    colnames(I) <- rownames(I)
  }
  else if(is.null(rownames(I)) & !is.null(colnames(I))){
    rownames(I) <- colnames(I)
  }
  else if(is.null(rownames(I)) & is.null(colnames(I))){
    colnames(I) <- rownames(I) <- 1:m
  }
  #function to find each partition of the full population P:
  partition <- function(x){
    if(length(x) == 1){
      return(TRUE)
    }
    else{
      return(all(apply(combn(x, 2),2,function(y){I[y[1],y[2]]})==1))
    }
  }
  pSet <- powerSet(x=1:m, m=m)[-1]
  names(pSet) <- character(length(pSet))
  names(pSet) <- unlist(lapply(pSet, function(x){paste0("P{", paste0(paste(x), collapse = ","), "}")}))
  reducedn <- unlist(lapply(pSet, partition))
  reducednames <- names(reducedn[reducedn])
  pi_partition <- numeric(length(reducednames))
  for(i in 1:length(reducednames)){
    pi_partition[i] <- readline(paste0("Population size of ", reducednames[i],":"))
  }
  pi_partition <- as.numeric(pi_partition)
  if(sum(pi_partition)!=1){
    return("ERROR: The sum of all pi_J must be equal to 1.")
  }
  pSetnames <- unlist(lapply(powerSet(reducednames, m=length(reducednames))[-1],
                      function(x){paste0(x, collapse = "v")}))
  test_in_P <- numeric(length(pSetnames))
  names(test_in_P) <- pSetnames
  for(i in 1:length(pSetnames)){
    test_in_P[i] <- readline(paste0("Is testing the hypothesis H_J: theta_J <= 0 planned in ",pSetnames[i], " (1 for yes, 2 for no)?"))
  }
  test_in_P <- names(test_in_P[test_in_P == 1])
  
}

#Correlation matrix function
corrpi <- function(object, wherecond = NULL, stagescond = NULL, knew = NULL, sigma2 = NULL){
  condition <- !is.null(wherecond) & !is.null(knew) & !is.null(stagescond)
  knew <- ifelse(!condition, 1, knew)
  testv <- testasvec(testmatrix = object$test, from = 1)
  #definition names of correlation matrix (cmnames)
  if(condition){
    cpops <- mapply(adds, wherecond, stagescond)
    cpops <- sapply(cpops, function(x){changetoconds(u=x, s=knew-1)})
    names(cpops) <- NULL
    cmnames <- c(names(testv[testv]), unlist(cpops))
  }
  else{
    cmnames <- names(testv[testv])
  }
  corrmat <- diag(length(cmnames))
  colnames(corrmat) <- rownames(corrmat) <- cmnames
  #filling the matrix
  cn <- combn(cmnames, 2)
  f <- function(x){
    iscondx <- c(is.cond(x[1]),is.cond(x[2]))
    w <- list()
    if(!condition){
      from <- c(1,1)
      for(i in 1:2){
        w[[i]] <- wUs(object=object, u=x[i], s=gets(x[i]), from=from[i])
      }
    }
    else{
      from <- numeric(2)
      for(i in 1:2){
        if(iscondx[i]){
          from[i] <- ifelse(gets(x[i], cond=T)<knew, 1, knew)
          w[[i]] <- wUs(object=object, u=rms2(x[i]), s=gets2(x[i]), from=from[i], remove.s = FALSE)
        }
        else{
          from[i] <- 1
          w[[i]] <- wUs(object=object, u=x[i], s=gets(x[i]), from=from[i])
        }
      }
      
    }
    rho_uv <- rhoUV(object=object, rms2(x[1]),rms2(x[2]))
    pos <- intersect(names(w[[1]]), names(w[[2]]))
    corrmat[x[1],x[2]] <<- rho_uv*sum(w[[1]][pos]*w[[2]][pos])
  }
  apply(cn, 2, f)
  corrmat <- corrmat + t(corrmat) - diag(length(cmnames))
  return(corrmat)
}

#function to compute critical values that satisfy PWER = alpha:
critpwer <- function(object, n, corr, WT, alpha, tau.type = "individual"){
  if(all(colnames(n)==colnames(object$prev))){
    #add n to object:
    object$n <- n
  }
  nU <- function(u){
    l <- object$n[,unlist(strsplit(u, split = "v"))]
    if(is.vector(l)){
      return(l)
    }
    else{
      return(rowSums(l))
    }
  }
  vecgrepl <- function(a,b){
    res <- logical(length(b))
    for(i in 1:length(b)){
      res[i] <- grepl(a,b[i], fixed = TRUE)
    }
    return(res)
  }
  vecgets <- function(b){
    gets <- function(u){
      as.numeric(substring(u, first = nchar(u)-1, last = nchar(u)-1))
    }
    res <- numeric(length(b))
    for(i in 1:length(b)){
      res[i] <- gets(b[i])
    }
    return(res)
  }
  #tau overall
  tau.o <- numeric(nrow(object$n))
  for(i in 1:nrow(object$n)){
    tau.o[i] <- sum(object$n[1:i,])/sum(object$n)
  }
  #tau individual
  tau.i <- n.u <- matrix(0, nr=nrow(object$test), nc=ncol(object$test))
  colnames(tau.i) <- colnames(n.u) <- colnames(object$test)
  rownames(tau.i) <- rownames(n.u) <- rownames(object$test)
  for(name in colnames(object$test)){
    n.u[,name] <- nU(name)
  }
  for(i in 1:nrow(object$test)){
    for(k in 1:ncol(object$test)){
      tau.i[i,k] <- sum(n.u[1:i,k])/sum(n.u[,k])
    }
  }
  pj <- function(j, c1){
    cnames <- colnames(corr)
    posj <- vecgrepl(a=j, b=cnames)
    if(tau.type=="overall"){
        crit <- c1*(tau.o/tau.o[1])^(WT-0.5)
        upper <- crit[vecgets(cnames[posj])]
        1-pmvnorm(upper=upper, corr = corr[posj,posj])[1]
    }
    else if(tau.type=="individual"){
      rms <- function(u){
        substring(u, first = 1, last = nchar(u)-3)
      }
      gets <- function(u){
        as.numeric(substring(u, first = nchar(u)-1, last = nchar(u)-1))
      }
      crit <- matrix(0, nr=nrow(object$test), nc=ncol(object$test))
      colnames(crit) <- colnames(object$test)
      rownames(crit) <- rownames(object$test)
      for(i in 1:nrow(object$test)){
        for(k in 1:ncol(object$test)){
          crit[i,k] <- c1*(tau.i[i,k]/tau.i[1,k])^(WT-0.5)
        }
      }
      upper <- numeric(0)
      for(name in cnames[posj]){
        upper <- c(upper, crit[gets(name),rms(name)])
      }
      1-pmvnorm(upper=upper, corr = corr[posj,posj])[1]
    }
  }
  pwer <- function(c1){
    pjvec <- numeric(length(object$prev))
    for(j in 1:length(pjvec)){
      pjvec[j] <- pj(j=names(object$prev)[j], c1=c1)
    }
    sum(object$prev*pjvec) - alpha
  }
  c1 <- uniroot(pwer, interval = c(qnorm(1-alpha), 5))$root
  if(tau.type=="overall"){
    return(c1*(tau.o/tau.o[1])^(WT-0.5))
  }
  else if(tau.type=="individual"){
    return(c1*sweep(tau.i, 2, tau.i[1,],'/')^(WT-0.5))
  }
}


#Example design II:
#oldD
prev <- c(.4,.4,.2); names(prev) = c("P{1}", "P{2}", "P{1,2}")
n <- matrix(c(80,80,40,0,0,40), nr = 2, byrow = TRUE)
colnames(n) <- names(prev)
rownames(n) <- 1:2
test <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE), nr = 2, byrow = TRUE)
colnames(test) <- c("P{1}vP{1,2}", "P{2}vP{1,2}", "P{1,2}")
rownames(test) <- 1:2
oldD <- list(prev = prev, n=n, test = test)

#newD
n <- matrix(c(80,80,40,80,0,40), nr = 2, byrow = TRUE)
colnames(n) <- names(prev)
rownames(n) <- 1:2
test <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE), nr = 2, byrow = TRUE)
colnames(test) <- c("P{1}vP{1,2}", "P{2}vP{1,2}", "P{1,2}", "P{1}")
rownames(test) <- 1:2
newD <- list(prev = prev, n=n, test = test)

knew = 2
oldCrit = c(2.082, 2.082)
theta_rel <- ThetaC(getPops(oldD), getPops(newD))
zint <- c(2, -0.7, 1.5, 2.5, 0.3)
decisionsold <- c(0,0,0,1,0); decisionsnew <- c(0,0,0,0,0)
names(zint) <- names(decisionsold) <- names(decisionsnew) <- colnames(theta_rel)
critAD(oldD, newD, knew, zint, decisionsold, decisionsnew, theta_rel, oldCrit, root.interval = c(-5,5))
#[[1]]
#[1] 0.01448802
#[[2]]
#[1] 2.727659


#EXAMPLE DESIGN I -> II
prev <- c(.4,.4,.2); names(prev) = c("P{1}", "P{2}", "P{1,2}")
n <- matrix(c(43,43,22,43,43,22), nr = 2, byrow = TRUE)
colnames(n) <- names(prev)
rownames(n) <- 1:2
test <- matrix(c(TRUE, TRUE, TRUE, TRUE), nr = 2, byrow = TRUE)
colnames(test) <- c("P{1}vP{1,2}", "P{2}vP{1,2}")
rownames(test) <- 1:2
oldD <- list(prev = prev, n=n, test = test)

#newD
n <- matrix(c(43, 43, 22, 0, 0, 108), nr = 2, byrow = TRUE)
colnames(n) <- names(prev)
rownames(n) <- 1:2
test <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE), nr = 2, byrow = TRUE)
colnames(test) <- c("P{1}vP{1,2}", "P{2}vP{1,2}", "P{1,2}")
rownames(test) <- 1:2
newD <- list(prev = prev, n=n, test = test)

knew = 2
oldCrit = c(2.246, 2.246)
theta_rel <- ThetaC(getPops(oldD), getPops(newD))
zint <- c(1.5, -0.7, 2, 2.38, 0.3)
decisionsold <- c(0,0,0,1,0)
names(zint) <- names(decisionsold) <- colnames(theta_rel)
critAD(oldD=oldD, newD=newD, knew=knew, zint=zint, decisionsold=decisionsold, 
       decisionsnew = decisionsold, oldCrit=oldCrit, root.interval = c(-5,5),
       theta_rel = theta_rel)
