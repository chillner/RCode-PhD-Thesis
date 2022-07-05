##################################
#CRP principle with PWER-control
##################################
source("GSDclassscript.R")
source("GSDDesignI") #to get Nroot function etc. (change to II or III to start with other design) 
library(gtools) #to use functions that compute permutations/combinations
library(doParallel) #to parallelize simulation runs
library(mvtnorm)

##############################################################################################################################################
#Some helping functions to characterize the population structure and help finding the relevant parameter configurations by means of ThetaC
##############################################################################################################################################
#Returns the power set of some set
powerset <- function(x) {
  sets <- lapply(1:(length(x)), function(i) combn(x, i, simplify = F))
  unlist(sets, recursive = F)
}
#Checks whether a certain function exists in the given population structure 
#(If Intersection(P1,P2) and Intersection(P2,P3) exist, but Intersection(P1,P3) doesn't, then P{1,3} shouldn't exist)
existspop <- function(intersections, x){
  cn <- combn(x,2)
  s <- 0
  for(j in 1:ncol(cn)){
    s <- s + intersections[cn[1,j], cn[2,j]]
  }
  return(ifelse(s == ncol(cn), TRUE, FALSE))
}
                 
#Find the partition of the overall population
popsToPartition <- function(intersections, setex = FALSE){
  
  f <- function(x){
    if(length(x) == 1){
      return(x)
    }
    else{
      if(existspop(intersections=intersections, x=x)){
        return(x)
      } else return(NULL)
    }
  }
  pset <- powerset(1:nrow(intersections))
  pset <- lapply(pset, f)
  pset <- pset[lengths(pset) != 0]
  if(!setex){
    return(unlist(lapply(pset, function(x){paste0("{",paste(x, collapse = ","),"}")})))
  }
  else{
    return(unlist(lapply(pset, function(x){paste(x, collapse = " & ")})))
  }
}
                 
#If populations are in form of disj. subpop., find the respective set of overall populations Pj                 
partitionToPops <- function(partitions){
  if(is.character(partitions)){
    gcomma <- grepl(",", partitions)
    gand <- grepl("&", partitions)
    if(any(gcomma)){
      splits <- strsplit(sapply(partitions[gcomma], function(x){substring(x, first = 2, last=nchar(x)-1)}), split = ",")
    }
    else if(any(gand)){
      splits <- strsplit(partitions[gand], split = " & ")
    }
    splits <- lapply(splits, as.numeric)
    pops <- unique(unlist(splits))
    intersections <- diag(length(pops))
    colnames(intersections) <- rownames(intersections) <- pops
    for(i in 1:length(splits)){
      cn <- combn(splits[[i]], 2)
      for(j in 1:ncol(cn)){
        intersections[cn[1,j], cn[2,j]] <- intersections[cn[2,j], cn[1,j]] <- 1
      }
    }
    return(list(mat = intersections, pops = pops))
  }
  
}


#######################################
#Finding adaptive critical value cnew
#######################################
#Function for finding Theta_combined:
ThetaC <- function(pop.old, pop.new){
  #combined string of population names:
  pC <- as.list(union(pop.old, pop.new))
  #populations with no intersections:
  opt <- c("= 0","< 0", "> 0")
  disj.pops <- pC[unlist(lapply(pC, function(x){noints(pop=x, popvec = pC[pC!=x])}))]
  if(length(disj.pops)!=0){
    d.pops.opt <- gtools::permutations(n = 2, r = length(disj.pops), v = opt[c(1,3)], 
                                       repeats.allowed = TRUE)
    d.pops.opt <- as.data.frame(d.pops.opt)
    colnames(d.pops.opt) <- unlist(disj.pops)
  }
  else{
    d.pops.opt <- list()
  }
  #Find all relevant intersections:
  new.pC <- pC[!(pC %in% disj.pops)] 
  repeat {
    tmp.pC <- list()
    co <- combn(x = new.pC, m = 2)
    for(i in 1:ncol(co)){
      intUV <- intersectUV(u = as.character(co[1,i]), v = as.character(co[2,i]))
      if(intUV %in% pC){
        tmp.pC[[i]] <- intUV
      }
    }
    tmp.pC <- as.list(unique(unlist(tmp.pC)))
    if(length(tmp.pC)==0){
      break
    }
    new.pC <- tmp.pC
  } 
  #remaining combined parameters:
  comb.pC <- setdiff(unlist(pC), union(new.pC, disj.pops))
  #Create minimal parameter set needed to find all configs:
  rel.opt<- gtools::permutations(n = 3, r = length(new.pC), v = opt, repeats.allowed = TRUE)
  rel.opt <- as.data.frame(rel.opt)
  colnames(rel.opt) <- new.pC
  rel.opt <- rel.opt[apply(rel.opt, 1, function(z){!all(z=="< 0"|z=="= 0") & !all(z=="> 0")}),]
  rel.opt <- rbind(rep("= 0", length(new.pC)), rel.opt)
  #all configs:
  f.comb <- function(z){
    counter <- matrix(FALSE, nr = 2, nc = length(comb.pC))
    colnames(counter) <- comb.pC
    rownames(counter) <- c("= 0", "> 0")
    L <- setNames(vector("list", ncol(counter)), comb.pC)
    for(j in 1:ncol(counter)){
      rfp <- rft(x=colnames(rel.opt), y=comb.pC[j])
      if(!all(z[rfp] == "= 0" | z[rfp] == "< 0")){
        counter["> 0", j] <- TRUE
        if(any(z[rfp] == "< 0")){
          counter["= 0",j] <- TRUE 
          #z[rfp][z[rfp] == "< 0"] <- "= 0"
        }
      }
      else{
        counter["= 0", j] <- TRUE
        #z[rfp][z[rfp] == "< 0"] <- "= 0"
      }
      L[[j]] <- rownames(counter)[counter[,j]]
    }
    #limit cases:
    z[z == "< 0"] <- "= 0"
    #all possibilities for comb.pC:
    comb.poss <- expand.grid(L)
    if(!is.null(nrow(d.pops.opt))){
      tmpres <- matrix(z, nr = nrow(d.pops.opt), nc = length(z), byrow = TRUE)
      tmpres <- as.data.frame(cbind(tmpres, d.pops.opt))
      colnames(tmpres) <- c(names(z), disj.pops)
      ntmp <- comb.poss[rep(1:nrow(comb.poss), 1, each = nrow(tmpres)),]
      ntmp <- ntmp[kronecker(1:nrow(comb.poss), c(0,nrow(comb.poss)), "+"),]
      res <- cbind(tmpres, ntmp)
    }
    else{
      res <- matrix(z, nr = nrow(comb.poss), nc = length(z), byrow = TRUE)
      res <- as.data.frame(cbind(res, comb.poss))
      colnames(res) <- c(names(z), colnames(comb.poss))
    }
    return(res)
  }
  tmplist <- list()
  for(i in 1:nrow(rel.opt)){
    tmplist[[i]] <- f.comb(rel.opt[i,])
    
  }
  allthetas <- unique(do.call(rbind, tmplist))[,unlist(pC)]
  rownames(allthetas) <- as.character(1:nrow(allthetas))
  return(allthetas)
}
#Function for finding cp.new
cp_new_loc <- function(d, corr, decisions, crit.adj, theta){
  #correlation matrix
  cnames <- colnames(corr)
  #all where PJ(...) = 1:
  countas1 <- logical(length(cnames))
  names(countas1) <- cnames
  #all test statistics needed for computation of P_theta_J(Z_U > c | zint) for each J
  IJ0theta <- function(J){
    posJ <- cnames[vecgrepl(a=J, b=cnames)]
    posJ <- posJ[sapply(posJ, function(x){is.cond(x) && (rms2(x) %in% colnames(d$test))})]
    if(length(posJ) == 0){character(0)}
    else{
      I0 <- sapply(posJ, function(x){which(rms2(x) == names(theta))})
      d0 <- sapply(posJ, function(x){which(rms2(x) == names(decisions))})
      f <- function(x){
        if(theta[I0[x]] == "= 0"){
          if(decisions[d0[x]] == 1){
            countas1[x] <<- TRUE
            TRUE
          }
          else{!(gets2(x) < knew)}
        }
        else{FALSE}
      }
      return(posJ[sapply(posJ, f)])
    }
  }
  #IJ0theta <- function(J){
    #posJ <- cnames[vecgrepl(a=J, b=cnames)]
    #posJ <- posJ[is.cond(posJ)]
    #if(length(posJ) == 0){return(character(0))}
    #else{
    #  I0 <- sapply(posJ, function(x){which(rms2(x) == names(theta))})
    #  d0 <- sapply(posJ, function(x){which(rms2(x) == names(decisions))})
    #  t0 <- sapply(posJ, function(x){d$test[gets2(x), rms2(x)]})
    #  return(posJ[theta[I0] == "= 0" & (decisions[d0] == 1 | t0)])
    #}
  #}
  #all J that need to be summed up over
  J0theta <- sapply(names(d$prev), function(J){length(IJ0theta(J)) != 0})
  if(sum(J0theta) == 0){return(0)}
  #J0theta <- sapply(names(d$prev), function(J){length(IJ0theta(J)) != 0})
  pJ <- function(J){
    pos <- IJ0theta(J)
    corrJ <- corr[pos, pos] #relevant submatrix of corr
    upperJ <- unlist(crit.adj[pos]) #"upper" component of pmvnorm(...)
    d$prev[J]*ifelse(any(countas1[pos]), 1, 1-pmvnorm(upper=upperJ, sigma = corrJ)[1]) 
  }
  return(sum(sapply(names(J0theta)[unname(J0theta)], pJ)))
}
#Function for finding cp.old
cp_old_loc <- function(d, corr, decisions, crit.adj, theta){
  #correlation matrix 
  cnames <- colnames(corr)
  #all where PJ(...) = 1:
  countas1 <- logical(length(cnames))
  names(countas1) <- cnames
  #all test statistics needed for computation of P_theta_J(Union_U Z_U > c | zint) for each J
  IJtheta <- function(J){
    posJ <- cnames[vecgrepl(a=J, b=cnames)]
    posJ <- posJ[sapply(posJ, function(x){rms2(x) %in% colnames(d$test)})]
    if(length(posJ) == 0){character(0)}
    else{
      I0 <- sapply(posJ, function(x){which(rms2(x) == names(theta))})
      d0 <- sapply(posJ, function(x){which(rms2(x) == names(decisions))})
      f <- function(x){
        if(theta[I0[x]] == "= 0"){
          if(decisions[d0[x]]){
            countas1[x] <<- TRUE
            TRUE
          }
          else{!(gets2(x) < knew)}
        }
        else if(theta[I0[x]] == "> 0"){
          !is.cond(x)
          #if(is.cond(x)){FALSE}
          #else{
          #  ifelse(gets2(x) < knew, decisions[d0[x]] == 1, TRUE)
          #}
        }
      }
      return(posJ[sapply(posJ, f)])
    }
  }
  #all J that need to be summed up over
  J0theta <- sapply(names(d$prev), function(J){length(IJtheta(J)) != 0})
  pJ <- function(J){
    pos <- IJtheta(J)
    if(length(pos) == 0){return(0)}
    else{
      corrJ <- corr[pos, pos] #relevant submatrix of corr
      upperJ <- unlist(crit.adj[pos]) #"upper" component of pmvnorm(...)
      return(d$prev[J]*ifelse(any(countas1[pos]), 1, 1-pmvnorm(upper=upperJ, sigma = corrJ)[1]))
    }
  }
  return(sum(sapply(names(J0theta)[unname(J0theta)], pJ)))
  ###################################
  #IJ0theta <- function(J){
  #  posJ <- cnames[vecgrepl(a=J, b=cnames)]
  #  posJ <- posJ[sapply(posJ, function(x){rms2(x) %in% colnames(d$test)})]
    #posJ <- posJ[is.cond(posJ)]
  #  if(length(posJ) == 0){return(character(0))}
  #  else{
   #   I0 <- sapply(posJ, function(x){which(rms2(x) == names(theta))})
  #    d0 <- sapply(posJ, function(x){which(rms2(x) == names(decisions))})
  #    t0 <- sapply(posJ, function(x){if(gets2(x) < knew)   else d$test[gets2(x), rms2(x)]})
   #   return(posJ[theta[I0] == "= 0" & (decisions[d0] == 1 | t0)])
    #}
  #}
  #IJ1theta <- function(J){
   # posJ <- cnames[vecgrepl(a=J, b=cnames)]
    #posJ <- posJ[!is.cond(posJ) & sapply(posJ, function(x){rms2(x) %in% colnames(d$test)})]
    #if(length(posJ) == 0){return(character(0))}
    #else{
     # I1 <- sapply(posJ, function(x){which(rms2(x) == names(theta))})
    #  d1 <- sapply(posJ, function(x){which(rms2(x) == names(decisions))})
    #  t1 <- sapply(posJ, function(x){d$test[gets2(x), rms2(x)]})
    #  posJ[theta[I1] == "> 0" & (t1 | decisions[d1] == 1)]
    #}
  #}
  #all J that need to be summed up over
  #J0theta <- sapply(names(d$prev), function(J){length(union(IJ0theta(J),IJ1theta(J))) != 0})
  #pJ <- function(J){
  #  pos <- union(IJ0theta(J),IJ1theta(J))
  #  corrJ <- corr[pos, pos] #relevant submatrix of corr
  #  upperJ <- unlist(crit.adj[pos]) #"upper" component of pmvnorm(...)
  #  d$prev[J]*ifelse(any(decisions[sapply(pos, rms2)] == 1), 1, 1-pmvnorm(upper=upperJ, sigma = corrJ)) 
  #}
  #sum(sapply(J0theta, pJ))
}
#Function for computing critical values by means of the CRP-principle for PWER-control
critAD <- function(oldD, newD, knew, zint, decisionsold, decisionsnew, oldCrit, theta_rel, solve.option = 1, root.interval){
    #Find ThetaC:
    U_old <-  getPops(oldD); U_new <- getPops(newD)
    #theta_rel <- ThetaC(U_old, U_new)
    K <- nstages(newD)
    #theta_rel <- ThetaC(pop.old = U_old, pop.new = U_new)
    ##correlation matrices:
    #new design
    tav <- testasvec(testmatrix = newD$test, from = knew)
    wherecond <- unique(rms(names(tav[tav])))
    stagescond <- setNames(object = vector("list", length(wherecond)), wherecond)
    for(w in wherecond){
      stagescond[[w]] <- as.numeric(names(newD$test[,w])[knew:K])
    }
    cmat <- corrpi(object = newD, knew = knew, wherecond = wherecond, stagescond = stagescond)
    #old design
    tav <- testasvec(testmatrix = oldD$test, from = knew)
    wherecond <- unique(rms(names(tav[tav])))
    stagescond <- setNames(object = vector("list", length(wherecond)), wherecond)
    for(w in wherecond){
      stagescond[[w]] <- as.numeric(names(oldD$test[,w])[knew:K])
    }
    cmat_old <- corrpi(object = oldD, knew = knew, wherecond = wherecond, stagescond = stagescond)
    #create cadj:
    cadj <- function(crit, design = "new"){
      if(design == "new"){ 
        d <- newD
        cnames <- colnames(cmat)
      }
      else if(design == "old"){
        d <- oldD
        cnames <- colnames(cmat_old)
      }
      cr <- function(x){
        if(!is.cond(x)){
          return(oldCrit[gets2(x)])
        }
        else{
          w <- wUs(object = d, u = x, s = gets2(x), from = 1)
          if(design == "new"){
            return((crit - w[1:(knew-1)]%*%zint[[rms2(x)]][as.numeric(names(w)[1:(knew-1)])])/sqrt(1-sum(w[1:(knew-1)]^2)))
          } 
          else if(design == "old"){
            return((oldCrit[gets2(x)] - w[1:(knew-1)]%*%zint[[rms2(x)]][as.numeric(names(w)[1:(knew-1)])])/sqrt(1-sum(w[1:(knew-1)]^2)))
          }
        }
      }
      crits <- sapply(cnames, cr)
      names(crits) <- cnames
      crits
    }
    #For each theta in thetaC find cp_old
    cpold <- list()
    #create a new oldCrit vector such that only "apply" has to be applied to compute each cp.old(theta):
    ocrit <- cadj(1, "old")
    cpoldtheta <- function(x){
      cp_old_loc(d=oldD, decisions=decisionsold, theta=x, corr=cmat_old, crit.adj = ocrit)
    }
    #cpold values
    cpold <- apply(theta_rel, 1, cpoldtheta)
    if(solve.option == 1){
      alpha_C <- min(cpold)
      fcpnew <- function(x){
        g <- function(y){
          cp_new_loc(d=newD, crit.adj=cadj(crit=x), theta=y, corr=cmat, decisions=decisionsnew)
        }
        max(apply(theta_rel, 1, g))
      }
      return(list(min(cpold), uniroot(function(z){fcpnew(z)-alpha_C}, interval = root.interval)$root))
    }
    else if(solve.option == 2){
      crnew <- numeric(nrow(theta_rel))
      for(i in 1:nrow(theta_rel)){
        cnew[i] <- uniroot(function(cr){fcpnew[[i]](cr)-cpoldtheta[i](cr)}, interval = root.interval)$root
      }
      return(max(cnew))
    }
}
############
#test:
#############
prev <- c(0.4,0.3,0.3); names(prev) <- c("P{1}", "P{2}", "P{1,2}")
n <- matrix(c(40,30,30,0,0,30), nr = 2, byrow = TRUE)
colnames(n) <- names(prev)
rownames(n) <- 1:2
test <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE), nr = 2, byrow = TRUE)
colnames(test) <- c("P{1}vP{1,2}", "P{2}vP{1,2}", "P{1,2}")
rownames(test) <- 1:2
oldD <- list(prev = prev, n = n, test = test)
wherecond <- c("P{1}","P{1,2}")
stagescond <- list(2,2)
#New design: H1 rejected, H2 retained => ignore rejection of H1 and test in P{1} and P{1,2} at stage 2
n.new <- n
n.new[2,1] <- 40
test.new <- matrix(c(T,T,F,F,F,F,F,T,T,F), nr = 2, byrow=T)
colnames(test.new) <- c("P{1}vP{1,2}", "P{2}vP{1,2}", "P{1,2}", "P{1}", "P{2}")
rownames(test.new) <- 1:2
newD <- list(prev = prev, n = n.new, test = test.new)
cmat <- corrpi(object = newD, wherecond = wherecond, stagescond=stagescond, knew = 2)
decisionsold <- c(1,0,0,0,0); decisionsnew <- c(0,0,0,0,0)
names(decisionsold) <- names(decisionsnew) <- c("P{1}vP{1,2}", "P{2}vP{1,2}", "P{1,2}", "P{1}", "P{2}")
oldCrit <- c(2.24,2.24)
root.interval = c(1,5)
zint <- list(3,1.5, 3,3); names(zint) <- names(decisions)
ThetaC(pop.old = getPops(oldD), pop.new= getPops(newD))

#####################################
#Simulation function
#####################################

#critical values for single stage design in Chapter 6.3
critSingle <- function(prev, alpha){
  piv <- prev[1:2]+prev[3]
  Corr <- matrix(c(1,prev[3]/sqrt(piv[1]*piv[2]), sqrt(prev[1]/piv[1]), 0, sqrt(prev[3]/piv[1]),
                   0, 1, 0, sqrt(prev[2]/piv[2]), sqrt(prev[3]/piv[2]),
                   0, 0, 1, 0, 0,
                   0, 0, 0, 1, 0,
                   0, 0, 0, 0, 1), nr = 5, byrow = T)
  Corr <- t(Corr) + Corr -diag(5)
  pwer <- function(crit){
    prev[1]*(1-pmvnorm(upper=rep(crit,2), corr = Corr[c(1,3),c(1,3)])) +
      prev[2]*(1-pmvnorm(upper=rep(crit,2), corr = Corr[c(2,4),c(2,4)])) +
      prev[3]*(1-pmvnorm(upper=rep(crit,3), corr= Corr[c(1,2,5),c(1,2,5)])) - alpha
  }
  uniroot(pwer, c(0,5))$root
}

#Function that runs simulations shown in Chapter 6.3
simAD <- function(Nsim, pops, prev, delta, knew = 2, seed = 321, theta_rel = theta_rel, N = 100,
                  newpops = NULL, lbound = qnorm(0.9), alpha = 0.025, beta=0.2, optparD1 = 0.5){
  
  set.seed(seed)
  p <- popsToPartition(pops)
  names(prev) <- p
  piv <- prev[1:2] + prev[3]
  names(delta) <- c("{1}","{2}","{1,2}")
  Delta <- c((delta[1]*prev[1]+delta[3]*prev[3])/(prev[1]+prev[3]),(delta[2]*prev[2]+delta[3]*prev[3])/(prev[2]+prev[3]))
  names(Delta) <- c("{1}v{1,2}","{2}v{1,2}")
  ########################
  ##INITIAL DESIGN
  ########################
  #Design I oldDesign object:
  D1 <- list(prev=prev, n = NULL, test = matrix(c(T,T,T,T), nr=2,
                                                dimnames = list(1:2, c("{1}v{1,2}","{2}v{1,2}"))))
  if(is.null(N)){
    ND1 <- Nroot(p=optparD1, prev=prev, delta=Delta, alpha = alpha,
                 beta = beta, spendingfct = NULL, tau.option="overall", powertype = "pwp",
                 search.interval = c(1,50000))
  } else{ND1 <- N}
  nD1 <- ND1*matrix(c(prev,prev), nr=2, byrow = T, dimnames = list(1:2, c("{1}","{2}","{1,2}")))
  D1$n <- nD1
  ##Simulate test statistics in disjoint subgroups P_J
  ZsimD1 <- list()
  for(k in 1:2){
    ZsimD1[[k]] <- list()
    for(i in 1:3){
      ZsimD1[[k]][[i]] <-  rnorm(Nsim, delta[i]*sqrt(D1$n[k,i]))
    }
    names(ZsimD1[[k]]) <- c("{1}","{2}", "{1,2}")
  }
  #combined test stats per stage:
  Zaccnames <- c("{1}v{1,2}(1)","{1}v{1,2}(2)", "{2}v{1,2}(1)", "{2}v{1,2}(2)")
  ZaccD1 <- as.data.frame(matrix(0, nr=Nsim, nc=4))
  colnames(ZaccD1) <- Zaccnames
  Z11comb <-sqrt(prev[1]/piv[1])*ZsimD1[[1]][[1]] + sqrt(prev[3]/piv[1])*ZsimD1[[1]][[3]]
  Z21comb <-sqrt(prev[2]/piv[2])*ZsimD1[[1]][[2]] + sqrt(prev[3]/piv[2])*ZsimD1[[1]][[3]]
  Z12comb <-sqrt(prev[1]/piv[1])*ZsimD1[[2]][[1]] + sqrt(prev[3]/piv[1])*ZsimD1[[2]][[3]]
  Z22comb <-sqrt(prev[2]/piv[2])*ZsimD1[[2]][[2]] + sqrt(prev[3]/piv[2])*ZsimD1[[2]][[3]]
  ##accrued stats:
  ##accrued stats for ZJ^(2)
  ZaccJ <- data.frame('{1}' = rep(0,Nsim), '{2}' = rep(0,Nsim), '{1,2}' = rep(0,Nsim))
  ZaccJ[,1] <- sqrt(.5)*ZsimD1[[1]][[1]] + sqrt(.5)*ZsimD1[[2]][[1]]
  ZaccJ[,2] <- sqrt(.5)*ZsimD1[[1]][[2]] + sqrt(.5)*ZsimD1[[2]][[2]]
  ZaccJ[,3] <- sqrt(.5)*ZsimD1[[1]][[3]] + sqrt(.5)*ZsimD1[[2]][[3]]
  #weights:
  w <- data.frame('{1}v{1,2}' = wUs(object=D1, u="{1}v{1,2}", s=2, remove.s = F),
                  '{2}v{1,2}' = wUs(object=D1, u="{2}v{1,2}", s=2, remove.s = F))
  ZaccD1[,'{1}v{1,2}(1)'] <- Z11comb
  ZaccD1[,'{2}v{1,2}(1)'] <- Z21comb
  ZaccD1[,'{1}v{1,2}(2)'] <- w[1,1]*Z11comb+w[2,1]*Z12comb
  ZaccD1[,'{2}v{1,2}(2)'] <- w[1,2]*Z21comb+w[2,2]*Z22comb
  #optimal critical values of old design:
  oldCritD1 <- critpwer(object=D1, n=D1$n, corr=corrpi(D1), WT=optparD1, alpha=alpha, tau.type = "overall")
  #stage 1/2 results (GSD):
  s1res <- ZaccD1[,c(1,3)] >= oldCritD1[1]
  s2res <- ZaccD1[,c(2,4)] >= oldCritD1[2]
  #interim data:
  z_int <- cbind(ZaccD1[,1], ZaccD1[,3], ZsimD1[[1]][[1]], ZsimD1[[1]][[2]], ZsimD1[[1]][[3]])
  #decisions:
  dec_old <- cbind(ifelse(s1res[,1],1,0),ifelse(s1res[,2],1,0), matrix(0, nr=Nsim, nc=3))
  #names z_int/dec_old
  colnames(z_int) <- colnames(dec_old) <- c("{1}v{1,2}","{2}v{1,2}", "{1}","{2}","{1,2}")
  #####################################
  ##RESULTS OF SIMULATIONS
  #####################################
  resultmatrix_uncond <- resultmatrix <- as.data.frame(matrix(0, nr = 6, nc = 8,
                                                              dimnames = list(c("GSD1","AD1","AD2","AD3", "AD4", "single"),
                                                                              c("H1", "H2", "H{1}", "H{2}", "H{1,2}",
                                                                                "H1*", "H2*", "E(N)"))))
  namendr <-  c("{1}v{1,2}", "{2}v{1,2}", "{1}", "{2}", "{1,2}","E(N)")
  #########SINGLE STAGE DESIGN##########
  cS <- critSingle(prev=prev, alpha = alpha)
  resultmatrix[6,] <- c(colMeans(cbind(ZaccD1[,c('{1}v{1,2}(2)','{2}v{1,2}(2)')] >= cS, ZaccJ >= cS,
                                       ZaccD1[,'{1}v{1,2}(2)'] >= cS|ZaccJ[,1] >= cS&ZaccJ[,3] >= cS,
                                       ZaccD1[,'{2}v{1,2}(2)'] >= cS|ZaccJ[,2] >= cS&ZaccJ[,3] >= cS)), 2*ND1)
  
  resultmatrix_uncond[6,] <- resultmatrix[6,]
  rm(ZaccD1, ZsimD1, Z11comb, Z12comb, Z21comb, Z22comb, ZaccJ, w, cS, Zaccnames)
  #########GSD1############
  colmeansnoNA <- function(x, na.rm =TRUE){
    whichAllNA <- apply(x, 2, function(y) all(is.na(y)))
    res <- numeric(ncol(x))
    res[!whichAllNA] <- colMeans(x[,!whichAllNA], na.rm=na.rm)
    return(res)
  }
  resultmatrix[1,1:2] <- resultmatrix_uncond[1,1:2] <- colMeans(cbind(s1res[,1] | s2res[,1], s1res[,2] | s2res[,2]))
  resultmatrix[1,8] <- resultmatrix_uncond[1,8] <- mean(apply(s1res, 1, function(x){ 
                                       if(x[1]&&x[2]) ND1 
                                       else if(!x[1]&&!x[2]) 2*ND1 
                                       else if(!x[1]&&x[2]) ND1*(1+piv[1]) 
                                       else if(x[1]&&!x[2]) ND1*(1+piv[2])}))
  resultmatrix_uncond[1,6:7] <- resultmatrix_uncond[1,1:2] 
  resultmatrix[1,6:7] <- resultmatrix[1,1:2]
  #########AD1#############
  AD1 <- function(j){
    res <- rep(NA,6)
    names(res) <- namendr
    if(s1res[j,1] && s1res[j,2]){
      return(c(1, 1, rep(NA,3), ND1))
    }
    else if(!s1res[j,1] & !s1res[j,2]){
      return(c(0, 0, rep(NA,3), ND1))
    }
    else{
      whichJ <- list()
      for(i in 1:2){
        if(!s1res[j,i]){
          whichJ[[i]] <- grep(as.character(i), names(res)[3:5], value = T)
        }
      }
      whichJ <- unique(unlist(whichJ))
      #ZwJ <- z_int[j,whichJ]
      if(all(z_int[j,whichJ]<lbound)){
        return(c(s1res[j,],rep(NA,3), ND1))
      }
      else if(all(z_int[j,whichJ] >= lbound)){
        #if(setequal(whichJ, names(res)[3:5])){
        #  return(c(s2res[j,], rep(NA,3)))
        #}
        if(setequal(whichJ, c("{1}","{1,2}"))){
          return(c(s2res[j,1], s1res[j,2], rep(NA,3), ND1+piv[1]*ND1))
        }
        else if(setequal(whichJ, c("{2}","{1,2}"))){
          return(c(s1res[j,1], s2res[j,2], rep(NA,3), ND1+piv[2]*ND1))
        }
      }
      else{
        w <- which(z_int[j,whichJ] >= lbound)
        res[1:2] <- s1res[j,]
        newDAD1 <- list(prev=prev, n=NULL, test=NULL)
        newDAD1$test <- cbind(D1$test, matrix(rep(F,6), nr=2, dimnames = list(c(1:2), c("{1}","{2}","{1,2}"))))
        newDAD1$test <- newDAD1$test[,c("{1}v{1,2}", "{2}v{1,2}", names(w))]
        newDAD1$test[2, names(w)] <- TRUE
        newDAD1$test[2, 1:2] <- FALSE
        newDAD1$n <- D1$n
        newDAD1$n[2,names(w)] <- sum(newDAD1$n[2,])*(newDAD1$prev[names(w)]/sum(newDAD1$prev[names(w)]))
        newDAD1$n[2, setdiff(colnames(newDAD1$n),names(w))] <- 0
        #expected N:
        #res[6:11] <- c(piv[1]*ND1,piv[2]*ND1, ifelse(newDAD1$n[2,] == 0, NA, colSums(newDAD1$n)), sum(newDAD1$n))
        res[6] <- sum(newDAD1$n)
        #generate new stage 2 z-scores in PJ[w]
        ZsimD1_2 <- numeric(length(w)); names(ZsimD1_2) <- names(w)
        wJ <- setNames(vector("list", length(w)), names(w))
        for(i in names(w)){
          wJ[[i]] <- wUs(object=newDAD1, u=i,s=2,remove.s = FALSE)
          ZsimD1_2[i] <- wJ[[i]][1]*z_int[j,i]  + wJ[[i]][2]* rnorm(1, delta[i]*sqrt(newDAD1$n[2,i]))
        }
        names(res) <- namendr
        res[names(w)] <- ZsimD1_2[names(w)] >= critAD(oldD = D1, newD = newDAD1, knew=knew, zint = z_int[j,], decisionsold = dec_old[j,],
                                             theta_rel=theta_rel, decisionsnew = dec_old[j,], oldCrit = oldCritD1, root.interval = c(-10,10))[[2]]
      }
    }
    return(res)
  }
  vnames <- c("AD1","namendr","s1res", "s2res", "z_int", "dec_old", "D1", "delta", "Delta",
              "Nsim", "knew", "theta_rel", "lbound", "alpha", "prev", "piv", "ND1", "oldCritD1")
  cl <- makeCluster(detectCores()-1, type = 'PSOCK')
  clusterExport(cl, as.list(c(ls_start,vnames)), envir = environment())
  clusterEvalQ(cl, library("mvtnorm", "gtools"))
  registerDoParallel(cl)
  #endres <- t(sapply(1:Nsim, FUN = AD1))
  endres <- do.call(rbind, clusterApply(cl, 1:Nsim, AD1))
  registerDoSEQ()
  closeAllConnections()
  rm(cl, AD1)
  resultmatrix[2,c(1:5,8)] <- colmeansnoNA(endres); resultmatrix[2,6:7] <- resultmatrix[2,1:2]
  endres[is.na(endres)] <- 0
  resultmatrix_uncond[2,c(1:5,8)] <- colMeans(endres); resultmatrix_uncond[2,6:7] <- resultmatrix_uncond[2,1:2]
  rm(endres)
  #########AD2#############
  AD2 <- function(j){
    res <- rep(NA,6)
    names(res) <- namendr
    if(s1res[j,1] && s1res[j,2]){
      return(c(1, 1, rep(NA,3),ND1))
    }
    else if(!s1res[j,1] & !s1res[j,2]){
      return(c(0, 0, rep(NA,3), ND1))
    }
    else{
      res[1:2] <- s1res[j,]
      newDAD2 <- list(prev=prev, n=NULL, test=NULL)
      newDAD2$test <- cbind(D1$test, c(F,T)); colnames(newDAD2$test) <- c("{1}v{1,2}", "{2}v{1,2}","{1,2}")
      newDAD2$test[2, 1:2] <- FALSE
      newDAD2$n <- D1$n
      newDAD2$n[2,] <- c(0,0,sum(D1$n[2,]))
      #expected N
      #res[6:11] <- c(piv[1]*ND1, piv[2]*ND1, rep(NA,2), sum(newDAD2$n[,3]), sum(newDAD1$n))
      res[6] <- sum(newDAD2$n)
      #generate new stage 2 z-score in P{1,2}
      w12 <- wUs(object = newDAD2, u = "{1,2}", s = 2, remove.s = FALSE)
      #ZsimD2_2 <- w12[1]*z_int[j,"{1,2}"] +  w12[2]*rnorm(1, delta[3]*sqrt(newDAD2$n[2,3]))
      names(res) <- namendr
      res["{1,2}"] <- w12[1]*z_int[j,"{1,2}"] + w12[2]*rnorm(1, delta[3]*sqrt(newDAD2$n[2,3])) >= 
                      critAD(oldD = D1, newD = newDAD2, knew=knew, zint = z_int[j,], decisionsold = dec_old[j,],
                      theta_rel=theta_rel, decisionsnew = dec_old[j,], oldCrit = oldCritD1, root.interval = c(-10,10))[[2]]
    }
    return(res)
  }
  vnames <- c("AD2","namendr","s1res", "s2res", "z_int", "dec_old", "D1", "delta", "Delta",
              "Nsim", "knew", "theta_rel", "alpha", "prev", "piv", "ND1", "oldCritD1")
  cl <- makeCluster(detectCores()-1, type = 'PSOCK')
  clusterExport(cl, as.list(c(ls_start,vnames)), envir = environment())
  clusterEvalQ(cl, library("mvtnorm", "gtools"))
  registerDoParallel(cl)
  endres <- do.call(rbind, clusterApply(cl, 1:Nsim, AD2))
  registerDoSEQ()
  closeAllConnections()
  rm(cl, AD2)
  resultmatrix[3,c(1:5,8)] <- colmeansnoNA(endres); resultmatrix[3,6:7] <- resultmatrix[3,1:2]
  endres[is.na(endres)] <- 0
  resultmatrix_uncond[3,c(1:5,8)] <- colMeans(endres); resultmatrix_uncond[3,6:7] <- resultmatrix_uncond[3,1:2]
  rm(endres)
  #########AD3#############
  AD3 <- function(j){
    res <- rep(NA,6)
    names(res) <- namendr
    if(s1res[j,1] && s1res[j,2]){
      return(c(1, 1, rep(NA,3),ND1))
    }
    else if(!s1res[j,1] & !s1res[j,2]){
      return(c(0, 0, rep(NA,3),ND1))
    }
    else{
      res[1:2] <- s1res[j,]
      whichJ <- list()
      for(i in 1:2){
        if(!s1res[j,i]){
          whichJ[[i]] <- grep(as.character(i), names(res)[3:5], value = T)
        }
      }
      whichJ <- unique(unlist(whichJ))
      newDAD3 <- list(prev=prev, n=NULL, test=NULL)
      newDAD3$test <- cbind(D1$test, matrix(rep(F,6), nr=2, dimnames = list(c(1:2), c("{1}","{2}","{1,2}"))))
      newDAD3$test <- newDAD3$test[,c("{1}v{1,2}", "{2}v{1,2}", whichJ)]
      newDAD3$test[2, whichJ] <- TRUE
      newDAD3$test[2, 1:2] <- FALSE
      newDAD3$n <- D1$n
      newDAD3$n[2,whichJ] <- sum(newDAD3$n[2,])*(newDAD3$prev[whichJ]/sum(newDAD3$prev[whichJ]))
      newDAD3$n[2, setdiff(colnames(newDAD3$n),whichJ)] <- 0
      #expected N:
      #res[6:11] <- c(piv[1]*ND1,piv[2]*ND1, ifelse(newDAD3$n[2,]==0, NA, colSums(newDAD3$n)), sum(newDAD3$n))
      res[6] <- sum(newDAD3$n)
      #generate new stage 2 z-scores in PJ[whichJ]
      ZsimD3_2 <- numeric(length(whichJ)); names(ZsimD3_2) <- whichJ
      wJ <- setNames(vector("list", length(whichJ)), whichJ)
      for(i in whichJ){
        wJ[[i]] <- wUs(object=newDAD3, u=i,s=2, remove.s = FALSE)
        ZsimD3_2[i] <- wJ[[i]][1]*z_int[j,i]  + wJ[[i]][2]*rnorm(1, delta[i]*sqrt(newDAD3$n[2,i]))
      }
      names(res) <- namendr
      res[whichJ] <- ZsimD3_2[whichJ] >= critAD(oldD = D1, newD = newDAD3, knew=knew, zint = z_int[j,], decisionsold = dec_old[j,],
                                                theta_rel=theta_rel, decisionsnew = dec_old[j,], oldCrit = oldCritD1, root.interval = c(-10,10))[[2]]
    }
    
    return(res)
  }
  vnames <- c("AD3","namendr","s1res", "s2res", "z_int", "dec_old", "D1", "delta", "Delta",
              "Nsim", "knew", "theta_rel", "alpha", "prev", "piv", "ND1", "oldCritD1")
  cl <- makeCluster(detectCores()-1, type = 'PSOCK')
  clusterExport(cl, as.list(c(ls_start,vnames)), envir = environment())
  clusterEvalQ(cl, library("mvtnorm", "gtools"))
  registerDoParallel(cl)
  endres <- do.call(rbind, clusterApply(cl, 1:Nsim, AD3))
  registerDoSEQ()
  closeAllConnections()
  rm(cl, AD3)
  resultmatrix[4,c(1:5,8)] <- colmeansnoNA(endres)
  #cond:
  ei1cond <- endres[!is.na(endres[,3])&!is.na(endres[,5]),]
  if(is.vector(ei1cond)){
    resultmatrix[4,6] <- as.numeric(endres[3] & endres[5])
    rm(ei1cond)
  }
  else{
    H1starcond <- mean(ei1cond[,1] | ei1cond[,3] & ei1cond[,5])
    resultmatrix[4,6] <- ifelse(is.nan(H1starcond), 0, H1starcond)
    rm(H1starcond, ei1cond)
  }
  ei2cond <- endres[!is.na(endres[,4])&!is.na(endres[,5]),]
  if(is.vector(ei2cond)){
    resultmatrix[4,7] <- as.numeric(endres[4] & endres[5])
    rm(ei2cond)
  }
  else{
    H2starcond <- mean(ei2cond[,2] | ei2cond[,4] & ei2cond[,5])
    resultmatrix[4,7] <- ifelse(is.nan(H2starcond), 0, H2starcond)
    rm(H2starcond, ei2cond)
  }
  #uncond:
  endres[is.na(endres)] <- 0
  resultmatrix_uncond[4,c(1:5,8)] <- colMeans(endres)
  resultmatrix_uncond[4,6:7] <- colMeans(cbind(endres[,1]|endres[,3]&endres[,5], endres[,2]|endres[,4]&endres[,5]))
  rm(endres)
  #########AD4#############
  AD4 <- function(j){
    res <- rep(NA,6)
    names(res) <- namendr
    if(s1res[j,1] && s1res[j,2]){
      return(c(1, 1, rep(NA,3), ND1))
    }
    else if(!s1res[j,1] & !s1res[j,2]){
      return(c(0, 0, rep(NA,3), ND1))
    }
    else{
      res[1:2] <- s1res[j,]
      whichJ <- list()
      for(i in 1:2){
        if(!s1res[j,i]){
          whichJ[[i]] <- grep(as.character(i), names(res)[3:5], value = T)
        }
      }
      whichJ <- unique(unlist(whichJ))
      newDAD4 <- list(prev=prev, n=NULL, test=NULL)
      newDAD4$n <- D1$n
      newDAD4$n[2,whichJ] <- sum(newDAD4$n[2,])*(newDAD4$prev[whichJ]/sum(newDAD4$prev[whichJ]))
      newDAD4$n[2, setdiff(colnames(newDAD4$n),whichJ)] <- 0
      newDAD4$test <- cbind(D1$test, matrix(rep(F,6), nr=2, dimnames = list(c(1:2), c("{1}","{2}","{1,2}"))))
      newDAD4$test <- newDAD4$test[,c("{1}v{1,2}", "{2}v{1,2}", whichJ)]
      newDAD4$test[2, whichJ] <- TRUE
      if(setequal(whichJ, c("{1}","{1,2}"))){
        newDAD4$test[2,1:2] <- c(TRUE,FALSE)
        #expected N:
        #res[6:11] <- c(piv[1]*ND1+ND1,piv[2]*ND1, ifelse(newDAD4$n[2,]==0, NA, colSums(newDAD4$n)), sum(newDAD4$n))
      } else if(setequal(whichJ, c("{2}","{1,2}"))){
        newDAD4$test[2,2:1] <- c(TRUE,FALSE)
        #expected N:ND1
        #res[6:11] <- c(piv[1]*ND1,piv[2]*ND1+ND1, ifelse(newDAD4$n[2,]==0, NA, colSums(newDAD4$n)), sum(newDAD4$n))
      }
      res[6] <- sum(newDAD4$n)
      #generate new stage 2 z-scores in PJ[whichJ]
      ZsimD4_2 <- numeric(length(whichJ)); names(ZsimD4_2) <- whichJ
      pop_ov <- paste0(whichJ, collapse="v")
      w_ov <- wUs(object=newDAD4, u=pop_ov, s=2, remove.s = FALSE)
      ZsimD4_2ov <- w_ov[1]*z_int[j,pop_ov]+w_ov[2]*rnorm(1, Delta[pop_ov]*sqrt(sum(newDAD4$n[2,whichJ]))) 
      wJ <- setNames(vector("list", length(whichJ)), whichJ)
      for(i in whichJ){
        wJ[[i]] <- wUs(object=newDAD4, u=i,s=2, remove.s = FALSE)
        ZsimD4_2[i] <- wJ[[i]][1]*z_int[j,i]  + wJ[[i]][2]*rnorm(1, delta[i]*sqrt(newDAD4$n[2,i]))
      }
      ZsimD4_2 <- c(ZsimD4_2ov, ZsimD4_2); names(ZsimD4_2) <- c(pop_ov, whichJ)
      names(res) <- namendr
      res[c(pop_ov, whichJ)] <- ZsimD4_2[c(pop_ov, whichJ)] >= 
                                critAD(oldD = D1, newD = newDAD4, knew=knew, zint = z_int[j,], decisionsold = dec_old[j,],
                                       theta_rel=theta_rel, decisionsnew = dec_old[j,], oldCrit = oldCritD1, root.interval = c(-10,10))[[2]]
    }
    return(res)
  }
  vnames <- c("AD4","namendr","s1res", "s2res", "z_int", "dec_old", "D1", "delta", "Delta",
              "Nsim", "knew", "theta_rel", "alpha", "prev", "piv", "ND1", "oldCritD1")
  cl <- makeCluster(detectCores()-1, type = 'PSOCK')
  clusterExport(cl, as.list(c(ls_start,vnames)), envir = environment())
  clusterEvalQ(cl, library("mvtnorm", "gtools"))
  registerDoParallel(cl)
  endres <- do.call(rbind, clusterApply(cl, 1:Nsim, AD4))
  registerDoSEQ()
  closeAllConnections()
  rm(cl, AD4, s1res, s2res)
  #cond:
  resultmatrix[5,c(1:5,8)] <- colmeansnoNA(endres)
  ei1cond <- endres[!is.na(endres[,3])&!is.na(endres[,5]),]
  if(is.vector(ei1cond)){
    resultmatrix[5,6] <- as.numeric(endres[3] & endres[5])
    rm(ei1cond)
  }
  else{
    H1starcond <- mean(ei1cond[,1] | ei1cond[,3] & ei1cond[,5])
    resultmatrix[5,6] <- ifelse(is.nan(H1starcond), 0, H1starcond)
    rm(H1starcond, ei1cond)
  }
  ei2cond <- endres[!is.na(endres[,4])&!is.na(endres[,5]),]
  if(is.vector(ei2cond)){
    resultmatrix[5,7] <- as.numeric(endres[4] & endres[5])
    rm(ei2cond)
  }
  else{
    H2starcond <- mean(ei2cond[,2] | ei2cond[,4] & ei2cond[,5])
    resultmatrix[5,7] <- ifelse(is.nan(H2starcond), 0, H2starcond)
    rm(H2starcond, ei2cond)
  }
  #uncond:
  endres[is.na(endres)] <- 0
  resultmatrix_uncond[5,c(1:5,8)] <- colMeans(endres)
  resultmatrix_uncond[5,6:7] <- colMeans(cbind(endres[,1]|endres[,3]&endres[,5], endres[,2]|endres[,4]&endres[,5]))
  rm(endres)
  
  ##################
  #RESULT LIST
  ##################
  rownames(resultmatrix) <- rownames(resultmatrix_uncond) <- c("GSD1","AD1","AD2","AD3", "AD4","single")
  colnames(resultmatrix) <- colnames(resultmatrix_uncond) <- c("H1", "H2", "H{1}", "H{2}", "H{1,2}",
                                                               "H1*", "H2*", "E(N)")
  scenario <- c(prev, delta, Delta, ND1, optparD1) 
  names(scenario) <- c(paste0("pi",names(delta)), paste0("d",names(delta)), paste0("D", names(Delta)), "N1", "WTpar")
  return(list(scenario = scenario, conditional = resultmatrix, unconditional = resultmatrix_uncond))
}

#Quick example for using simAD
#Nsim = 1000 #values > 10^4 might take a while...
#pops = matrix(1, nr=2, nc=2) #two intersecting populations
#delta=.3*rep(1,3) #same effect in each population P{1}, P{2}, P{1,2}
#knew=2 #mid-trial change happens at stage 2
#seed=321  #some seed
#lbound = qnorm(.9) decision boundary for AD1
#alpha=.025 #2.5% type I error rate for initial GSD
#beta=.2 #80 Power for the initial design
#optparD1 = .5 #Pocock design for initial design
#prev = c(.4,.4,.2) #pi{1}=pi{2} = 0.4, p{1,2} = 0.2
#N = 100 #overall sample size of stage 1 (can also be computed by the function if set to NULL)
