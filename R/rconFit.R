
##
## Dette er for iterativ fitting og for user defined fitting...
##

rconFitIterative <- function(m,control=rcox.control()){
  if (control$trace>=3) cat("...rconFitIterative\n")  
  S <- m$S;  K<-m$K;  n<-m$n

  converged <- FALSE
  logLvec   <- prevlogL <- ellK(K,S,n)
  varIndex  <- 1:nrow(S)

  switch(control$representation,
         'mixrep'={    gc <- m$mixrep},
         'stdrep'={    gc <- m$stdrep})  
  
  CS   <- gc$CS
  ECCU <- makeIncMatList(gc$ECCU,S,type='ecc',full=FALSE)
  ECCR <- makeIncMatList(gc$ECCR,S,type='ecc',full=FALSE)
  VCCU <- makeIncMatList(gc$VCCU,S,type='vcc',full=FALSE)
  VCCR <- makeIncMatList(gc$VCCR,S,type='vcc',full=FALSE)
  
  val <- list(K=K)

  while(!converged){
    Kprev <-K

    if (length(CS)>0){
      val <- fitIPSset(CS, K,S,n,varIndex, control)
      K <- val$K;
    }

    if (length(VCCU)>0){
       switch(control$VCCU,
              'IPS'={val <- fitIPSset(VCCU, K,S,n,varIndex,control)},
              'NR' ={val <- fitNR    (VCCU, K,S,n,varIndex,type='vertices',control)}
              )
       K <- val$K;
     }

    if (length(VCCR)>0){
      switch(control$VCCR,
             'NR' ={val <- fitNR    (VCCR, K,S,n,varIndex,type='vertices',control)}
             )
      K <- val$K
    }

    if (length(ECCU)>0){
      switch(control$ECCU,
             'IPS'={val <- fitIPSedge(ECCU, K,S,n,varIndex,type='edges',control)},
             'NR' ={val <- fitNR     (ECCU, K,S,n,varIndex,type='edges',control)}
             )
      K <- val$K;
    }

    if (length(ECCR)>0){
      switch(control$ECCR,
             'NR' ={val <- fitNR (ECCR, K,S,n,varIndex,type='edges',control)}
             )
      K <- val$K;
    }
    
    logL    <- ellK(K,S,n);
    logLvec <- c(logLvec, logL)
    if ((logL-prevlogL) < control$logLepsilon)
      converged <- TRUE
    else {
      prevlogL <- logL; Kprev <- K
    }
    if (control$trace>=3)
      cat("...End of NIPS loop,  logL:", logL,"\n")
  }  
  return(list(K=K, logL=logL, logLvec=logLvec))
}


fitNR <- function(x, K, S, n, varIndex,type,control){
  if (control$trace>=4){
    cat("....Modified Newton:",type, ":", paste(x),"\n");
  }
  for (i in 1:length(x)){
    cl      <- x[[i]]; #print (cl)
    idx     <- sort(unique(unlist(cl)))
    cidx    <- setdiff(1:nrow(S),idx)
    incmat2 <- attr(cl,'info')$incmat2
    trSS2   <- attr(cl,'info')$trSS2
    if (length(cidx)>0)
      subt <- K[idx,cidx,drop=F]%*%cholSolve(K[cidx,cidx,drop=F])%*%K[cidx,idx,drop=F]
    else
      subt <- 0
    val<-modifiedNewton (K, S, n, idx, subt, incmat2, cl, trSS2, type, control)
    K <- val$K
    if (control$trace>=5)
      cat(".....it:", val$NRi, "logL:", val$logL, "Fit:", paste(cl), "\n")    
  }
  return(list(K=K))
}


modifiedNewton <- function(K, S, n, idx, subt, incmat2, cl, trSS2, type, control){
  ts <- proc.time()
  ll <- ellK(K,S,n); iii <- 0; f <- n-1; itmax <- control$maxit
  prev.adj2 <- 0

  if (type=='vertices'){
    aidx <- aidx2 <- matrix(rep(cl,2),ncol=2)
  } else {
    aidx <- matrix(unlist(cl),ncol=2,byrow=TRUE)
    aidx2 <- aidx[,2:1,drop=FALSE]
  }

  ## K <- external(K, subt, incmat2, f, idx, control$IPMepsilon)
  repeat{ 
    Sigma2   <-  cholSolve(K[idx,idx]-subt)
    
    ## Metode 1
    trIS     <-  trXY(incmat2,Sigma2)
    trISIS   <-  trXYXY(incmat2,Sigma2)

    Delta2   <-  trIS - trSS2
    ##adj2     <-  Delta2 /(trISIS + 0.5*f*Delta2^2 )
    adj2     <-  Delta2 /(trISIS + 0.5*Delta2^2 )
    K[aidx] <- K[aidx2] <- K[aidx] + adj2

    iii <- iii+1
    dadj2 <- (adj2 - prev.adj2)
    prev.adj2 <- adj2
    if ( (abs(dadj2)< control$IPMepsilon) | (iii>=itmax) )
      break
  }
  ll2   <- ellK(K,S,n)
  value <- list(K=K,NRi=iii,cl=c(cl),logL=ll2)
  return(value)
}


fitIPSedge <- function(x, K, S, n, varIndex,type,control){
  if (control$trace>=4){
    cat("\nFitting edge with modified IPS:",type," : ");
    cat(paste(x),"\n")
    #lapply(x,function(a)print(c(a)))  
  }
  for (i in 1:length(x)){
    cl      <- x[[i]]; #print (cl)
    idx     <- sort(unique(unlist(cl)))
    cidx    <- setdiff(1:nrow(S),idx)
    if (length(cidx)>0)
      subt <- K[idx,cidx,drop=F]%*%cholSolve(K[cidx,cidx,drop=F])%*%K[cidx,idx,drop=F]
    else
      subt <- 0
    val<-modifiedIPS (K, S, n, idx, subt, cl)
    K <- val$K
  }
  if (control$trace>=4){
    cat("logL:", val$logL,"Fit:",paste(unlist(x)),"\n")
  }
  return(list(K=K))
}

modifiedIPS <- function(K,S,n,idx,subt,cl){

  ts  <- proc.time()
  K11 <- K[idx,idx]  
  KK  <- K11-subt
  a12 <- subt[1,2]
  s   <- S[idx[1],idx[2]]
  sqrtD <- sqrt((2*s*a12+1)^2 + 4*s*(-s*a12^2-a12+s*prod(diag(KK))))
  x1  <- (-(2*s*a12+1) + sqrtD)/(-2*s)

  K[idx[1],idx[2]] <-  K[idx[2],idx[1]] <-x1
  
  ll2 <- ellK(K,S,n)
  return(list(K=K, NRi=0, logL=ll2))
}

fitIPSset <- function(term,K,S,n,varIndex,control){
  if (control$trace>=4){
    cat("\nFitting set with classical IPS:\n");
    lapply(term,function(a)print(c(a)))
  }
  cliques <- term
  my.complement <- function(C) return(setdiff(varIndex,C))
  cliques.complements <- lapply(cliques, my.complement)

  if(length(cliques.complements[[1]])==0){
    return(list(K=solve(S)))
  }

  for(j in 1:length(cliques)){
    C <- cliques[[j]]
    notC <- cliques.complements[[j]]
    K[C,C] <- solve( S[C,C,drop=FALSE] ) +

      K[C,notC,drop=FALSE]%*%solve(K[notC,notC,drop=FALSE])%*%K[notC,C,drop=FALSE]
  }

  logL <- ellK(K,S,n)
  if (control$trace>=4)
    cat("logL:", logL,"\n")
  return(list(K=K))
}



