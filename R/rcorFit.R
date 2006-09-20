
rcorFitIterative <- function(m,control=rcox.control()){

  if (control$trace>=3)
    cat("...rcorFitIterative\n")
  
  tstart <- proc.time()
  conv   <- control$logLepsilon
  itmax  <- control$maxit
  rconModel     <- m$call
  method <- rconModel$method
  
  VCC       <- c(m$stdrep$VCCU, m$stdrep$VCCR)
  nam       <- m$varnames
  n         <- m$n
  S         <- m$S
  C.curr    <- diag(1, nrow(S));  #print("C.curr (Start)"); print(C.curr)
  ll        <- d.ll <- NULL
  prev.logL <- iii<- 0
  a         <- rescaleC(S,C.curr,VCC); #print("a"); print(a)
  
  controlNew                <- control
  controlNew$VCCU           <- controlNew$VCCR <- FALSE
  controlNew$representation <- 'stdrep'
  controlNew$method         <- 'user'

  rconModel$type    <- 'rcon'
  rconModel$method  <- 'user'
  rconModel$data    <- NULL
  rconModel$n       <- n
  rconModel$control <- controlNew
  rconModel$fit     <- TRUE
  ##print(rconModel)
  
  repeat{    
    iii <- iii + 1
   
    S.curr           <- diag(a) %*% S %*% diag(a)
    dimnames(S.curr) <- dimnames(S)
    ##print("S.curr"); print(S.curr)
    
    ## Estimate rho
    rconModel$Kstart <- C.curr
    rconModel$S      <- S.curr    
    rconModelnew     <- eval(rconModel)
    C.curr           <- rconModelnew$fit$K;  ##print("C.curr"); print(C.curr)
    ## Calculate K
    K.curr    <- diag(a) %*% C.curr %*% diag(a)

    ## Estimate a
    a         <- rescaleC(S, C.curr, VCC)
    
    ## Monitoring convergence
    curr.logL <- ellK(K.curr, S, n)
    ll        <- c(ll, curr.logL)
    diff.logL <- curr.logL-prev.logL
    prev.logL <- curr.logL
    d.ll      <- c(d.ll, diff.logL)

    if (control$trace>=3)
      cat ("...rcorFit (iteration)! logL:", curr.logL, "\n");

    if ((abs(diff.logL)/abs(prev.logL) < conv) || iii>=itmax){
      break
    }
  }
  dimnames(K.curr) <- dimnames(C.curr)
  usedTime <- (proc.time()-tstart)[3]
  value <- list(K=K.curr, logL=curr.logL, logLvec=ll, time=usedTime)

  if (control$trace>=3)
    cat ("...rcorFit done! logL:", curr.logL, "\n");
  return(value)
}


rescaleC <- function(S,K,NSi,cstart=NULL,itmax=100){
  #print(S); print(K); print(NSi)
  
  d <- nrow(S)
  if (is.null(cstart))
    cstart <- rep(1,d)
  cprev <- cstart
  iii   <- 0

  all <- unique(unlist(NSi))
  all <- all[order(all)]
  repeat{
    iii <- iii + 1
    for (i in 1:length(NSi)){
      cns   <- NSi[[i]]
      compl <- setdiff(all,cns)
      Ai <- sum(S[cns,cns] * K[cns,cns])
      Bi <- sum(S[compl,cns] * K[compl,cns] *cstart[compl])
      Ci <- length(cns)
      Di <- (Bi^2 + 4*Ai*Ci)
      xi <- (-Bi + sqrt(Di))/(2*Ai)
      cstart[cns] <- xi
      ##cat("Ai:", Ai, "Bi:", Bi, "Ci:", Ci, "\n")
    }
    if ((sum( (cprev-cstart)^2) < 1e-9) | iii>=itmax){
      #cat("Iterations:",iii,"\n\n")
      break
    }
    cprev <- cstart
  }
  attr(cstart,"iterations") <- iii
  return(cstart)
}


