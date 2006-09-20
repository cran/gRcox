getInf   <- function(m,K=m$fit$K,scale='original') UseMethod('getInf')

getInf.rcon <- function(m,K=m$fit$K,scale='original'){
  ## cat("getInf.rcon\n")
  S<-m$S;   f <- m$n-1;   gc <- m$stdrep
  theta <- getThetaFromK(m,K=K)
  dimnames(K)<- dimnames(S)
  Sigma <- solve(K)
  J <- matrix(0,nrow=length(theta), ncol=length(theta))

  if (is.null(m$vccTerms)){
    vccTerms  <- makeIncMatList(c(gc$VCCU,gc$VCCR),S,Sinv=NULL,type='vcc',full=TRUE)
    eccTerms  <- makeIncMatList(c(gc$ECCU,gc$ECCR),S,Sinv=NULL,type='ecc',full=TRUE)
  } else {
    vccTerms <- m$vccTerms
    eccTerms <- m$eccTerms
  }

  lvcc <- length(vccTerms)


  for (u in 1:length(vccTerms)){
    ##print("V,V")
    term.u <- vccTerms[[u]]
    Ku <- attr(term.u,"info")$incmat2;    ##print(c(term.u))
    for (v in 1:length(vccTerms)){
      term.v <- vccTerms[[v]]
      Kv <- attr(term.v,"info")$incmat2;  ##print(c(term.v))

      ##val <- (f/2)* tr(Ku %*% Sigma %*% Kv %*% Sigma)
      val <- (f/2)* trXYZY(Ku,Sigma, Kv)      
      if (scale=='free')
        J[u,v] <- val * (K[term.u[1],term.u[1]]*K[term.v[1],term.v[1]])        
      else
        J[u,v] <- val 
    }    
  }

  if (length(eccTerms)>0){
    ##print("V,E")    
    for (u in 1:length(vccTerms)){
      term.u <- vccTerms[[u]]
      Ku <- attr(term.u,"info")$incmat2;    ##  print(c(term.u))
      for (v in 1:length(eccTerms)){
        term.v <- eccTerms[[v]]
        Kv <- attr(term.v,"info")$incmat2;  ##  print(c(term.v))
        ##val <- (f/2)*tr(Ku %*% Sigma %*% Kv %*% Sigma)
        val <- (f/2)* trXYZY(Ku, Sigma, Kv)      
        if (scale=='free')
          J[u,v+lvcc] <- J[v+lvcc,u] <- val * K[term.u[1],term.u[1]]
        else
          J[u,v+lvcc] <- J[v+lvcc,u] <- val 
      }    
    }
    ##print("E,E")    
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]]
      Ku <- attr(term.u,"info")$incmat2
      for (v in 1:length(eccTerms)){
        term.v <- eccTerms[[v]]
        Kv <- attr(term.v,"info")$incmat2
        ##J[u+lvcc,v+lvcc] <- (f/2)* tr(Ku %*% Sigma %*% Kv %*% Sigma);
        J[u+lvcc,v+lvcc] <- (f/2)* trXYZY(Ku, Sigma, Kv)      
      }    
    }
  }

  return(J)
}

getInf.rcor <- function(m,K=m$fit$K,scale='original'){
  ## cat("getInf.rcor\n")
  S<-m$S;   f <- m$n-1;   gc <- m$stdrep
  theta <- getThetaFromK(m,K=K)
  C       <- cov2cor(K)
  Cinv    <- solve(C)
  A       <- diag(sqrt(diag(K)))
  Lam     <- A

  Sigma   <- solve(K)
  J <- matrix(0,nrow=length(theta), ncol=length(theta))


  if (is.null(m$vccTerms)){
    vccTerms  <- makeIncMatList(c(gc$VCCU,gc$VCCR),S,Sinv=NULL,type='vcc',full=TRUE)
    eccTerms  <- makeIncMatList(c(gc$ECCU,gc$ECCR),S,Sinv=NULL,type='ecc',full=TRUE)
  } else {
    vccTerms <- m$vccTerms
    eccTerms <- m$eccTerms
  }
            

  lvcc <- length(vccTerms)

  for (u in 1:length(vccTerms)){
    ##print(u)
    term.u <- vccTerms[[u]]
    Ku <- attr(term.u,"info")$incmat2;    ##print(c(term.u))
    for (v in 1:length(vccTerms)){
      term.v <- vccTerms[[v]]
      Kv <- attr(term.v,"info")$incmat2;  ##print(c(term.v))
      ##val <- (2*f)* tr(Ku %*% A %*% S %*% A %*% Kv %*% C)
      val <- (2*f)* trX(Ku %*% Cinv %*% Kv %*% C)      
      if (scale=='free')
        J[u,v] <- val        
      else
        J[u,v] <- val / (A[term.u[1],term.u[1]]*A[term.v[1],term.v[1]])
    }    
  }

  if (length(eccTerms)>0){
    ##print("V,E")    
    for (u in 1:length(vccTerms)){
      term.u <- vccTerms[[u]]
      Ku <- attr(term.u,"info")$incmat2;    ##  print(c(term.u))
      for (v in 1:length(eccTerms)){
        term.v <- eccTerms[[v]]
        Kv <- attr(term.v,"info")$incmat2;  ##  print(c(term.v))
        ##val <- f*tr(Ku %*% A %*% S %*% A %*% Kv)
        val <- f*trX(Ku %*% Cinv %*% Kv)
        if (scale=='free')
          J[u,v+lvcc] <- J[v+lvcc,u] <- val
        else
          J[u,v+lvcc] <- J[v+lvcc,u] <- val /(A[term.u[1],term.u[1]]);
      }    
    }

    ##print("E,E")    
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]]
      Ku <- attr(term.u,"info")$incmat2
      for (v in 1:length(eccTerms)){
        term.v <- eccTerms[[v]]
        Kv <- attr(term.v,"info")$incmat2
        J[u+lvcc,v+lvcc] <- (f/2)* trX(Ku %*% Cinv %*% Kv %*% Cinv);
      }    
    }
  }
  return(J)
}



########################################################################

getScore   <- function(m,K=m$fit$K,scale='original') UseMethod('getScore')


getScore.rcon <- function(m,K=m$K,scale='original'){
  S<-m$S;   f <- m$n-1;   gc <- m$stdrep
  Sigma <- solve(K)
  
  if (is.null(m$vccTerms)){
    vccTerms  <- makeIncMatList(c(gc$VCCU,gc$VCCR),S,Sinv=NULL,type='vcc',full=TRUE)
    eccTerms  <- makeIncMatList(c(gc$ECCU,gc$ECCR),S,Sinv=NULL,type='ecc',full=TRUE)
  } else {
    vccTerms <- m$vccTerms
    eccTerms <- m$eccTerms
  }
  
  lvcc <- length(vccTerms)

  score <- rep(NA, length(vccTerms)+length(eccTerms))
  for (u in 1:length(vccTerms)){
    term.u <- vccTerms[[u]]
    ##print(term.u)
    Ku <- attr(term.u,"info")$incmat2
    ##val <-  (f/2)*(tr(Ku%*%Sigma) - tr(Ku%*%S))
    ##val <-  (f/2)*tr(Ku%*%( Sigma-S))
    val <-  (f/2)*trXY(Ku, (Sigma-S))
    if (scale=='free')
      score[u] <- val * K[term.u[1],term.u[1]]
    else
      score[u] <- val 
  }    

  if (length(eccTerms)>0){
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]]
      Ku <- attr(term.u,"info")$incmat2
      ##val <- (f/2)*(tr(Ku%*%Sigma) - tr(Ku%*%S))
      ##val <-  (f/2)*tr(Ku%*%(Sigma-S))
      val <-  (f/2)*trXY(Ku, (Sigma-S))
      score[u+lvcc] <- val
    }    
  }  

  return(score)
}


getScore.rcor <- function(m,K=m$fit$K, scale='original'){
  S<-m$S;   f <- m$n-1;   gc <- m$stdrep
  theta   <- getThetaFromK(m,K=K) ## Lav om, bruges kun til dim...
  C       <- cov2cor(K); 
  Cinv    <- solve(C)
  A       <- diag(sqrt(diag(K)))
  Lam     <- A
  ## print(A%*%C%*%A);  print(K)
  score <- rep(NA, length(theta))

  if (is.null(m$vccTerms)){
    vccTerms  <- makeIncMatList(c(gc$VCCU,gc$VCCR),S,Sinv=NULL,type='vcc',full=TRUE)
    eccTerms  <- makeIncMatList(c(gc$ECCU,gc$ECCR),S,Sinv=NULL,type='ecc',full=TRUE)
  } else {
    vccTerms <- m$vccTerms
    eccTerms <- m$eccTerms
  }
            

  lvcc <- length(vccTerms)

  for (u in 1:length(vccTerms)){
    term.u <- vccTerms[[u]]
    ##print(term.u)
    Ku <- attr(term.u,"info")$incmat2
    val <-  f*(sum(diag(Ku)) - sum(diag(Ku %*% C %*% A %*% S %*% A)))
    if (scale=='free')
      score[u] <- val
    else
      score[u] <- val / (A[term.u[1],term.u[1]])
  }    

  if (length(eccTerms)>0){
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]]
      Ku <- attr(term.u,"info")$incmat2
      score[u+lvcc] <- f*(trX(Ku%*%Cinv) - trX(Ku %*% A %*% S %*% A))
    }    
  }  
  return(score)
}



