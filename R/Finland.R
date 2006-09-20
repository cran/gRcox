rconFitHYD <- function(m,control=rcox.control()){
  
  if (control$trace>=3) cat("...rconFitHYD\n")
  S<-m$S;  n<-m$n;  K<-m$K; 
  p <- nrow(S);  f <- n-1;   gc <- m$stdrep

  vccTerms  <- makeIncMatList(m$canrep$vcc,S,Sinv=NULL,type='vcc',full=TRUE)
  eccTerms  <- makeIncMatList(m$canrep$ecc,S,Sinv=NULL,type='ecc',full=TRUE)
  allTerms  <- c(vccTerms, eccTerms)

  A <- matrix(0, ncol=length(allTerms), nrow=length(allTerms))
  B <- rep(0,length(allTerms))

  for (u in 1:length(allTerms)){
    term.u <- allTerms[[u]];     ###print("term.u");print(c(term.u)); 
    Ku     <- attr(term.u,"info")$incmat2
    idxu   <- attr(term.u,"info")$useidx
    bu     <- trX(Ku)*1
    B[u]   <- bu
    for (v in 1:length(allTerms)){
      term.v <- allTerms[[v]]      ##print("term.v");print(c(term.v))
      Kv     <- attr(term.v,"info")$incmat2
      idxv   <- attr(term.v,"info")$useidx
      used   <- unique(c(idxu,idxv))

      ##auv <- sum(diag(Ku%*%S%*%Kv))
      ## opt 0
      auv <- trX(Ku%*%S%*%Kv)

      ## opt 1
      ##auv<- trXYZ(Ku,S,Kv)

      ## opt 2
      #auv <- trX(Ku[used,used,drop=FALSE]%*%S[used,used,drop=FALSE]%*%Kv[used,used,drop=FALSE]) 
      ## opt 3
      #auv <- trXYZ(Ku[used,used,drop=FALSE],S[used,used,drop=FALSE],
      #             Kv[used,used,drop=FALSE]) 
      ##cat("auv:", auv, "auv2:",auv2,"diff", auv-auv2,"\n")
      A[u,v] <- auv
    }
  }

  theta <- solve(A,B)
  K2 <- S; K2[,] <- 0; ## print(K2)
 
  for (u in 1:length(vccTerms)){
    term.u <- vccTerms[[u]]
    for (e in term.u){
      K2[e,e] <- theta[u]
    }
  }

  lvcc <- length(vccTerms)
  if (length(eccTerms)>0){
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]]
      theta.cur <- theta[u+lvcc]
      for (e in term.u){
        K2[e[1],e[2]] <- K2[e[2],e[1]] <- theta.cur
      }
    }
  }

  K <- K2

  if (det(K)<0){
    cat("K not positive definite, regularising...\n")
    KKK <- regulariseS(K)
    K <- findKinModelPrimitive(m,KKK,type='rcon')
  }

  ll2 <- ellK(K,S,n)
  return(list(K=K, logL=ll2, vccTerms=vccTerms, eccTerms=eccTerms))
}

rcorFitHYD <- function(s0,control=rcox.control()){

  conv <- 0.001
  if(control$trace>=3) cat("rcorFitHYD\n")
  n             <- s0$n
  S             <- s0$S
  rconModel     <- s0$call

  controlNew        <- s0$control
  controlNew$method <- 'hyd'
  rconModel$method  <- 'hyd'
  rconModel$Kstart  <- diag(diag(S))

  rconModel$type       <- 'rcon'
  rconModel$fit        <- TRUE
  rconModel$control    <- controlNew

  fit <- eval(rconModel)$fit
  rconModel$data    <- NULL
  rconModel$n       <- n

  curr.S    <- s0$S
  C.curr    <- fit$K;

  Lambda  <- diag(sqrt(diag(C.curr)));  
  curr.S  <- Lambda %*% curr.S %*% Lambda
  dimnames(curr.S) <- dimnames(S)
  
  rconModel$S    <- curr.S
  rconModelOut   <- eval(rconModel)
  C.curr         <- rconModelOut$fit$K

  K.curr    <- Lambda %*% C.curr %*% Lambda


  
  if (det(K.curr)<0){
    cat("K not positive definite, regularising...\n")
    KKK <- regulariseS(K.curr)
    K.curr <- findKinModelPrimitive(m,KKK,type='rcon')
  }
  
  logL <- ellK(K.curr, S, n)

  value <- list(K=K.curr, logL=logL)
  return(value)
}
















