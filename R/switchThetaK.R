getThetaFromK   <- function(m,K=m$fit$K,scale='original') UseMethod('getThetaFromK')

getThetaFromK.rcon <- function(m,K=m$fit$K,scale='original'){
  ##cat("getThetaFromK.rcon\n")
  
  S<-m$S
  gc <- m$stdrep
  dimnames(K)<- dimnames(S)
  

  if (is.null(m$vccTerms)){
    vccTerms  <- makeIncMatList(c(gc$VCCU,gc$VCCR),S,Sinv=NULL,type='vcc',full=TRUE)
    eccTerms  <- makeIncMatList(c(gc$ECCU,gc$ECCR),S,Sinv=NULL,type='ecc',full=TRUE)
  } else {
    vccTerms <- m$vccTerms
    eccTerms <- m$eccTerms
  }

  lvcc <- length(vccTerms)
  
  theta <- rep(NA, length(vccTerms)+length(eccTerms))
  
  for (u in 1:length(vccTerms)){
    term.u <- vccTerms[[u]]
    term.u <- term.u[[1]]
    term.u <- rep(term.u,2)[1:2]
    if (scale=='free')
      val <- log(K[term.u[1],term.u[2]]) # lambda
    else
      val <- K[term.u[1],term.u[2]]      # eta
    theta[u] <- val
  }
  if (length(eccTerms)>0){
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]]
      term.u <- term.u[[1]]
      term.u <- rep(term.u,2)[1:2]
      val <- K[term.u[1],term.u[2]] 
      theta[u+lvcc] <- val 
    }
  }
  return(theta)
}

getThetaFromK.rcor <- function(m,K=m$fit$K,scale='original'){
  ##cat("getThetaFromK.rcor\n")
  S<-m$S
  gc <- m$stdrep
  dimnames(K)<- dimnames(S)
  C <- cov2cor(K)


  if (is.null(m$vccTerms)){
    vccTerms  <- makeIncMatList(c(gc$VCCU,gc$VCCR),S,Sinv=NULL,type='vcc',full=TRUE)
    eccTerms  <- makeIncMatList(c(gc$ECCU,gc$ECCR),S,Sinv=NULL,type='ecc',full=TRUE)
  } else {
    vccTerms <- m$vccTerms
    eccTerms <- m$eccTerms
  }

  lvcc <- length(vccTerms)

  theta <- rep(NA, length(vccTerms)+length(eccTerms))
  
  for (u in 1:length(vccTerms)){
    term.u <- vccTerms[[u]]
    term.u <- term.u[[1]]
    term.u <- rep(term.u,2)[1:2]
    if (scale=='free')
      val <- log(sqrt(K[term.u[1],term.u[2]])) # lambda
    else
      val <- sqrt(K[term.u[1],term.u[2]])      # eta
    theta[u] <- val
  }
  if (length(eccTerms)>0){
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]]
      term.u <- term.u[[1]]
      term.u <- rep(term.u,2)[1:2]
      #val <- K[term.u[1],term.u[2]] / sqrt(K[term.u[1],term.u[1]]*K[term.u[2],term.u[2]])
      val <- C[term.u[1],term.u[2]] 
      theta[u+lvcc] <- val 
    }
  }
  return(theta)
}



#####################################################################################

getKFromTheta <- function(m,theta=m$theta,scale='original') UseMethod('getKFromTheta')

getKFromTheta.rcon <- function(m,theta=m$theta,scale='original'){
  ##cat("getKFromTheta.rcon\n")
  S  <- m$S
  gc <- m$stdrep

  if (is.null(m$vccTerms)){
    vccTerms  <- makeIncMatList(c(gc$VCCU,gc$VCCR),S,Sinv=NULL,type='vcc',full=TRUE)
    eccTerms  <- makeIncMatList(c(gc$ECCU,gc$ECCR),S,Sinv=NULL,type='ecc',full=TRUE)
  } else {
    vccTerms <- m$vccTerms
    eccTerms <- m$eccTerms
  }

  lvcc <- length(vccTerms)

  K <- S; K[,] <- 0
  for (u in 1:length(vccTerms)){
    term.u <- vccTerms[[u]]
    for (j in 1:length(term.u)){
      term.uj <- term.u[[j]]
      term.uj <- rep(term.uj,2)[1:2]
      ##if (scale=='free')
      ##  val <- exp(theta[u])
      ##else
      val <- theta[u]
      K[term.uj[1],term.uj[2]] <- val 
    }  
  }

  if (scale=='free')
    diag(K) <- exp(diag(K))
  
  if (length(eccTerms)>0){
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]];       
      val <- theta[u+lvcc]
      for (j in 1:length(term.u)){
        term.uj <- term.u[[j]]
        term.uj <- rep(term.uj,2)[1:2]
        K[term.uj[1],term.uj[2]] <- K[term.uj[2],term.uj[1]] <- val
        ##*(sqrt(K[term.uj[1],term.uj[1]] * K[term.uj[2],term.uj[2]]))        
      }  
    }
  }


  ##   time2<- proc.time()
#   cat("Time1:", (time2-time1)[3], "\n")

#   allTerms <- c(vccTerms, eccTerms)
#   K2 <- S; K2[,]<-0
#   for (u in 1:length(allTerms)){
#     term.u <- allTerms[[u]]
#     K2<-K2+attr(term.u,"info")$incmat2 * theta[u]
#   }
#   diag(K2) <- exp(diag(K2))

#   time3<- proc.time()
# #  cat("Time12", (time3-time2)[3], "\n")

#  print(K2); 
  return(K)
}

getKFromTheta.rcor <- function(m,theta=m$theta,scale='original'){
  ##cat("getKFromTheta.rcor\n")
  
  S<-m$S
  K <- S; K[,] <- 0;   dimnames(K)<- dimnames(S)

  gc <- m$stdrep


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
    for (j in 1:length(term.u)){
      term.uj <- term.u[[j]]
      term.uj <- rep(term.uj,2)[1:2]
      if (scale=='free')
        val <- exp(theta[u]*2)
      else
        val <- theta[u]^2
      K[term.uj[1],term.uj[2]] <- val 
    }  
  }

  if (length(eccTerms)>0){
    for (u in 1:length(eccTerms)){
      term.u <- eccTerms[[u]];       
      val <- theta[u+lvcc]
      for (j in 1:length(term.u)){
        term.uj <- term.u[[j]]
        term.uj <- rep(term.uj,2)[1:2]
        K[term.uj[1],term.uj[2]] <- K[term.uj[2],term.uj[1]] <- val*(sqrt(K[term.uj[1],term.uj[1]] * K[term.uj[2],term.uj[2]]))        
      }  
    }
  }
  return(K)
}













