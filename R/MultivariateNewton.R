
multivariateNewton <- function(m, K = m$Kstart, control=m$control){

  stephalf <- function(disc1, theta){
    ##cat("stephalf: logL (start):", logL, "prevlogL:", prevlogL, "\n")
    itcount  <- 1
    stepsize <- 1
    repeat{ 
      
      thetaNew  <- theta+disc1
      Knew      <- getKFromTheta(m, thetaNew, scale='free')
      logL      <- ellK(Knew, m$S, m$n)
      difflogL  <- logL - prevlogL      
      if (control$trace>=5) 
        if (stepsize<1)
          cat(".....iteration: stepsize:",stepsize,
              "logL:", logL, "prevlogL:", prevlogL,
              "difflogL:", difflogL, "\n")      
        else
          cat(".....iteration:",
              "logL:", logL, "prevlogL:", prevlogL,
              "difflogL:", difflogL, "\n")      
      
      if (difflogL < -1e-6){
        stepsize <- stepsize / 2
        if (control$trace>=5) 
          cat(".....Halving stepsize; stepsize:", stepsize, "\n")
        disc1   <- disc1 * stepsize
        itcount <- itcount + 1
      }
      else
        return(list(theta=thetaNew,logL=logL,K=Knew, difflogL=difflogL))
    }
  }

  S  <- m$S
  gc <- m$stdrep
  
  if (is.null(m$vccTerms)){
    vccTerms  <- makeIncMatList(c(gc$VCCU,gc$VCCR),S,Sinv=NULL,type='vcc',full=TRUE)
    eccTerms  <- makeIncMatList(c(gc$ECCU,gc$ECCR),S,Sinv=NULL,type='ecc',full=TRUE)
  } else {
    vccTerms <- m$vccTerms
    eccTerms <- m$eccTerms
  }

  logL  <- prevlogL <- ellK(K,m$S,m$n)
  if (control$trace>=3)
    cat ("...MultivariateNewton, logL (start):", logL, "\n")
  f <- m$n-1
  
  theta <- getThetaFromK(m,K=K,scale='free')
      
  tstart <- proc.time()
  if (control$maxit>0){
    for(i in 1:control$maxit){    
      S         <- getScore(m,K=K,scale='free')
      J         <- getInf  (m,K=K,scale='free')

      ##disc1     <- qr.solve(J+S%*%t(S)/f, S)
      if (i==1)
        disc1     <- try(qr.solve(J+S%*%t(S)/(f), S))
      else
        disc1     <- try(qr.solve(J+S%*%t(S)/f, S))
      if (class(disc1)=='try-error'){
        print(solve(J))
        print(round(J,2));
        print(eigen(J)$values)
        print(round(S,2))
        print(eigen(J+S%*%t(S))$values)
      }
      tmp       <- stephalf(disc1, theta)
      thetaNew  <- tmp$theta
      logL      <- tmp$logL
      K         <- tmp$K
      difflogL  <- tmp$difflogL
      if (control$trace>=4)
        cat("....Iteration:",i,"logL:",logL,  "difflogL:", difflogL,"\n")
      
      if( difflogL<control$logLepsilon ){
        break
      } else {
        prevlogL <- logL; theta <-thetaNew
      }
    }      
  }
  logL <- ellK(K,m$S,m$n)
  return(list(K=K, logL=logL))
}

#   if (is.na(logL)){
#     K <- regularizeS(K)
#     logL  <- prevlogL <- ellK(K,m$S,m$n)
#     if (control$trace>=3) cat ("...MultivariateNewton, logL (start):", logL, "\n")
#     if (control$trace>=3) cat ("...MultivariateNewton, logL (start):", logL, "\n")
#   }
