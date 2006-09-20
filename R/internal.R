
getAtomCC <- function(x){  x[(lapply(x,length)==1)]  }
getCompCC <- function(x){  x[!(lapply(x,length)==1)] }

getParmVecEntry <- function(ccui, CC){ ## Sæt et andet sted....
  which(unlist(lapply(lapply(CC, is.element, list(ccui)),any)))
}

str2formula <- function(str){
  as.formula(paste("~",as.character(str)))
}

cliqueToEdges <- function (x) {
  if (length(x) == 1)
    return(x)
  else {
    value <- NULL
    lnx <- length(x)
    for (i in 1:(lnx-1)){
      for (j in (i+1):lnx){
        value <- c(value,list(c(x[i],x[j])))
      }
    }
    return(value)
  }
}

ellK <- function(K, S, n){
  ##print("det(K)"); print(det(K))
  ##print(trXY(K,S)); print(sum(diag(K%*%S)))
  value <- n/2*(log(det(K)) - trXY(K,S))
  ##value <- n/2*(log(det(K)) - nrow(S))
  return(value)
}

cholSolve <- function(ma)
  chol2inv(chol(  ma  ))

listOrder         <- function(x) UseMethod('listOrder')
listOrder.numeric <- function(x){ x[order(x)] }
listOrder.list    <- function(x){ lapply(x,function(v)listOrder(v)) }
listOrder.default <- function(x){ x }

getIndex           <- function(x,vn)UseMethod('getIndex')
getIndex.list      <- function(x,vn){ lapply(x, getIndex, vn) }
getIndex.character <- function(x,vn){ match(x,vn) }
getIndex.integer   <- function(x,vn){ x }
getIndex.numeric   <- function(x,vn){ x }

getNames         <- function(x,vn) UseMethod('getNames')
getNames.numeric <- function(x,vn){ vn[x] }
getNames.list    <- function(x,vn){ lapply(x, getNames, vn) }


getFormula <- function(x,vn){
  if (is.null(x))
    return(NULL)
  ecc<-x
  a<-getNames(ecc,vn)

  value<-lapply(a,function(b){
    v<- lapply(b, function(d){paste(d,collapse=':')})
    v <- paste(v,collapse='+')
    as.formula(paste("~",v))
  })
  return(value)
}


splitForm <- function(x){
  val <- switch(class(x),
                'formula'={
                  mf<-paste(x)[2]
                  val<-strsplit(unlist(strsplit(mf," \\+ ")),":")
                  val
                },
                'list'={
                  x }
                )
  return(val)
}


joinForm <- function(x){
  val <-switch(class(x),
         'list'={
           val <- as.formula(paste("~",paste(unlist(lapply(x, paste,collapse=':')),
                                             collapse='+')))
           val
         },
         'formula'={
           x }
         )
  return(val)
}


getSlot<-function(m,s){
  s <- tolower(s);
  switch(s,
    'vcc'={c(m$stdrep$VCCU,m$stdrep$VCCR)},
    'ecc'={c(m$stdrep$ECCU,m$stdrep$ECCR)}
    )
}


getdf <- function(m1){
  nr <- nrow(m1$Kstart)
  np<-length(getCC(m1,type='ecc'))+length(getCC(m1,type='vcc'))
  df <- nr*(nr+1)/2 - np
  return(df)
}








# mparsemodel <- function(m2){
#   a<-mparse(NULL, c(m2$stdrep$ECCU,m2$stdrep$ECCR), c(m2$stdrep$VCCU,m2$stdrep$VCCR))
#   m2$stdrep <- a$stdrep
#   m2$mixrep <- a$mixrep
#   m2$fit <- NULL
#   return(m2)
# }



# getdf <- function(K){
#   Kup <- K[upper.tri(K,diag=TRUE)]
#   d   <- unique(as.numeric(Kup))
#   np <- length(d[d!=0])
#   nr <- nrow(K)
#   df <- nr*(nr+1)/2 - np
#   return(df)
# }



# ##
# ## Get variance of K[i,j]-K[u,v] - RCON
# ##
# varDiffConc <- function(ij,uv,K,n){
#  getCovConc(ij,ij,K,n)+getCovConc(uv,uv,K,n)-2*getCovConc(ij,uv,K,n)
# }

# ##
# ## Get variance of C[i,j]-C[u,v] - RCOR
# ##
# varDiffCorr <- function(ij,uv,K,n){
#  getCovCorr(ij,ij,K,n)+getCovCorr(uv,uv,K,n)-2*getCovCorr(ij,uv,K,n)
# }

# ##
# ## Get cov(K[i,j], K[u,v])
# ##
# getCovConc <- function(ij,uv,K,n){
#   p <- nrow(K)
#   if (length(ij)==1) ij <- rep(ij,2)
#   if (length(uv)==1) uv <- rep(uv,2)
#   i<-ij[1]; j<-ij[2]; u<-uv[1]; v<-uv[2]
#   value <- (K[i,u]*K[j,v]+K[i,v]*K[j,u] +
#             (2*K[i,j]*K[u,v])/(n-p-1))/((n-p)*(n-p-1)*(n-p-3))
#   return(value)
# }




# regularizeS <- function(SSS){
#   ev<-length(which(eigen(SSS)$value>1e-1))
#   if (ev == nrow(SSS))
#     return(SSS)

#   ## Else do bisection
#   a.left <- 0; a.right <- 1
#   dSSS <- diag(diag(SSS))
#   alpha <- alpha.prev <- 1
  
#   repeat {
# #    print(c(a.left,a.right))
#     a.cand <- (a.left+a.right)/2
#     SSScand<-SSS*(1-a.cand)+ dSSS*a.cand
#     ev<-length(which(eigen(SSScand)$value>1e-1))
#     if (ev == nrow(SSS)){
#       a.left <- a.cand
#     } else {
#       a.right <- a.cand
#     }
#     if (a.right - a.left < 1e-4)
#       break()
#   }
#   #print(a.cand)
#   SSScand<-SSS*(1-a.cand)+ dSSS*a.cand
#   return(SSScand)
# }


regulariseS <- function(SSS){
  eig<- eigen(SSS)
  ev <- eig$values>0
  SSScand <- eig$vectors[,ev] %*% diag(eig$values[ev]) %*% t(eig$vectors[,ev])
  return(SSScand)
}
