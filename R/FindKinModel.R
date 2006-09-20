###
### Modifies KS such that it is in the model m by averaging elements of KS
### which are restricted to being identical
###

findKinModel <- function(m) UseMethod("findKinModel")

findKinModel.rcon <- function(m){
  rconFitHYD(m)
}

findKinModel.rcor <- function(m){
  n             <- m$n
  S             <- m$S
  rconModel     <- m$call

  rconModel$method     <- 'hyd'
  rconModel$type       <- 'rcon'
  rconModel$fit        <- FALSE
  rconModel$S          <- S
  tmp    <- eval(rconModel)
  C.curr <- tmp$Kstart
  
  val    <- findKinModelPrimitive(m, C.curr, type="rcor")
  
  val<- val
  return(list(K=val, logL=ellK(val, S, n), vccTerms=tmp$vccTerms, eccTerms=tmp$eccTerms))
}


#######################################################################

findKinModelPrimitive <- function(m,KS,type='rcon'){
  if (is.null(KS))
    return(NULL)

  switch(type,
         'rcor'={
           C <- findKinModelPrimitive(m,cov2cor(KS),type='rcon')
           A <- findKinModelPrimitive(m,diag(diag(KS)),type='rcon')
           A <- diag(sqrt(diag(A)))
           KKmod <- A%*%C%*%A
         },
         
         'rcon'={
           gc <- m$stdrep
           KK <- K <- matrix(0, ncol=ncol(KS), nrow=nrow(KS))
           
           VCC <- c(gc$VCCU,gc$VCCR)
           ECC <- c(gc$ECCU,gc$ECCR)
           
           for (i in seq(along=ECC)){
             x<-ECC[[i]]
             id<-matrix(unlist(x),ncol=2,byrow=TRUE)
             KK[id] <- mean(KS[id])
           }
           KK <- KK + t(KK)
           for (i in seq(along=VCC)){
             x<-VCC[[i]]
             diag(KK)[x] <- mean(diag(KS)[x])
           }
           KKdiag <- diag(diag(KK))           
           alpha  <- 0
           repeat{
             KKmod <- KK + alpha*KKdiag
             if (min(eigen(KKmod)$values) > 0)
               break()
             alpha <- alpha+0.1
           }                    
         }
         )
  return(KKmod)
}

