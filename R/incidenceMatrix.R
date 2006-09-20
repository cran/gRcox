###
### Creating incidence matrices for ECC/VCC
###

makeIncMatList <- function(x,S=NULL,Sinv=NULL,type,full=FALSE){
  lapply(x, function(a){
    v<-makeIncMat(a,S,Sinv,type,full);
    attr(a,'info')<-v;
    return(a)})
}

makeIncMat <- function(x,S=NULL,Sinv=NULL,type,full=FALSE){
  switch(type,
         'ecc'={makeIncMatECC(x,S,Sinv,full)},
         'vcc'={makeIncMatVCC(x,S,Sinv,full)})
}


###
### Functions below here are for *internal use only*
###

makeIncMatVCC <- function(x,S=NULL,Sinv=NULL,full=FALSE){

  if (is.null(x)) return(NULL)
  ccmatrix <- cbind(v1=x,v2=x)
  incmat <- matrix(0,nrow=nrow(S),ncol=nrow(S))
  for (e in x){
    incmat[e[1],e[1]] <- 1
  }
 
  idx     <- sort(unique(unlist(x)))
  incmat2 <- diag(1,length(idx))

  storage.mode(incmat) <- storage.mode(incmat2) <- 'double'

  value         <- list(ccmatrix=ccmatrix,useidx=idx)
  value$incmat2 <-  if(full) incmat else incmat2
  
  if (!is.null(S)){
    value$trSS2 <- trXY(incmat2,S[idx,idx,drop=FALSE])
  }
  
  if (!is.null(Sinv)){
    value$trSS2inv       <- trXY(incmat2,Sinv[idx,idx,drop=FALSE])
  }
  return(value)
}

makeIncMatECC <- function(x,S=NULL,Sinv=NULL,full=FALSE){

  if (is.null(x)) return(NULL)
  ccmatrix <- do.call("rbind",x)
  incmat <- matrix(0,nrow=nrow(S),ncol=nrow(S))
  for (e in x){
    incmat[e[1],e[2]] <- incmat[e[2],e[1]] <- 1
  }
  
  idx    <- sort(unique(unlist(x)))  
  incmat2 <- incmat[idx,idx]
  S2      <- S[idx,idx]
  storage.mode(incmat) <- storage.mode(incmat2) <- 'double'
  value <- list(ccmatrix=ccmatrix,useidx=idx)
  value$incmat2 <-  if(full) incmat else incmat2

  if (!is.null(S)){
    value$trSS2 <- trXY(incmat2,S[idx,idx,drop=FALSE])
  }

  if (!is.null(Sinv)){
    trSS2inv <- trXY(incmat2,Sinv[idx,idx])
    value$trSS2inv <- trSS2inv
  }
  return(value)
}


