joinVCC <- function(cc,m,fit=TRUE){
  joinCC(cc,m,type='vcc',fit)
}

joinECC <- function(cc,m,fit=TRUE){
  joinCC(cc,m,type='ecc',fit)
}

joinCC <- function(cc,m,type='ecc',fit=TRUE){
  .primitiveJoinSplitCC(cc, m, type, fit, operation='join')
}

splitVCC <- function(cc,m,fit=TRUE){
  splitCC(cc,m,type='vcc',fit)
}

splitECC <- function(cc,m,fit=TRUE){
  splitCC(cc,m,type='ecc',fit)
}

splitCC <- function(cc,m,type='ecc',fit=TRUE){
  .primitiveJoinSplitCC(cc, m, type, fit, operation='split')
}

addECC <- function(cc,m,fit=TRUE){
  .primitiveAddDropECC(cc,m,fit,operation='add')
}

dropECC <- function(cc,m,fit=TRUE){
  .primitiveAddDropECC(cc,m,fit,operation='drop')
}



########################################################################

.primitiveJoinSplitCC <- function(cc, m, type, fit=TRUE, operation){
  
  mcall <- m$call  
  switch(type,
         'ecc'={   cclist <- mcall$ecc },
         'vcc'={   cclist <- mcall$vcc })

  switch(operation,
         'join'={
           cc <- lapply(cc, joinForm)
           ccnum<- unlist(lapply(cc, getIdxOfCC, cclist))                      
           tmp <- lapply(cclist[ccnum],paste)
           tmp <- lapply(tmp,function(x)x[2])
           tmp <- as.formula(paste("~",
                                   paste(unlist(lapply(tmp,
                                                       strsplit, "\\+")),collapse="+"),
                                   sep=''))
           cclist  <- c(list(tmp),cclist[-ccnum])           
         },
         'split'={
           cc <- joinForm(cc)
           cc<- getIdxOfCC(cc, cclist)           
           newi<-unlist(strsplit(paste(cclist[[cc]])[2],".\\+."))
           newi<-lapply(paste("~",newi),as.formula)
           cclist <- c(cclist[-cc],newi)
         })

  switch(type,
         'ecc'={   mcall$ecc <- cclist },
         'vcc'={   mcall$vcc <- cclist })
  mcall$fit <- fit;  mcall <- eval(mcall);  return(mcall)
}

.primitiveAddDropECC <- function(cc,m,fit=TRUE,operation){
  mcall <- m$call
  cclist <- mcall$ecc

  switch(operation,
         'add'={
           cc <- joinForm(cc)
           mcall$ecc <- c(cc,cclist)
         },
         'drop'={
           cc<- getIdxOfCC(cc, cclist)
           mcall$ecc <- cclist[-cc]           
         })
  mcall$fit <- fit;  mcall <- eval(mcall);  return(mcall)
}










