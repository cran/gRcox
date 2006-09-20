
###
### Pairwise comparison of ECC/VCC based on wald/lr tests
###


compareVCC <- function(cc1, cc2, mf, statistic=c('wald','dev'),
                       criterion=c('test','aic','bic')){
  compareCC(cc1,cc2,mf,type='vcc',statistic=statistic)
}

compareECC <- function(cc1, cc2, mf, statistic=c('wald','dev'),
                       criterion=c('test','aic','bic')){
  compareCC(cc1,cc2,mf,type='ecc',statistic=statistic)
}

compareAllVCC <- function(mf, order=TRUE, statistic=c('wald','dev'),
                       criterion=c('test','aic','bic')){
  compareAllCC(mf,type='vcc',order=order, statistic=statistic)
}

compareAllECC <- function(mf, order=TRUE, statistic=c('wald','dev'),
                       criterion=c('test','aic','bic')){
  compareAllCC(mf,type='ecc',order=order, statistic=statistic)
}


compareAllCC <- function(mf, type='ecc', order=TRUE,statistic=c('wald','dev'),
                       criterion=c('test','aic','bic')){
  
  if (mf$control$trace>=2)
    cat("..compareAllCC\n")

  statistic <- match.arg(statistic)
  criterion <- match.arg(criterion)
  
  mcall <- mf$call
  switch(type,
         'ecc'={   cclist <- mcall$ecc },
         'vcc'={   cclist <- mcall$vcc })
  result <- list()

  if (length(cclist)<=1)
    return()

  for (u in 1:(length(cclist)-1)){
    for (v in (u+1):length(cclist)){
      ccu    <- cclist[[u]]
      ccv    <- cclist[[v]]
      tmp    <- compareCC(ccu, ccv, mf, type=type, statistic=statistic,
                          criterion=criterion)
      result <- c(result,list(tmp)) 
    }
  }

  
  val <- cclist2df(result)

  if (order==TRUE){
    switch(criterion,
           'test'={
             val <- val[order(val$p, decreasing=TRUE),]},
           'aic'={
             val <- val[order(val$aic, decreasing=FALSE),]},
           'bic'={
             val <- val[order(val$bic, decreasing=FALSE),]})
  }
  rownames(val) <- 1:nrow(val)
  return(val)
}


  
compareCC <- function(cc1,cc2,mf,type='ecc',statistic=c('wald','dev'),
                       criterion=c('test','aic','bic')){
  statistic <- match.arg(statistic)
  criterion <- match.arg(criterion)
  
  if (mf$isfit==FALSE){
    mf <- rcoxFit(mf);   ##  cat("Fitting model...\n")
  }

  modelType <- mf$type
  logLmf    <- mf$fit$logL
  ###dfmf      <- getdf(mf)

  
  vn        <- mf$varnames
  n         <- mf$n
  theta     <- mf$fit$theta
  Jinv      <- mf$fit$Jinv
  
  cc1i <- getIndex(splitForm(cc1),vn)[[1]]
  cc2i <- getIndex(splitForm(cc2),vn)[[1]]

  ## Calculate contrast, sdev etc..
  offsetvcc <- ifelse (type=='vcc', 0, length(getCC(mf,type='vcc', form='numeric')))
  CC        <- getCC(mf,type=type, form='numeric')
  iu        <- getParmVecEntry(cc1i, CC)+offsetvcc
  iv        <- getParmVecEntry(cc2i, CC)+offsetvcc
  Diff      <- theta[iu]   - theta[iv]
  vDiff     <- Jinv[iu,iu] + Jinv[iv,iv] - 2*Jinv[iu,iv]

  if (statistic=='dev'){
    mtmp     <- joinCC(list(cc1,cc2), mf, type=type)
    logL     <- mtmp$fit$logL
    statisticstat <- 2*(logLmf - logL)
    statname <- 'DEVstat'
  } else {
    statisticstat  <- Diff^2/vDiff
    statname <- 'Waldstat'
  }
  aic <- statisticstat - 2
  bic <- statisticstat - log(mf$n)
  p.value  <- 1-pchisq(statisticstat, df=1)

  cc1v          <- c(cc1=paste(cc1)[2], cc2=paste(cc2)[2])
  stat          <- c(Diff, sqrt(vDiff), statisticstat, p.value, aic, bic)
  names(stat)   <- c("diff", "se.diff",statname,"p", "aic", "bic")
  result        <- list('cc'= cc1v, stat=stat)
  class(result) <- 'cctest'
  return(result)
}




###
### Test if ECC parameter is zero (i.e.\ deletion of an ECC)
###

zeroECC <- function(ecc, mf, statistic=c('wald','dev'),
                       criterion=c('test','aic','bic')){

  statistic <- match.arg(statistic)
  criterion <- match.arg(criterion)
  if (!mf$isfit)
    mf <- rcoxFit(mf)
  logLmf <- mf$fit$logL

  vn <- mf$varnames
  K <- mf$fit$K

  theta     <- mf$fit$theta
  Jinv      <- mf$fit$Jinv
  
  ccui<-getIndex(splitForm(ecc),vn)[[1]]
  ccui <- rep(ccui,2)[1:2]

  offsetvcc <- length(getCC(mf,type='vcc', form='numeric'))
  CC        <- getCC(mf,type='ecc', form='numeric')
  iu        <- getParmVecEntry(ccui, CC)+offsetvcc
  Diff      <- theta[iu]
  vDiff     <- Jinv[iu,iu] 

  est <- K[ccui[1],ccui[2]]
  if (statistic=='dev'){
    mtmp<-dropECC(ecc, mf)
    logL     <- mtmp$fit$logL
    statisticstat <- 2*(logLmf - logL)
    statname <- 'DEVstat'
  } else {
    statisticstat  <- Diff^2/vDiff
    statname <- 'Waldstat'
  }
  
  aic <- statisticstat - 2
  bic <- statisticstat - log(mf$n)
  p.value  <- 1-pchisq(statisticstat, df=1)
  stat          <- c(Diff, sqrt(vDiff), statisticstat, p.value, aic, bic)
  names(stat)   <- c("diff", "se.diff",statname,"p", "aic", "bic")

  #stat          <- c(est=Diff, se.est=sqrt(vDiff), stat=statisticstat, p=p.value)
  #names(stat)   <- c("est", "se.est",statname,"p")

  result        <- list(cc=list(cc1=paste(ecc)[2]), stat=stat)
  class(result) <- 'cctest'  
  return(result)
}

zeroAllECC <- function(mf, order=TRUE, statistic=c('wald','dev'),
                       criterion=c('test','aic','bic')){
  statistic <- match.arg(statistic)
  criterion <- match.arg(criterion)
  if (!mf$isfit)
    mf <- rcoxFit(mf)

  mcall  <- mf$call
  cclist <- mcall$ecc
  result <- NULL
  for (u in 1:(length(cclist))){ 
    tmp <- zeroECC(cclist[[u]], mf,statistic=statistic)
    result <- c(result,list(tmp)) 
  }

  result <- cclist2df(result)  
  if (order==TRUE){
    result <- result[order(result$p,decreasing=TRUE),]
    rownames(result) <- 1:nrow(result)
  }
  return(result)
}


cclist2df <- function(b){
  ccframe<-data.frame(do.call("rbind", lapply(b,function(x)x$cc)))
  if (ncol(ccframe)>1)
    names(ccframe) <- paste("cc",1:ncol(ccframe),sep='')
  else
    names(ccframe) <- "cc"
  cbind(ccframe, do.call("rbind", lapply(b,function(x)x$stat)))
}

print.cctest <- function(x,...){
  b<-paste(paste(names(x$cc),x$cc, sep=': '),"\n")
  lapply(b,cat)
  print(x$stat)
  return(x)
}




















