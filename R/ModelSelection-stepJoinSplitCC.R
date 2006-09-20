
###
###  Stepwise join of ECC/VCC
###

showIfDisplay <- function(mf, display){
  if (display)
    display(mf)
}

stepJoinVCC <-function(m ,alpha=0.05, details=TRUE,
                       criterion=c('test','aic','bic'),
                       statistic=c('wald','dev'), display=FALSE){
  stepJoinCC(m ,type='vcc',
             alpha=alpha, details=details, statistic=statistic,
             criterion=criterion,display=display)
}

stepJoinECC <-function(m ,alpha=0.05,details=TRUE,
                       criterion=c('test','aic','bic'),
                       statistic=c('wald','dev'), display=FALSE){
  stepJoinCC(m ,type='ecc', alpha=alpha, details=details,
             statistic=statistic, criterion=criterion,display=display )
}


stepSplitVCC <-function(m ,alpha=0.05, details=TRUE,
                       criterion=c('test','aic','bic'),
                       display=FALSE){
  stepSplitCC(m ,type='vcc',  alpha=alpha,
              details=details, criterion=criterion,display=display)
}


stepSplitECC <-function(m ,alpha=0.05, details=TRUE,
                       criterion=c('test','aic','bic'),
                       display=FALSE){
  stepSplitCC(m ,type='ecc',
              alpha=alpha, details=details, 
              criterion=criterion,display=display)
}




stepJoinCC <-function(m ,type='ecc', alpha=0.05,details=TRUE,
                      criterion=c('test','aic','bic'),
                      statistic=c('wald','dev'), display=FALSE){
  mf <- m
  if (mf$control$trace>=2)
    cat("..stepJoinCC\n")
  
  statistic <- match.arg(statistic)
  criterion <- match.arg(criterion)
  
  if (!mf$isfit)
    mf <- rcoxFit(mf)

  showIfDisplay(mf, display)
  
  tr <- compareAllCC(mf,type=type,order=TRUE,statistic=statistic, criterion=criterion)

  xxxx <- switch(criterion,
                 'test'={
                   tr <- tr[order(tr$p, decreasing=TRUE),]
                   tr[1,"p"]>alpha},
                 'aic'={
                   tr <- tr[order(tr$aic, decreasing=FALSE),]
                   tr[1,"aic"]<0},
                 'bic'={
                   tr <- tr[order(tr$bic, decreasing=FALSE),]
                   tr[1,"bic"]<0})
  rownames(tr) <- 1:nrow(tr)
  if (!xxxx)
    return()

  repeat{
    if (xxxx){
      idx<-tr[1,c('cc1','cc2')]
      printIfDetails(tr, details)
      
      cat("Joining", type,"\n");
      cat(paste("cc1:", idx[[1]],"\n"))
      cat(paste("cc2:", idx[[2]],"\n"))
      
      aa <- as.list(tr[1,-c(1,2)])
      aa2<-as.numeric(aa);  names(aa2) <-names(aa); print(aa2)
      
      idx <- lapply(unclass(idx), function(a) paste("~",a))
      idx <- lapply(idx, as.formula)
      
      m  <- joinCC(idx,mf,type=type,fit=FALSE)
      mf <- rcoxFit(m)
      showIfDisplay(mf, display)
    
      tr <- compareAllCC(mf,type=type,order=TRUE,statistic=statistic, criterion=criterion)

      if(is.null(tr))
        break()
      xxxx <- switch(criterion,
                     'test'={
                       tr <- tr[order(tr$p, decreasing=TRUE),]
                       tr[1,"p"]>alpha},
                     'aic'={
                       tr <- tr[order(tr$aic, decreasing=FALSE),]
                       tr[1,"aic"]<0},
                     'bic'={
                       tr <- tr[order(tr$bic, decreasing=FALSE),]
                       tr[1,"bic"]<0})
      rownames(tr) <- 1:nrow(tr)

    } else{
      break()
    }
  }
  return(mf)
}



testSplitAllCC <- function(bf3, type='ecc'){

  gettest <- function(mf){
    if (!is.null(mf$fit))
      c(logL=mf$fit$logL, df=mf$fit$df)
  }
  
  cclist2 <- getCC(bf3, type=type, form='list')
  cclist1 <- getCC(bf3, type=type)
  composite<-cclist1[lapply(cclist2,length) > 1]
  
  result <- NULL
  
  for (u in 1:length(composite)){
    ccu <- composite[[u]]
    tmp <- splitCC(ccu, bf3, type=type)
    t1<-gettest(bf3)
    t2<-gettest(tmp)
    Q <- 2*(t2['logL']-t1['logL'])
    names(Q)<-NULL
    df <- t1['df']-t2['df']
    aic <- Q - 2*df
    bic <- Q - log(bf3$n)
    names(df)<-NULL
    p.value <- 1-pchisq(Q,df)
    val <- data.frame(u=u,cc=paste(ccu)[2],DEVstat=Q, df=df, p=p.value, aic=aic, bic=bic)
    result <- rbind(result,val)
  }
  result <- result[order(result$p,decreasing=TRUE),]
  result[,"u"]<-NULL
  return(result)
}

stepSplitCC <-function(m ,type='ecc', alpha=0.05, details=TRUE,
                       criterion=c('test','aic','bic'),
                       display=FALSE){

  mf <- m
  criterion <- match.arg(criterion)
  if (!mf$isfit)
    mf <- rcoxFit(mf)

  showIfDisplay(mf, display)
  
  cclist2 <- getCC(mf, type=type, form='list')
  cclist1 <- getCC(mf, type=type)
  composite<-cclist1[lapply(cclist2,length) > 1]
  res <- testSplitAllCC(mf,type)

  res$aic <- -res$aic
  res$bic <- -res$bic

  xxxx <- switch(criterion,
                 'test'={
                   res <- res[order(res$p, decreasing=FALSE),]
                   res[1,"p"]<alpha},
                 'aic'={
                   res <- res[order(res$aic, decreasing=FALSE),]
                   res[1,"aic"]<0},
                 'bic'={
                   res <- res[order(res$bic, decreasing=FALSE),]
                   res[1,"bic"]<0})
  rownames(res) <- 1:nrow(res)

  if (!xxxx)
    return()
  
  
  repeat{
    if(xxxx){
      
      printIfDetails(res, details)

      cat("Splitting", type, ":", paste(res[1,"cc"]), "\n")
      aa <- as.list(res[1,-c(1)])
      aa2<-as.numeric(aa); names(aa2) <-names(aa);  print(aa2)

      mf <- splitCC(composite[[res[1,'cc']]],
                    mf, type=type)
      
      showIfDisplay(mf, display)
      
      cclist2 <- getCC(mf, type=type, form='list')
      cclist1 <- getCC(mf, type=type)
      composite<-cclist1[lapply(cclist2,length) > 1]

      if (length(composite)>0){
        res <- testSplitAllCC(mf,type)
        
        res$aic <- -res$aic
        res$bic <- -res$bic

        xxxx <- switch(criterion,
                       'test'={
                         res <- res[order(res$p, decreasing=FALSE),]
                         res[1,"p"]<alpha},
                       'aic'={
                         res <- res[order(res$aic, decreasing=FALSE),]
                         res[1,"aic"]<0},
                       'bic'={
                         res <- res[order(res$bic, decreasing=FALSE),]
                         res[1,"bic"]<0})
        rownames(res) <- 1:nrow(res)
      } else {
        break()
      }
    } else {
      break()
    }
  }
  return(mf)
}


stepAddECC <- function(m, alpha=0.05,
                       criterion=c('test','aic','bic'),details=TRUE,display=FALSE){
  repeat{
    v <- addECCTest(m,alpha,criterion=criterion,
                     details=details)
    if (!is.null(v)){
      m <- v
      showIfDisplay(m,display)
    }
    else
      break()
  }
  return(m)
}


stepDropECC <- function(m, alpha=0.05, statistic=c('wald','dev'),
                       criterion=c('test','aic','bic'),details=TRUE, display=FALSE){
  statistic <- match.arg(statistic)
  repeat{
    v <- dropECCTest(m,alpha,statistic=statistic,criterion=criterion,
                     details=details)
    if (!is.null(v)){
      m <- v
      showIfDisplay(m,display)
    }
    else
      break()
  }
  return(m)
}



