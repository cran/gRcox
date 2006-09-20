##
## Forwards/backwards selection among edge colour classes
##

joinECCTest <- function(m, alpha=0.05,statistic=c('wald','dev'),
                       criterion=c('test','aic','bic'), details=TRUE){
  joinCCTest(m, type='ecc',alpha,statistic, criterion,details)
}

joinVCCTest <- function(m, alpha=0.05,statistic=c('wald','dev'),
                       criterion=c('test','aic','bic'), details=TRUE){
  joinCCTest(m, type='vcc',alpha,statistic, criterion, details)
}

splitECCTest <- function(m, alpha=0.05,
                       criterion=c('test','aic','bic'),  details=TRUE){
  splitCCTest(m, type='ecc',alpha, criterion, details)
}

splitVCCTest <- function(m, alpha=0.05,
                       criterion=c('test','aic','bic'),details=TRUE){
  splitCCTest(m, type='vcc',alpha, criterion, details)
}





printIfDetails <- function(val, details){
  if (details==TRUE & !is.null(val)){
    if (nrow(val)>3){
      print(val[1:3,])
      cat("....\n")
    } else {
      print(val)
    }
  }
  
}

joinCCTest <- function(m, type='ecc',alpha=0.05,statistic=c('wald','dev'),
                       criterion=c('test','aic','bic'), details=TRUE){
  mf <- m
  statistic <- match.arg(statistic)
  criterion <- match.arg(criterion)
  
  val <- compareAllCC(mf, type=type, order=TRUE,statistic=statistic)

  printIfDetails(val, details)
  
  if (!is.null(val)){
    xxxx <- switch(criterion,
                   'test'={
                     val <- val[order(val$p, decreasing=TRUE),]
                     val[1,"p"]>alpha},
                   'aic'={
                     val <- val[order(val$aic, decreasing=FALSE),]
                     val[1,"aic"]<0},
                   'bic'={
                     val <- val[order(val$bic, decreasing=FALSE),]
                     val[1,"bic"]<0})
    rownames(val) <- 1:nrow(val)
    
    if (xxxx){
      cat("Joining", type,"\n");
      cat(paste("cc1:", val[1,"cc1"], "\n"))
      cat(paste("cc2:", val[1,"cc2"], "\n"))
      
      ##cat(paste(lapply(val[1,c("cc1","cc2")], paste),
      ##  collapse="\n"),"\n")
      
      aa <- as.list(val[1,-c(1,2)])
      aa2<-as.numeric(aa) 
      names(aa2) <-names(aa)
      print(aa2)
      edge<- list(str2formula(val[1,"cc1"]),str2formula(val[1,"cc2"]))
      mf<-joinCC(edge,mf)
      return(mf)
    }
  }
}


splitCCTest <- function(m, type='ecc',alpha=0.05,
                       criterion=c('test','aic','bic'), details=TRUE){
  mf <- m
  criterion <- match.arg(criterion)
  logLstart <- mf$fit$logL
  dfstart   <- mf$fit$df

  cc1  <- getCC(mf,type=type, form='list')
  cc2  <- getCC(mf,type=type, form='formula')

  cc1<-cc2[lapply(cc1,length)>1]; ## Only composite classes can be split

  foo2 <- function(x){
    tmp<-splitCC(x,mf,type=type)  
    df <- dfstart-tmp$fit$df
    statisticstat <- 2*(tmp$fit$logL-logLstart) 
    val <- data.frame(cc=paste(x)[2], DEVstat=statisticstat,
                      df=df,  p=1-pchisq(statisticstat,df), 
                      aic= statisticstat-2*df, bic=statisticstat-log(mf$n)*df)
    return(val)
  }

  dev2 <- lapply(cc1, foo2)
  dev2 <- do.call("rbind",dev2)

  p    <- dev2$p;  idx  <- which.min(p)

  ## Change sign of AIC, BIC since we add to the model
  dev2$aic <- -dev2$aic
  dev2$bic <- -dev2$bic
  xxxx <- switch(criterion,
                 'test'={
                   dev2 <- dev2[order(dev2$p, decreasing=TRUE),]
                   dev2[1,"p"]<alpha},
                 'aic'={
                   dev2 <- dev2[order(dev2$aic, decreasing=FALSE),]
                   dev2[1,"aic"]<0},
                 'bic'={
                   dev2 <- dev2[order(dev2$bic, decreasing=FALSE),]
                   dev2[1,"bic"]<0})

  rownames(dev2) <- 1:nrow(dev2)

  ### RET DETTE...
  printIfDetails(dev2, details)
  if (xxxx){
    m <- splitCC(cc1[[idx]], mf, type=type)
    cat("Splitting:", type, paste(cc1[[idx]])[2],"\n")
    aa <- as.list(dev2[1,-c(1)])
    aa2<-as.numeric(aa) 
    names(aa2) <-names(aa)
    print(aa2)
    attr(m,"p") <- p[idx]
    return(m)
  }
}


# addECCTest <- function(m, alpha=0.05,
#                        criterion=c('test','aic','bic'), details=TRUE){

#   criterion     <- match.arg(criterion)
#   currentEdges  <- unlist(getCC(m,form='list'),recursive=FALSE)
#   allPossibleEdges <- namesToPairs(m$varnames)
  
#   candidates <- setdiff(allPossibleEdges, currentEdges)

#   if (m$control$trace>=3)
#     cat("...addECCTest - edges to test:", length(candidates),"\n")

#   logLstart <- m$fit$logL
    
#   foo2 <- function(x) {
#     cc=paste(x,collapse=':')
#     if (m$control$trace>=3)
#       cat("...Edge:", cc, "\n")
#     tmp <- addECC(list(x),m)
#     DEVstat <- 2*(tmp$fit$logL-logLstart)
#     p <- 1-pchisq(DEVstat,1)
#     aic <- DEVstat - 2
#     bic <- DEVstat - log(m$n)
#     vvv <- data.frame(cc=cc, DEVstat=DEVstat, p=p, aic=aic, bic=bic)
#     return(vvv)
#   }
#   timexxx <- proc.time()
#   dev2    <- lapply(candidates, foo2)
#   if (m$control$trace>=3)
#     cat("...Time for testing:", (proc.time()-timexxx)[3],"\n")
  
#   val <- do.call("rbind", dev2)

#   ## Change sign of AIC, BIC since we add to the model
#   val$aic <- -val$aic
#   val$bic <- -val$bic

#   xxxx <- switch(criterion,
#                  'test'={
#                    val <- val[order(val$p,   decreasing=FALSE),]
#                    val[1,"p"]<alpha},
#                  'aic'={
#                    val <- val[order(val$aic, decreasing=FALSE),]
#                    val[1,"aic"]<0},
#                  'bic'={
#                    val <- val[order(val$bic, decreasing=FALSE),]
#                    val[1,"bic"]<0})
#   rownames(val) <- 1:nrow(val)

#   printIfDetails(val,details)

#   if (xxxx){
#     edge  <- val[1,'cc']
#     edgef <- as.formula(paste("~",edge))
#     m     <- addECC(edgef, m)
#     cat(paste("Adding edge:\n"))
#     cat(paste(edge, "\n"))

#     aa <- as.list(val[1,-c(1)])
#     aa2<-as.numeric(aa);  names(aa2) <-names(aa);  print(aa2)

#     return(m)
#   }
# }




dropECCTest <- function(m, alpha=0.05,statistic=c('wald','dev'),
                        criterion=c('test','aic','bic'), details=TRUE){
  statistic <- match.arg(statistic)
  criterion <- match.arg(criterion)
  
  val <- zeroAllECC(m,order=TRUE,statistic=statistic)
  printIfDetails(val, details)

  xxxx <- switch(criterion,
                 'test'={
                   val <- val[order(val$p, decreasing=TRUE),]
                   val[1,"p"]>alpha},
                 'aic'={
                   val <- val[order(val$aic, decreasing=FALSE),]
                   val[1,"aic"]<0},
                 'bic'={
                   val <- val[order(val$bic, decreasing=FALSE),]
                   val[1,"bic"]<0})
  rownames(val) <- 1:nrow(val)
  
  if (xxxx){
    ecc<-str2formula(val[1,"cc"])
    cat(paste("Dropping edge; criterion:",criterion,"\n"))
    cat(paste(unlist(val[1,"cc"]),collapse=":"),"\n")
    aa <- as.list(val[1,-c(1)])
    aa2<-as.numeric(aa) 
    names(aa2) <-names(aa)
    print(aa2)
    
    m<-dropECC(ecc,m)
    return(m)
  }
}

    ##attr(m,"droppedECC") <- val[1,"cc"]
    ##attr(m,"p") <- val[1,"p"]




testAddEdges <- function(m, edges=NULL, var=NULL, alpha=0.05,
                   criterion=c('test','aic','bic')){

  criterion        <- match.arg(criterion)
  currentEdges     <- unlist(getCC(m,form='list'),recursive=FALSE)

  if (is.null(edges) & is.null(var)){
    allEdges <- namesToPairs(m$varnames)
  } else {
    if (!is.null(var)){
      allEdges <- namesToPairs(var)
    } else {
      allEdges <- edges
    }
  }
  candidates       <- setdiff(allEdges, currentEdges)
  if (length(candidates)==0)
    return()
  
  logLstart        <- m$fit$logL
    
  testAddSingleEdge <- function(x) {
    tmp    <- addECC(list(x),m)
    DEVstat <- 2*(tmp$fit$logL-logLstart)
    cc     <- paste(x,collapse=':')
    vvv    <- data.frame(cc=cc, DEVstat=DEVstat)
    return(vvv)
  }
  
  dev2 <- lapply(candidates, testAddSingleEdge)
  val <- do.call("rbind", dev2)

  ## Sign of AIC, BIC are changed since we add to the model
  val <- transform(val, p=1-pchisq(DEVstat,1), aic=-(DEVstat-2), bic=-(DEVstat-log(m$n)))

  switch(criterion,
         'test'={
           val <- val[order(val$p,   decreasing=FALSE),]
         },
         'aic'={
           val <- val[order(val$aic, decreasing=FALSE),]
         },
         'bic'={
           val <- val[order(val$bic, decreasing=FALSE),]
         })
  rownames(val) <- 1:nrow(val)
  return(val)
}


addECCTest <- function(m, alpha=0.05,
                       criterion=c('test','aic','bic'), details=TRUE){

  criterion     <- match.arg(criterion)
  val <- testAddEdges(m,criterion=criterion)
  if (is.null(val))
    return()

  printIfDetails(val,details)
  
  xxxx <- switch(criterion,
                 'test'={  val[1,"p"]<alpha},
                 'aic'={   val[1,"aic"]<0},
                 'bic'={   val[1,"bic"]<0})

  if (xxxx){
    edge  <- val[1,'cc']
    edgef <- as.formula(paste("~",edge))
    m     <- addECC(edgef, m)
    cat(paste("Adding edge; criterion:",criterion,"\n"))
    cat(paste(edge, "\n"))

    aa <- as.list(val[1,-c(1)])
    aa2<-as.numeric(aa);  names(aa2) <-names(aa);  print(aa2)

    return(m)
  }
}



