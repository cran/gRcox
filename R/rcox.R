rcox <- function(gm=NULL, vcc=NULL, ecc=NULL, 
                 type=c('rcon','rcor'),
                 method =
                 c("scoring",
                   "ipm",
                   "hyd",
                   "user"                   
                   ),
                 fit=TRUE, data=NULL, S=NULL, n=NULL,
                 Kstart,
                 control=rcox.control()){
  
  mcall   <- match.call()
  type    <- match.arg(type)
  method  <- match.arg(method)

  ## Get (S,n) if not supplied
  if (missing(data) && missing(S)){
    stop("Either data or S matrix must be given...")
  } else {
    if (missing(S)){
      if (control$trace >=1)
        cat(".Calculating S\n")
      S <- cov(data)
      n <- nrow(data)
    }
  }
  
  vn <- colnames(S)
  if (is.null(vn))
    stop("S must have column names\n")

  mcall$S <- S
  mcall$n <- n
  
  ## Find the used names and extract corresponding subset of data
  usedNames <- unique(unlist(c(lapply(vcc,splitForm),lapply(ecc, splitForm),
                               splitForm(gm))))

  usedIdx<-match(usedNames,vn)
  usedIdx<-usedIdx[order(usedIdx)]

  S    <- S[usedIdx,usedIdx]
  Sinv <- try(cholSolve(S))
  if (class(Sinv)=="try-error"){
    cat("S is not invertible")
    Sinv <- NULL
  }
  vn <- colnames(S)

  if (!missing(data))
    data <- data[,usedIdx]
  
  ## Turn variable names into indices and order these

  formToIndex <- function(x,vn){
    if(!is.null(x)){
      z <-listOrder(getIndex(x,vn))
      return( if (length(z)>0) z)
    }}

  gmN  <- splitForm(gm)
  gmI  <- formToIndex(gmN, vn)
  eccN <- lapply(ecc, splitForm)
  eccI <- formToIndex(eccN, vn)
  vccN <- lapply(vcc, splitForm)
  vccI <- formToIndex(vccN, vn)

  canrep <- findCanonicalRepresentation(gmI, eccI, vccI) 
  value  <- findFitRepresentations(canrep$vcc, canrep$ecc)
  
  ## Delete specification in gm
  mcall$gm  <- NULL

  ## Make proper forumlae for vcc, ecc
  mcall$ecc <- getFormula(getSlot(value,'ecc'),vn)
  mcall$vcc <- getFormula(getSlot(value,'vcc'),vn)
  
  ## Handle fitting control according to model type and fitting method.
  control <- setControl(control, method, type)

  value$call    <- mcall
  value$canrep  <- canrep
  value$data    <- data 
  value$S       <- S
  value$Sinv    <- Sinv
  value$control <- control
  value$n       <- n
  value$varnames<- vn
  value$type    <- type
  value$isfit   <- FALSE
  class(value)  <- c('rcox', type)

  ## Create Kstart
  ##
  if (missing(Kstart)){
    timexxx <- proc.time()
    xxx <- findKinModel(value)
    value$Kstart   <- xxx$K
    value$vccTerms <- xxx$vccTerms;
    value$eccTerms <- xxx$eccTerms
    if (is.na(xxx$logL)){
      cat("WARNING: Kstart not pos def; regularising\n")
      KKK <- regulariseS(xxx$K)
      KKK <- findKinModelPrimitive(value,KKK,type)
      print(ellK(KKK,S,n))
      value$Kstart <- KKK
    }
    value$initlogL <- xxx$logL
    if (control$trace>=3)      
      cat("...Finding Kstart in model - time:", (proc.time()-timexxx)[3],"\n");
  } else {
    value$Kstart <- Kstart
    if (control$trace>=3)      
      cat("...Kstart is specified\n");
  }
  
  if (fit==TRUE){  
    value <- rcoxFit(value)
  }
  return(value)
}


setControl <- function(control,method,type){
  switch(type,
         'rcon'={
           switch(method,
                  "scoring"={
                    control$representation='stdrep';
                    control$VCCU='NR'; control$VCCR='NR';
                    control$ECCU='NR'; control$ECCR='NR'},
                  "ipm"={ # Was "mixed"
                    control$representation='mixrep';
                    control$VCCU='IPS'; control$VCCR='NR';
                    control$ECCU='IPS'; control$ECCR='NR'}
                )           
         },
         'rcor'={
           switch(method,
                  "scoring"={
                    control$representation='stdrep';
                    control$VCCU=NULL;  control$VCCR=NULL;
                    control$ECCU='NR'; control$ECCR='NR'},                  
                  "ipm"={ # Was "mixed"
                    control$representation='stdrep';
                    control$VCCU=NULL; control$VCCR=NULL;
                    control$ECCU='NR'; control$ECCR='NR'}
                  )
         } )

  control$method <- method  ## Not a nice way, but...
  if (control$trace>=4){
    cat("....setControl type:", type, "method:", method,"\n")
  } 
  return(control)
}


