getCC <- function(mf, type='ecc', form=c('formula','list','numeric')){
  form <- match.arg(form)
  v <- switch(type,
    'ecc'={mf$call$ecc},
    'vcc'={mf$call$vcc})
  value <- switch(form,
    'formula'={v},
    'list'={lapply(v, splitForm)},
    'numeric'={getIndex(lapply(v, splitForm), mf$varnames)})
  return(value)
}

getIdxOfCC <- function(term,cclist){

  splitcc <- function(term){
                                        #cat("splitcc:\n")
                                        #print(class(term))
    if (class(term)=='formula')
      value <- lapply(strsplit(paste(text=term[[2]]),":"),sort)
    else
      value <- as.list(c("",term))
                                        #print(value)
    return(value)
  }

  a<- splitcc(term)
  x<-lapply(cclist, splitcc)
  idx<-which(unlist(lapply(x,setequal,a)))
  value <-
    if (length(idx)>0)
      idx
    else
      NULL
  return(value)
  }
