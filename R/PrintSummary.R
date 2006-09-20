
shprint <- function(l){
  lapply(l,function(x) cat(paste(" ",x, collapse='.'),"\n"))
  invisible(return(l))
}

print.rcox <- function(x,...){
  m <- x$stdrep
  vcc<-c(m$VCCU,m$VCCR)  
  ecc<-c(m$ECCU,m$ECCR)
  
  cat("\nModel type:", toupper(x$type), "\n")  
  cat("vcc:", deparse(x$call$vcc),"\n")
  cat("ecc:", deparse(x$call$ecc),"\n")

  if (!is.null(x$fit)){
    cat("logL:", x$fit$logL, "df:", x$fit$df, "dimension:", x$fit$dim,
        "time taken:", x$fit$timetaken, "\n")
  }
}

summary.rcox <- function(object,...){
  x <- object
  print(x)
  
  if (!is.null(x$fit)){
    cat("\n")
    cat("C\\K (partial correlations \\ concentrations):\n")
    print(round(x$fit$KC,8))
  }
}
