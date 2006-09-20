
namesToPairs <- function(x){
  idx <- 1:length(x)
  val <- NULL
  for (i in 1:length(idx)){
    val <- c(val, lapply(idx[-(1:i)],c,i))
  }
  val <- lapply(val,sort)
  value <- lapply(val,function(y){x[y]})
  return(value)
}
