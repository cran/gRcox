display <- function(x) UseMethod("display")
display.rcox <- function(x){

  m2 <- x
  eccList <- getCC(m2,'ecc','list')
  vccList <- getCC(m2,'vcc','list')
  
  vccColors<-heat.colors(length(vccList))
  V <- unlist(vccList)
  V <- V[order(V)]
  vertexColors <- NULL
  for (i in 1:length(vccList)){
    tmp <- vccList[[i]]
    #print(tmp)
    if (length(tmp)==1){
       vcolor <- "white"    
    } else {
       vcolor <- vccColors[i] 
    }
    d <- c(vstr = rep(vcolor, length(tmp)))
    names(d) <- tmp
    vertexColors <- c(vertexColors, d)    
   }

   nAttrs <- list()
   nAttrs$fillcolor <- vertexColors

  G <- new("graphNEL", nodes=V)
  plot(G, nodeAttrs=nAttrs)
  
  eccColors<-topo.colors(length(eccList))

  edgeColors <- NULL
  if (length(eccList)>0){    
    for (i in 1:length(eccList)){
      tmp <- eccList[[i]]; ltmp <- length(tmp)
      for (j in 1:ltmp){
        ee <- tmp[[j]]
        ee <- ee[order(ee)]
                                        #print(ee)
        G <- addEdge(ee[1], ee[2], G, weight=1)
        estr <- paste(ee[1],"~",ee[2],sep='')
        if (ltmp > 1){
          ecolor <- eccColors[i]
          d <- c(estr = ecolor)
          names(d) <- estr
          edgeColors <- c(edgeColors, d)                  
      }
      }
    }
  }
    
  if (!is.null(edgeColors))
    eAttrs <- list(color=edgeColors)
  else
    eAttrs <- list()

  plot(G, "neato", nodeAttrs = nAttrs, edgeAttrs = eAttrs)


}
