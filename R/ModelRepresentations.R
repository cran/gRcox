## Turn cliques into edges. Note that a singleton variable i is represented as (i NA).
## Organize edges in data frame

##lprint <- function(x)print(paste(x))

findCanonicalRepresentation<- function(gm=NULL, ecc=NULL, vcc=NULL){

  isMember <- function(x,y){ # Is x (a vector) a member of y (a list of vectors)
    unlist(lapply(y,function(yy){all(yy==x)}))
  }

  varnames <- unique(unlist(c(gm,ecc,vcc)))

  gmEdges  <<- unlist(lapply(gm, cliqueToEdges),recursive=FALSE)
  ##print("gmEdges");  lprint(gmEdges);
  gmVertices <<- unique(unlist(gm))
  ##print("gmVertices"); print(gmVertices)

  ## Remove edges whic are really only vertices
  gmEdges <- gmEdges[lapply(gmEdges,length)>1]
  
  ## Edges in gm which are also in ecc must be identified in gm
  eccEdges <<- unlist(ecc,recursive=FALSE)
  ##print("eccEdges");  lprint(eccEdges)
  gmEdgesIneccEdges <- unlist(lapply(lapply(gmEdges, isMember, eccEdges),any))
  ##print("gmEdgesIneccEdges");   print(gmEdgesIneccEdges)
  ##print(gmEdges[gmEdgesIneccEdges])
  if (length(gmEdges[gmEdgesIneccEdges])>0)
    gmEdges[gmEdgesIneccEdges] <- NULL
  ##print("gmEdges (after removal)");  lprint(gmEdges);
  ## Turn remaining edges in gm into atomic ecc's; add these to ecc
  ##print("ecc (as given)"); lprint(ecc)
  ecc <- c(ecc, lapply(gmEdges,list))
  ##print("ecc (after adding edges from gm)"); lprint(ecc)

  ## Vertices given by gm which are also in vcc must be identified in gm
  vccVertices <- unlist(vcc)
  ##print("vccVertices (as given)"); lprint(vccVertices)
  gmVertices <- setdiff(gmVertices, vccVertices)
  ## Turn remaining vertives in gm into atomic vcc's; add these to vcc
  vcc <- c(vcc, lapply(gmVertices,list))
  vcc <- lapply(vcc, unlist, recursive=FALSE)
  ##print("vcc (after adding vertices from gm)"); lprint(vcc)
  
  ## We now have a standard representation of the model:
  ##cat("Standard representation:\n")
  ##cat("vcc:", paste(vcc),"\n")
  ##cat("ecc:", paste(ecc),"\n")
  return(list(vcc=vcc, ecc=ecc))  
}




findFitRepresentations <-function(vccI,eccI){

  vcc <- vccI; ecc <- eccI
  atomECC     <- getAtomCC(ecc)
  compECC     <- getCompCC(ecc)
  atomVCC     <- getAtomCC(vcc)
  compVCC     <- getCompCC(vcc)


  CS            <- NULL
  atomVCCvec  <- unlist(atomVCC)
  dropidx       <- NULL

  for (e in atomECC){
    e <- unlist(e);
    if (!any(is.na(match(e,atomVCCvec)))){ ## if T, the edge has atomic endpoints
      CS <- c(CS, list(e))
      dropidx <- c(dropidx, TRUE)
    } else {
      dropidx <- c(dropidx, FALSE)  
    }
  }

  vert    <- unique(unlist(CS))
  vertidx <- match(vert,atomVCCvec)
  
  atomVCC[vertidx] <- NULL  ## Now neutral sets are removed
  atomECC[dropidx] <- NULL  ## Now neutral sets are removed

  atomVCC <- unlist(atomVCC,recursive=FALSE)
  compVCC <- lapply(compVCC,unlist, recursive=FALSE)

  atomVCC <- if(length(atomVCC)>0) atomVCC
  compVCC <- if(length(compVCC)>0) compVCC
  atomECC <- if(length(atomECC)>0) atomECC
  compECC <- if(length(compECC)>0) compECC
  
  
  #lprint(m1$mixrep$CS); lprint(CS)
  #lprint(m1$mixrep$VCCU); lprint(atomVCC)
  #lprint(m1$mixrep$VCCR); lprint(compVCC)
  
  #lprint(m1$mixrep$ECCU); lprint(atomECC)
  #lprint(m1$mixrep$ECCR); lprint(compECC)
  
  mixrep <- list(CS=CS,VCCU=atomVCC, VCCR=compVCC, ECCU=atomECC, ECCR=compECC)
  #lprint(mixrep)
  #lprint(m1$mixrep)
  
  atomECC     <- getAtomCC(ecc)
  compECC     <- getCompCC(ecc)
  atomVCC     <- getAtomCC(vcc)
  compVCC     <- getCompCC(vcc)
  
  stdrep<-list(VCCU=lapply(getAtomCC(vcc),unlist,recursive=FALSE), 
               VCCR=lapply(getCompCC(vcc),unlist,recursive=FALSE), 
               ECCU=getAtomCC(ecc), ECCR=getCompCC(ecc))
  #lprint(stdrep)
  #lprint(m1$stdrep)
  
  return(list(stdrep=stdrep,mixrep=mixrep))  
}











# mparse <- function(gm=NULL, ecc=NULL, vcc=NULL){

#   varnames <- unique(unlist(c(gm,ecc,vcc)))
  
#   gme <- unlist(lapply(c(gm,unlist(ecc,recursive=FALSE)), cliqueToEdges),recursive=FALSE)

  
#   if (!is.null(gme)){
#     gme <- listOrder(gme)
#     gme <- lapply(gme, function(x) if (length(x)==1) c(x,NA) else x)
#     imat <- matrix(unlist(gme),ncol=2,byrow=TRUE)
#   } else {
#     x<- as.list(unlist(vcc))
#       imat <- cbind(x,x)
#   }
  
#   extravar <- setdiff(unlist(vcc),as.numeric(as.matrix(imat)))
#   if (length(extravar)>0)
#     imat <- rbind(imat,cbind(extravar,NA))
#   imat <- data.frame(unique(imat))
#   names(imat) <- c('a','b')
#   imat$a  <- as.numeric(imat$a) ## DIRTY HACK
#   imat$b  <- as.numeric(imat$b) ## DIRTY HACK

#   #print(imat)
#   if (all(is.na(imat$b)) | all(imat$a==imat$b))
#     imat <- data.frame('a'=NA,'b'=NA)

  
#   ## Add 'colour number' to each vertex colour class
#   ## Organize the classes in two data frames

#   if (!is.null(vcc)){
#     for (i in seq(along=vcc)){
#       x <- vcc[[i]]
#       if (length(x)==1)
#         vcc[[i]] <- c(x,NA)
#       else
#         vcc[[i]] <- as.numeric(t(cbind(x,i)))
#     }
#     vccmat  <- matrix(unlist(vcc), ncol=2,byrow=TRUE)
#   } else {
#     vccmat <- cbind(varnames,NA)
#   }
#   vccamat <- data.frame(vccmat); names(vccamat) <- c('a','a.col')
#   vccbmat <- data.frame(vccmat); names(vccbmat) <- c('b','b.col')

#   ## Add 'colour number' to each edge colour class
#   ## Organize the classes in a data frame
#   if (!is.null(ecc)){
#     for (i in seq(along=ecc)){
#       x<- ecc[[i]]
#       for (j in seq(along=x)){
#         if (length(x)==1){
#           x[[j]] <- c(x[[j]],NA)
#         } else {
#           x[[j]] <- c(x[[j]],i)
#         }
#         ecc[[i]] <- x
#       }
#     }
#     eccmat <- matrix(unlist(ecc), ncol=3,byrow=TRUE)
#     eccmat <- data.frame(eccmat);
#     names(eccmat) <- c('a','b','ab.col')
#     eccmat$a  <- as.numeric(eccmat$a) ## DIRTY HACK
#     eccmat$b  <- as.numeric(eccmat$b) ## DIRTY HACK

#   } else {
#     eccmat <- data.frame('a'=NA,'b'=NA, 'ab.col'=NA)
#   }

# #   if (!is.null(imat))
# #     imat.eccmat <- merge(imat,eccmat,all=TRUE)
# #   else
# #     imat.eccmat <- NULL

# #   if (!is.null(imat.eccmat)){
# #     aaaa<-merge(merge(imat.eccmat,vccamat,all=TRUE),vccbmat,all=TRUE)
# #   } else {
# #     aaaa<-merge(vccamat,vccbmat,all=TRUE)
# #   }
# #   print("aaaa");  print(aaaa)

#   ## Join all information into data frame with 5 columns
#   allmat <- merge(merge(merge(imat,eccmat,all=TRUE),vccamat,all=TRUE),vccbmat,all=TRUE)
#   #allmat <-aaaa
#   #print("allmat")
#   #print(allmat)
  
#   ## Extract rows with single variables which are not connected to other variables
#   ## All other rows contains edges
#   singlevarB<-!complete.cases(allmat[,1:2])
#   singlevar <- allmat[singlevarB,]
#   edgesMat  <- allmat[!singlevarB,]

#   ## Find edges with no restrictions and with no restrictions on endpoints
#   ## Find also the complement where NA is replaced by 0
#   zero.edgesB<-apply(is.na(edgesMat[,3:5]),1,all)
#   allZero <- edgesMat[zero.edgesB,]
#   allNZero <- edgesMat[!zero.edgesB,]
#   allNZero[is.na(allNZero)] <- 0

#   ## Find complete set
#   az <- as.matrix(allZero[,1:2])
#   dimnames(az) <- NULL
#   CS <- unlist(apply(az,1, list),recursive=FALSE)
#   CS <- listOrder(CS)
  
#   ## Consider edge colour classes
#   ed <- allNZero[,1:3]
#   elist <- split(ed, ed$ab.col)
  
#   ## 1) Singleton edges
#   e0 <- elist$'0'
#   if (!is.null(e0)){
#     e0 <- as.matrix(e0[,c('a','b')])
#     dimnames(e0)<-NULL
#     e0 <- apply(e0,1,list)
#   } else {
#     e0 <- NULL
#   }
  
#   ## 2) Non-singleton edges
#   e1 <- elist
#   e1$'0' <- NULL
#   names(e1)<-NULL
#   e2<- lapply(e1, function(x) {
#     x<-as.matrix(x[,1:2]); dimnames(x)<-NULL;
#     x<-unlist(apply(x,1,list),recursive=FALSE);
#   })
  
#   ##ECCU <- listOrder(e0) # Singletons
#   ##ECCR <- listOrder(e2) # Manytons
  
#   ECCU <- listOrder(if(length(e0)!=0) e0)
#   ECCR <- listOrder(if(length(e2)!=0) e2)
  
#   ## Consider vertex colour classes
  
#   allNZero <- rbind(allNZero, singlevar)
#   allNZero[,3:5][is.na(allNZero[,3:5])] <- 0
  
#   a <- allNZero[,c('a','a.col')]
#   b <- allNZero[,c('b','b.col')]
  
#   names(a) <- names(b) <- c('a','a.col')
#   ab <- unique(rbind(a,b))
#   ab <- ab[!is.na(ab[,1]),]
#   vlist <- split(ab,ab$a.col)
  
#   ## 1) Singleton vertices
#   v0 <- vlist$'0'
#   if (!is.null(v0)){
#     v0 <- (v0)[,'a']
#     v0 <- as.list (v0)
#   } else {
#     v0 <- NULL
#   }
  
  
#   ## 2) Non-singleton vertices
#   v1 <- vlist
#   v1$'0' <- NULL
#   names(v1)<- NULL
#   v1<-lapply(v1,function(x)x[,'a'])
  
  
#   VCCU <- listOrder(if(length(v0)!=0) v0)
#   VCCR <- listOrder(if(length(v1)!=0) v1)
  
#   ## Remove elements from VCCU which are already contained in CS
  
#   VCCU <- as.list(setdiff(unlist(VCCU), unlist(CS)))
#   VCCU <- listOrder(if(length(VCCU)!=0) VCCU)
  
#   mixrep <- list(CS=CS,VCCU=VCCU,VCCR=VCCR,ECCU=ECCU,ECCR=ECCR)
  
#   VCCU <- c(VCCU,as.list(unique(unlist(CS))))
#   VCCU <- listOrder(if(length(VCCU)!=0) VCCU)
  
   
#   stdrep <- list(VCCU=VCCU,  VCCR=VCCR,  ECCU=c(ECCU,lapply(CS,list)), ECCR=ECCR)
  
#   return(list(stdrep=stdrep,mixrep=mixrep))
# }
