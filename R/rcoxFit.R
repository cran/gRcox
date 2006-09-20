
rcoxFit  <- function(m,control=m$control){

  ts <- proc.time()
  
  val       <- rcoxFitPrimitive(m, control=control)
  val$theta <- getThetaFromK(m,K=val$K)
  val$J     <- getInf       (m,K=val$K)
  
  val$timetaken <- (proc.time()-ts)[3]
  val$df    <- getdf(m)
  val$dim   <- length(getCC(m,type='ecc'))+length(getCC(m,type='vcc'))

  val$C     <- CorMatrix <- -cov2cor(val$K)
  KC        <- val$K
  KC[lower.tri(KC)] <- CorMatrix[lower.tri(CorMatrix)]
  val$KC    <- KC
  val$Jinv  <- solve(val$J)
  m$isfit   <- TRUE
  m$fit     <- val
  return(m)  
}

rcoxFitPrimitive <- function(m,control=m$control) UseMethod('rcoxFitPrimitive')

rcoxFitPrimitive.rcon <- function(m,control=m$control){
  method <- control$method
  if (control$trace>=3) cat("...rcoxFitPrimitive.rcon\n")  
  switch(method,
         'scoring'={       multivariateNewton(m,K=m$K,control=m$control)         },
         'ipm' ={          rconFitIterative(m,control)         },
         'user'={          rconFitIterative(m,control)         },                  
         'hyd'={           copyInitial(m)         }         
         ##'newton'={        rconFitIterative(m,control)        },
         )
} 

rcoxFitPrimitive.rcor <- function(m,control=m$control){
  method <- control$method
  ##print(method)
  switch(method,
         'scoring'={        multivariateNewton(m,K=m$K,control=m$control)        },
         'ipm'    ={        rcorFitIterative(m,control)         },
         'hyd'    ={        copyInitial(m)         }
         ##'newton'={         rcorFitIterative(m,control)         },           
         )
} 


copyInitial <- function(m){
  return(list(K=m$Kstart, logL=m$initlogL, vccTerms=m$vccTerms, eccTerms=m$eccTerms))
}
