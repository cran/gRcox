
rcox.control <- function (logLepsilon = 1e-03,
                          IPMepsilon=1e-6,
                          maxit = 25,
                          method=NULL,
                          representation='mixrep',
                          VCCU='IPS',
                          VCCR='NR',
                          ECCU='IPS',
                          ECCR='NR',
                          trace = 0) 
{
  if (!is.numeric(logLepsilon) || logLepsilon <= 0) 
    stop("value of logLepsilon must be > 0")
  if (!is.numeric(maxit) || maxit <= 0) 
    stop("maximum number of iterations must be > 0")
  list(logLepsilon   = logLepsilon,
       IPMepsilon = IPMepsilon,
       maxit     = maxit,
       method    = method,
       representation=representation,
       VCCU=VCCU,
       VCCR=VCCR,
       ECCU=ECCU,
       ECCR=ECCR,
       trace = trace)
}
