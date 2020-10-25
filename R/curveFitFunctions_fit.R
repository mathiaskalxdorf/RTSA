fitSigmoidTR <- function(xVec, yVec, startPars, maxAttempts, fixT0=TRUE){
  ## Fit sigmoidal model to a vector of TPP-TR measurements
  strSigm <- fctSigmoidTR(deriv=0)
  fitFct <- as.formula(paste("y ~", strSigm))
  varyPars <- 0
  attempts <- 0
  repeatLoop <- TRUE
  
  ## Check if number of non-missing values is sufficient
  ## (NLS can only handle data with at least three non-missing values)
  validValues <- !is.na(yVec)
  if (sum(validValues) <=2){
    m <- NA
    class(m) <- "try-error"
  } else{  
    ## Perform fit
    while(repeatLoop & attempts < maxAttempts){
      parTmp <- startPars * (1 + varyPars*(runif(1, -0.2, 0.2)))      
      m <- try(nls(formula=fitFct, start=parTmp, data=list(x=xVec, y=yVec), 
                   na.action=na.exclude, algorithm="port", 
                   lower=c(0.0, 1e-5, 1e-5), upper=c(1.5, 15000, 250)), 
               silent=TRUE)
      attempts <- attempts + 1
      varyPars <- 1
      if (class(m)!="try-error") {
        repeatLoop <- FALSE
      }
    }
  }
  return(m)
}

fitSigmoidTR2 <- function(xVec, yVec, startPars, maxAttempts, T0){
  ## Fit sigmoidal model to a vector of TPP-TR measurements
  strSigm <- fctSigmoidTR2(deriv=0,T0)
  fitFct <- as.formula(paste("y ~", strSigm))
  varyPars <- 0
  attempts <- 0
  repeatLoop <- TRUE
  
  ## Check if number of non-missing values is sufficient
  ## (NLS can only handle data with at least three non-missing values)
  validValues <- !is.na(yVec)
  if (sum(validValues) <=2){
    m <- NA
    class(m) <- "try-error"
  } else{  
    ## Perform fit
    while(repeatLoop & attempts < maxAttempts){
      parTmp <- startPars * (1 + varyPars*(runif(1, -0.2, 0.2)))      
      m <- try(nls(formula=fitFct, start=parTmp, data=list(x=xVec, y=yVec), 
                   na.action=na.exclude, algorithm="port", 
                   lower=c(0.0, 1e-5, 1e-5), upper=c(1.5, 15000, 250)), 
               silent=TRUE)
      attempts <- attempts + 1
      varyPars <- 1
      if (class(m)!="try-error") {
        repeatLoop <- FALSE
      }
    }
  }
  return(m)
}

fitSigmoidCCR <- function(xVec, yVec, hill_init, pec50_init, slopeBounds, 
                          concBounds){
  ## Fit sigmoidal model to a vector of TPP-CCR measurements
  
  ## Prepare model fit:
  strSigm <- fctSigmoidCCR()
  fitFct <- as.formula(paste("y ~", strSigm))
  
  ## Attempt model fit by numerical optimization with nls:
  lower <- c(slopeBounds[1], concBounds[1])
  upper <- c(slopeBounds[2], concBounds[2])
  startPars <- list(hill=hill_init, infl=pec50_init)
  m <- try(nls(formula=fitFct, algorithm="port", data=list(x=xVec, y=yVec), 
               start=startPars, lower=lower, upper=upper, na.action=na.exclude),
           silent=TRUE)
  
  ## Check if fit was successful and if estimated parameters have sufficient quality:
  retry <- FALSE
  if(class(m) == "try-error") {
    retry <- TRUE
  } else {
    ## If fit was successful extract pEC50 and Hill slope for quality check
    coeffsTmp <-coef(m) 
    hill  <- coeffsTmp["hill"]
    pec50 <- coeffsTmp["infl"]
    if (!(pec50 >= concBounds[1] & pec50 <= concBounds[2] & sign(hill)==sign(hill_init))){
      retry <- TRUE
    }
  }
  
  ## If fit was not successful, or did not yield satisfactory curve parameters, 
  ## repeat by 'naive' grid search algorithm with nls2:
  if (retry==TRUE){
    startNLS2 <- list(hill=slopeBounds, infl=concBounds)
    m <- try(nls2(formula=fitFct, algorithm="grid-search", 
                  data  = list(x=xVec, y=yVec), 
                  start = startNLS2, na.action=na.exclude),
             silent=TRUE)  
  }
  
  return(m)
}

fctSigmoidTR <- function(deriv){
  ## Return sigmoidal function or its derivatives used for curve fit 
  ## in TPP-TR experiments
  if (deriv == 0){
    fctStr <- "(1 - Pl) * 1 / (1+exp(-(a/x-b))) + Pl"
  } else if (deriv == 1){
    fctStr <- "-((1 - Pl) * (exp(-(a/x - b)) * (a/x^2))/(1 + exp(-(a/x - b)))^2)"
  } else if (deriv == 2){
    fctStr <- "-((1 - Pl) * 1 * (exp(-(a/x - b)) * (a/x^2) * (a/x^2) - exp(-(a/x - b)) * (a * (2 * x)/(x^2)^2))/
             (1 + exp(-(a/x - b)))^2 - (1 - Pl) * 1 * (exp(-(a/x - b)) * (a/x^2)) *
             (2 * (exp(-(a/x - b)) * (a/x^2) * (1 + exp(-(a/x - b)))))/((1 + exp(-(a/x - b)))^2)^2)"
  }
  return(fctStr)
}
fctSigmoidTR2 <- function(deriv,T0){
  ## Return sigmoidal function or its derivatives used for curve fit 
  ## in TPP-TR experiments
  if (deriv == 0){
    fctStr <- paste("(",as.character(T0)," - Pl) * 1 / (1+exp(-(a/x-b))) + Pl",sep="")
  } else if (deriv == 1){
    fctStr <- "-((1 - Pl) * (exp(-(a/x - b)) * (a/x^2))/(1 + exp(-(a/x - b)))^2)"
  } else if (deriv == 2){
    fctStr <- "-((1 - Pl) * 1 * (exp(-(a/x - b)) * (a/x^2) * (a/x^2) - exp(-(a/x - b)) * (a * (2 * x)/(x^2)^2))/
             (1 + exp(-(a/x - b)))^2 - (1 - Pl) * 1 * (exp(-(a/x - b)) * (a/x^2)) *
             (2 * (exp(-(a/x - b)) * (a/x^2) * (1 + exp(-(a/x - b)))))/((1 + exp(-(a/x - b)))^2)^2)"
  }
  return(fctStr)
}

fctSigmoidCCR <- function(){
  ## Return sigmoidal function for curve fit in TPP-CCR experiments
  fctStr <- "1 / (1 + exp((infl - x) * hill))"
  return(fctStr)
}
