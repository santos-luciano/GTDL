#'@name MaxGTDL
#'@title Maximum probability estimate of the GTDL package
#'
#'@description This is the main interface of the maxLik package and the function that performs the maximum probability estimate of the GTDL package
#'
#'@param start vector of parameters to obtaind maximum likelihood.
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'@references
#'
#' Aarset,M. V. How to identify bathtub hazard rate. IEEE transactions on reliability. Two full days, Washington, DC area, v. r-36, n. 1 1987.
#' Mackenzie, G. Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family. Journal of the Royal Statistical Society. Series D (The Statistician),  v. 45, n. 1, p. 21-34. 1996.
#'
#'@examples
#'
#' # times data (from Aarset, 1987))
#'
#'data(artset1987)
#'mod <- MaxGTDL(c(1,-0.05,-1))
#'


#'@rdname MaxGTDL
#'@export
MaxGTDL <- function(start,t,...){
  
  op <- suppressWarnings(optim(par = start,fn = likeGTDL,method = "BFGS",t1 = t,hessian = TRUE))
  se <- sqrt(diag(solve(op$hessian)))
  z <- op$par/se
  pvalue <- 2 * (1 - stats::pnorm(abs(z)))
  TAB <- cbind(Estimate = op$par, Std.Error = se,
               z.value = z, `Pr(>|z|)` = pvalue)
  mTab <- list( Lik = op$value,
                Converged = op$convergence, Coefficients = TAB)
  return(mTab)
  
}

likeGTDL <- function(param,t1){ 
  f1 <- sum(dGTDL(param = param,t = t1,log = TRUE))
  return(-f1)
}


