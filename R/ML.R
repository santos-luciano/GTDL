like2 <- function(t,censur,para){
  l <- (hGTDL(t = t,param = para)^censur)*sGTDL(t = t,param = para)
  ll <- sum(log(l))
  return(-ll)  
}

#'@name max.GTDL
#'@title Maximum probability estimate of the GTDL package
#'
#'@param start vector of parameters to obtaind maximum likelihood.
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param censur non-negative random variable that represents whether the sample is censored or not.
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'@references
#'
#'Colosimo, E. A and  Giolo, S. R. Análise de sobrevivência aplicada.  Edgard Blucher: São Paulo. 2006.
#'
#'
#'@examples
#'
#'data(hepatis)
#'maxparam <- max.GTDL(start = c(1,-0.05,-1),t = hepatitis$t, censur = hepatitis$censured)
#'maxparam


#'@rdname max.GTDL
#'@export
maxGTDL <- function(start,t,censur,method = "BFGS" ){
  
  op <- suppressWarnings(optim(par = start,fn = like2,
                               method = method,t = t,censur = censur,hessian = TRUE))
  se <- sqrt(diag(solve(op$hessian)))
  z <- op$par/se
  pvalue <- 2 * (1 - stats::pnorm(abs(z)))
  TAB <- cbind(Estimate = op$par, Std.Error = se,
               z.value = z, `Pr(>|z|)` = pvalue)
  mTab <- list( Lik = op$value,
                Converged = op$convergence, Coefficients = TAB)
  return(mTab)
  
}



