#'@name max.GTDL
#'@title 
#'
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'
#'
#'Colosimo, E. A.,() Análise de sobrevivência aplicada

like2 <- function(t,censur,para){
  l <- (hGTDL(t = t,param = para)^censur)*sGTDL(t = t,param = para)
  ll <- sum(log(l))
  return(-ll)  
}


#'@rdname max.GTDL
#'@export

max.GTDL <- function(start,t,censur){
  
  op <- suppressWarnings(optim(par = start,fn = like2,
                               method = "BFGS",t = t,censur = censur,hessian = TRUE))
  se <- sqrt(diag(solve(op$hessian)))
  z <- op$par/se
  pvalue <- 2 * (1 - stats::pnorm(abs(z)))
  TAB <- cbind(Estimate = op$par, Std.Error = se,
               z.value = z, `Pr(>|z|)` = pvalue)
  mTab <- list( Lik = op$value,
                Converged = op$convergence, Coefficients = TAB)
  return(mTab)
  
}


