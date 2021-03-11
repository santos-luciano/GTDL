#'@name MaxGTDL2
#'@title 
#'
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'
#'
#'Colosimo, E. A.,() Análise de sobrevivência aplicada

likeGTDL2 <- function(t1,censura,para){
  l <- (hGTDL(t = t1,param = para)^censura)*sGTDL(t = t1,param = para)
  ll <- prod(l)
  return(ll)  
}


#'@rdname MaxGTDL2
#'@export

MaxGTDL2 <- function(start,t,censur,...){
  
  op <- suppressWarnings(optim(par = start,fn = likeGTDL2,method = "BFGS",t1 = t,censura = censur,hessian = TRUE))
  se <- sqrt(diag(solve(op$hessian)))
  z <- op$par/se
  pvalue <- 2 * (1 - stats::pnorm(abs(z)))
  TAB <- cbind(Estimate = op$par, Std.Error = se,
               z.value = z, `Pr(>|z|)` = pvalue)
  mTab <- list( Lik = op$value,
                Converged = op$convergence, Coefficients = TAB)
  return(mTab)
  
}

