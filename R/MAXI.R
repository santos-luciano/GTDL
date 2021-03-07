#'@name MaxGTDL
#'@title Maximum probability estimate of the GTDL package
#'
#'@description This is the main interface of the maxLik package and the function that performs the maximum probability estimate of the GTDL package
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'@references
#'
#' Aarset,M. V.,(1987). How to identify bathtub hazard rate. IEEE transactions on reliability.
#' Mackenzie,G.,(1996). Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family. Journal of the Royal Statistical Society. Series D (The Statistician), (45), 21-34.
#'
#'@examples
#'
#' ## times data (from Arset, 1987))
#'
#'data(artset1987)
#'t<-artset1987
#'max<-MaxGTDL(c(1,-0.05,-1))
#'

likeGTDL<-function(param,t1){ 
  f1<-sum(dGTDL(param = param,t = t1,log = TRUE))
  return(-f1)
}



#'@rdname MaxGTDL
#'@export
#'

MaxGTDL<-function(start,t,...){
  
  op<-suppressWarnings(optim(par = start,fn = likeGTDL,method = "BFGS",t1 = t,hessian = TRUE))
  se <- sqrt(diag(solve(op$hessian)))
  z <- op$par/se
  pvalue <- 2 * (1 - stats::pnorm(abs(z)))
  TAB <- cbind(Estimate = op$par, Std.Error = se,
               z.value = z, `Pr(>|z|)` = pvalue)
  mTab <- list( Lik = op$value,
                Converged = op$convergence, Coefficients = TAB)
  return(mTab)
  
}
