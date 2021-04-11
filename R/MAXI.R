likeGTDL <- function(param,t){ 
  f1 <- sum(dGTDL(param = param,t = t,log = TRUE))
  return(-f1)
}

#'@name MaxGTDL
#'@title Maximum probability estimate of the GTDL package
#'
#'@description The maximum likelihood estimation of the GTDL distribution
#'
#'@param start vector of parameters to obtaind maximum likelihood.
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'@references
#'
#' Aarset, M. V. (1987). How to Identify a Bathtub Hazard Rate. IEEE Transactions on Reliability, 36, 106–108.
#' Mackenzie, G. (1996) Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family. Journal of the Royal Statistical Society. 
#' Series D (The Statistician). (45). 21-34.
#'
#'@examples
#' # times data (from Aarset, 1987))
#'data(artset1987)
#'mod <- MaxGTDL(c(1,-0.05,-1),t = artset1987)
#'

#'@rdname MaxGTDL
#'@export
MaxGTDL <- function(start,t,method = 'BFGS'){
  op <- suppressWarnings(optim(par = start,fn = likeGTDL,method = method,t = t,hessian = TRUE))
  se <- sqrt(diag(solve(op$hessian)))
  z <- op$par/se
  pvalue <- 2 * (1 - stats::pnorm(abs(z)))
  TAB <- cbind(Estimate = op$par, Std.Error = se,
               z.value = z, `Pr(>|z|)` = pvalue)
  mTab <- list( Lik = op$value,
                Converged = op$convergence, Coefficients = TAB)
  rownames(mTab$Coefficients) <- c("lambda","alpha","beta")
  return(mTab)
}

