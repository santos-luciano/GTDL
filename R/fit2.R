like2 <- function(t,formula,censur,para){
  x_aux <- model.matrix(formula)
  x <- matrix(x_aux[,-1],ncol = (ncol(x_aux)-1))
  p <- ncol(data.matrix(x))
  ll <- NULL
  for(i in 1:dim(x)[1]){
    gama_aux <- x[i,]%*%matrix(para[3:(p+2)])
    para_aux <- c(para[1:2],gama_aux)  
    l <- (hGTDL(t = t[i],param = para_aux)^censur[i])*sGTDL(t = t[i],param = para_aux)
    ll[i] <- log(l)
  }
  return(-sum(ll))  
}

#'@title Maximum probability estimate of the GTDL package
#'
#'@param start vector of parameters to obtaind maximum likelihood.
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param censur non-negative random variable that represents whether the sample is censored or not.
#'@param formula The structure matrix of covariates of dimension n x p
#'@param method The method to be used
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'@references
#'
#'Colosimo, E. A and  Giolo, S. R. Análise de sobrevivência aplicada.  Edgard Blucher: São Paulo. 2006.
#'
#'@examples
#'
#'require(survival)
#'lung <- lung[-14,]
#'lung$ph.ecog[lung$ph.ecog==3]<-2
#'formula <- ~lung$sex+factor(lung$ph.ecog)+lung$age
#'censur <- ifelse(lung$status==1,0,1)
#'start <- c(0.03,0.05,-1,0.7,2,-0.1)
#'fit.model <- mle2.GTDL(t = lung$time,start = start,
#'                      formula = formula,
#'                      censur = censur)
#'fit.model
#'
#'@rdname fit2
#'@export
#'@import stats
mle2.GTDL <- function(t,start,formula,censur,method = "BFGS"){
   op <- suppressWarnings(optim(par = start,fn = like2,
                               t = t,
                               method = method,
                               formula = formula,
                               censur = censur,
                               hessian = TRUE))
  se <- sqrt(diag(solve(op$hessian)))
  z <- op$par/se
  pvalue <- 2 * (1 - stats::pnorm(abs(z)))
  TAB <- cbind(Estimate = op$par, Std.Error = se,
               z.value = z, `Pr(>|z|)` = pvalue)
  mTab <- list( Lik = op$value,
                Converged = op$convergence, Coefficients = TAB)
  rownames(mTab$Coefficients) <- c("lambda","alpha",paste0("beta ",c(1:(dim(mTab$Coefficients)[1]-2))))
  mTab$Lik <- round(mTab$Lik,4)
  mTab$Coefficients <- round(mTab$Coefficients,4)
  return(mTab)
  
}
