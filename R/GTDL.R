#'@name GTDL
#'@title The Distribution GTDL
#'
#'@description Density, survival function, failure function and random generation for the GTDL distribution.
#'
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param alpha,gamma scalars.
#'@param lambda non-negative.
#'@param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'
#'
#'@references
#'
#' Mackenzie,G. ,(2016),Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family
#'
#'@examples
#'t <- seq(0,20,by = 0.01)
#'lambda <- 1
#'alpha <- -0.05
#'gamma <- -1
#'denGTDL<- dGTDL(t,lambda,alpha,gamma,log = FALSE)
#'plot(x = t,y = denGTDL)

#'@rdname GTDL 
#'@export
dGTDL<-function(t,lambda,alpha,gamma,log = FALSE){
  if(log == FALSE)
  ((lambda*exp(t*alpha+gamma))/(1+exp(t*alpha+gamma)))*((1+exp(t*alpha+gamma))/(1+exp(gamma)))^(-lambda/alpha)
  else log((lambda*exp(t*alpha+gamma)/(1+exp(t*alpha+gamma)))*((1+exp(t*alpha+gamma))/(1+exp(gamma)))^(-lambda/alpha))
  }
#'@rdname GTDL 
#'@export
hGTDL<-function(t,lambda,alpha,gamma){
  (lambda*exp(t*alpha+gamma))/(1+exp(t*alpha+gamma))
}
#'@rdname GTDL 
#'@export
sGTDL<-function(t,lambda,alpha,gamma){
  dGTDL(t,lambda,alpha,gamma)/hGTDL(t,lambda,alpha,gamma)
}

