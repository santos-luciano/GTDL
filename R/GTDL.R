#'@name GTDL
#'@title The Distribution GTDL
#'
#'@description Density, survival function, failure function and random generation for the GTDL distribution.
#'
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param alpha,gamma scalars.
#'@param lambda non-negative.
#'@param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'@param pcensura sample censorship rate
#'@param log logical; if TRUE, probabilities p are given as log(p).
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'
#'
#'@references
#'
#' Mackenzie,G.,(2016).Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family. Journal of the Royal Statistical Society. Series D (The Statistician), Vol. 45, No. 1 (1996), pp. 21-34
#'
#'@examples
#' 
#' library(GTDL)
#' t <- seq(0,20,by = 0.1)
#' lambda <- 1.00
#' alpha <- -0.05
#' gamma <- -1.00
#' y1 <- hGTDL(t,lambda,alpha,gamma)
#' y2 <- sGTDL(t,lambda,alpha,gamma)
#' y3 <- dGTDL(t,lambda,alpha,gamma,log = FALSE)
#' tt <- as.matrix(cbind(t,t,t))
#' yy <- as.matrix(cbind(y1,y2,y3))
#' matplot(tt,yy,type="l",xlab="time",ylab="",lty = 1:3,col=1:3,lwd=2)
#' 
#' 
#' y1 <- hGTDL(t,1,0.5,-1.0)
#' y2 <- hGTDL(t,1,0.25,-1.0)
#' y3 <- hGTDL(t,1,-0.25,1.0)
#' y4 <- hGTDL(t,1,-0.50,1.0)
#' y5 <- hGTDL(t,1,-0.06,-1.6)
#' tt <- as.matrix(cbind(t,t,t,t,t))
#' yy <- as.matrix(cbind(y1,y2,y3,y4,y5))
#' matplot(tt,yy,type="l",xlab="time",ylab="Hazard function",lty = 1:3,col=1:3,lwd=2)
#' 
NULL


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




#'@rdname GTDL
#'@export

rGTDL<-function(n,lambda,alpha,gamma){
  u<-runif(n)
  t<-(1/alpha)*(log((1+exp(gamma))*(1-u)^(-alpha/lambda)-1)-gamma)
  return(t)
}
