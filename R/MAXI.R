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
#' Aarset,M. V.,(1987).How to identify bathtub hazard rate.Ieee transactions on reliability, vol. r-36, no. 1.
#' Mackenzie,G.,(2016).Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family. Journal of the Royal Statistical Society. Series D (The Statistician), Vol. 45, No. 1 (1996), pp. 21-34
#'
#'@examples
#'
#' ## times data (from Arset, 1987))
#'
#'data(artset1987)
#'t<-artset1987
#'max<-MaxGTDL(c(1,-0.05,-1))
#'
#'@rdname MaxGTDL
#'@export

likeGTDL<-function(param){
  lambda<-param[1]
  alpha<-param[2]
  gamma<-param[3]
  f<-sum(dGTDL(t,lambda,alpha,gamma,log = TRUE))
  return(f)
}

#'@rdname MaxGTDL
#'@export
#'@import maxLik
#'

MaxGTDL<-function(start,...){
  aux<-suppressWarnings(maxLik(likeGTDL,start = start
              ,grad = NULL,hess = NULL))
  return(aux)
}
