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



#'@rdname MaxGTDL
#'@export
#'

MaxGTDL<-function(start,t,...){
  likeGTDL<-function(param,t){ 
    f1<-sum(dGTDL(param = param,t = t,log = TRUE))
    return(-f1)
  }
  aux<-suppressWarnings(optim(par = start,fn = likeGTDL,method = "BFGS",t = t,hessian = TRUE))
              
  return(aux)
}
