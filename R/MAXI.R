#'@name MaxGTDL
#'@title Maximum probability estimate of the GTDL package
#'
#'@description This is the main interface of the maxLik package and the function that performs the maximum probability estimate of the GTDL package
#'
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'
#'
#'@references
#'
#' How to identify bathtub hazard rate (Aarset,1987 )
#'
#'@examples
#'t<-c(0.1,0.2,1,1,1,1,1,2,3,6,7,11,12,18,18,18,18,18,21,32,36,40,45,46,47,50,55,60,63,63,67,67,67,67,72,75,79,82,82,83,84,84,84,85,85,85,85,86,86)
#'max<-MaxGTDL(c())




MaxGTDL<-function(start,...){
  likeGTDL<-function(param){
    lambda<-param[1]
    alpha<-param[2]
    gamma<-param[3]
    f<-sum(dGTDL(t,lambda,alpha,gamma,log = TRUE))
    return(f)
  }
  aux<-maxLik(likeGTDL,start = start
              ,grad = NULL,hess = NULL)
  return(aux)
}

