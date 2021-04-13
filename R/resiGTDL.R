#'@name resi.GTDL
#'@title Residual value of the GTDL distribution 
#'
#'
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param param These are the lambda, alpha and gamma parameters of the GTDL distribution.
#'@param censura absence of the occurrence of the event at the time of analysis.

#'
#'@describe
#'
#'
#'
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'
#'@examples
#'
#'data(lung)
#'lung <- lung[-14,]
#'lung$ph.ecog[lung$ph.ecog==3]<-2
#'formula <- ~lung$sex+factor(lung$ph.ecog)+lung$age
#'censur <- ifelse(lung$status==1,0,1)
#'start <- c(0.03,0.05,-1,0.7,2,-0.1)
#'max_para <- maxGTDL(t = lung$time,start = start,
#'            formula = formula,
#'            censur = censur)
#'x <- quantile.GTDL(t = lung$time,param = max_para$Coefficients[,1],
#'              censura = censur)
#'envole.GTDL(x)

#'@rdname resi.GTDL
#'@export

quantileGTDL <- function(t,formula,param,censura){
  x.aux <- model.matrix(formula)
  x <- x.aux[,-1]
  p <- ncol(data.matrix(x))
  conf <- NULL
  for(i in 1:length(t)){
    gamma.aux <-  x[i,]%*%matrix(param[3:(p+2)])
    param.aux <- c(param[1:2],gamma.aux)
    conf[i] <- reability.GTDL(t=t[i],param=param.aux)  
  }
  
  qr <- qnorm(censura* (1 - conf) + (1-censura)*runif(length(t),1-conf))
  return(qr)
}


#'@rdname resi.GTDL
#'@export

envelope.GTDL <- function(x){
  U	         <- x
  n	         <- length(x)
  d2s 	     <- sort(U)
  xq2 	     <- qnorm(ppoints(n))
  Xsim 	     <- matrix(0, 100, n)
  for(i in 1:100){
    u2       <- rnorm(n)
    Xsim[i,] <- u2
  }
  Xsim2      <- apply(Xsim, 1, sort)
  d21        <- matrix(0, n, 1)
  d22        <- matrix(0, n, 1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,], 0.025)
    d22[i]  <- quantile(Xsim2[i,], 0.975)
  }
  d2med      <- apply(Xsim2, 1, mean)
  fy         <- range(d2s, d21, d22)
  plot(xq2, d2s, xlab = quote("Theoretical Quantiles"),
       ylab = quote("Quantile Residuals"), 
       pch = 20, ylim = fy)
  par(new = T)
  plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
  par(new = T)
  plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "", lwd=2)
  par(new = T)
  plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
}

