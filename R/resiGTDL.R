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
       pch = 20, ylim = fy,cex.axis = 1.2,
       cex.lab = 1.2, cex = 0.6, bg = 5)
  par(new = T)
  plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
  par(new = T)
  plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "", lwd=2)
  par(new = T)
  plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "", lwd=1.2)
}

#'@title Residual value of the GTDL distribution 
#'
#'
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param param These are the lambda, alpha and gamma parameters of the GTDL distribution.
#'@param censur absence of the occurrence of the event at the time of analysis.
#'@param formula The structure matrix of covariates of dimension n x p
#'@describe TESTE
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'
#'@examples
#'
#'### Example 1
#'
#'require(survival)
#'data(lung)
#'lung <- lung[-14,]
#'lung$ph.ecog[lung$ph.ecog==3]<-2
#'t1 <- lung$time
#'formula1 <- ~lung$sex+factor(lung$ph.ecog)+lung$age
#'censur1 <- ifelse(lung$status==1,0,1)
#'start1 <- c(0.03,0.05,-1,0.7,2,-0.1)
#'fit.model1 <- mle2.GTDL(t = t1,start = start1,
#'            formula = formula1,
#'            censur = censur1)
#'r1 <- q.GTDL(t = t1,formula = formula1 ,param = fit.model1$Coefficients[,1],
#'              censur = censur1)
#'r1
#'
#'### Example 2
#'
#'data(tumor)
#'t2 <- tumor$time
#'formula2 <- ~tumor$group
#'censur2 <- tumor$censured
#'start2 <- c(1,-0.05,1.7)
#'fit.model2 <- mle2.GTDL(t = t2,start = start2,
#'                        formula = formula2,
#'                        censur = censur2)
#'r2 <- q.GTDL(t = t2,formula = formula2, param = fit.model2$Coefficients[,1],
#'             censur = censur2)
#'r2


#'@rdname resiGTDL
#'@export
#'@import stats
#'@import graphics

q.GTDL <- function(t,formula,param,censur){
  x.aux <- model.matrix(formula)
  x <- matrix(x.aux[,-1],ncol = (ncol(x.aux)-1))
  p <- ncol(data.matrix(x))
  conf <- NULL
  for(i in 1:length(t)){
    gamma.aux <-  x[i,]%*%matrix(param[3:(p+2)])
    param.aux <- c(param[1:2],gamma.aux)
    conf[i] <- reability.GTDL(t=t[i],param=param.aux)  
  }
  
  qr <- qnorm(censur* (1 - conf) + (1-censur)*runif(length(t),1-conf))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  graphics::par(mfrow = c(1, 2), pty = "s", col = "royalblue")
  plot(qr, xlab = "Index", ylab = "Quantile residuals",
       pch = 15, main = "", cex.axis = 1.2,
       cex.lab = 1.2, cex = 0.6, bg = 5)
  envelope.GTDL(qr)
  return(qr)
}


