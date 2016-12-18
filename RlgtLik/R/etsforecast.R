forecast.etsLGT <- function(object, h=ifelse(object$m>1, 2*object$m, 10),
    level=c(80,95), npaths=5000, PI=FALSE, ...)
{
  # Check inputs
  #if(h>2000 | h<=0)
  if(h <= 0)
    stop("Forecast horizon out of bounds")
  
  if(min(level) > 0 & max(level) < 1)
    level <- 100*level
  else if(min(level) < 0 | max(level) > 99.99)
    stop("Confidence limit out of range")
  
  # Order levels
  level <- sort(level)
  
  n <- length(object$x)
  #damped <- as.logical(object$components[4])
  
  f <- pegelsfcast.C(h,object,level=level,bootstrap=bootstrap,npaths=npaths)
  
  tsp.x <- tsp(object$x)
  if(!is.null(tsp.x))
    start.f <- tsp(object$x)[2] + 1/object$m
  else
    start.f <- length(object$x)+1
  out <- list(model=object,mean=ts(f$mu,frequency=object$m,start=start.f),level=level,x=object$x)

  #  if(PI)
#  {
#    if(!is.null(f$var))
#    {
#      out$lower <- out$upper <- ts(matrix(NA,ncol=length(level),nrow=h))
#      for(i in 1:length(level))
#      {
#        marg.error <- sqrt(f$var) * abs(qnorm((100-level[i])/200))
#        out$lower[,i] <- out$mean - marg.error
#        out$upper[,i] <- out$mean + marg.error
#      }
#      tsp(out$lower) <- tsp(out$upper) <- tsp(out$mean)
#    }
#    else if(!is.null(f$lower))
#    {
#      out$lower <- ts(f$lower)
#      out$upper <- ts(f$upper)
#      tsp(out$lower) <- tsp(out$upper) <- tsp(out$mean)
#    }
#    else if(PI)
#      warning("No prediction intervals for this model")
#    else if(biasadj)
#      warning("No bias adjustment possible")
#  }
  
  out$fitted <- fitted(object)
  out$method <- object$method
  out$residuals <- residuals(object)
  
#  if(!is.null(lambda))
#  {
#    #out$x <- InvBoxCox(object$x,lambda)
#    #out$fitted <- InvBoxCox(out$fitted,lambda)
#    out$mean <- InvBoxCox(out$mean,lambda)
#    if(biasadj){
#      out$mean <- InvBoxCoxf(x = out, lambda = lambda)
#    }
#    if(PI)  # PI = TRUE
#    {
#      out$lower <- InvBoxCox(out$lower,lambda)
#      out$upper <- InvBoxCox(out$upper,lambda)
#    }
#  }
  if(!PI)
    out$lower <- out$upper <- out$level <- NULL
  
  return(structure(out,class="forecast"))
}

pegelsfcast.C <- function(h,obj,npaths,level,bootstrap)
{
#  y.paths <- matrix(NA,nrow=npaths,ncol=h)
#  obj$lambda <- NULL # No need to transform these here as we do it later.
#  for(i in 1:npaths)
#    y.paths[i,] <- simulate.ets(obj, h, future=TRUE, bootstrap=bootstrap)
  
  print("pegelsfcast.C")
  
  
  
  y.f <- .C("etsforecast",
      as.double(obj$state[length(obj$x)+1,]),
      as.integer(obj$m),
      as.integer(switch(obj$components[2],"N"=0,"A"=1,"M"=2)),
      as.integer(switch(obj$components[3],"N"=0,"A"=1,"M"=2)),
      as.double(obj$par["alpha"]),
      as.double(obj$par["beta"]),
      as.double(ifelse(obj$components[4]=="FALSE",1,obj$par["phi"])),
      as.double(obj$par["lambda"]),
      as.double(obj$par["rho"]),
      as.integer(h),
      as.double(numeric(h)),
      PACKAGE="RlgtLik")[[11]]
  
#browser()
  
  print(y.f)
  
  
  if(abs(y.f[1]+99999) < 1e-7)
    stop("Problem with multiplicative damped trend")
  
#  lower <- apply(y.paths,2,quantile,0.5 - level/200, type=8, na.rm=TRUE)
#  upper <- apply(y.paths,2,quantile,0.5 + level/200, type=8, na.rm=TRUE)
#  if(length(level)>1)
#  {
#    lower <- t(lower)
#    upper <- t(upper)
#  }
  return(list(mu=y.f,lower=y.f,upper=y.f))
}

