blgt.sMAPE <- function(yp, yt)
{
  mean( abs(yp - yt)/(yp + yt) )*200
}

blgt.Bayes.Loss <- function(z.hat, rv) {
  sum(apply(rv, 1, function(s){ blgt.sMAPE(s, z.hat)}))
}

yf <- readRDS("./experiments/extra/m3.other.homo.rds")
yf <- yf$forecast
H <- 8

bayes.forecasts <- lapply(yf, function(series){
  # start from the median
  z.init <- series$yf.med
  res <- optim(z.init, blgt.Bayes.Loss, rv = series$yf)
  res$par
})


