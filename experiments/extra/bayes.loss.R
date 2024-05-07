library(Mcomp)

blgt.sMAPE <- function(yp, yt)
{
  mean( abs(yp - yt)/(yp + yt) )*200
}

blgt.Bayes.Loss <- function(z.hat, rv) {
  sum(apply(rv, 1, function(s){ blgt.sMAPE(s, z.hat)}))
}

print(Sys.time())

for (category in c("yearly", "monthly", "quarterly", "other")) {
  # yf <- readRDS("./experiments/extra/m3.other.homo.rds")
  yf <- readRDS(paste0("/data/sherilyn_bak/M3nu2/m3.", category, ".homo.rds"))
  yf <- yf$forecast
  # yf <- yf[1:2]
  # yf[[1]]$yf <- yf[[1]]$yf[1:2,]
  # yf[[2]]$yf <- yf[[2]]$yf[1:2,]
  H <- ncol(yf[[1]]$yf)
  
  bayes.forecasts <- lapply(yf, function(series){
    # start from the median
    z.init <- series$yf.med
    res <- optim(z.init, blgt.Bayes.Loss, rv = series$yf)
    res$par
  })
  
  M3.data <- subset(M3,category) 
  
  sMAPE <- lapply(seq_along(bayes.forecasts), function(i){
    blgt.sMAPE(bayes.forecasts[[i]], M3.data[[i]]$xx)
  })
  
  sMAPE <- unlist(sMAPE)
  
  print(paste(category, "series, sMAPE = ", mean(sMAPE)))
}

print(Sys.time())

