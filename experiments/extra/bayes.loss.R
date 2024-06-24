library(Mcomp)
library(optimx)

blgt.sMAPE <- function(yp, yt)
{
  mean( abs(yp - yt)/(yp + yt) )*200
}

sMAPE.mtx <- function(yf,Y,Y.sum)
{
  yf.sum = sum(yf)
  num = rowSums(abs(sweep(Y,2,yf)))
  sMAPE = sum(num/(yf.sum+Y.sum))
  sMAPE
}

print(Sys.time())

for (category in c("yearly", "quarterly", "other")) {
  for (variant in c("homo", "hetero")) {
    # yf <- readRDS("./experiments/extra/m3.other.homo.rds")
    yf <- readRDS(paste0("/data/sherilyn_bak/M3nu2/m3.", category, ".", variant, ".rds"))
    yf <- yf$forecast
    # yf <- yf[1:2]
    # yf[[1]]$yf <- yf[[1]]$yf[1:2,]
    # yf[[2]]$yf <- yf[[2]]$yf[1:2,]
    
    
    bayes.forecasts <- lapply(yf, function(series){
      # random sample 1e4
      # idx <- sample(1:nrow(series$yf), 1e4)
      # Y <- series$yf[idx,]
      Y <- series$yf
      
      # horizon
      H <- ncol(Y)
  
      Y.sum = rowSums(Y)
      
      rv <- optimx(colMeans(Y),sMAPE.mtx,Y=Y, Y.sum=Y.sum, lower=rep(0,H), method="L-BFGS-B")
      
      res <- rep(0, H)
      for (i in 1:H) {
        res[i] <- rv[[paste0("p", i)]]
      }
      
      res
    })
    
    M3.data <- subset(M3,category) 
    
    sMAPE <- lapply(seq_along(bayes.forecasts), function(i){
      blgt.sMAPE(bayes.forecasts[[i]], M3.data[[i]]$xx)
    })
    
    sMAPE <- unlist(sMAPE)
    
    print(paste(category, variant, "series, sMAPE = ", mean(sMAPE)))
  }
}

print(Sys.time())

