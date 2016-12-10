require(Mcomp)
require(forecast)

#options(error=recover)
#
#M3.data <- append(subset(M3,"yearly"), append(subset(M3,"quarterly"), subset(M3,"monthly")))

max_forecast_horizon <- max(unlist(lapply(M3, function(x) {length(x$xx)})))
#nseries <- length(M3.data)


complete_forecasts <- list()

appendList <- function (x, val) 
{
  stopifnot(is.list(x), is.list(val))
  xnames <- if(is.null(names(x))) 1:length(x) else names(x)
  
  
  for (v in if(is.null(names(val))) 1:length(val) else names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
          appendList(x[[v]], val[[v]])
        else {
          #c(x[[v]], val[[v]])
          #browser()
          res <- if(is.vector(x[[v]])) c(x[[v]], rep(NA, max_forecast_horizon - length(x[[v]])))
              else x[[v]]
          rbind(res, 
              c(val[[v]], rep(NA, max_forecast_horizon - length(val[[v]]))))
        }
  }
  x
}

complete_forecasts <- appendList(forecasts[[1]],forecasts[[2]])

for(i in 3:length(forecasts)) {
  
  print(paste0(i, " : ", length(complete_forecasts)))
  
  complete_forecasts <- appendList(complete_forecasts,forecasts[[i]])

}

complete_forecasts <- append(complete_forecasts, lapply(M3Forecast, function(x) {as.matrix(x)}))

complete_forecasts[["AAM1"]] <- NULL
complete_forecasts[["AAM2"]] <- NULL

scale <- c(unlist(lapply(M3.data,function(x){
              
              l <- 1
              if(x$period == "QUARTERLY") l <- 4
              if(x$period == "MONTHLY") l <- 12
              
              mean(abs(diff(x$x,l)))
            })))

calc_error_measures <- function(res_method) {
  
  if(is.data.frame(res_method) || !is.list(res_method)) {
    
    res <- NULL
    
    #646
    for(x in 1:nseries) {
      #if(nrow(res_method))
      curr <- accuracy(res_method[x,], M3.data[[x]]$xx)
      mase <- curr[,"MAE"]/scale[x]
      
      smape <- mean(200 * abs(res_method[x,1:length(M3.data[[x]]$xx)] - M3.data[[x]]$xx) / (abs(res_method[x,1:length(M3.data[[x]]$xx)]) + abs(M3.data[[x]]$xx)))
      
      res <- rbind(res, cbind(curr, MASE=mase, sMAPE=smape))
    }
    
    #browser()
    #rownames(res) <- NULL #rownames(res_method)
    res  
    
  }
  else lapply(res_method, calc_error_measures)
}

error_measures <- calc_error_measures(complete_forecasts)

flat_error_measures <- list()

flattenErrorMeasures <- function(name, x) {
  
  if(is.data.frame(x) || !is.list(x)) {
    flat_error_measures[[name]] <<- x
  } else {
    
    for(n in names(x)) {
      
      flattenErrorMeasures(if(name=="") n else paste0(name,".",n), x[[n]])
    }
  }
  
}

flattenErrorMeasures("",error_measures)

names(flat_error_measures)

errors <- list()
for (measureName in colnames(flat_error_measures[[1]])) {
  errors[[measureName]] <- as.data.frame(lapply(flat_error_measures, function(x) { x[,measureName] }))   
}

errors_filtered <- errors


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#load("/home/bergmeir/ownCloud/estancia_Australia/experiments/bagging_M3_20141118_161100_eval_errors.RData")

#library(forecast)
#library(Mcomp)
#
#M3.data <- append(subset(M3,"yearly"), append(subset(M3,"quarterly"), subset(M3,"monthly")))
#nseries <- length(M3.data)


#colMeans(errors_filtered[["sMAPE"]])

ranks.mase <- t(apply(errors_filtered[["MASE"]],1,function(x) {rank(unlist(x))}))
ranks.smape <- t(apply(errors_filtered[["sMAPE"]],1,function(x) {rank(unlist(x))}))

ind <- list()

ind[["yearly"]] <- which(unlist(lapply(M3.data, function(x) {x$period == "YEARLY"})))
ind[["quarterly"]] <- which(unlist(lapply(M3.data, function(x) {x$period == "QUARTERLY"})))
ind[["monthly"]] <- which(unlist(lapply(M3.data, function(x) {x$period == "MONTHLY"})))

rtable <- lapply(ind, function (c_ind) {
      cbind(rank.mase=colMeans(ranks.mase[c_ind,]), mean.mase=colMeans(errors_filtered[["MASE"]][c_ind,]),
          rank.smape=colMeans(ranks.smape[c_ind,]), mean.smape=colMeans(errors_filtered[["sMAPE"]][c_ind,]))
    })



lapply(ind, function (c_ind) {
      cbind(colMeans(errors_filtered[["RMSE"]][c_ind,]))
    })

print(rtable)

rtable$yearly[order(rtable$yearly[,"rank.mase"]),]

lapply(rtable, function(x) {
      x[order(x[,"rank.mase"]),]
    })


