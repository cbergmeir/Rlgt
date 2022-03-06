# The script was used to calculate LGT accuracy with and without the Daniel's bug fixed (Daniel found it)

library("Rlgt")
library("Mcomp")

getMeasures = function(ground, forecasted, time, training) {
  result = list()
  result$MAPE = mean(abs(ground - forecasted)/ground)*100
  result$MASE = mean(abs(ground - forecasted))/mean(abs(training[2:length(training)] - training[1:(length(training)-1)]))
  result$sMAPE = mean(abs(ground - forecasted)/abs(ground + forecasted))*200
  result$TIME = time
  return(result)
}

forecastAndMeasure = function(ts, type="") {
  start_time = proc.time()
  if(type == "forecast_ets") {
    model = ets(ts$x, lgt=FALSE)
    forecast = forecast(model, h=h)$mean
  } else if(type == "forecast_aan") {
    model = ets(ts$x, model="AAN", lgt=FALSE)
    forecast = forecast(model, h=h, simulate=TRUE, npaths=0)$mean
  } else if(type == "forecast_lgt") {
    model = ets(ts$x, lgt=TRUE)
    forecast = forecast(model, h=h, simulate=TRUE, npaths=0)$mean
  } else if(type %in% c("lgt", "nostudent", "nohet", "noglobal", "ets")){
    rlgt_model <- rlgt(ts$x,
                       control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=3000, NUM_OF_CHAINS=4, NUM_OF_CORES=4),
                       verbose=F,
                       experimental=type)
    forecast = forecast(rlgt_model, h=h)$mean
  } else {
    stop("Unknown type...")
  }
  time = (proc.time() - start_time)["elapsed"]
  result = getMeasures(ts$xx, forecast, time, ts$x)
}

appendToResults = function(results, result, s) {
  for(name in names(result)) {
    results[[name]][[s]] = c(results[[name]][[s]], result[[name]])
  }
  return(results)
}

data <- subset(M3,"yearly")

results = list()
results$MAPE = results$MASE = results$sMAPE = results$TIME = list()

for(i in 1:length(data)) {
# for(i in 201:220) {
# for(i in 1:3) {
  print(paste("************ Timeseries:", i, "/", length(data)))
  start_time_ = proc.time()
  ts = data[[i]]
  h = length(ts$xx)
  
  result = forecastAndMeasure(ts, type="lgt")
  results =  appendToResults(results, result, "lgt")
}

results
saveRDS(results, file = "results_bug_fixed_another_way_plus_other_changes.rds")


df = NULL
for(name in c("sMAPE", "MASE", "MAPE", "TIME")) {
  # print(name)
  row = cn = NULL
  for(method in names(results[[name]])) {
    # print(method)
    m = mean(results[[name]][[method]])
    # print(m)
    row = c(row, m)
    cn = c(cn, method)
  }
  names(row) = cn
  df_ = data.frame(row)
  colnames(df_) = name
  df = rbind(df, data.frame(t(df_)))
}
library(plyr)
ren = c(
  "lgt"="LGT (RLGT pkg)"
  )

df = data.frame(t(plyr::rename(df, ren)))

print(df)
# with the bug
# sMAPE     MASE     MAPE     TIME
# LGT (RLGT pkg) 15.21259 2.498124 19.60923 37.93929

# with bug fixed
# sMAPE     MASE     MAPE     TIME
# LGT (RLGT pkg) 15.29431 2.503083 19.75226 37.26857

# with the bug fixed in LGT.stan
# // b[t] = bSm*(l[t]-l[t-1]) + (1-bSm)*b[t-1] ;
# b[t] = bSm*(l[t]-l[t-1]) + (1-bSm)*locTrendFract*b[t-1] ;
# sMAPE     MASE     MAPE     TIME
# LGT (RLGT pkg) 15.15063 2.480465 19.56943 42.94946

# With the change:
# l[t] = levSm*(y[t]-r[t]) + (1-levSm)*(expVal[t]-r[t]) ;
# sMAPE     MASE     MAPE    TIME
# LGT (RLGT pkg) 15.61546 2.564389 20.16579 41.0796


# With the change:
#	l[t] = levSm*(y[t]-r[t]) + (1-levSm)*(l[t-1]+locTrendFract*b[t-1]) ;
# sMAPE    MASE    MAPE     TIME
# LGT (RLGT pkg) 15.55677 2.56444 20.0642 40.72133
