library("Rlgt")
library("Mcomp")
library("forecast")

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
  } else if(type == "lgt_daniel") {
    train = list(as.numeric(ts$x))
    fut = list(as.numeric(ts$xx))
    rv = Rlgt::blgt.multi.forecast(train, fut, n.samples=1e4, parallel=F)
    forecast = rv$forecast[[1]]$yf.med
  } else if(type %in% c("lgt", "nostudent", "nohet", "noglobal", "ets")) {
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


start_time = Sys.time()

for(i in 1:length(data)) {
# for(i in 201:220) {
# for(i in 1:32) {
# for(i in 1:3) {
  print(paste("************ Timeseries:", i, "/", length(data)))
  start_time_ = proc.time()
  ts = data[[i]]
  h = length(ts$xx)
  
  result = forecastAndMeasure(ts, type="forecast_ets")
  results =  appendToResults(results, result, "forecast_ets")
  
  result = forecastAndMeasure(ts, type="forecast_aan")
  results =  appendToResults(results, result, "forecast_aan")
  
  result = forecastAndMeasure(ts, type="forecast_lgt")
  results =  appendToResults(results, result, "forecast_lgt")
  
  result = forecastAndMeasure(ts, type="lgt")
  results =  appendToResults(results, result, "lgt")

  result = forecastAndMeasure(ts, type="lgt_daniel")
  results =  appendToResults(results, result, "lgt_daniel")
  
  result = forecastAndMeasure(ts, type="nostudent")
  results =  appendToResults(results, result, "lgt_no_student")

  result = forecastAndMeasure(ts, type="nohet")
  results =  appendToResults(results, result, "nohet")

  result = forecastAndMeasure(ts, type="noglobal")
  results =  appendToResults(results, result, "noglobal")

  result = forecastAndMeasure(ts, type="ets")
  results =  appendToResults(results, result, "bayesian_ets")
}
end_time = Sys.time()
print(end_time - start_time)

results
# saveRDS(results, file = "results_M3_yearly_with_daniels_sampler.rds")
# results = readRDS(file = "results_M3_yearly_with_daniels_sampler.rds")

df = NULL
# for(name in c("sMAPE", "MASE", "MAPE", "TIME")) {
for(name in c("sMAPE", "MASE", "TIME")) {
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
ren = c("forecast_ets"="ETS/ZZZ (forecast pkg)",
        "forecast_aan"="ETS/AAN (forecast pkg)",
        "forecast_lgt"="LGT (forecast pkg)",
        "lgt"="LGT (RLGT pkg, STAN sampler)",
        "lgt_daniel"="LGT (RLGT pkg, new sampler)",
        "lgt_no_student"="LGT w/o student dist",
        "nohet"="LGT w/o heteroscedasticity",
        "noglobal"="LGT w/o global trend",
        "bayesian_ets"="ETS/AAN (bayesian)")

df = data.frame(t(plyr::rename(df, ren)))

print(df)

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
xdf = xtable(df)
xdf

# With new sampler and bugs fixed
# sMAPE     MASE         TIME
# ETS/ZZZ (forecast pkg)       17.24400 2.839472  0.012589147
# ETS/AAN (forecast pkg)       19.06704 3.057809  0.005596899
# LGT (forecast pkg)           17.10807 2.908374  0.005441860
# LGT (RLGT pkg, STAN sampler) 15.17822 2.488689 41.072992248
# LGT (RLGT pkg, new sampler)  14.99662 2.534608  4.506682171
# LGT w/o student dist         15.33969 2.500897 27.891162791
# LGT w/o heteroscedasticity   15.06496 2.581932 35.890542636
# LGT w/o global trend         16.21114 2.715059 23.522992248
# ETS/AAN (bayesian)           15.93463 2.815171 16.958682171

# Old results
# sMAPE     MASE     MAPE         TIME
# ETS/ZZZ (forecast pkg)     17.24400 2.839472 20.94054  0.012325581
# ETS/AAN (forecast pkg)     19.06704 3.057809 23.96589  0.005534884
# LGT (forecast pkg)         17.10807 2.908374 23.18540  0.006620155
# LGT (RLGT pkg)             15.18116 2.492097 19.65700 40.297472868
# LGT w/o student dist       15.39737 2.505321 19.75991 28.451007752
# LGT w/o heteroscedasticity 15.03690 2.546494 19.83144 36.632279070
# LGT w/o global trend       16.19220 2.708261 19.66141 25.280387597
# ETS/AAN (bayesian)         19.03227 3.504452 21.83422 17.736403101