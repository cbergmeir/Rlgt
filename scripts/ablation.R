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


start_time = Sys.time()

for(i in 1:length(data)) {
  # for(i in 201:220) {
  # for(i in 1:32) {
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
saveRDS(results, file = "results_M3_yearly_with_bug_fixed.rds")


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
ren = c("forecast_ets"="ETS/ZZZ (forecast pkg)",
        "forecast_aan"="ETS/AAN (forecast pkg)",
        "forecast_lgt"="LGT (forecast pkg)",
        "lgt"="LGT (RLGT pkg)",
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



# print(c(mean(results$RMSE$forecast_ets), mean(results$RMSE$forecast_lgt), mean(results$RMSE$bayesian_ets)))
# print(c(mean(results$MAE$forecast_ets), mean(results$MAE$forecast_lgt), mean(results$MAE$bayesian_ets)))
# print(c(mean(results$SMAPE$forecast_ets), mean(results$SMAPE$forecast_lgt), mean(results$SMAPE$bayesian_ets)))

# All ts ets vs lgt-optim_c vs ets-bayesian 
# [1] 1180.765 1334.218 1336.680
# [1] 1023.947 1161.244 1159.948
# [1] 17.24400 17.10807 19.06207


# 20 ts (1:20), MAX_NUM_OF_REPEATS=1
# [1] 1532.047 1834.627 1468.438 1276.862 1378.471 1471.407
# [1] 1303.563 1520.330 1244.418 1078.893 1167.447 1239.461
# [1] 19.48781 21.14873 18.25671 15.79901 17.14307 17.62256

# 20 ts (21:40), MAX_NUM_OF_REPEATS=1
# [1] 1795.730 2253.457 1570.703 1745.291 1494.157 1541.634
# [1] 1510.224 1944.384 1305.791 1449.456 1238.000 1278.945
# [1] 22.10774 38.02309 20.08168 25.20936 19.58681 20.11632

# 20 ts (221:240), MAX_NUM_OF_REPEATS=1
# [1] 684.5877 944.2440 674.7359 678.6479 646.6793 672.1985
# [1] 609.6987 842.6257 600.2941 599.4939 570.3095 597.7290
# [1] 13.59033 22.19283 13.85329 13.21668 13.03935 13.60897

# 20 ts (201:220), MAX_NUM_OF_REPEATS=1
# [1] 1095.0732 1429.4380 1031.1539 1047.4615 1103.3432  998.4857
# [1]  932.3666 1267.3319  869.4085  880.9428  921.5788  841.1711
# [1] 17.03815 36.22884 17.49029 16.91098 18.23641 16.88800

# All time series
# [1] 1180.765 1449.000 1076.315 1103.771 1105.007 1062.305
# [1] 1023.9474 1273.0357  924.6725  951.3847  951.480  912.3506
# [1] 17.24400 25.64161 15.48017 16.24636 16.04185 15.20522
