library(Rlgt)
#set.seed(12)
options(width=180)

data("iclaims.example")

# predict initial unemployment claims of US based on google search data
curr_series <- "iclaims"
data.train  <- iclaims.example$claims
sizeTestSet <- 10

mod <- list()
forecasts <- list()
#--------------------------------
#Fit LGT model
debug(rlgt)
mod[["LGT"]] <- rlgt(data.train, model.type ="LGT", nCores = 4, 
                     nChains = 4, control = lgt.control(MAX_NUM_OF_REPEATS = 10, 
                                                        NUM_OF_ITER = 2000), 
                     verbose=TRUE)
# print the model details
print(mod[["LGT"]])

# With Regressions
x_mat <- as.matrix(iclaims.example[, c('trend.unemploy', 'trend.filling',
                                       'trend.job')])
# debug(rlgt)
mod[["LGT_Reg"]] <- rlgt(data.train, model.type ="LGT", nCores = 4, 
                         xreg = x_mat,
                         nChains = 4, control = lgt.control(
                           MAX_NUM_OF_REPEATS = 10, 
                           NUM_OF_ITER = 2000), 
                         verbose=TRUE)
# print the model details
print(mod[["LGT_Reg"]])

# print the interval for all vars
posterior_interval(mod[["LGT_Reg"]])
debug(forecast)
forecasts[["LGT_Reg"]] <- forecast(mod[["LGT_Reg"]], h = sizeTestSet/2, 
                               level=c(80, 95, 98))
plot(forecasts[["LGT_Reg"]], main=paste(curr_series,'by LGT_Reg'))

data.train  <- ts(iclaims.example$claims, start = 1, frequency = 52)
sizeTestSet <- frequency(data.train)
#--------------------------------
#Fit SGT model
mod[["SGT"]] <- rlgt(data.train, model.type="SGT", nCores=4, nChains=4,
                     control=lgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000), 
                     verbose=TRUE)
# print the model details
print(mod[["SGT"]])

# print the interval for all vars
posterior_interval(mod[["SGT"]])
forecasts[["SGT"]] <- forecast(mod[["SGT"]], h = sizeTestSet, 
                               level=c(80, 95, 98))
plot(forecasts[["SGT"]],main=paste(curr_series,'by SGT'))


