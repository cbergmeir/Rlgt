library(Rlgt)
library(rstan)
#set.seed(12)
options(width=180)
data("iclaims.example")

# Data setup --------------------------------
# predict initial unemployment claims of US based on google search data
curr_series <- "iclaims"
y  <- iclaims.example$claims
y  <- ts(y, start = 1, frequency = 52)
sizeTestSet <- frequency(y)
y.train <- y[1:(length(y) - sizeTestSet)]
y.train  <- ts(y.train, start = 1, frequency = 52)

# Regression Matrix
x.mat <- as.matrix(
  iclaims.example[, c('trend.unemploy', 'trend.filling', 'trend.job')])
x.mat.train <- x.mat[1:(nrow(x.mat) - sizeTestSet),]
x.mat.test <- x.mat[(nrow(x.mat) - sizeTestSet + 1):nrow(x.mat),]
mod <- list()
forecasts <- list()
#--------------------------------
#Fit SGT model without Regression 
mod[["SGT"]] <- rlgt(y.train, model.type="SGT", nCores=4, nChains=4,
                     control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=1000), 
                     verbose=TRUE)
# print the model details
print(mod[["SGT"]])

# print the interval for all vars
posterior_interval(mod[["SGT"]])
forecasts[["SGT"]] <- forecast(mod[["SGT"]],  
                               h = sizeTestSet,
                               level=c(80, 95, 98))
plot(forecasts[["SGT"]], main=paste(curr_series,'by SGT'))
lines(y, col = 'red')

#--------------------------------
#Fit SGT model with Regression 
mod[["SGT_REG"]] <- rlgt(y.train, model.type="SGT", nCores=4, nChains=4,
                         xreg = x.mat.train,
                         control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=1000),
                         verbose=TRUE)
# print the model details
print(mod[["SGT_REG"]])

# print the interval for all vars
posterior_interval(mod[["SGT_REG"]])
forecasts[["SGT_REG"]] <- forecast(mod[["SGT_REG"]],
                                   x.mat.test,
                                   level=c(80, 95, 98))
plot(forecasts[["SGT_REG"]], main=paste(curr_series,'by SGT with Regression'))
lines(y, col = 'red')

