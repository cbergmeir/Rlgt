rm(list=ls())
library(Rlgt)
library(dplyr)
options(width=180)
data("umcsent.example")

mod <- list()
forecasts <- list()
#--------------------------------
#Fit LGT Regression model 
curr_series <- "consumer.sent"
sizeTestSet <- 24
y  <- umcsent.example$consumer.sent
train_index <- 1:(length(y) - sizeTestSet)
test_index  <- (1:length(y))[-train_index]
y.train <- y[train_index]
y.test  <- y[test_index]

# Regression Matrix
x.mat <- as.matrix(umcsent.example[, c("search.engine",
                                       "financial.planning", "bus.news",      
                                       "investing", "energy.utilities")])

x.mat.train <- x.mat[train_index,]
x.mat.test <- x.mat[test_index,]
# Without Regression 
mod[["LGT"]] <- rlgt(y.train, 
                     control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=5000), 
                     verbose=TRUE)
# With Regression 
mod[["LGT_REG"]] <- rlgt(y.train, 
                         xreg = x.mat.train,
												 control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=10000), 
                         verbose=TRUE)

# Test 1 regressor
mod[["LGT_REG1"]] <- rlgt(y.train, 
                         xreg = x.mat.train[,1],
                         control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=8000), 
                         verbose=TRUE)


forecasts[["LGT"]] <- forecast(mod[["LGT"]],
                               h = sizeTestSet,
                               level=c(80, 95, 98))
plot(forecasts[["LGT"]], main=paste(curr_series,'by LGT'))
points(test_index, y.test, col=1, type='b',lwd=2)

# print the model details
print(mod[["LGT_REG"]])
# print the interval for all vars
posterior_interval(mod[["LGT_REG"]])
forecasts[["LGT_REG"]] <- forecast(mod[["LGT_REG"]],                                   
                                   x.mat.test,
																	 h = sizeTestSet,
                                   level=c(80, 95, 98))
plot(forecasts[["LGT_REG"]], main=paste(curr_series,'by LGT with Regression'))
points(y, col = 'red')

#--------------------------------
# Fit SGT Regression model 
data("iclaims.example")
curr_series <- "iclaims"
y  <- iclaims.example$claims
acf(y, lag.max=55)  #this is weekly data
length(y)
SEASONALITY=365.25/7
sizeTestSet <- frequency(y)
train_index <- 1:(length(y) - sizeTestSet)
test_index  <- (1:length(y))[-train_index]
y.train <- y[train_index]
y.test  <- y[test_index]
# y.train  <- ts(y.train, start = 1, frequency = 52)

# Regression Matrix
x.mat <- as.matrix(iclaims.example[, c('trend.unemploy', 'trend.filling', 'trend.job')])
x.mat.train <- x.mat[train_index,]
x.mat.test  <- x.mat[test_index,]

# Without Regression,  
mod[["SGT"]] <- rlgt(y.train, 
                     seasonality = SEASONALITY,
                     control=rlgt.control(NUM_OF_ITER=10000), 
                     verbose=TRUE)
# With Regression 
mod[["SGT_REG"]] <- rlgt(y.train, 
                         seasonality = SEASONALITY,
                         xreg = x.mat.train,
                         control=rlgt.control(NUM_OF_ITER=10000),
                         verbose=TRUE)
# print the model details

# print the interval for all vars
forecasts[["SGT"]] <- forecast(mod[["SGT"]],  
                               h = sizeTestSet,
                               level=c(80, 95, 98))
plot(forecasts[["SGT"]], main=paste(curr_series,'by SGT'))
#points(y, col = 'red')

# print the interval for all vars
posterior_interval(mod[["SGT_REG"]])
forecasts[["SGT_REG"]] <- forecast(mod[["SGT_REG"]],
                                   xreg = matrix(x.mat.test,nrow=1,ncol=3),
                                   level=c(80, 95, 98))
plot(forecasts[["SGT_REG"]], main=paste(curr_series,'by SGT with Regression'))
#points(y, col = 'red')

# Run with 1 Regressor
# mod[["SGT_REG"]] <- rlgt(y.train,
#                          seasonality = frequency(y),
#                          xreg = x.mat.train[,1,drop = FALSE],
#                          control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000),
#                          verbose=TRUE)

