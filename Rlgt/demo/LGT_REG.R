library(Rlgt)
library(rstan)
#set.seed(12)
options(width=180)
data("umcsent.example")

# Data setup --------------------------------
# predict initial unemployment claims of US based on google search data
curr_series <- "consumer.sent"
umcsent.example2 <- umcsent.example %>%
  mutate_at(.vars = c("consumer.sent", "search.engine",
                      "financial.planning", "bus.news",      
                      "investing", "energy.utilities"),
            .funs = log)
sizeTestSet <- 24
y  <- umcsent.example2$consumer.sent
y.train <- y[1:(length(y) - sizeTestSet)]

# Regression Matrix
x.mat <- as.matrix(
  umcsent.example2[, c("search.engine",
                       "financial.planning", "bus.news",      
                       "investing", "energy.utilities")])
x.mat.train <- x.mat[1:(nrow(x.mat) - sizeTestSet),]
x.mat.test <- x.mat[(nrow(x.mat) - sizeTestSet + 1):nrow(x.mat),]
mod <- list()
forecasts <- list()
#--------------------------------
#Fit LGT model without Regression 
mod[["LGT"]] <- rlgt(y.train, model.type="LGT", nCores=4, nChains=4,
                     control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=1000), 
                     verbose=TRUE)
# print the model details
print(mod[["LGT"]])

# print the interval for all vars
posterior_interval(mod[["LGT"]])
forecasts[["LGT"]] <- forecast(mod[["LGT"]],  
                               h = sizeTestSet,
                               level=c(80, 95, 98))
plot(forecasts[["LGT"]], main=paste(curr_series,'by LGT'))
lines(y, col = 'red')

#--------------------------------
#Fit LGT model with Regression 
mod[["LGT_REG"]] <- rlgt(y.train, model.type="LGT", nCores=4, nChains=4,
                         xreg = x.mat.train,
                         control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=1000), 
                         verbose=TRUE)
# print the model details
print(mod[["LGT_REG"]])

# print the interval for all vars
posterior_interval(mod[["LGT_REG"]])
forecasts[["LGT_REG"]] <- forecast(mod[["LGT_REG"]],                                   
                                   x.mat.test,
                                   level=c(80, 95, 98))
plot(forecasts[["LGT_REG"]], main=paste(curr_series,'by LGT with Regression'))
lines(y, col = 'red')