library(Rlgt)
library(dplyr)
library(tidyr)
library(lubridate)
data_raw <- read.csv('../data-raw/kaggle_retail/sales data-set 2.csv', stringsAsFactors = FALSE) 
data_fix1 <- data_raw %>%
  mutate(Date = as.Date(Date, tryFormats = c("%d/%m/%Y", "%Y-%m-%d")))

kaggle.retail.sales <-  subset(data_fix1, Store == 1 & Dept == 1)

plot(x = kaggle.retail.sales$Date, y = kaggle.retail.sales$Weekly_Sales)
lines(x = kaggle.retail.sales$Date, y = kaggle.retail.sales$Weekly_Sales, col = 'red')
points(x = kaggle.retail.sales$Date, y = kaggle.retail.sales$IsHoliday * 50000, col = 'blue')

sizeTestSet <- 12
y <- kaggle.retail.sales$Weekly_Sales
# Regression Matrix
x.mat <- as.matrix(kaggle.retail.sales[, c('IsHoliday')])
x.mat.train <- x.mat[1:(nrow(x.mat) - sizeTestSet),]
x.mat.train[,2] <- 0

x.mat.test <- x.mat[(nrow(x.mat) - sizeTestSet + 1):nrow(x.mat),]

y.train <- y[1:(length(y) - sizeTestSet)]
y.train  <- ts(y.train, start = 1, frequency = 52)


mod[["SGT_REG"]] <- rlgt(y.train, 
                         xreg = x.mat.train,
                         control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=4000),
                         verbose=TRUE)
# print the model details
print(mod[["SGT"]])