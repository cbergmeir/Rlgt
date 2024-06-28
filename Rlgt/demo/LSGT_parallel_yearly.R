## Build data and test with yearly series
library(Mcomp)
library(Rlgt)
M3.data <- subset(M3,"yearly")

train.data = list()
future.data = list()

for (i in 1:645) {
  train.data[[i]] = as.numeric(M3.data[[i]]$x)
  future.data[[i]] = as.numeric(M3.data[[i]]$xx) 
}
## Test -- change below to test more series
w.series = 1:20
# w.series = 1:645        # uncomment to test all series

# run in parallel by default
s = system.time({rv=blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4, homoscedastic = T)})

s                         # overall timing info
s[[3]] / length(w.series) # per series time

mean(rv$sMAPE)            # performance in terms of mean sMAPE
mean(rv$InCI)/6           # coverage of prediction intervals -- should be close to 95%