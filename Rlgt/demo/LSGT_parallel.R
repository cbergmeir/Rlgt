## Build data and test with quarterly series
library(Mcomp)
library(Rlgt)
M3.data <- subset(M3,"quarterly")
 
train.data = list()
future.data = list()
 
for (i in 1:756) {
   train.data[[i]] = as.numeric(M3.data[[i]]$x)
   future.data[[i]] = as.numeric(M3.data[[i]]$xx) 
}
## Test -- change below to test more series
w.series = 1:20
# w.series = 1:756        # uncomment to test all series

# run in parallel by default
s = system.time({rv=blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4, m = 4)})

s                         # overall timing info
s[[3]] / length(w.series) # per series time
 
mean(rv$sMAPE)            # performance in terms of mean sMAPE
mean(rv$InCI)/8           # coverage of prediction intervals -- should be close to 95%

# can also specify not run in parallel
s = system.time({rv=blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4, parallel = F, m = 4)})

