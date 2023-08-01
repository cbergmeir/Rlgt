data <- read.csv("C:/Users/alexa/OneDrive/Desktop/RA/data/tesla-google-trends.csv")

dd = ts(data[73:nrow(data),2], start=c(2010,1), frequency=12)

plot(dd, type="l", ylab="", xlab="")

d = data[73:nrow(data),2]
plot(d)

library("Rlgt")


rlgt_model <- rlgt(d,
   control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=3000, NUM_OF_CHAINS=4, NUM_OF_CORES=4),
   verbose=T)

debugonce(Rlgt:::forecast.rlgtfit)
forec = Rlgt:::forecast.rlgtfit(rlgt_model, h=48, NUM_OF_TRIALS = 5000, level=90)
plot(forec)

plot(dd, type="l", ylab="", xlab="", ylim=c(0,160), xlim=c(2010,2025.1))
for(i in 1:10) {
  ddd = ts(forec$yf[i,], start = c(2021,7), frequency = 12)
  lines(ddd, col='lightgrey')
}
ddd = ts(forec$mean, start = c(2021,7), frequency = 12)
ddd_l = ts(as.vector(forec$lower), start = c(2021,7), frequency = 12)
ddd_u = ts(as.vector(forec$upper), start = c(2021,7), frequency = 12)
lines(ddd, col='black', lwd=2)
lines(ddd_l, col='black', lwd=0.2)
lines(ddd_u, col='black', lwd=0.2)

################################################################################

# library(forecast)
# 
# model = ets(d, model="AAN", lgt=FALSE)
# f = forecast(model, h=48)
# plot(f)
