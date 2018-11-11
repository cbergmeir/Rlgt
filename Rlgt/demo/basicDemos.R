# TODO: Add comment
# 
# Author: slawek
###############################################################################

library(forecast)
library(Rlgt)

options(width=180)
if (.Platform$OS.type=="windows")  memory.limit(6000)


################ LGT
theDataSet=BJsales 

tsdisplay(theDataSet)
frequency(theDataSet)

horizon=10
train=theDataSet[1:(length(theDataSet)-horizon)]
actuals=theDataSet[(length(theDataSet)+1-horizon):length(theDataSet)]

model <- rlgt(train, 
		control=rlgt.control(NUM_OF_ITER=5000,MAX_NUM_OF_REPEATS=1))   

forec= forecast(model, h = horizon)

plot(forec, main="BJsales by LGT")
xs=seq(from=length(train)+1,to=length(train)+ length(actuals))
lines(xs,actuals, col=1, type='b',lwd=2)	

sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
print(paste("sMAPE:",signif(sMAPE,3)))


################ LGT with regression
regDataSet=BJsales.lead 

regTrain=regDataSet[1:(length(theDataSet)-horizon)]
regTest=regDataSet[(length(theDataSet)+1-horizon):length(theDataSet)]

regModel <- rlgt(train, xreg = regTrain, 
		control=rlgt.control(NUM_OF_ITER=5000,MAX_NUM_OF_REPEATS=1), 
		verbose=TRUE)

forec= forecast(regModel, regTest, h = horizon)

plot(forec, main="BJsales with lead regressor by LGT")
xs=seq(from=length(train)+1,to=length(train)+ length(actuals))
lines(xs,actuals, col=1, type='b',lwd=2)	

sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
print(paste("sMAPE:",signif(sMAPE,3)))


################### SGT
theDataSet=AirPassengers

tsdisplay(theDataSet)
frequency(theDataSet)

horizon=2*frequency(theDataSet)
train=theDataSet[1:(length(theDataSet)-horizon)]
actuals=theDataSet[(length(theDataSet)+1-horizon):length(theDataSet)]

rstanmodel <- rlgt(train, seasonality=frequency(theDataSet), 
		seasonality.type="generalized", level.method="seasAvg",
		control=rlgt.control(NUM_OF_ITER=10000, MAX_NUM_OF_REPEATS=1))   

forec= forecast(rstanmodel, h = length(actuals))

plot(forec, main="AirPassengers by SGT")
xs=seq(from=length(train)+1,to=length(train)+ length(actuals))
lines(xs,actuals, col=1, type='b',lwd=2)	

sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
print(paste("sMAPE:",signif(sMAPE,3)))


#########################