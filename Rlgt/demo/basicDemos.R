# TODO: Add comment
# 
# Author: slawek
###############################################################################

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


################### SGT, time series input and forecast
theDataSet=AirPassengers

tsdisplay(theDataSet)
frequency(theDataSet)
tspOrg <- tsp(theDataSet)

horizon=2*frequency(theDataSet)
train=theDataSet[1:(length(theDataSet)-horizon)]  #forecast the last horizon of the series
#class(train) -> numeric
actuals=theDataSet[(length(theDataSet)+1-horizon):length(theDataSet)]

#making them back time series
train=ts(train, start=tspOrg[1], frequency=tspOrg[3])
tspt=tsp(train)
actuals=ts(actuals, start=tspt[2]+1/tspt[3], frequency=tspt[3])

rstanmodel <- rlgt(train,  
		seasonality.type="generalized", level.method="seasAvg",
		control=rlgt.control(NUM_OF_ITER=10000, MAX_NUM_OF_REPEATS=1))   

forec= forecast(rstanmodel, h = length(actuals))

plot(forec, main="AirPassengers by SGT")
lines(actuals, lwd=2)

sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
print(paste("sMAPE:",signif(sMAPE,3),"%"))


################### SGT, numeric input and forecast
theDataSet=AirPassengers

tsdisplay(theDataSet)
frequency(theDataSet)

horizon=2*frequency(theDataSet)
train=theDataSet[1:(length(theDataSet)-2*horizon)] #forestast second last horizon from the series
actuals=theDataSet[(length(theDataSet)+1-2*horizon):(length(theDataSet)-horizon)]

rstanmodel <- rlgt(train, seasonality=frequency(theDataSet), 
		seasonality.type="generalized", level.method="seasAvg",
		control=rlgt.control(NUM_OF_ITER=10000, MAX_NUM_OF_REPEATS=1))   

forec= forecast(rstanmodel, h = length(actuals))

plot(forec, main="AirPassengers by SGT")
xs=seq(from=length(train)+1,to=length(train)+ length(actuals))
lines(xs,actuals, col=1, type='b',lwd=2)	

sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
print(paste("sMAPE:",signif(sMAPE,3),"%"))


#########################  SGT on dual-seasonality time series, msts input and forecast
theDataSet=taylor

tsdisplay(theDataSet)

#taylor is a dual seasonality time series (48,336), but is treated as a single seasonality series of the larger frequency (336)
seasonality=frequency(theDataSet)  #larger seasonality
horizon=seasonality
train=theDataSet[1:(2*seasonality)]  #using first two weeks
actuals=theDataSet[(2*seasonality+1):(2*seasonality+horizon)] #to forecast the third one
#class(taylor) -> "msts" "ts"  
#class(train) -> 	"msts" "integer"
#making them back time series
train=msts(train, seasonal.periods=attributes(taylor)$msts)
tspx <- tsp(train)
actuals=msts(actuals, seasonal.periods=attributes(taylor)$msts, start=tspx[2] + 1/seasonality)

rstanmodel <- rlgt(train,  #because seaonality2 is not specified, a single seasonality model (of seasonality equal to the largest seasonality, 336)  will be used 
		level.method="seasAvg",
		control=rlgt.control(NUM_OF_ITER=5000, MAX_NUM_OF_REPEATS=1))   

forec= forecast(rstanmodel, h = length(actuals))

plot(forec, main="Taylor by SGT")
lines(actuals, lwd=2)

sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
print(paste("sMAPE:",signif(sMAPE,3),"%"))


#########################  SGT on dual-seasonality time series, numeric input and forecast
theDataSet=taylor

tsdisplay(theDataSet)

#taylor is a dual seasonality time series (48,336), but is treated as a single seasonality series of the larger frequency (336)
seasonality=frequency(theDataSet)  #larger seasonality
horizon=seasonality
train=theDataSet[(seasonality+1):(3*seasonality)]  #using weeks 2 and 3
actuals=theDataSet[(3*seasonality+1):(3*seasonality+horizon)] #to forecast the fourth one
#class(taylor) ->  "msts" "ts"  
#class(train) -> "msts"    "integer"

rstanmodel <- rlgt(train,  #because seaonality2 is not specified, a single seasonality model (of seasonality equal to the largest seasonality, 336)  will be used 
		#level.method="seasAvg",
		control=rlgt.control(NUM_OF_ITER=5000, MAX_NUM_OF_REPEATS=1))   

forec= forecast(rstanmodel, h = length(actuals))

plot(forec, main="Taylor by SGT")
xs=seq(from=length(train)+1,to=length(train)+ length(actuals))
lines(xs,actuals, col=1, type='b',lwd=2)	

sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
print(paste("sMAPE:",signif(sMAPE,3),"%"))
