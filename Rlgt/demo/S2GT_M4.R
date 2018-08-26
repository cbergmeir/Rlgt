# Testing dual seasonalty method S2GT on hourly subset of the M4 Forecasting Competition set
# We are showing two ways of passing data as either vector of numbers or msts object

options(width=180)

library(Rlgt)
# install.packages("devtools")
# devtools::install_github("carlanetto/M4comp2018")
library(M4comp2018)

H=48
SEASONALITY=24
SEASONALITY2=168

#names(M4[[1]])
hourly=Filter(function(l) l$period == "Hourly", M4)
hourly=sample(hourly) #shuffle

NUM_OF_CASES=min(length(hourly),NUM_OF_CASES)  # If you let it run its full course (comment out next line), you should see a good result :-)
NUM_OF_CASES=3

quantileLoss<-function(forec, actual, tau) {
	diff=actual-forec
	pinBallL=pmax(diff*tau,diff*(tau-1))
	mean(pinBallL/actual)*200
}

#i=1
sumSMAPE=0; sumQ99Loss=0; sumQ95Loss=0; sumQ5Loss=0;
numOfCases95pExceeded=0; numOfCases5pExceeded=0;
for (i in 1:NUM_OF_CASES) {
	seriesName=hourly[[i]]$st
	
	if (i%%3==0) {  #just for demo and testing. In your code stick to one of the alternatives
		trainData = as.numeric(hourly[[i]]$x) #"naked" vector, so both seasonalities need to be specified in control
		actuals = as.numeric(hourly[[i]]$xx) # actuals have to be matching trainData; both are of numeric class
		rstanmodel <- rlgt(trainData, model="S2GT", nCores=4, nChains=4,
			control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000, #longer time series, say several hundred points-long, require smaller number of iterations
			SEASONALITY=SEASONALITY, SEASONALITY2=SEASONALITY2, MAX_TREE_DEPTH = 12 ), 
			verbose=TRUE)	
	}	else if (i%%3==1) {
		trainData = hourly[[i]]$x # trainData is of ts class, so the SEASONALITY will be extracted from it. SEASONALITY2 has to be specified,  
		actuals = hourly[[i]]$xx  # class of actuals has to be the same as one of trainData; both are of ts class
		rstanmodel <- rlgt(trainData, model="S2GT", nCores=4, nChains=4,
			control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000, 
			SEASONALITY2=SEASONALITY2, MAX_TREE_DEPTH = 12 ), 
			verbose=TRUE)			
	}	else {
		#we convert input and actuals from ts to msts class
		trainData=msts(hourly[[i]]$x, seasonal.periods=c(SEASONALITY,SEASONALITY2), ts.frequency =SEASONALITY, start=start(hourly[[i]]$x))
		actuals=msts(hourly[[i]]$xx, seasonal.periods=c(SEASONALITY,SEASONALITY2), ts.frequency =SEASONALITY, start=start(hourly[[i]]$xx))
		rstanmodel <- rlgt(trainData, model="S2GT", nCores=4, nChains=4,
			control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000, 
			MAX_TREE_DEPTH = 12 ),
			verbose=TRUE)
	}
	
	forec= forecast(rstanmodel, h = H, level=c(90,98))
	# str(forec, max.level=1)
	plot(forec, main=seriesName)
	
	sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
	sumSMAPE=sumSMAPE+sMAPE
	
	numOfCases95pExceeded=numOfCases95pExceeded+sum(actuals>forec$upper[,1])
	numOfCases5pExceeded=numOfCases5pExceeded+sum(actuals<forec$lower[,1])
	
	q95Loss=quantileLoss(forec$upper[,1], actuals, 0.95)
	sumQ95Loss=sumQ95Loss+q95Loss
	q99Loss=quantileLoss(forec$upper[,2], actuals, 0.99)
	sumQ99Loss=sumQ99Loss+q99Loss
	q5Loss=quantileLoss(forec$lower[,1], actuals, 0.05)
	sumQ5Loss=sumQ5Loss+q5Loss
	print(paste0(seriesName," sMAPE:",signif(sMAPE,3) ,' q5Loss:',signif(q5Loss,3),' q95Loss:',signif(q95Loss,3),' q99Loss:',signif(q99Loss,3) ))
}
sMAPE=sumSMAPE/NUM_OF_CASES
q95Loss=sumQ95Loss/NUM_OF_CASES
q99Loss=sumQ99Loss/NUM_OF_CASES
q5Loss=sumQ5Loss/NUM_OF_CASES
exceed95=numOfCases95pExceeded/(NUM_OF_CASES*H)*100
exceed5=numOfCases5pExceeded/(NUM_OF_CASES*H)*100
print(paste0("SUMMARY: Num of cases:", NUM_OF_CASES, ", sMAPE:",signif(sMAPE,4),
  ', % of time 95p exceeded:',signif(exceed95,4), ', % of time 5p exceeded:',signif(exceed5,4), 
	', q5Loss:',signif(q5Loss,4),', q95Loss:',signif(q95Loss,4),', q99Loss:',signif(q99Loss,4) ))

