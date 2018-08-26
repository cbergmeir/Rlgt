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

NUM_OF_CASES=5
NUM_OF_CASES=min(length(hourly),NUM_OF_CASES)


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
	actuals = as.vector(hourly[[i]]$xx)
	
	if (i%%4==0) {
		trainData = as.vector(hourly[[i]]$x) #"naked" vector, so both seasonalities need to be specified in control
		rstanmodel <- rlgt(trainData, model="S2GT", nCores=4, nChains=4,
			control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000, #longer time series, say several hundred points-long, require smaller number of iterations
			SEASONALITY=SEASONALITY, SEASONALITY2=SEASONALITY2, MAX_TREE_DEPTH = 12 ), 
			verbose=TRUE)	
	}	else if (i%%4==1) {
		trainData = hourly[[i]]$x # Although data is of class ts, and it "knows" its single seasonality, SEASONALITY specified in control has priority. SEASONALITY2 has to be specified.   
		rstanmodel <- rlgt(trainData, model="S2GT", nCores=4, nChains=4,
			control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000, 
			SEASONALITY=SEASONALITY, SEASONALITY2=SEASONALITY2, MAX_TREE_DEPTH = 12 ), 
			verbose=TRUE)	
	}	else if (i%%4==2) {
		trainData = hourly[[i]]$x # Here class(trainData)=="ts". If SEASONALITY not specified, we extract it as frequency(trainData). SEASONALITY2 has to be specified.   
		rstanmodel <- rlgt(trainData, model="S2GT", nCores=4, nChains=4,
			control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000, 
			SEASONALITY2=SEASONALITY2, MAX_TREE_DEPTH = 12 ), 
			verbose=TRUE)	
  } else if (i%%4==3) {
		trainData=msts(hourly[[i]]$x, seasonal.periods=c(SEASONALITY,SEASONALITY2))   #both seasonalities will be extracted from trainData 
		rstanmodel <- rlgt(trainData, model="S2GT", nCores=4, nChains=4,
			control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000, 
			MAX_TREE_DEPTH = 12), 
				verbose=TRUE)	
	}	else {
		trainData=msts(hourly[[i]]$x, seasonal.periods=c(SEASONALITY,SEASONALITY2)) # Although data is of class msts, seasonalities specified in rlgt.control have priority (in this example specified unnecessarily)     
		rstanmodel <- rlgt(trainData, model="S2GT", nCores=4, nChains=4,
			control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000, 
			SEASONALITY=SEASONALITY, SEASONALITY2=SEASONALITY2, MAX_TREE_DEPTH = 12 ),
			verbose=TRUE)
	}
	
	forec= forecast(rstanmodel, h = H, level=c(90,98))
	
	sMAPE=mean(abs(forec$median-actuals)/(forec$median+actuals))*200
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

