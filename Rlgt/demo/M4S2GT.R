#testing dual seasonalty method S2GT on hourly subset of the M4 competition set

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
#NUM_OF_CASES=min(length(hourly),NUM_OF_CASES)
NUM_OF_CASES=5

quantileLoss<-function(forec, actual, tau) {
	diff=actual-forec
	pinBallL=pmax(diff*tau,diff*(tau-1))
	mean(pinBallL/actual)*200
}

#i=1
sumSMAPE=0; sumQ99Loss=0; sumQ95Loss=0; sumQ5Loss=0;
numOfCases95pExceeded=0; numOfCases5pExceeded=0;
for (i in 1:NUM_OF_CASES) {
	series=hourly[[i]]$st
	data.train = as.vector(hourly[[i]]$x)  #class(hourly[[i]]$x) -> ts
	MAX_SERIES_LEN_IN_SEAS=7  #we can truncate long series, say longer than 5-7 largest seasonalities
	if (length(data.train)>MAX_SERIES_LEN_IN_SEAS*SEASONALITY2) {
		data.train=data.train[(data.train-MAX_SERIES_LEN_IN_SEAS*SEASONALITY2+1):length(data.train)]
	}
	data.test = as.vector(hourly[[i]]$xx)
	
	rstanmodel <- rlgt(data.train, model="S2GT", nCores=4, nChains=4,
		control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=500, #longer time series, say several hundred points-long, require smaller number of iterations
		SEASONALITY=SEASONALITY, SEASONALITY2=SEASONALITY2 ), #dual seasonality can't be inferred here
		verbose=TRUE)
	forec= forecast(rstanmodel, h = H, level=c(90,98))
	
	sMAPE=mean(abs(forec$median-data.test)/(forec$median+data.test))*200
	sumSMAPE=sumSMAPE+sMAPE
	
	numOfCases95pExceeded=numOfCases95pExceeded+sum(data.test>forec$upper[,1])
	numOfCases5pExceeded=numOfCases5pExceeded+sum(data.test<forec$lower[,1])
	
	q95Loss=quantileLoss(forec$upper[,1], data.test, 0.95)
	sumQ95Loss=sumQ95Loss+q95Loss
	q99Loss=quantileLoss(forec$upper[,2], data.test, 0.99)
	sumQ99Loss=sumQ99Loss+q99Loss
	q5Loss=quantileLoss(forec$lower[,1], data.test, 0.05)
	sumQ5Loss=sumQ5Loss+q5Loss
	print(paste0(series," sMAPE:",signif(sMAPE,3) ,' q5Loss:',signif(q5Loss,3),' q95Loss:',signif(q95Loss,3),' q99Loss:',signif(q99Loss,3) ))
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

