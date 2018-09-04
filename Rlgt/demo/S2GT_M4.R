# Testing dual seasonalty method S2GT on hourly subset of the M4 Forecasting Competition set
# We are showing three ways of passing data: as a vector of numbers, ts, or msts object. 
# Finally, generalized seasonality version is used. 
#These are long series and dual seasonality models are more complicated to fit, 
#each step may take from a few minutes to over 1 hour.

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

NUM_OF_CASES=length(hourly) #running over 400 cases would take quite a few days :-)
NUM_OF_CASES=10

quantileLoss<-function(forec, actual, tau) {
	diff=actual-forec
	pinBallL=pmax(diff*tau,diff*(tau-1))
	mean(pinBallL/actual)*200
}

legend_str_vect=NULL; legend_cols_vect=NULL; legend_char_vect=NULL
legend_str_vect=c(legend_str_vect,"forecast")  
legend_cols_vect=c(legend_cols_vect, 'blue')
legend_char_vect=c(legend_char_vect,'-')

legend_str_vect=c(legend_str_vect,"actuals")  #used for short displays
legend_cols_vect=c(legend_cols_vect, 'black')
legend_char_vect=c(legend_char_vect,'-')


i<-1
sumSMAPE=0; sumQ99Loss=0; sumQ95Loss=0; sumQ5Loss=0;
numOfCases95pExceeded=0; numOfCases5pExceeded=0;
for (i in 1:NUM_OF_CASES) {
	seriesName=hourly[[i]]$st
	
	if (i==1) {  #just for demo and testing. In your code stick to one of the alternatives
		trainData = as.numeric(hourly[[i]]$x) #"naked" vector, so both seasonalities need to be specified in control
		actuals = as.numeric(hourly[[i]]$xx) # actuals have to be matching trainData; both are of numeric class
		rstanmodel <- rlgt(trainData, seasonality=SEASONALITY, seasonality2=SEASONALITY2,
			control=rlgt.control(MAX_NUM_OF_REPEATS=2, NUM_OF_ITER=1000, #longer time series, say several hundred points-long, require smaller number of iterations
			MAX_TREE_DEPTH = 12), 
			verbose=TRUE)	
	}	else if (i==2) {
		trainData = hourly[[i]]$x # trainData is of ts class, so the SEASONALITY will be extracted from it. SEASONALITY2 has to be specified,  
		actuals = hourly[[i]]$xx  # class of actuals has to be the same as one of trainData; both are of ts class
		rstanmodel <- rlgt(trainData, seasonality2=SEASONALITY2,
			control=rlgt.control(MAX_NUM_OF_REPEATS=2, NUM_OF_ITER=1000), 
			verbose=TRUE)			
	} else  if (i==3)  {
		trainData = hourly[[i]]$x # trainData is of ts class, so the SEASONALITY will be extracted from it. SEASONALITY2 has to be specified,  
		actuals = hourly[[i]]$xx  # class of actuals has to be the same as one of trainData; both are of ts class
		rstanmodel <- rlgt(trainData, seasonality.type="generalized",  #also try generalized seasonality
				seasonality2=SEASONALITY2,
				control=rlgt.control(MAX_NUM_OF_REPEATS=2, NUM_OF_ITER=1000), 
				verbose=TRUE)					
	} else {
		trainData = hourly[[i]]$x 
		actuals = hourly[[i]]$xx  
		rstanmodel <- rlgt(trainData, seasonality2=SEASONALITY2,
				control=rlgt.control(NUM_OF_ITER=1000), 
				verbose=TRUE)			
		
	}
	
	forec= forecast(rstanmodel, h = H, level=c(90,98))
	# str(forec, max.level=1)
	plot(forec, main=seriesName)
	
	if (inherits(trainData,"ts")) {
		lines(actuals, col=1, lwd=2)	
	} else {
		xs=seq(from=length(trainData)+1,to=length(trainData)+ length(actuals))
		lines(xs,actuals, col=1, type='b',lwd=2)	
	}
	legend("topleft", legend_str_vect,
			pch=legend_char_vect, 
			col=legend_cols_vect, cex=1)
	
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

