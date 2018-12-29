# Testing dual seasonalty method S2GT and SGT on hourly subset of the M4 Forecasting Competition set
# We are showing three ways of passing data: as a vector of numbers, ts, or msts object.
#Then, we are tryin all possible combinations of seasonality types, error and level methods.
#These are long series and dual seasonality models are more complicated to fit, 
#each step may take from a few minutes to over 1 hour (especially slow can be models with the "innov" method for error size)


options(width=180)
if (.Platform$OS.type=="windows")  memory.limit(10000)

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
#NUM_OF_CASES=10

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


i<-1; forecasts=list()
sumSMAPE=0; sumQ99Loss=0; sumQ95Loss=0; sumQ5Loss=0;
numOfCases95pExceeded=0; numOfCases5pExceeded=0;


# ways to pass seasonal data
for (i in 1:6) {
	ii=sample(NUM_OF_CASES,1)
	seriesName=hourly[[ii]]$st
	print(paste("starting",seriesName))
	
	if (i==1) {  #just for demo and testing. In your code stick to one of the alternatives
		trainData = as.numeric(hourly[[ii]]$x) #"naked" vector, so both seasonalities need to be specified in control
		actuals = as.numeric(hourly[[ii]]$xx) # actuals have to be matching trainData; both are of numeric class
		rstanmodel <- rlgt(trainData, seasonality=SEASONALITY, seasonality2=SEASONALITY2, #specifying seasonality2>1 means we are running S2GT 
			control=rlgt.control(MAX_TREE_DEPTH = 12), 
			verbose=TRUE)	
	}	else if (i==2) {
		trainData = hourly[[ii]]$x # trainData is of ts class, so the SEASONALITY will be extracted from it.   
		actuals = hourly[[ii]]$xx  # class of actuals has to be the same as one of trainData; both are of ts class
		rstanmodel <- rlgt(trainData, seasonality2=SEASONALITY2,  #specifying seasonality2>1 means we are running S2GT
			control=rlgt.control(MAX_RHAT_ALLOWED=1.01), 
			verbose=TRUE)			
	}	else if (i==3){
		#we convert input and actuals from ts to msts class. Although msts has both seasonalities
		trainData=msts(hourly[[ii]]$x, seasonal.periods=c(SEASONALITY,SEASONALITY2), ts.frequency =SEASONALITY, start=start(hourly[[ii]]$x))
		actuals=msts(hourly[[ii]]$xx, seasonal.periods=c(SEASONALITY,SEASONALITY2), ts.frequency =SEASONALITY, start=start(hourly[[ii]]$xx))
 		rstanmodel <- rlgt(trainData, seasonality2=SEASONALITY2,  #specifying seasonality2>1 means we are running S2GT
			control=rlgt.control(NUM_OF_ITER=5000),
			verbose=TRUE)
	} else if (i==4) {  #just for demo and testing. In your code stick to one of the alternatives
		trainData = as.numeric(hourly[[ii]]$x) #"naked" vector, so both seasonalities need to be specified in control
		actuals = as.numeric(hourly[[ii]]$xx) # actuals have to be matching trainData; both are of numeric class
		rstanmodel <- rlgt(trainData, seasonality=SEASONALITY, #no seasonality2 specified -> we are running SGT  
				verbose=TRUE)	
	}	else if (i==5) {
		trainData = hourly[[ii]]$x # trainData is of ts class, so the SEASONALITY will be extracted from it.   
		actuals = hourly[[ii]]$xx  # class of actuals has to be the same as one of trainData; both are of ts class
		rstanmodel <- rlgt(trainData, #no seasonality2 specified -> we are running SGT
				control=rlgt.control(ADAPT_DELTA=0.95, MAX_TREE_DEPTH = 12),
				verbose=TRUE)			
	}	else if (i==6){
		#we convert input and actuals from ts to msts class. Although msts has both seasonalities
		trainData=msts(hourly[[ii]]$x, seasonal.periods=c(SEASONALITY,SEASONALITY2), ts.frequency =SEASONALITY, start=start(hourly[[ii]]$x))
		actuals=msts(hourly[[ii]]$xx, seasonal.periods=c(SEASONALITY,SEASONALITY2), ts.frequency =SEASONALITY, start=start(hourly[[ii]]$xx))
		rstanmodel <- rlgt(trainData, #no seasonality2 specified -> we are running SGT
				control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=5000),
				verbose=TRUE)
	}
	
	forec= forecast(rstanmodel, h = H, level=c(90,98))
	forecasts[[seriesName]]<-list(mean=forec$mean, lower=forec$lower, upper=forec$upper)
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
	

#just exercising every possible combination of options
i=7
for (seasonality.type in c("multiplicative","generalized")){
	for (level.method in c("HW", "seasAvg","HW_sAvg")) {
		for (error.size.method in c("std","innov")) {
			ii=sample(NUM_OF_CASES,1)
			seriesName=hourly[[ii]]$st
			print(paste("starting",seriesName))
			ii=sample(NUM_OF_CASES,1)
			trainData = hourly[[ii]]$x 
			actuals = hourly[[ii]]$xx  
			rstanmodel <- rlgt(trainData,  #seasonality2=SEASONALITY2, seasonality2 not specified -> we are running SGT   
					seasonality.type=seasonality.type, level.method=level.method, error.size.method=error.size.method, 
					verbose=TRUE)					
			forec= forecast(rstanmodel, h = H, level=c(90,98))
			# str(forec, max.level=1)
			
			forecasts[[seriesName]]<-list(mean=forec$mean, lower=forec$lower, upper=forec$upper)
			plot(forec, main=paste(seriesName,seasonality.type,level.method,error.size.method))
			
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
			i=i+1    #needed for summaries
		}
	}
}

sMAPE=sumSMAPE/i
q95Loss=sumQ95Loss/i
q99Loss=sumQ99Loss/i
q5Loss=sumQ5Loss/i
exceed95=numOfCases95pExceeded/(i*H)*100
exceed5=numOfCases5pExceeded/(i*H)*100
print(paste0("SUMMARY: Num of cases:", i, ", sMAPE:",signif(sMAPE,4),
  ', % of time 95p exceeded:',signif(exceed95,4), ', % of time 5p exceeded:',signif(exceed5,4), 
	', q5Loss:',signif(q5Loss,4),', q95Loss:',signif(q95Loss,4),', q99Loss:',signif(q99Loss,4) ))

