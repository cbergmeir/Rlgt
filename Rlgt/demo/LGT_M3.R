# Test lgt - M3 yearly data

#install.packages("Mcomp")
library(Mcomp)
library(Rlgt)
#set.seed(12)

options(width=180)

M3.data <- subset(M3,"yearly")
M3.data <- sample(M3.data) #shuffle
NUM_OF_CASES=length(M3.data)  # If you let it run its full course (comment out next line), you should see a very good result :-)
#NUM_OF_CASES=5

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

#i<-1; str(M3.data[[i]])
H=length(M3.data[[1]]$xx)
sumSMAPE=0; sumQ99Loss=0; sumQ95Loss=0; sumQ5Loss=0;
numOfCases95pExceeded=0; numOfCases5pExceeded=0;
for (i in 1:NUM_OF_CASES) {
	series=M3.data[[i]]$sn
	if(i==1) { #just for demo and testing. In your code stick to one of the two alternatives. Plotting, etc. is easier with ts inputs.
		trainData <- as.numeric(M3.data[[i]]$x)  #"naked", numeric vector
		actuals <- as.numeric(M3.data[[i]]$xx)   # class of actuals has to be the same
		rstanmodel <- rlgt(trainData, 
				control=rlgt.control(NUM_OF_ITER=4000), 
				verbose=FALSE)
	} else {
		trainData <- M3.data[[i]]$x   # trainData is of ts class
		actuals <- M3.data[[i]]$xx    # class of actuals has to be the same
		rstanmodel <- rlgt(trainData,
				control=rlgt.control(NUM_OF_ITER=4000), 
				verbose=FALSE)
	}
	forec= forecast(rstanmodel, h = H, level=c(90,98))
	
	plot(forec, main=series)
	
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
	print(paste0(series," sMAPE:",signif(sMAPE,3) ,' q5Loss:',signif(q5Loss,3),' q95Loss:',signif(q95Loss,3),' q99Loss:',signif(q99Loss,3) ))
}
sMAPE=sumSMAPE/i
q95Loss=sumQ95Loss/i
q99Loss=sumQ99Loss/i
q5Loss=sumQ5Loss/i
exceed95=numOfCases95pExceeded/(i*H)*100
exceed5=numOfCases5pExceeded/(i*H)*100
print(paste0("SUMMARY: Num of cases:", i, ", sMAPE:",signif(sMAPE,3),
	', % of time 95p exceeded:',signif(exceed95,3), ', % of time 5p exceeded:',signif(exceed5,3), 
	', q5Loss:',signif(q5Loss,3),', q95Loss:',signif(q95Loss,3),', q99Loss:',signif(q99Loss,3) ))

