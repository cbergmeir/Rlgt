# Test SGT - M3 quarterly data
# If you let it run  full course, you should see a good result :-)

library("Mcomp")
library("Rlgt")
#set.seed(12)

options(width=180)

M3.data <- subset(M3,"quarterly")
M3.data <- sample(M3.data) #shuffle

quantileLoss<-function(forec, actual, tau) {
	diff=actual-forec
	pinBallL=pmax(diff*tau,diff*(tau-1))
	mean(pinBallL/actual)*200
}

#iter<-1; str(M3.data[[iter]])
sumSMAPE=0; sumQ99Loss=0; sumQ95Loss=0; sumQ5Loss=0;
for (iter in 1:length(M3.data)) {
	series=M3.data[[iter]]$sn
  data.train <- M3.data[[iter]]$x
  data.test <- M3.data[[iter]]$xx
  rstanmodel <- fit.lgt(data.train, model="SGT", nCores=4, nChains=4,
                          control=lgt.control(SEASONALITY=4, MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=5000), 
                          verbose=TRUE)
	forec= forecast(rstanmodel, h = length(data.test), level=c(90,98))
	
	sMAPE=mean(abs(forec$median-data.test)/(forec$median+data.test))*200
	sumSMAPE=sumSMAPE+sMAPE
	q95Loss=quantileLoss(forec$upper[,1], data.test, 0.95)
	sumQ95Loss=sumQ95Loss+q95Loss
	q99Loss=quantileLoss(forec$upper[,2], data.test, 0.99)
	sumQ99Loss=sumQ99Loss+q99Loss
	q5Loss=quantileLoss(forec$lower[,1], data.test, 0.05)
	sumQ5Loss=sumQ5Loss+q5Loss
	print(paste0(series," sMAPE:",signif(sMAPE,3) ,' q5Loss:',signif(q5Loss,3),' q95Loss:',signif(q95Loss,3),' q99Loss:',signif(q99Loss,3) ))
}
sMAPE=sumSMAPE/length(M3.data)
q95Loss=sumQ95Loss/length(M3.data)
q99Loss=sumQ99Loss/length(M3.data)
q5Loss=sumQ5Loss/length(M3.data)
print(paste0("SUMMARY: Num of cases:", length(M3.data), " sMAPE:",signif(sMAPE,4),' q5Loss:',signif(q5Loss,3),' q95Loss:',signif(q95Loss,3),' q99Loss:',signif(q99Loss,3) ))

