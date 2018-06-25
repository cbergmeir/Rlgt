# Test SGT - M3 monthly data
# This is a most computationally demanding subset of M3, so also a good opportunity to showcase parallel execution
# We will use all cores on the server - the more of them the faster the whole process will be.
# If you let it run full course, you should see a good result :-)

library(Mcomp)
library(Rlgt)
library(doParallel)
#set.seed(12)
options(width=180)

numOfCores=parallel:::detectCores()
CLUSTER_SIZE=as.integer(numOfCores/3)  #Each cluster will use 4 cores, so it looks like oversubscription, 
#but processes in one cluster tend not to finish at the same time.



M3.data <- subset(M3,"monthly")
SEASONALITY=12
M3.data <- sample(M3.data) #shuffle

quantileLoss<-function(forec, actual, tau) {
	diff=actual-forec
	pinBallL=pmax(diff*tau,diff*(tau-1))
	mean(pinBallL/actual)*200
}


#stopCluster(cl) # this is useful if rerunning the script in the same session. First time it will report error - it is OK
#Sys.sleep(2)
cl = makeCluster(CLUSTER_SIZE, outfile="rclust.out")
sink(); graphics.off();
stanOutputPath='stanOutput.txt' #in working directory
unlink(stanOutputPath)

firstTime=TRUE; i=1
sumSMAPE=0; sumQ99Loss=0; sumQ95Loss=0; sumQ5Loss=0;
#ret_df=foreach(i=1:length(M3.data), .combine=rbind, .inorder=FALSE, .packages=c("sn","rstan","Rlgt")) %dopar% { # the long loop :-)
for (i in 1:length(M3.data)) {	
	if (firstTime==TRUE) { #needs to be done inside every new slave process
		sink(file=stanOutputPath, append=TRUE, split=TRUE)
		firstTime=FALSE
	}
	series=M3.data[[iter]]$sn
	data.train <- M3.data[[iter]]$x
	data.test <- M3.data[[iter]]$xx
	rstanmodel <- fit.lgt(data.train, model="gSGT", nCores=4, nChains=4,
			control=lgt.control(SEASONALITY=SEASONALITY, MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=5000), 
			verbose=TRUE)
	forec= forecast(rstanmodel, h = length(data.test), level=c(90,98))
	sMAPE=mean(abs(forec$median-data.test)/(forec$median+data.test))*200
	q95Loss=quantileLoss(forec$upper[,1], data.test, 0.95)	
	q99Loss=quantileLoss(forec$upper[,2], data.test, 0.99)
	q5Loss=quantileLoss(forec$lower[,1], data.test, 0.05)
	data.frame(series=series, sMAPE=sMAPE, q5Loss=q5Loss, q95Loss=q95Loss, q99Loss=q99Loss)
}	




sumSMAPE=sumSMAPE+sMAPE
sumQ95Loss=sumQ95Loss+q95Loss
sumQ99Loss=sumQ99Loss+q99Loss
sumQ5Loss=sumQ5Loss+q5Loss	

print(paste0(series," sMAPE:",signif(sMAPE,3) ,' q5Loss:',signif(q5Loss,3),' q95Loss:',signif(q95Loss,3),' q99Loss:',signif(q99Loss,3) ))





#iter<-1; str(M3.data[[iter]])

for (iter in 1:length(M3.data)) {
	
}
sMAPE=sumSMAPE/length(M3.data)
q95Loss=sumQ95Loss/length(M3.data)
q99Loss=sumQ99Loss/length(M3.data)
q5Loss=sumQ5Loss/length(M3.data)
print(paste0("SUMMARY: Num of cases:", length(M3.data), " sMAPE:",signif(sMAPE,4),' q5Loss:',signif(q5Loss,3),' q95Loss:',signif(q95Loss,3),' q99Loss:',signif(q99Loss,3) ))

