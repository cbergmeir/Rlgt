# Test SGT - M3 monthly data
# This is a most computationally demanding subset of M3, so also a good opportunity to showcase parallel execution
# We will use all cores on the server - the more of them the faster the whole process will be.
# After a few days (on a 6-core computer) you will see something like this:
# Num of cases:1428, sMAPE:13.53, % of time 95p exceeded:6.26, % of time 5p exceeded:5.155, q5Loss:3.913, q95Loss:6.645, q99Loss:2.125"

library(Mcomp)
library(Rlgt)
library(doParallel)
#set.seed(12)
options(width=180)

numOfCores=parallel:::detectCores()
CLUSTER_SIZE=as.integer(numOfCores/4)  #Each cluster will use 4 cores

M3.data <- subset(M3,"monthly")
SEASONALITY=12
M3.data <- sample(M3.data) #shuffle
NUM_OF_CASES=length(M3.data)
#NUM_OF_CASES=20
H=length(M3.data[[1]]$xx)

quantileLoss<-function(forec, actual, tau) {
	diff=actual-forec
	pinBallL=pmax(diff*tau,diff*(tau-1))
	mean(pinBallL/actual)*200
}


#stopCluster(cl) # this is useful if rerunning the script in the same session. First time it will report error - it is OK
#Sys.sleep(2)
cl = makeCluster(CLUSTER_SIZE, outfile="rclust.out")
registerDoParallel(cl)
sink(); graphics.off();
stanOutputPath='stanOutput.txt' #in working directory
unlink(stanOutputPath)

i=1
ret_df=foreach(i=1:NUM_OF_CASES, .combine=rbind, .inorder=FALSE, .packages=c("Rlgt")) %dopar% { # the long loop :-)	
	if (i==1) { #needs to be done inside every new slave process
		sink(file=stanOutputPath, append=TRUE, split=TRUE)
	}
	series=M3.data[[i]]$sn
	data.train <- M3.data[[i]]$x
	data.test <- M3.data[[i]]$xx
	rstanmodel <- rlgt(data.train, model="gSGT", nCores=4, nChains=4,
		control=rlgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=5000), #we do not need to specify SEASONALITY in rlgt.control, as it is extracted from M3.data[[i]]$x
		verbose=FALSE)
	forec= forecast(rstanmodel, h = H, level=c(90,98))
	sMAPE=mean(abs(forec$mean-data.test)/(forec$mean+data.test))*200
	q95Loss=quantileLoss(forec$upper[,1], data.test, 0.95)	
	q99Loss=quantileLoss(forec$upper[,2], data.test, 0.99)
	q5Loss=quantileLoss(forec$lower[,1], data.test, 0.05)
	data.frame(series=series, sMAPE=sMAPE, 
		q5Loss=q5Loss, q95Loss=q95Loss, q99Loss=q99Loss,
		numOfCases95pExceeded=sum(data.test>forec$upper[,1]),
		numOfCases5pExceeded=sum(data.test<forec$lower[,1]))
}	
sMAPE=mean(ret_df$sMAPE)
q95Loss=mean(ret_df$q95Loss)
q99Loss=mean(ret_df$q99Loss)
q5Loss=mean(ret_df$q5Loss)
exceed95=mean(ret_df$numOfCases95pExceeded)/H*100
exceed5=mean(ret_df$numOfCases5pExceeded)/H*100
print(paste0("SUMMARY: Num of cases:", NUM_OF_CASES, ", sMAPE:",signif(sMAPE,4),
  ', % of time 95p exceeded:',signif(exceed95,4), ', % of time 5p exceeded:',signif(exceed5,4), 
	', q5Loss:',signif(q5Loss,4),', q95Loss:',signif(q95Loss,4),', q99Loss:',signif(q99Loss,4) ))
