# LGT with custom Gibbs sampler
library(Mcomp)
library(Rlgt)

quantileLoss<-function(forec, actual, tau) {
  diff=actual-forec
  pinBallL=pmax(diff*tau,diff*(tau-1))
  mean(pinBallL/actual)*200
}

M3.data <- subset(M3,"yearly")
M3.data <- sample(M3.data) #shuffle
w.series <- length(M3.data) 

H <- length(M3.data[[1]]$xx)

sumSMAPE=0; sumQ99Loss=0; sumQ95Loss=0; sumQ5Loss=0;
numOfCases95pExceeded=0; numOfCases5pExceeded=0;

for (i in 1:w.series) {
  
  series <- M3.data[[i]]$sn
  
  trainData <- as.numeric(M3.data[[i]]$x) #"naked", numeric vector
  actuals <- as.numeric(M3.data[[i]]$xx)   # class of actuals has to be the same
  model <- rlgt(trainData, 
                control=rlgt.control(NUM_OF_ITER=4000),
                method = "Custom_Gibbs",
                verbose=FALSE)
  
  forec <- forecast(model, h = H, level=c(90,98))
  
  plot(forec, main=series)
  xs <- seq(from=length(trainData)+1,to=length(trainData)+ length(actuals))
  lines(xs,actuals, col=1, type='b',lwd=2)
  
  sMAPE <- mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
  sumSMAPE=sumSMAPE+sMAPE
  
  numOfCases95pExceeded=numOfCases95pExceeded+sum(actuals>forec$upper[,1])
  numOfCases5pExceeded=numOfCases5pExceeded+sum(actuals<forec$lower[,1])
  
  q95Loss <- quantileLoss(forec$upper[,1], actuals, 0.95)
  q99Loss <- quantileLoss(forec$upper[,2], actuals, 0.99)
  q5Loss <- quantileLoss(forec$lower[,1], actuals, 0.05)
  
  sumQ95Loss=sumQ95Loss+q95Loss
  sumQ99Loss=sumQ99Loss+q99Loss
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

