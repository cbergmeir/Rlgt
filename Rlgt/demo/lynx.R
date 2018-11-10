# see https://andrewgelman.com/2017/03/01/facebooks-prophet-uses-stan  :-)
# 

library(Rlgt)

str(lynx)

SEASONALITY=38  #by looking at the graph and running acf()

train=lynx[1:80]
actuals=lynx[81:length(lynx)]

rstanmodel <- rlgt(train, seasonality=SEASONALITY,
		control=rlgt.control(NUM_OF_ITER=10000),   
		verbose=TRUE)   

forec= forecast(rstanmodel, h = length(actuals))

plot(forec)
xs=seq(from=length(train)+1,to=length(train)+ length(actuals))
lines(xs,actuals, col=1, type='b',lwd=2)	

sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
msqrt=sqrt(mean(forec$mean-actuals)^2)
print(paste("sMAPE:",signif(sMAPE,3), "mse", signif(msqrt)))

