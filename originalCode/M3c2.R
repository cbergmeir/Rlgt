# Author: Slawek Smyl 
# Feb 2016
# Added my GT algorithm to prof. R Hyndman's script 
# attached to his blog http://robjhyndman.com/hyndsight/show-me-the-evidence/
# Please note: 
# 1. Depending on your computer hardware (number and speed of the cores), 
#   it will take between 1 to 7 days to run this script
# 2. This is an example code, accompanying a paper, not a production-quality code.
###############################################################################
# install.packages("Mcomp")
# install.packages("rstan") #but first install RTools

USE_DB=FALSE  #if using the database, 
# 1) create the destination table first (look at the end of the file for an example script)
# 2) create ODBC 64-bit system name "slawek" pointing to the destination database. For SQL Server use the latest ODBC Driver

require(Mcomp)
library(rstan)
library(sn)
library(doParallel)
library(R2HTML)
if (USE_DB) {
	library(RODBC)
}

CAUCHY_SD_DIV=200
runName='std200' #if using a database, make it different for every run

numOfCores=parallel:::detectCores()
if (numOfCores==8) {
	CLUSTER_SIZE=3
} else if (numOfCores==16) {
	CLUSTER_SIZE=9
} else {
	CLUSTER_SIZE=as.integer(numOfCores/2)
} 

OUTPUT_DIR="c:/reports/m3"
if (USE_DB) {
  FORECAST_TABLE='ForecastM3C'
  ODBC_SOURCE_NAME='slawek'
	fchannel<- odbcConnect(ODBC_SOURCE_NAME) #ODBC 64-bit system name needs to be created. For SQL Server use latest ODBC Driver  
	odbcGetInfo(fchannel)
}

STAN_CLUSTER_SIZE=4
rstan_options(auto_write = TRUE)
options(mc.cores = STAN_CLUSTER_SIZE, width=180)


MAX_RHAT_ALLOWED=1.005 #see Stan's manual
NUM_OF_ITER=5000
MAX_NUM_OF_REPEATS=3
MIN_SIGMA=0.001 # for numerical stability
MIN_VAL=0.001
MAX_VAL=1e38 # max for real type of SQL Server, hopefully never hit :-) 
MIN_NU=1.5
MAX_NU=20
MIN_POW=-0.5;
MAX_POW=1

Sys.setenv(TZ = "GMT");
Sys.timezone()
today=Sys.Date();
now=Sys.time()
imageWidth=1000; imageHeight=400
NUM_OF_TRIALS=5000  
NUM_OF_TRIALS_TO_SHOW=100
reportTitle='M3'
allSave_df=NULL


#non-seasonal model
parametersX <- c("l", "b", "nu", "sigma", "levSm",  "bSm", 
		"powx", "coefTrend",  "powTrend", "offsetSigma", "locTrendFract")
modelX <- '
	data {  
		real<lower=0> CAUCHY_SD;
		real MIN_POW;  real MAX_POW;
		real<lower=0> MIN_SIGMA;
		real<lower=1> MIN_NU; real<lower=1> MAX_NU;
		int<lower=1> N;
		vector<lower=0>[N] y;
		//real<lower=0> alphaBETA; real<lower=0> betaBETA; 
	}
	parameters {
		real<lower=MIN_NU,upper=MAX_NU> nu; 
		real<lower=0> sigma;
		real <lower=0,upper=1>levSm;
		real <lower=0,upper=1> bSm;
		real <lower=0,upper=1> powx;
		real <lower=0,upper=1> powTrendBeta;
		real coefTrend;
		real <lower=MIN_SIGMA> offsetSigma;
		real <lower=-1,upper=1> locTrendFract;
	} 
	transformed parameters {
		real <lower=MIN_POW,upper=MAX_POW>powTrend;
		vector[N] l; vector[N] b;
		
		l[1] <- y[1]; b[1] <- 0;
		powTrend<- (MAX_POW-MIN_POW)*powTrendBeta+MIN_POW;
		
		for (t in 2:N) {
			l[t]  <- levSm*y[t] + (1-levSm)*l[t-1] ;
			b[t]  <- bSm*(l[t]-l[t-1]) + (1-bSm)*b[t-1] ;
		}
	}
	model {
		sigma ~ cauchy(0,CAUCHY_SD) T[0,];
		offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
		coefTrend ~ cauchy(0,CAUCHY_SD);
		//powTrendBeta ~ beta(alphaBETA, betaBETA);
		
		for (t in 2:N) {
		y[t] ~ student_t(nu, l[t-1]+coefTrend*fabs(l[t-1])^powTrend+locTrendFract*b[t-1], 
		  sigma*fabs(l[t-1])^powx+ offsetSigma);
	}
}
'  
stanModelX <- stan_model(model_code=modelX)
#str(stanModel)


#seasonal
parametersS <- c("l", "s", "sSm","nu", "sigma", "levSm", 
		"powx", "coefTrend", "powTrend", "offsetSigma")
modelS <- '
	data {  
		int<lower=2> SEASONALITY;
		real<lower=0> CAUCHY_SD;
		real MIN_POW;  real MAX_POW;
		real<lower=0> MIN_SIGMA;
		real<lower=1> MIN_NU; real<lower=1> MAX_NU;
		int<lower=1> N;
		vector<lower=0>[N] y;
		//real<lower=0> alphaBETA; real<lower=0> betaBETA; 
	}
	parameters {
		real<lower=MIN_NU,upper=MAX_NU> nu; 
		real<lower=0> sigma;
		real <lower=0,upper=1>levSm;
		real <lower=0,upper=1>sSm;
		real <lower=0,upper=1>powx;
		real <lower=0,upper=1> powTrendBeta;
		real coefTrend;
		real <lower=MIN_SIGMA> offsetSigma;
		vector[SEASONALITY] initSu; #unnormalized
	} 
	transformed parameters {
		real <lower=MIN_POW,upper=MAX_POW>powTrend;
		vector[N] l;
		vector[SEASONALITY] inits;
		vector[N+SEASONALITY] s;
		real sumsu;
		
		sumsu <- 0;
		for (i in 1:SEASONALITY) 
		  sumsu <- sumsu+ initSu[i];
		for (i in 1:SEASONALITY) 
		  inits[i] <- initSu[i]*SEASONALITY/sumsu;
		
		for (i in 1:SEASONALITY) {
		  s[i] <- inits[i];
		}
		s[SEASONALITY+1] <- inits[1];
		
		l[1] <- y[1]/(s[1]);
		powTrend<- (MAX_POW-MIN_POW)*powTrendBeta+MIN_POW;
		
		for (t in 2:N) {
		  l[t]  <- levSm*y[t]/(s[t]) + (1-levSm)*l[t-1] ;  
		  s[t+SEASONALITY] <- sSm*y[t]/l[t]+(1-sSm)*s[t];
		}
	}
	model {
		real expVal;
		
		//powx ~ beta(alphaBETA,betaBETA);
		sigma ~ cauchy(0,CAUCHY_SD) T[0,];
		offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
		coefTrend ~ cauchy(0, CAUCHY_SD);
		for (t in 1:SEASONALITY) {
		  initSu[t] ~ normal (1, 0.3) T[0.01,];
		}
		
		for (t in 2:N) {
		  expVal <- (l[t-1]+ coefTrend*fabs(l[t-1])^powTrend)*s[t];
		  y[t] ~ student_t(nu, expVal, sigma*fabs(expVal)^powx+ offsetSigma);
		}
	}
' 
stanModelS <- stan_model(model_code=modelS)
#str(stanModelS)



BeginReport = function(title) {
	cat(paste("<html><head><title>",
					title, "</title></head>",
					"<body bgcolor=#D0D0D0>", sep = ""),
			file = htmlFilePath, append = FALSE)
	cat(paste("<html><head><title>",
					title, "</title></head>",
					"<body bgcolor=#D0D0D0>", sep = ""),
			file = htmlFilePath2, append = FALSE)
}
EndReport = function() {
	cat("<hr size=1></body></html> \n ",
			file = htmlFilePath, append = TRUE)
	cat("<hr size=1></body></html> \n ",
			file = htmlFilePath2, append = TRUE)
}
WriteHead = function(head) {
	cat("</p><br><br><h2>",head,"</h2>\n ",
			file = htmlFilePath, append = TRUE)
	cat("</p><br><br><h2>",head,"</h2>\n ",
			file = htmlFilePath2, append = TRUE)
}
InsertGraph=function(imageFileName) {
	relPath2=paste('images\\',imageFileName, sep='')
	relPath=paste(today,'_',runName,'\\images\\',imageFileName, sep='')
	HTMLInsertGraph(relPath, Caption = '', WidthHTML=imageWidth, HeightHTML=imageHeight, file=htmlFilePath)
	HTMLInsertGraph(relPath2, Caption = '', WidthHTML=imageWidth, HeightHTML=imageHeight, file=htmlFilePath2)
}
WriteHTML=function(html) {
	cat(html," \n", file = htmlFilePath, append = TRUE)
	cat(html," \n", file = htmlFilePath2, append = TRUE)
}


####
nseries <- length(M3)
theta <- as.matrix(M3Forecast$THETA)
fpro <- as.matrix(M3Forecast$ForecastPro)
fcx <- as.matrix(M3Forecast$ForcX)
bjauto <- as.matrix(M3Forecast$`B-J auto`)
ab1 <- as.matrix(M3Forecast$AutoBox1)
ab2 <- as.matrix(M3Forecast$AutoBox2)
ab3 <- as.matrix(M3Forecast$AutoBox3)


DAY_DIR=paste(OUTPUT_DIR,paste(today,runName,sep='_'),sep='/')
stanOutputPath=paste(DAY_DIR, 'stanOutput.txt',sep='/')
IMAGES_DIR=paste(DAY_DIR,'images',sep='/')

dir.create(IMAGES_DIR, recursive=T)
htmlFilePath=paste(OUTPUT_DIR, '/',reportTitle,'.html',sep='')
htmlFilePath2=paste(DAY_DIR, '/',reportTitle,'.html',sep='')
unlink(htmlFilePath)
unlink(htmlFilePath2)

title_=paste(reportTitle,'created on',Sys.Date())
BeginReport(title_)
WriteHTML(title_)
WriteHead(reportTitle)


stopCluster(cl) # this is useful if rerunning the script in the same session. First time it will report error - it is OK
Sys.sleep(2)
CLUST_OUT_PATH=paste(OUTPUT_DIR,"rclust.txt",sep='/')
file.remove(CLUST_OUT_PATH)
cl <- makeCluster(CLUSTER_SIZE, outfile=CLUST_OUT_PATH)
registerDoParallel(cl)

sink(); graphics.off();
unlink(stanOutputPath)

firstTime=TRUE; idr=225
ret_df<-foreach(idr=1:nseries, .combine=rbind, .packages=c("sn","rstan","RODBC","R2HTML","Mcomp")) %dopar% { # the long loop :-)
	if (firstTime==TRUE) { #needs to be done inside every new slave process
		fchannel<- odbcConnect(ODBC_SOURCE_NAME)
		odbcGetInfo(fchannel)
		sink(file=stanOutputPath, append=TRUE, split=TRUE)
		firstTime=FALSE
	}
	# str(M3[[1]]); str(M3[[700]]); str(M3[[1700]]);   
	
	series=M3[[idr]]$sn
	ets_pred <- as.numeric(forecast(ets(M3[[idr]]$x),h=18,PI=FALSE)$mean)
	aarima_pred <- as.numeric(forecast(auto.arima(M3[[idr]]$x),h=18)$mean)
	hybrid_pred <- 0.5*(ets_pred + aarima_pred)
	pred_df=data.frame(series=series,ets=ets_pred, aarima=aarima_pred, hybrid=hybrid_pred, lgtso=NA)
	# str(pred_df)
	 
	proceed=TRUE
	if (USE_DB) {
		doneAlready_query=paste("select predQ50 from ",FORECAST_TABLE,
				" where [run]='",runName,  
				"' and series='",series,
				"' order by [predictionHorizon]",  sep='')
		doneAlready_df=sqlQuery(fchannel,doneAlready_query, stringsAsFactors =F)
		if (nrow(doneAlready_df)>0) {
			pred_df$lgtso[1:length(doneAlready_df$predQ50)]=doneAlready_df$predQ50
			proceed=FALSE
		}
	}
	
	if (proceed) {
		x=as.numeric(M3[[idr]]$x)
		n=M3[[idr]]$n
		y=x+rnorm(n,0,sd=abs(min(x))*0.0001) # I found that adding a bit of jitter is helping Stan in case of some flat series		
		# summary(y)		
		yy=as.numeric(M3[[idr]]$xx)
		maxPredictionHorizon=M3[[idr]]$h
		
		plot(c(y,yy), main=series,  col=c(rep(1,n),rep(2,maxPredictionHorizon)), ylab='',xlab='years') #just on screen for debugging
		CauchySd=max(y)/CAUCHY_SD_DIV
		print(paste("CauchySd",round(CauchySd,2)))
		print(paste(Sys.time(),"starting ", series))  
		
		variable=M3[[idr]]$period
		if (variable=="QUARTERLY") {
			seasonality=4
			stanModel=stanModelS
			parameters=parametersS
		} else if (variable=="MONTHLY") {
			seasonality=12
			stanModel=stanModelS
			parameters=parametersS
		} else {
			seasonality=1
			stanModel=stanModelX
			parameters=parametersX
		}
	
		data <- list(SEASONALITY=seasonality, CAUCHY_SD=CauchySd, 
				MIN_POW=MIN_POW, MAX_POW=MAX_POW, 
				MIN_SIGMA=MIN_SIGMA, 
				MIN_NU=MIN_NU,  MAX_NU=MAX_NU, 
				#alphaBETA=alphaBETA, betaBETA=betaBETA,
				y=y, N=n) # to be passed on to Stan
		
		
		avgRHat=1e20; irep=1
		for (irep in 1:MAX_NUM_OF_REPEATS) {
			initializations <- list();
			for (irr in 1:STAN_CLUSTER_SIZE) {
					initializations[[irr]]=list( 
							initSu=runif(seasonality,min=0.8,max=1.2)#this initialization will be updated inside Stan, not too important
					)
			}
			samples1=
					sampling(#control=list(adapt_delta = 0.9, max_treedepth=11),
							stanModel,   
							data=data, 
							init=initializations,
							pars=parameters,
							iter=NUM_OF_ITER*2^(irep-1),
							chains=STAN_CLUSTER_SIZE,
							cores=STAN_CLUSTER_SIZE,
							open_progress=F,
							refresh = 1000)
			samples1
			
			ainfo=summary(samples1)
			RHats=ainfo$summary[,10]
			RHats=as.numeric(RHats[is.finite(RHats)])
			currRHat=mean(RHats, na.rm=T)
			if (currRHat<=MAX_RHAT_ALLOWED) {
				samples=samples1
				avgRHat=currRHat
				print(samples)
				print(paste("avgRHat",avgRHat))
				break
			} else {
				if (currRHat<avgRHat) {#in the worst case this is at least once executed, because avgRHat is initialized high
					samples=samples1
					avgRHat=currRHat
					print(samples)
					print(paste("avgRHat",avgRHat))
				} else {
					print ("worse...")
					print(paste("currRHat",currRHat))
				}
				print (paste("trying to do better...",series))
			}
			#str(samples, max.level =4)
		}#repeat if needed
		print(summary(do.call(rbind, args = get_sampler_params(samples1, inc_warmup = F)), digits = 2)) #diagnostics
		
		l=extract(samples)$l
		lM=apply(l,2,mean)
		#plot(lM, type='l')
		mc_size=dim(l)[1]
		mc_range=1:mc_size
		lastDataIndex=dim(l)[2]
		lastLevel=l[,lastDataIndex]
		#str(lastLevel); summary(lastLevel)
		lastLevelM=mean(lastLevel)
		paste("lastLevel",round(lastLevelM,1))
		
		levSm <- extract(samples)$levSm
		summary(levSm)
		levSmM=mean(levSm)
		#print(paste("levSm",levSmM))
	
		nu <- extract(samples)$nu
		summary(nu)
		nuM=mean(nu)
		#print(paste("nu",nuM))
	
		coefTrend <- extract(samples)$coefTrend
		summary(coefTrend)
		coefTrendM=mean(coefTrend)
		#print(paste("coefTrend",coefTrendM))
	
		powTrend <- extract(samples)$powTrend
		summary(powTrend)
		powTrendM=mean(powTrend)
		#print(paste("powTrend",powTrendM))
		
		sigma <- extract(samples)$sigma
		summary(sigma)
		sigmaM=mean(sigma)
		#print(paste("sigma",sigmaM))
	
		powx <- extract(samples)$powx
		summary(powx)
		powxM=mean(powx)
		#print(paste("powx",powxM))
	
		offsetSigma=extract(samples)$offsetSigma
		summary(offsetSigma)
		offsetSigmaM=mean(offsetSigma)
		paste("offsetSigma",round(offsetSigmaM,2))
		
		if (seasonality>1) {
			sSm <- extract(samples)$sSm
			#summary(sSm)
			#hist(sSm,breaks=200)
			sSmM=mean(sSm)
			#print(paste("sSm",sSmM))
			
			s=extract(samples)$s
			sM=apply(s,2,mean)
			cat("s:"); print(sM)
			print(summary(sM))
			#plot(sM,type='l')
			
		  plot(c(y,yy), main=series,  col=c(rep(1,n),rep(2,maxPredictionHorizon)), ylab='',xlab='years') #just on screen for debugging
			lines(lM*sM[1:length(lM)],col='red')
			lines(lM,col='blue')
			
			lastDetails=paste("coefT=",round(coefTrendM,2), ", powT=",round(powTrendM,2),
					", nu=",round(nuM,2), 
					", sigma=",round(sigmaM,2),", powx=",round(powxM,2), ", Sigma0=",round(offsetSigmaM,2), 
					", lSm=",round(levSmM,2),	", sSm=",round(sSmM,2), sep='')
		} else {
			b=extract(samples)$b
			bM=apply(b,2,mean)
			lastB=b[,lastDataIndex]
			lastBM=mean(lastB)
			paste("lastB",round(lastBM,2))
			
			bSm <- extract(samples)$bSm
			summary(bSm)
			bSmM=mean(bSm)
			paste("bSm",round(bSmM,2))
			
			locTrendFract <- extract(samples)$locTrendFract
			print(summary(locTrendFract))
			locTrendFractM=mean(locTrendFract)
			paste("locTrendFract",round(locTrendFractM,2))
			
			lastDetails=paste("coefT=",round(coefTrendM,2), ", powT=",round(powTrendM,2),
					", sigma=",round(sigmaM,2),", powx=",round(powxM,2),", offsetS=",round(offsetSigmaM,2), 
					", nu=",round(nuM,2), ", lSm=",round(levSmM,2), 
					", lastB=",round(lastBM,1), ", bSm=",round(bSmM,2), 
					", lTFract=", round(locTrendFractM,2), sep='')
		}
		
		print(paste(Sys.time(),series))
		print(lastDetails)
		
		
		t=1; irun=1
		yf=matrix(lastLevelM,nrow=NUM_OF_TRIALS, ncol=maxPredictionHorizon)
		if (seasonality>1) {
			sM=rep(1,seasonality+maxPredictionHorizon)
			for (irun in 1:NUM_OF_TRIALS) {
				indx=sample(mc_range,1)
				prevLevel=mean(lastLevel[indx], na.rm=T)
				powxM=mean(powx[indx], na.rm=T)
				powTrendM=mean(powTrend[indx], na.rm=T)
				coefTrendM=mean(coefTrend[indx], na.rm=T)
				levSmM=mean(levSm[indx], na.rm=T)
				sSmM=mean(sSm[indx], na.rm=T)
				nuM=mean(nu[indx], na.rm=T)
				sigmaM=mean(sigma[indx], na.rm=T)
				offsetSigmaM=mean(offsetSigma[indx], na.rm=T)
				lastsM=s[indx,]
		
				for (iw in 1:seasonality) {
					sM[iw]=lastsM[length(lastsM)-seasonality+iw]
				}
				
				for (t in 1:maxPredictionHorizon) {
					seasonA=sM[t]
					expVal <- (prevLevel+ coefTrendM*abs(prevLevel)^powTrendM)*seasonA;
					error=rst(n=1, xi=0, 
							omega=sigmaM*abs(expVal)^powxM+offsetSigmaM, 
							alpha=0, nu=nuM)
					yf[irun,t]=max(MIN_VAL,expVal+error)
					currLevel=max(MIN_VAL,levSmM*yf[irun,t]/seasonA + (1-levSmM)*prevLevel) ;
					#s[t+seasonality] <- sSm*y[t]/l[t]+(1-sSm)*s[t];
					sM[t+seasonality] <- sM[t];
					
					if (currLevel>MIN_VAL) {
						prevLevel=currLevel
					} 
				}
			}
		} else {
			for (irun in 1:NUM_OF_TRIALS) {
				indx=sample(mc_range,1)
				prevLevel=mean(lastLevel[indx], na.rm=T)
				powxM=mean(powx[indx], na.rm=T)
				powTrendM=mean(powTrend[indx], na.rm=T)
				coefTrendM=mean(coefTrend[indx], na.rm=T)
				bM=mean(lastB[indx], na.rm=T)
				levSmM=mean(levSm[indx], na.rm=T)
				bSmM=mean(bSm[indx], na.rm=T)
				nuM=mean(nu[indx], na.rm=T)
				sigmaM=mean(sigma[indx], na.rm=T)
				offsetSigmaM=mean(offsetSigma[indx], na.rm=T)
				locTrendFractM=mean(locTrendFract[indx], na.rm=T)
				
				for (t in 1:maxPredictionHorizon) { 
					error=rst(n=1, 
							xi=coefTrendM*(abs(prevLevel))^powTrendM + locTrendFractM*bM , 
							omega=sigmaM*(abs(prevLevel))^powxM+offsetSigmaM, alpha=0, nu=nuM)
					yf[irun,t]=min(MAX_VAL,max(MIN_VAL,prevLevel+error))
					currLevel=max(MIN_VAL,levSmM*yf[irun,t] + (1-levSmM)*prevLevel) ;
					
					if (currLevel>MIN_VAL) {
						bM= bSmM*(currLevel-prevLevel)+(1-bSmM)*bM
						prevLevel=currLevel
					} 
				}
				# yf[irun,]
			}
		}
		avgYfs=apply(yf,2,quantile,probs=c(0.05,0.5,0.95,0.99,0.01))
		pred_df$lgtso[1:(dim(avgYfs)[2])] <-avgYfs[2,]
		
		save_df=NULL; predictionHorizon=1
		for (predictionHorizon in 1:maxPredictionHorizon) {
			if (is.null(save_df))
				save_df=data.frame(run=runName, series=series,
						predictionHorizon=predictionHorizon, 
						dateTimeOfPrediction=now, trueValue=yy[predictionHorizon], predQ1=avgYfs[5,predictionHorizon],
						predQ5=avgYfs[1,predictionHorizon], predQ50=avgYfs[2,predictionHorizon], 
						predQ95=avgYfs[3,predictionHorizon], predQ99=avgYfs[4,predictionHorizon], 
						variable=variable, algDetails=lastDetails)
			else 
				save_df=rbind(save_df, data.frame(run=runName, series=series,
								predictionHorizon=predictionHorizon, 
								dateTimeOfPrediction=now, trueValue=yy[predictionHorizon], predQ1=avgYfs[5,predictionHorizon],
								predQ5=avgYfs[1,predictionHorizon], predQ50=avgYfs[2,predictionHorizon], 
								predQ95=avgYfs[3,predictionHorizon], predQ99=avgYfs[4,predictionHorizon], 
								variable=variable, algDetails=lastDetails))
		}
		if (USE_DB) {
		  sqlSave(fchannel, save_df, tablename=FORECAST_TABLE, append=T, rownames=F,verbose=F)
		} 
		
		imageFileName=paste(series,'.png',sep='')
		filePath=paste(IMAGES_DIR, imageFileName,sep='/')
		unlink(filePath)
		png(file=filePath, bg="white", width = imageWidth, height = imageHeight, pointsize=16)
		
		ymax=max(c(y,yy),max(avgYfs[3,]))
		ymin=min(min(c(y,yy)),min(avgYfs[1,]))
		plot(c(y,yy), main=series,  col=c(rep(1,n),rep(2,maxPredictionHorizon)), ylab='',xlab=lastDetails,ylim=c(ymin,ymax) )
		
		for (irun in 1:NUM_OF_TRIALS_TO_SHOW) 
			lines((n+1):(n+maxPredictionHorizon),yf[irun,], col='gray')
		
		lines((n+1):(n+maxPredictionHorizon),avgYfs[1,], col='pink',lwd=2)
		lines((n+1):(n+maxPredictionHorizon),avgYfs[2,], col='blue',lwd=2)
		lines((n+1):(n+maxPredictionHorizon),avgYfs[3,], col='pink',lwd=2)
		#lines((n+1):(n+maxPredictionHorizon),avgYfs[5,], col='red',lwd=2)
		#lines((n+1):(n+maxPredictionHorizon),avgYfs[4,], col='red',lwd=2)
		
    lines(c(y,yy), main=series,  col=c(rep(1,n),rep(2,maxPredictionHorizon)),type='b')
		
		dev.off() ;
		InsertGraph(imageFileName)
	} #proceed
	pred_df
}
ret_df$series=as.character(ret_df$series)
# str(ret_df)
write.csv(ret_df[!is.na(ret_df$lgtso),c("series","lgtso")], file=paste(DAY_DIR, '/','lgtso.csv',sep=''),row.names = F)


#repack; i=1
lgtso<-ets1 <- aarima <- hybrid <- matrix(NA,nrow=nseries,ncol=18)
for(i in 1:nseries) {
	if (i<10) {
		series=paste('N000',i,sep='')
	} else if (i<100) {
		series=paste('N00',i,sep='')
	} else if (i<1000) {
		series=paste('N0',i,sep='')
	} else {
		series=paste('N',i,sep='')
	} 
	oneSeries_df=ret_df[ret_df$series==series,]
	lgtso[i,]=oneSeries_df$lgtso
	ets1[i,]=oneSeries_df$ets
	aarima[i,]=oneSeries_df$aarima
	hybrid[i,]=oneSeries_df$hybrid
}



# Compute accuracy
mase <- mape <- smape <- wmape <- matrix(NA,nrow=11,ncol=nseries)
f <- matrix(NA, nrow=11, ncol=18)
for(i in 1:nseries)
{
	x <- M3[[i]]$xx
	n <- length(x)
	f[1,1:n] <- theta[i,1:n]
	f[2,1:n] <- fpro[i,1:n]
	f[3,1:n] <- fcx[i,1:n]
	f[4,1:n] <- bjauto[i,1:n]
	f[5,1:n] <- ab1[i,1:n]
	f[6,1:n] <- ab2[i,1:n]
	f[7,1:n] <- ab3[i,1:n]
	f[8,1:n] <- ets1[i,1:n]
	f[9,1:n] <- aarima[i,1:n]
	f[10,1:n] <- hybrid[i,1:n]
	f[11,1:n] <- lgtso[i,1:n]
	scale <- mean(abs(diff(M3[[i]]$x, lag=frequency(x))))
	for(j in 1:11)
	{
		mape[j,i] <- mean(abs((x-f[j,1:n])/x))*100
		wmape[j,i] <- sum(abs(x-f[j,1:n]))/sum(abs(x))*100
		smape[j,i] <- mean(abs(x-f[j,1:n])/(abs(x)+abs(f[j,1:n])))*200
		mase[j,i] <- mean(abs(x-f[j,1:n])/scale)
	}
}

# All series
NUM_OF_METRICS=4
NUM_OF_CATEGORIES=4
m3table <- matrix(NA, nrow=11, ncol=NUM_OF_METRICS+NUM_OF_METRICS*NUM_OF_CATEGORIES)
m3table[,1] <- rowMeans(mape,na.rm=TRUE)
m3table[,2] <- rowMeans(smape)
m3table[,3] <- rowMeans(wmape)
m3table[,3+1] <- rowMeans(mase)
m3table[,4+1] <- rowMeans(mape[,1:645],na.rm=TRUE)
m3table[,5+1] <- rowMeans(smape[,1:645])
m3table[,6+1] <- rowMeans(wmape[,1:645])
m3table[,6+2] <- rowMeans(mase[,1:645])
m3table[,7+2] <- rowMeans(mape[,(1+645):(645+756)],na.rm=TRUE)
m3table[,8+2] <- rowMeans(smape[,(1+645):(645+756)])
m3table[,9+2] <- rowMeans(wmape[,(1+645):(645+756)])
m3table[,9+3] <- rowMeans(mase[,(1+645):(645+756)])
m3table[,10+3] <- rowMeans(mape[,(1+645+756):(645+756+1428)],na.rm=TRUE)
m3table[,11+3] <- rowMeans(smape[,(1+645+756):(645+756+1428)])
m3table[,12+3] <- rowMeans(wmape[,(1+645+756):(645+756+1428)])
m3table[,12+4] <- rowMeans(mase[,(1+645+756):(645+756+1428)])
m3table[,13+4] <- rowMeans(mape[,(1+645+756+1428):3003],na.rm=TRUE)
m3table[,14+4] <- rowMeans(smape[,(1+645+756+1428):3003])
m3table[,15+4] <- rowMeans(wmape[,(1+645+756+1428):3003])
m3table[,15+5] <- rowMeans(mase[,(1+645+756+1428):3003])
rownames(m3table) <- c("Theta","ForecastPro","ForecastX","BJauto",
		"Autobox1","Autobox2","Autobox3",
		"ETS","AutoARIMA","Hybrid",'GT')
colnames(m3table) <- c(
		"MAPE","sMAPE","wMAPE","MASE",
		"MAPEy","sMAPEy","wMAPEy","MASEy",
		"MAPEq","sMAPEq","wMAPEq","MASEq",
		"MAPEm","sMAPEm","wMAPEm","MASEm",
		"MAPEo","sMAPEo","wMAPEo","MASEo")
j <- order(m3table[,3])
round(m3table[j,],2)




#output table creation script for SQL Server
#CREATE TABLE [dbo].[ForecastM3C](
#		[run] [varchar](64) NOT NULL,
#		[series] [varchar](20) NOT NULL,
#		[predictionHorizon] [int] NOT NULL,
#		[dateTimeOfPrediction] [datetime] NOT NULL,
#		[trueValue] [real] NULL,
#		[predQ1] [real] NULL,
#		[predQ5] [real] NULL,
#		[predQ50] [real] NULL,
#		[predQ95] [real] NULL,
#		[predQ99] [real] NULL,
#		[variable] [varchar](20) NOT NULL,
#		[algDetails] [varchar](150) NULL,
#		CONSTRAINT [ForecastM3C_pk] PRIMARY KEY CLUSTERED (
#				[run] ASC,
#				[series] ASC,
#				[predictionHorizon] ASC
#		)WITH (PAD_INDEX = OFF, STATISTICS_NORECOMPUTE = OFF, IGNORE_DUP_KEY = OFF, ALLOW_ROW_LOCKS = ON, ALLOW_PAGE_LOCKS = ON) ON [PRIMARY]
#) ON [PRIMARY]




