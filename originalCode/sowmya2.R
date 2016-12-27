# GT. weekly+weekly special(partial)+yearly seasonality
# separate sigma for s2
# Dec 2016

# install.packages(c("forecast"))
# install.packages(c("rstan","doParallel"))
# install.packages(c("sqldf","R2HTML"))
# install.packages(c("sn","RODBC"))

#library(forecast)
library(sn)
library(rstan)
library(RODBC)
library(sqldf)
library(doParallel)
library(R2HTML)

######################################################################

OUTPUT_DIR="e:/reports/Store"
FORECAST_TABLE='Store'
reportTitle='Store'
inputPath="E:\\progs\\sowmya\\sowmya.csv"
NUM_OF_PAST_POINTS_TO_SHOW=as.integer(365*6/12)
################################

PRODUCTION_MODE=FALSE

numOfCores=parallel:::detectCores()
if (numOfCores==8) {
	CLUSTER_SIZE=3
} else if (numOfCores>=16) {
	CLUSTER_SIZE=as.integer(numOfCores/2)
} else {
	CLUSTER_SIZE=as.integer(numOfCores/3)
} 
STAN_CLUSTER_SIZE=4
rstan_options(auto_write = TRUE)
options(mc.cores = STAN_CLUSTER_SIZE, width=180)

MIN_MEMORY_MB=3999
if (memory.limit()<MIN_MEMORY_MB) {
	memory.limit(MIN_MEMORY_MB)
}	
memory.limit()


MAX_RHAT_ALLOWED=1.005
NUM_OF_ITER=1000
MAX_NUM_OF_REPEATS=3
RUNMED_SIZE=0
STEP_BACK_WINDOW=0
MIN_DAYS_TO_MAKE_FORECAST=2*365
SHOW_ALG_DETAILS=T

MAX_NU=15 # if changing, change also inside models
MIN_NU=1.4
MIN_SIGMA=0.05
MIN_VAL=0.1
MAX_VAL=1e38 # roughly what SQL Server real can take
MIN_POW=-0.5
MAX_POW=1
LOC_TREND_FACT_ALPHA=1
LOC_TREND_FACT_BETA=1

POW_TREND_BETA=1
POWX_ALPHA=1

######################
CAUCHY_SD_DIV=200
MAX_SIZE_OF_LOOK_BACK_WINDOW=365*5;
mainPreProcessingAlgorithm=""

POWX_BETA=3
POW_TREND_ALPHA=2
Stan_predictionAlgorithm=paste('SSGT2_powTA',POW_TREND_ALPHA,'_powXB',POWX_BETA, sep='')

Sys.setenv(TZ = "GMT");
Sys.timezone()
today=Sys.Date();
START_DATE=as.Date('2013-7-1')
ANCHOR_STEP=30
MAX_FORECAST_HORIZON=30
now=Sys.time()
imageWidth=1300; imageHeight=600
NUM_OF_TRIALS=2000  
FREQ=1
NUM_OF_TRIALS_TO_SHOW=50
if (NUM_OF_TRIALS<NUM_OF_TRIALS_TO_SHOW)
  NUM_OF_TRIALS_TO_SHOW=NUM_OF_TRIALS
S2_SEASONALITY=7
QUANTS=c(5,50,95, 99)/100 #but you can't change them without also changeing code and the detination table structure 

parametersS = c("l", "s", "sSm", "s2", "s2Sm", "sy", "sySm",
		"nu", "levSm", "coefTrend", "powTrend",
		"sigma", "powx",  "offsetSigma",
		"sigma2", "powx2",  "offsetSigma2")
modelS = '
	data {  
    int<lower=1> NUM_OF_S2_DAYS;
		real<lower=0> CAUCHY_SD;
		real MIN_POW;  real MAX_POW;
		real<lower=0> MIN_SIGMA;
		real<lower=1> MIN_NU; real<lower=1> MAX_NU;
		int<lower=1> N;
		vector<lower=0>[N] y;
		int<lower=0,upper=1> S2_INDEX[N];
		real<lower=0> TREND_ALPHA; real<lower=0> TREND_BETA;
		real<lower=0> POWX_ALPHA;  real<lower=0> POWX_BETA;
		int<lower=2> S2_SEASONALITY; //must be multiple of 7
		vector[S2_SEASONALITY] INIT_S2U; #unnormalized
	}
  transformed data {
		int<lower=2> YEARLY_SEASONALITY;
		int<lower=1> NWS; #length of vector of yearly seasonalities

		YEARLY_SEASONALITY=52;
		NWS=N/7+1+YEARLY_SEASONALITY;
  }
	parameters {
		real<lower=MIN_NU,upper=MAX_NU> nu; 
		real<lower=0> sigma;
		real <lower=0,upper=1>powx;
		real <lower=MIN_SIGMA> offsetSigma;
		real<lower=0> sigma2;
		real <lower=0,upper=1>powx2;
		real <lower=MIN_SIGMA> offsetSigma2;
		real <lower=0,upper=1>levSm;
		real <lower=0,upper=1>sSm;
		real <lower=0,upper=1>s2Sm;
		real <lower=0,upper=1>sySm;
		real <lower=0,upper=1> powTrendBeta;
		real coefTrend;
		vector[7] initSu; #unnormalized
		vector[S2_SEASONALITY] initS2u; #unnormalized
		vector[YEARLY_SEASONALITY] initSyu; #unnormalized
	} 
	transformed parameters {
		real <lower=MIN_POW,upper=MAX_POW>powTrend;
		vector[N] l;

		vector[7] inits; 
		vector[N+7+1] s;
		
		vector[S2_SEASONALITY] inits2;
		vector[NUM_OF_S2_DAYS+S2_SEASONALITY+1] s2;
		
		vector[YEARLY_SEASONALITY] initSy;
		vector[NWS] sy;

		real sumsu; real sumy;

		sumsu = 0;
		for (i in 1:7) 
		  sumsu = sumsu+ initSu[i];
		for (i in 1:7) 
		  inits[i] = initSu[i]*7/sumsu;
		for (i in 1:7) {
		  s[i] = inits[i];
		}
		s[7+1] = inits[1];


		sumsu = 0;
		for (i in 1:S2_SEASONALITY) 
		  sumsu = sumsu+ initS2u[i];
		for (i in 1:S2_SEASONALITY) 
		  inits2[i] = initS2u[i]*S2_SEASONALITY/sumsu;
		for (i in 1:S2_SEASONALITY) {
		  s2[i] = inits2[i]*INIT_S2U[i];
		}
		
		sumsu = 0;
		for (i in 1:YEARLY_SEASONALITY) 
		  sumsu = sumsu+ initSyu[i];
		for (i in 1:YEARLY_SEASONALITY) 
		  initSy[i] = initSyu[i]*YEARLY_SEASONALITY/sumsu;
		for (i in 1:YEARLY_SEASONALITY) {
		  sy[i] = initSy[i];
		}


		if (S2_INDEX[1]==0) {
		  l[1] = y[1]/(s[1]*sy[1]);
    } else {
			l[1] = y[1]/(s2[1]*sy[1]);
		}
		powTrend= (MAX_POW-MIN_POW)*powTrendBeta+MIN_POW;

		{
			int is2; int iy;
			sumy=0; is2=0; iy=1;
			for (t in 2:N) {
				if (S2_INDEX[t]==0) {
			  	l[t]  = levSm*y[t]/(s[t]*sy[iy]) + (1-levSm)*l[t-1] ;  
			  	s[t+7] = sSm*y[t]/(l[t]*sy[iy])+(1-sSm)*s[t];
					sumy=sumy+y[t]/(l[t]*s[t]);
				} else {//this code will not work if the first day of the series (in-sample) is in inside the Thanksgiving week. This possibility needs be removed in the calling R code. 
					is2=is2+1;
					l[t]  = levSm*y[t]/(s2[is2]*sy[iy]) + (1-levSm)*l[t-1] ;  
			  	s2[is2+S2_SEASONALITY] = s2Sm*y[t]/(l[t]*sy[iy])+(1-s2Sm)*s2[is2];
					s[t+7] = s[t];
	        sumy=sumy+y[t]/(l[t]*s2[is2]);
				}
	      if ((t-1)%7==0) {
	        sy[iy+YEARLY_SEASONALITY] = sySm*sumy/7+(1-sySm)*sy[iy];
				  sumy=0;
				  iy=iy+1;
	      }
			}
		}
	}
	model {
		real expVal; int is2; int iy;

		powTrendBeta ~ beta(TREND_ALPHA, TREND_BETA);
		powx ~  beta(POWX_ALPHA, POWX_BETA);
		sigma ~ cauchy(0,CAUCHY_SD) T[0,];
		offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
		sigma2 ~ cauchy(0,5*CAUCHY_SD) T[0,];
		offsetSigma2 ~ cauchy(MIN_SIGMA,5*CAUCHY_SD) T[MIN_SIGMA,];
		coefTrend ~ cauchy(0, CAUCHY_SD);

		for (t in 1:7) 
		  initSu[t] ~ normal (1, 0.2) T[0.05,];
    for (t in 1:S2_SEASONALITY) 
		  initS2u[t] ~ normal (1, 0.2) T[0.1,];//this is just a disturbance on top of externally-provided initialization
		for (t in 1:YEARLY_SEASONALITY) 
		  initSyu[t]~ normal (1, 0.3) T[0.05,];
 
		is2=0;
		for (t in 2:N) {
			iy=(t-1)/7+1;
			if (S2_INDEX[t]==0) {
		  	expVal = (l[t-1]+ coefTrend*fabs(l[t-1])^powTrend)*s[t]*sy[iy];
				y[t] ~ student_t(nu, expVal, sigma*fabs(expVal)^powx+ offsetSigma);	
			} else {
				is2=is2+1;
				expVal = (l[t-1]+ coefTrend*fabs(l[t-1])^powTrend)*s2[is2]*sy[iy];
				y[t] ~ student_t(nu, expVal, sigma2*fabs(expVal)^powx2+ offsetSigma2);
			}
		}
	}
' 
stanModelS = stan_model(model_code=modelS)
#str(stanModelS)

##############################################
comp_df=read.csv(inputPath)
comp_df$OrderDate=as.Date(comp_df$OrderDate, format="%m/%d/%Y")
str(comp_df)

plot(comp_df$OrderDate, comp_df$TotalOrders, type='l')
which(diff(comp_df$OrderDate)!=1)
#973
comp_df$OrderDate[972:974]
add_df=comp_df[comp_df$OrderDate==comp_df$OrderDate[973-6],]
add_df$OrderDate=add_df$OrderDate+7
comp_df=rbind(comp_df,add_df)
comp_df=comp_df[order (comp_df$OrderDate),]
which(diff(comp_df$OrderDate)!=1)
plot(comp_df$OrderDate, comp_df$TotalOrders, type='l')

summary(comp_df$OrderDate)
dayOfYear_=as.numeric(format(comp_df$OrderDate, format = "%j"))
weekOfYear=(dayOfYear_-1)%/%7+1

LENGTH_OF_THANKS_GIVING_PERIOD=7
startingDatesOfThanksgivingPeriod=c(as.Date('2013-11-27'),as.Date('2014-11-26'), 
		as.Date('2015-11-25'), as.Date('2016-11-23'),as.Date('2017-11-22'))# UPDATE YEARLY
pastDates=sort(unique(comp_df$OrderDate))
futureDates=seq(from=max(pastDates)+1, to=as.Date('2018-11-1'), by='day') # UPDATE YEARLY

allDates=c(pastDates,futureDates)
thanksGivingWeeksIndex=allDates>=startingDatesOfThanksgivingPeriod[1] & allDates<startingDatesOfThanksgivingPeriod[1]+LENGTH_OF_THANKS_GIVING_PERIOD
for (i in 2:length(startingDatesOfThanksgivingPeriod)) {
	thanksGivingWeeksIndex=thanksGivingWeeksIndex |
			allDates>=startingDatesOfThanksgivingPeriod[i] & allDates<startingDatesOfThanksgivingPeriod[i]+LENGTH_OF_THANKS_GIVING_PERIOD
}			

#which(thanksGivingWeeksIndex)
thanksGivingIndexNumeric=as.numeric(thanksGivingWeeksIndex)
numOfThanksgivingPeriodDays=sum(thanksGivingIndexNumeric)
maxDateInDb=max(comp_df$OrderDate)

comp_df$zone='total'
locations=unique(comp_df$zone)

#estimate Thanksgiving seasonality 
indx1=comp_df$OrderDate<startingDatesOfThanksgivingPeriod[1] & comp_df$OrderDate>=startingDatesOfThanksgivingPeriod[1]-7 
indx2=comp_df$OrderDate>=startingDatesOfThanksgivingPeriod[1]+LENGTH_OF_THANKS_GIVING_PERIOD & comp_df$OrderDate<startingDatesOfThanksgivingPeriod[1]+7+LENGTH_OF_THANKS_GIVING_PERIOD
which(indx1 | indx2)
mean1=mean(comp_df$TotalOrders[indx1 | indx2])
tgIndx=comp_df$OrderDate>=startingDatesOfThanksgivingPeriod[1] & comp_df$OrderDate<startingDatesOfThanksgivingPeriod[1]+LENGTH_OF_THANKS_GIVING_PERIOD
which(tgIndx)
tgCoefs1=comp_df$TotalOrders[tgIndx]/mean1

indx1=comp_df$OrderDate<startingDatesOfThanksgivingPeriod[2] & comp_df$OrderDate>=startingDatesOfThanksgivingPeriod[2]-7 
indx2=comp_df$OrderDate>=startingDatesOfThanksgivingPeriod[2]+LENGTH_OF_THANKS_GIVING_PERIOD & comp_df$OrderDate<startingDatesOfThanksgivingPeriod[2]+7+LENGTH_OF_THANKS_GIVING_PERIOD
which(indx1 | indx2)
mean2=mean(comp_df$TotalOrders[indx1 | indx2])
tgIndx=comp_df$OrderDate>=startingDatesOfThanksgivingPeriod[2] & comp_df$OrderDate<startingDatesOfThanksgivingPeriod[2]+LENGTH_OF_THANKS_GIVING_PERIOD
which(tgIndx)
tgCoefs2=comp_df$TotalOrders[tgIndx]/mean2
tgCoefs=rowMeans(cbind(tgCoefs1, tgCoefs2))



################################
BeginReport = function(title) {
	cat(paste("<html><head><title>",
					title, "</title></head><link rel=stylesheet href=\"R2HTML.css\" type=text/css>",
					"<body bgcolor=#D0D0D0>", sep = ""),
			file = htmlFilePath, append = FALSE)
	cat(paste("<html><head><title>",
					title, "</title></head><link rel=stylesheet href=\"..\\R2HTML.css\" type=text/css>",
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
	relPath=paste(today,'_',predictionAlgorithm,'_',preProcessingAlgorithm,'\\images\\',imageFileName, sep='')
	HTMLInsertGraph(relPath, Caption = '', WidthHTML=imageWidth, HeightHTML=imageHeight, file=htmlFilePath)
	HTMLInsertGraph(relPath2, Caption = '', WidthHTML=imageWidth, HeightHTML=imageHeight, file=htmlFilePath2)
}
WriteHTML=function(html) {
	cat(html," \n", file = htmlFilePath, append = TRUE)
	cat(html," \n", file = htmlFilePath2, append = TRUE)
}


plotQ=function (item, itemKeys,
		pastDates, yPast, futureDates, forecastVect, knownFutureValues=NULL,
		quants=NULL, quantileMatrix=NULL, pathsMatrix=NULL,
		indxOfQuantToDisp=2,  xlabel='', ylabel='', logAxis='') {
	
	numOfColors=length(quants)
	#kolorki=heat.colors(numOfColors, alpha = 0.2)
	kolorki=rainbow(numOfColors, alpha = 0.2, start=0, end=0.4)
	if (quants[1]<quants[length(quants)]) {
		kolorki=rev(kolorki)
	}
	
	quantToDisp=NULL
	maxy=max(yPast, forecastVect, na.rm=T)
	if (!is.null(quantileMatrix)) {
		if (nrow(quantileMatrix)<indxOfQuantToDisp) {
			indxOfQuantToDisp=1
		}
		quantToDisp=quants[indxOfQuantToDisp]
		maxy=max(maxy, quantileMatrix[indxOfQuantToDisp,], na.rm=T)
	}
	miny=min(yPast, forecastVect, na.rm=T)
	
	mtitle=paste(max(pastDates),": ", item,", ",itemKeys,
			', expected:',as.integer(forecastVect[length(forecastVect)]),sep='')
	if (!is.null(quantToDisp)) {
		mtitle=paste(mtitle,', ',quantToDisp,'p:',
				as.integer(quantileMatrix[indxOfQuantToDisp,length(quantileMatrix[indxOfQuantToDisp, ])]),sep='')
	}
	
	plot(pastDates, yPast, 
			xlim=c(min(pastDates),max(futureDates)),ylim=c(miny,maxy),
			main=mtitle,
			xlab=xlabel, ylab=ylabel, type='l' ,lwd=2, log=logAxis)
	
	if (!is.null(pathsMatrix)) {
		for (iq in 1:nrow(pathsMatrix)) {
			lines(futureDates,pathsMatrix[iq,], col='gray')	
		}
	}
	
	lines(futureDates,forecastVect, col='blue',lwd=2) #pred
	if (!is.null(knownFutureValues) && length(knownFutureValues)>0) {
		lines(futureDates[1:length(knownFutureValues)],knownFutureValues,lwd=2) #known future values
	}
	
	if (!is.null(quantileMatrix)) {
		for (iq in 1:nrow(quantileMatrix)) {
			lines(futureDates,quantileMatrix[iq,], col=kolorki[iq],lwd=2)
		}
	}
	if (quants[1]<quants[length(quants)]) {
		legend("topleft", col = c('black', 'blue', rev(kolorki)), bty="n",
				legend = c("Time series history", "Expected", paste(rev(quants)*100,'p',sep='')), 
				lty = c(1,1, rep(NA,length(quants))), 
				lwd = c(2,2, rep(NA,length(quants))), 
				pch = c(NA, NA, rep(15,length(quants))), pt.cex = 2.5)
		
	} else {
	  legend("topleft", col = c('black', 'blue', kolorki), bty="n",
			legend = c("Time series history", "Expected", paste(quants*100,'p',sep='')), 
			lty = c(1,1, rep(NA,length(quants))), 
			lwd = c(2,2, rep(NA,length(quants))), 
			pch = c(NA, NA, rep(15,length(quants))), pt.cex = 2.5)
  }
}


############################################################################
preProcessingAlgorithm=mainPreProcessingAlgorithm
predictionAlgorithm=Stan_predictionAlgorithm
DAY_DIR=paste(OUTPUT_DIR,paste(today,predictionAlgorithm,sep='_'),sep='/')
stanOutputPath=paste(DAY_DIR, 'stanOutput.txt',sep='/')
IMAGES_DIR=paste(DAY_DIR,'images',sep='/')

dir.create(IMAGES_DIR, recursive=T)
htmlFilePath=paste(OUTPUT_DIR, '/',reportTitle,'.html',sep='')
htmlFilePath2=paste(DAY_DIR, '/',reportTitle,'.html',sep='')
unlink(htmlFilePath)
unlink(htmlFilePath2)

title_=paste(reportTitle,'created on',Sys.Date(), '.  Max date in db=',maxDateInDb)
BeginReport(title_)
WriteHTML(title_)
WriteHead(reportTitle)

###################################################################################################################
fchannel= odbcConnect('slawek')
odbcGetInfo(fchannel)

print(Stan_predictionAlgorithm)
stopCluster(cl)
Sys.sleep(2)
CLUST_OUT_PATH="c:/temp/rclust.txt"
file.remove(CLUST_OUT_PATH)
cl = makeCluster(CLUSTER_SIZE, outfile=CLUST_OUT_PATH)
registerDoParallel(cl)


sink(); graphics.off();
unlink(stanOutputPath)

firstTime=TRUE; idr=1
#foreach(idr=1:length(locations), .packages=c("sn","rstan","RODBC","R2HTML"), .inorder=F) %dopar% {
#if (firstTime==TRUE) {
#	fchannel= odbcConnect('slawek')
#	odbcGetInfo(fchannel)
#	sink(file=stanOutputPath, append=TRUE, split=TRUE)
#	firstTime=FALSE
#}
for(idr in 1:length(locations)) {
  location=as.character(locations[idr])  
  oneLocation_df=comp_df[comp_df$zone==location,]

	if (nrow(oneLocation_df)>MIN_DAYS_TO_MAKE_FORECAST/FREQ) {
    maxAnchorDate=max(oneLocation_df$OrderDate)
		firstAnchorDate=min(oneLocation_df$OrderDate)+MIN_DAYS_TO_MAKE_FORECAST
		allPastDates=sort(unique(oneLocation_df$OrderDate))
		
		prevPredMaxDate_query=paste("select max(anchorDate) maxAnchorDate from ",FORECAST_TABLE,
				" where zone='",location,
				"' and predictionAlgorithm='",predictionAlgorithm,
				"'",  sep='')
		prevPredMaxDate_df=sqlQuery(fchannel,prevPredMaxDate_query, stringsAsFactors =F)
		maxDateInDB=prevPredMaxDate_df$maxAnchorDate[1]
		if (!is.na(maxDateInDB)) {
			maxDateInDB=as.Date(maxDateInDB)
		}
    
    
		if (PRODUCTION_MODE) {
			anchorDates=maxAnchorDate
		} else {
			if (!is.na(maxDateInDB)) {
				if (maxAnchorDate>maxDateInDB) {
					anchorDates=rev(seq(from=maxAnchorDate, to=maxDateInDB, by=-ANCHOR_STEP))	
					if (maxDateInDB %in% anchorDates) {
						anchorDates=anchorDates[-1]
					}
				} else {
					anchorDates=NULL
				}
			} else {#nothing in db
				if (maxAnchorDate>firstAnchorDate) {
					anchorDates=rev(seq(from=maxAnchorDate, to=firstAnchorDate, by=-ANCHOR_STEP))	
				} else {
					anchorDates=NULL
				}
			}
		}
			
    if (!is.null(anchorDates)) {
			ysmooth=oneLocation_df$TotalOrders
			if (min(ysmooth)<1) {
				ysmooth=ysmooth+1+min(ysmooth)
			}
			summary(ysmooth)
			plot(oneLocation_df$OrderDate,ysmooth, main=paste(location), type='l')
			
			ianchDate=1
			#foreach(ianchDate = 1:length(anchorDates), .packages=c("sn","rstan","RODBC","R2HTML"), .inorder=F) %dopar% {
	    for (ianchDate in 1:length(anchorDates)) {
				if (firstTime==TRUE) {
					fchannel= odbcConnect('slawek')
					odbcGetInfo(fchannel)
					#sink(file=stanOutputPath, append=TRUE, split=TRUE)
					firstTime=FALSE
				}
				# anchorDate=as.Date('2016-11-30')
	      anchorDate=anchorDates[ianchDate]
        print("")
        print(paste(Sys.time(),"starting ", location, "on", anchorDate))
         
				pastDatesIndx=oneLocation_df$OrderDate<=anchorDate & oneLocation_df$OrderDate>anchorDate-MAX_SIZE_OF_LOOK_BACK_WINDOW 
				pastDates=oneLocation_df$OrderDate[pastDatesIndx]
				
				futureDatesIndex=oneLocation_df$OrderDate>anchorDate & oneLocation_df$OrderDate<=anchorDate+MAX_FORECAST_HORIZON
				yFuture=ysmooth[futureDatesIndex]
				futureDates=oneLocation_df$OrderDate[futureDatesIndex]
				if (length(futureDates)<MAX_FORECAST_HORIZON) {#futureDates are used for plotting
					futureDates=allDates[allDates>anchorDate & allDates<=anchorDate+MAX_FORECAST_HORIZON]
				}
				
				y=ysmooth[pastDatesIndx]
				#y=y+rnorm(length(y),mean=0, sd=min(abs(y))*0.0001) #jitter added
				n = length(y)
				s2Index=thanksGivingIndexNumeric[1:n]
				s2FutureIndex=thanksGivingIndexNumeric[(n+1):(n+MAX_FORECAST_HORIZON)]
				#lines(anchorDate,y[n], type='p',col='red') #on screen
        
	      CauchySd=max(y)/CAUCHY_SD_DIV
				print(paste("CauchySd",round(CauchySd,2),sep=': '))
	      data = list(NUM_OF_S2_DAYS=numOfThanksgivingPeriodDays,
					CAUCHY_SD=CauchySd, S2_INDEX=s2Index,
	      	MIN_POW=MIN_POW, MAX_POW=MAX_POW, 
  				MIN_SIGMA=MIN_SIGMA, MIN_VAL=MIN_VAL, 
    			MIN_NU=MIN_NU,  MAX_NU=MAX_NU, 
					POWX_ALPHA=POWX_ALPHA, POWX_BETA=POWX_BETA,
					TREND_ALPHA=POW_TREND_ALPHA, TREND_BETA=POW_TREND_BETA,
					S2_SEASONALITY=S2_SEASONALITY, INIT_S2U=tgCoefs,
	        y=y, N=n) # to be passed on to Stan
        
				avgRHat=1e20; irep=1
				for (irep in 1:MAX_NUM_OF_REPEATS) {
					initializations = list();
					for (irr in 1:CLUSTER_SIZE) {
						initializations[[irr]]=list(
							initSu=rnorm(7,1,0.1),
						  initS2u=rnorm(S2_SEASONALITY,1,0.1),	
						  initSyu=rnorm(52,1,0.1)
						)
					}
					
					samples1=
							sampling(stanModelS,   
									data=data, 
									init=initializations,
									pars=parametersS,
									iter=NUM_OF_ITER*2^(irep-1),
									chains=STAN_CLUSTER_SIZE,
									cores=STAN_CLUSTER_SIZE,
									open_progress=F,
									refresh = 100)
					# print(samples1)
					
					ainfo=summary(samples1)
					RHats=ainfo$summary[,10]
					RHats=as.numeric(RHats[is.finite(RHats)])
					currRHat=mean(RHats, na.rm=T)
					print(samples1)
					print(paste("currRHat",currRHat))
					if (currRHat<=MAX_RHAT_ALLOWED) {
						samples=samples1
						avgRHat=currRHat
						break
					} else {
						if (currRHat<avgRHat) {#in the worst case this is at least once executed, because avgRHat is initialized high
							samples=samples1
							avgRHat=currRHat
						} else {
							print ("worse...")
							print(paste("currRHat",currRHat))
						}
						print (paste("trying to do better...",location))
					}
					#str(samples, max.level =4)
				}#repeat if needed
				
      
        l=extract(samples)$l
        mc_size=dim(l)[1]
        mc_range=1:mc_size
        lastDataIndex=dim(l)[2]
        lastLevel=l[,lastDataIndex]
        #str(lastLevel); summary(lastLevel)
        lastLevelM=mean(lastLevel)
	      paste("lastLevel",lastLevelM)
	      
				sSm <- extract(samples)$sSm
				#summary(sSm)
				sSmM=mean(sSm)
				paste("sSm",sSmM)
				
				s=extract(samples)$s
				sM=apply(s,2,mean)
				#summary(sM)
        lastValidPositionIns=max(which(sM>0.01))
				cat("s:"); print(sM[1:lastValidPositionIns])

				s2Sm <- extract(samples)$s2Sm
				#summary(s2Sm)
				s2SmM=mean(s2Sm)
				paste("s2Sm",s2SmM)
				
				s2=extract(samples)$s2
				s2M=apply(s2,2,mean)
				#summary(sM)
				lastValidPositionIns2=max(which(s2M>0.01))
				cat("s2:"); print(s2M[1:lastValidPositionIns2])
				
				sySm <- extract(samples)$sySm
				#summary(sySm)
				sySmM=mean(sySm)
				paste("sySm",sySmM)
				
				sy=extract(samples)$sy
				syM=apply(sy,2,mean)
				#summary(syM)
				lastValidPositionInsy=max(which(syM>0.01))
				cat("sy:"); print(syM[1:lastValidPositionInsy])
				
	      sigma = extract(samples)$sigma
        #summary(sigma)
        sigmaM=mean(sigma)
        paste("sigma",sigmaM)

        powx = extract(samples)$powx
        #summary(powx)
        powxM=mean(powx)
        paste("powx",powxM)

  			offsetSigma=extract(samples)$offsetSigma
				#summary(offsetSigma)
				offsetSigmaM=mean(offsetSigma)
				paste("offsetSigma",round(offsetSigmaM,2))

				sigma2 = extract(samples)$sigma2
				#summary(sigma)
				sigma2M=mean(sigma2)
				paste("sigma2",sigma2M)
				
				powx2 = extract(samples)$powx2
				#summary(powx)
				powx2M=mean(powx2)
				paste("powx2",powx2M)
				
				offsetSigma2=extract(samples)$offsetSigma2
				#summary(offsetSigma)
				offsetSigma2M=mean(offsetSigma2)
				paste("offsetSigma2",round(offsetSigma2M,2))
				
        powTrend = extract(samples)$powTrend
				#summary(powTrend)
				powTrendM=mean(powTrend)
				paste("powTrend",round(powTrendM,2))

				coefTrend = extract(samples)$coefTrend
				#summary(coefTrend)
				coefTrendM=mean(coefTrend)
				paste("coefTrend",round(coefTrendM,2))
	
        muM=0
        #print(paste("mu",muM))
        
        levSm = extract(samples)$levSm
        #summary(levSm)
        levSmM=mean(levSm)
        paste("levSm",levSmM)

        nu = extract(samples)$nu
        #summary(nu)
        nuM=mean(nu)
        paste("nu",nuM)
				
				print(paste(Sys.time(),location,anchorDate))
				lastDetails=paste("coefT=",round(coefTrendM,2), ", powT=",round(powTrendM,2),
  				", sigma=",round(sigmaM,2),", powx=",round(powxM,2),", offsetS=",round(offsetSigmaM,2),
					", sigma2=",round(sigma2M,2),", powx2=",round(powx2M,2),", offsetS2=",round(offsetSigma2M,2), 
					", nu=",round(nuM,2), ", lSm=",round(levSmM,2), 
    			", sSm=",round(sSmM,2), ", s2Sm=",round(s2SmM,2),
					", sySm=",round(sySmM,2), sep='')
				print(lastDetails)
        
				t=1; irun=1
				sSRecent=rep(1,7+MAX_FORECAST_HORIZON)
				yf=matrix(lastLevelM,nrow=NUM_OF_TRIALS, ncol=MAX_FORECAST_HORIZON)
				for (irun in 1:NUM_OF_TRIALS) {
					index=sample(mc_range,1)
					prevLevel=lastLevel[index]
					powTrendS=powTrend[index]
					coefTrendS=coefTrend[index]
					levSmS=levSm[index]
					sSRecent[1:7]=s[index,(lastValidPositionIns-7+1):lastValidPositionIns]
					s2SRecent=s2[index,(lastValidPositionIns2-S2_SEASONALITY+1):lastValidPositionIns2]# assumption here is that the forecast horizon<1year
					sySRecent=sy[index,(lastValidPositionInsy-52+1):lastValidPositionInsy]  # assumption here is that the forecast horizon<1year
					nuS=nu[index]
					sigmaS=sigma[index]
					powxS=powx[index]
					offsetSigmaS=offsetSigma[index]
					sigma2S=sigma2[index]
					powx2S=powx2[index]
					offsetSigma2S=offsetSigma2[index]
					
					is2=0;t=1
					for (t in 1:MAX_FORECAST_HORIZON) {
						iy=t%/%7+1;
						if (s2FutureIndex[t]==0) {
							seasonA=sSRecent[t]*sySRecent[iy]
							expVal <- (prevLevel+ coefTrendS*abs(prevLevel)^powTrendS)*seasonA;
							omega=sigmaS*abs(expVal)^powxS+offsetSigmaS
						} else {
							is2=is2+1;
							seasonA=s2SRecent[is2]*sySRecent[iy]
							expVal <- (prevLevel+ coefTrendS*abs(prevLevel)^powTrendS)*seasonA;
							omega=sigma2S*abs(expVal)^powx2S+offsetSigma2S
						}
						error=rst(n=1, xi=0, omega= omega, alpha=0, nu=nuM)
						yf[irun,t]=max(MIN_VAL,expVal+error)
						currLevel=max(MIN_VAL,levSmM*yf[irun,t]/seasonA + (1-levSmM)*prevLevel) ;
		
						sSRecent[t+7] <- sSRecent[t];
						if (currLevel>MIN_VAL) {
							prevLevel=currLevel
						} 
					}
					# yf[irun,]
				}
				#	yf
				#}
        
				
        avgYfs=apply(yf,2,quantile,probs=QUANTS)	
				forecastVect=pmax(MIN_VAL,pmin(MAX_VAL,avgYfs[2,]))	
				
				
				#versionNumber=paste(gsub('-','',as.character(anchorDate),fixed=T),Version_Number,sep='')
				save_df=NULL; predictionHorizon=1; currMonth=-1
				for (predictionHorizon in (1:MAX_FORECAST_HORIZON)) {
					#predictionAlgorithm, zone, anchorDate, targetDate, trueValue, predQ5, predQ50, predQ95, predQ99, dateTimeOfPrediction, algDetails
					targetDate=anchorDate+predictionHorizon*FREQ
					newMonth=as.POSIXlt(targetDate)$mon+1
					#if (newMonth!=currMonth) {
	          if (is.null(save_df)){
	            save_df=data.frame(predictionAlgorithm=predictionAlgorithm, zone=location,
									anchorDate=anchorDate, targetDate=targetDate, trueValue=NA, 
									predQ5=avgYfs[1,predictionHorizon], predQ50=avgYfs[2,predictionHorizon],  
									predQ95=avgYfs[3,predictionHorizon],predQ99=avgYfs[4,predictionHorizon],
									dateTimeOfPrediction=now, algDetails=lastDetails)
	          } else {
	            save_df=rbind(save_df, data.frame(predictionAlgorithm=predictionAlgorithm, zone=location,
											anchorDate=anchorDate, targetDate=targetDate, trueValue=NA, 
											predQ5=avgYfs[1,predictionHorizon], predQ50=avgYfs[2,predictionHorizon],  
											predQ95=avgYfs[3,predictionHorizon],predQ99=avgYfs[4,predictionHorizon],
											dateTimeOfPrediction=now, algDetails=lastDetails))
						}
						currMonth=newMonth
	        #}
				}
        sqlSave(fchannel, save_df, tablename=FORECAST_TABLE, append=T, rownames=F,verbose=F)  
      
      
				#intermediate plots
				imageFileNam=paste(gsub(' ','_',location, fixed=T),anchorDate, sep='_')
				imageFileName=paste(imageFileNam,'.png',sep='')
				filePath=paste(IMAGES_DIR, imageFileName,sep='/')
				unlink(filePath)
				png(file=filePath, bg="white", width = imageWidth, height = imageHeight, pointsize=16)
				
				plotQ(reportTitle, location, 
						pastDates[(length(pastDates)-NUM_OF_PAST_POINTS_TO_SHOW+1):length(pastDates)], 
						y[(length(pastDates)-NUM_OF_PAST_POINTS_TO_SHOW+1):length(pastDates)], 
						futureDates, forecastVect, 
						knownFutureValues=yFuture,
						quants=QUANTS[-2], quantileMatrix=avgYfs[-2,], 
						pathsMatrix=yf[sample(nrow(yf),NUM_OF_TRIALS_TO_SHOW),],
						indxOfQuantToDisp=3, xlabel=lastDetails, ylabel='TotalOrders', logAxis='y')
				
        dev.off() ;
        InsertGraph(imageFileName)
			}
      
			#update trueValue 
			irow=1
			trueValueNotFilledQuery=paste("select distinct targetDate from ",FORECAST_TABLE,
					" where zone='",location,"'  and targetDate<='",max(oneLocation_df$OrderDate),
					"' and predictionAlgorithm='",predictionAlgorithm,        
					"' and trueValue is null 
							group by targetDate order by targetDate", sep='')
			trueValueNotFilled_df=sqlQuery(fchannel,trueValueNotFilledQuery,rows_at_time = 1000, stringsAsFactors =F )
			for (irow in 1:nrow(trueValueNotFilled_df)) {
				targetDate=trueValueNotFilled_df[irow,1]
				matchingDate=oneLocation_df$OrderDate[oneLocation_df$OrderDate<=targetDate & targetDate<(oneLocation_df$OrderDate+FREQ)]
				if (!is.na(matchingDate) && length(matchingDate)>0) {
					trueValue=oneLocation_df[oneLocation_df$OrderDate==matchingDate,"TotalOrders"]
					trueValueUpdateQuery=paste("update ",FORECAST_TABLE," set trueValue=",trueValue,
							" where zone='",location,
							"' and targetDate='",targetDate,
							"' and predictionAlgorithm='",predictionAlgorithm, 
							"'", sep='')
					sqlQuery(fchannel,trueValueUpdateQuery)
				}
			}
		}#anchor dates not null 
	} #enough days

  gc()
} #through locations
EndReport()
sink()
print(Stan_predictionAlgorithm)



#CREATE TABLE [dbo].[Store](
#		[predictionAlgorithm] [varchar](64) NOT NULL,
#		[zone] [varchar](40) NOT NULL,
#		[anchorDate] [date] NOT NULL,
#		[targetDate] [date] NOT NULL,
#		[trueValue] [real] NULL,
#		[predQ5] [real] NULL,
#		[predQ50] [real] NULL,
#		[predQ95] [real] NULL,
#		[predQ99] [real] NULL,
#		[dateTimeOfPrediction] [datetime] NOT NULL,
#		[algDetails] [varchar](640) NULL,
#		CONSTRAINT [store_pk] PRIMARY KEY CLUSTERED 
#				(
#				[predictionAlgorithm] ASC,
#				[zone] ASC,
#				[anchorDate] ASC,
#				[targetDate] ASC
#		))
