# prod Driver, unified
# Dec 2016
# Slawek
#This contains a few help functions related to plotting and the HTML report creation
#and the main script that does either production mode run (use latest date) or 
#full backtest, depending on PRODUCTION_MODE variable

library(sn)
library(rstan)
library(RODBC)
library(sqldf)
library(doParallel)
library(R2HTML)
library(csdfutils)
library(forecast)

######################################################################

## Parameter loading from Database
#configDF <- Fetch.Configuration("TPS_Midterm")
#if(file.exists("constants.R")){
#	file.remove("constants.R")
#}
#Save.Configuration(configDF, "constants.R")
#if(file.exists("constants.R")){
#	source("constants.R")
#}else{
#	Log.Error("Parameter is not loaded successfully!")
#	# TODO (yopyon) : terminate the R process
#}

############################################################################
#HTML report functions
BeginReport = function(title, htmlFilePath, htmlFilePath2) {
	cat(paste("<html><head><title>",
					title, "</title></head><link rel=stylesheet href=\"R2HTML.css\" type=text/css>",
					"<body bgcolor=#D0D0D0>", sep = ""),
			file = htmlFilePath, append = FALSE)
	cat(paste("<html><head><title>",
					title, "</title></head><link rel=stylesheet href=\"..\\R2HTML.css\" type=text/css>",
					"<body bgcolor=#D0D0D0>", sep = ""),
			file = htmlFilePath2, append = FALSE)
}
EndReport = function(htmlFilePath, htmlFilePath2) {
	cat("<hr size=1></body></html> \n ",
			file = htmlFilePath, append = TRUE)
	cat("<hr size=1></body></html> \n ",
			file = htmlFilePath2, append = TRUE)
}
WriteHead = function(head, htmlFilePath, htmlFilePath2) {
	cat("</p><br><br><h2>",head,"</h2>\n ",
			file = htmlFilePath, append = TRUE)
	cat("</p><br><br><h2>",head,"</h2>\n ",
			file = htmlFilePath2, append = TRUE)
}
InsertGraph=function(imageFileName, htmlFilePath, htmlFilePath2) {
	relPath2=paste('images\\',imageFileName, sep='')
	relPath=paste(today,'_',predictionAlgorithm,'\\images\\',imageFileName, sep='')
	HTMLInsertGraph(relPath, Caption = '', WidthHTML=imageWidth, HeightHTML=imageHeight, file=htmlFilePath)
	HTMLInsertGraph(relPath2, Caption = '', WidthHTML=imageWidth, HeightHTML=imageHeight, file=htmlFilePath2)
}
WriteHTML=function(html, htmlFilePath, htmlFilePath2) {
	cat(html," \n", file = htmlFilePath, append = TRUE)
	cat(html," \n", file = htmlFilePath2, append = TRUE)
}

#main plot
#plotQ(Forecast_Name, item, pastDates, y, futureDates, forecastVect, 
#  						knownFutureValues=yFuture,
#  						quants=QUANTS[-indexOfMedian], quantileMatrix=avgYfs[-indexOfMedian,], 
#  						pathsMatrix=pathsMatrix,
#  						indxOfQuantToDisp=indexOfHigh, xlabel=lastDetails, ylabel=VAR, logAxis=LOG_FOR_DISPLAY)
#item=Forecast_Name;  itemKeys=item; 	yPast=y; knownFutureValues=yFuture;
#quants=QUANTS[-indexOfMedian]; quantileMatrix=avgYfs[-indexOfMedian,]; pathsMatrix=pathsMatrix
#indxOfQuantToDisp=indexOfHigh; 
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
				round(quantileMatrix[indxOfQuantToDisp,length(quantileMatrix[indxOfQuantToDisp, ])]),sep='')
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

#####################
Log.Info(paste("starting",predictionAlgorithm,"forecast"))

indexOfMedian=which(PERCENTILES==50)
indexOfHigh=which(PERCENTILES==95)
if (length(indexOfHigh)==0) {
	indexOfHigh=length(PERCENTILES)
}
#the code below assumes that the high percentiles are more numerous and importan. Also that low percentiles are subset of mirrors to high percentiles
coverageLevels=PERCENTILES[PERCENTILES>50]*2-100
highPercentiles=PERCENTILES[PERCENTILES>50]
lowPercentiles=PERCENTILES[PERCENTILES<50]
lowerTailsInForecast=50-coverageLevels/2

positionsOfLowerTailsInForecastLower=list()
for (lpc in lowPercentiles) {
	position=which(lpc == lowerTailsInForecast)
	if (length(position)>0 && position>0) {
		positionsOfLowerTailsInForecastLower[[length(positionsOfLowerTailsInForecastLower)+1]] =position
	}	
}
positionsOfLowerTailsInForecastLower=as.numeric(positionsOfLowerTailsInForecastLower)

maxDateInDb=max(bigSum_df$date)
DAY_DIR=paste(OUTPUT_DIR,paste(today,predictionAlgorithm,sep='_'),sep='/')
stanOutputPath=paste(DAY_DIR, 'stanOutput.txt',sep='/')
IMAGES_DIR=paste(DAY_DIR,'images',sep='/')

dir.create(IMAGES_DIR, recursive=T)
htmlFilePath=paste(OUTPUT_DIR, '/',reportTitle,'.html',sep='')
htmlFilePath2=paste(DAY_DIR, '/',reportTitle,'.html',sep='')
unlink(htmlFilePath)
unlink(htmlFilePath2)

title_=paste(reportTitle,'created on',Sys.Date(), '.  Max date in db=',maxDateInDb)
BeginReport(title_, htmlFilePath, htmlFilePath2)
WriteHTML(title_, htmlFilePath, htmlFilePath2)
WriteHead(reportTitle, htmlFilePath, htmlFilePath2)


tryCatch(
	stopCluster(cl) ,
	error = function (e) print("OK, cluster not initialized")
)
Sys.sleep(2)
if(file.exists(CLUST_OUT_PATH)){
	file.remove(CLUST_OUT_PATH)
}
cl <- makeCluster(CLUSTER_SIZE, outfile=CLUST_OUT_PATH)
registerDoParallel(cl)

sink(); graphics.off();
if(file.exists(stanOutputPath)){
	unlink(stanOutputPath)
}

items=unique(sort(bigSum_df$item))
firstTime=TRUE; idr=1
foreach(idr=1:length(items), .packages=c("sn","rstan","RODBC","R2HTML", "csdfutils","forecast"), .inorder=FALSE) %dopar% {
  #for(idr in 1:length(items)) { #disable the above and enable this line for debugging
	if (firstTime==TRUE) {
		fchannel<- Fetch.DBConnection(SQL.SERVER.NAME, DATABASE.NAME, Fetch.UID(), Fetch.Passwd())
		odbcGetInfo(fchannel)
		#sink(file=stanOutputPath, append=TRUE, split=TRUE) # TODO(Slawek): Isn't this redundant?
		firstTime=FALSE
	}
	item=as.character(items[idr])  
	oneItem_df=bigSum_df[bigSum_df$item==item,]
	
	if (nrow(oneItem_df)>=MIN_NUMBER_OF_DATA_POINTS_MAKE_FORECAST) {
		maxAnchorDate=max(oneItem_df$date)
		firstAnchorDate=min(oneItem_df$date)+MIN_NUMBER_OF_DATA_POINTS_MAKE_FORECAST*FREQ
		allPastDates=sort(unique(oneItem_df$date))
		maxDateInDB=getMaxDateInDB(oneItem_df) #in reader
		
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
			yall=oneItem_df$value
			plot(oneItem_df$date,yall, main=paste(item), type='l')
			# ianchDate=1
			for (ianchDate in 1:length(anchorDates)) {
				anchorDate=anchorDates[ianchDate]
				print("")
				print(paste(Sys.time(),"starting ", item,"on", anchorDate))
				
				pastDatesIndx=oneItem_df$date<=anchorDate & oneItem_df$date>anchorDate-MAX_SIZE_OF_LOOK_BACK_WINDOW 
				pastDates=oneItem_df$date[pastDatesIndx]
				
				maxDateBeforeAnchorDate=max(oneItem_df$date[oneItem_df$date<=anchorDate])
				futureDates=seq(from=(maxDateBeforeAnchorDate+FREQ), to=(maxDateBeforeAnchorDate+MAX_FORECAST_HORIZON*FREQ),by=FREQ)
				yFuture=yall[oneItem_df$date%in%futureDates]
				
				if (length(pastDates)>=MIN_NUMBER_OF_DATA_POINTS_MAKE_FORECAST) {
					y=yall[pastDatesIndx]
					y=y+rnorm(length(y),mean=0, sd=min(abs(y))*0.0001) #jitter added
					n <- length(y)
					lines(anchorDate,y[n], type='p',col='red') #on screen
  				if (length(pastDates)>=MIN_DAYS_TO_MAKE_FORECAST/FREQ) {
  					CauchySd=max(y)/CAUCHY_SD_DIV
  					Log.Info(paste("CauchySd",round(CauchySd,2),sep=': '))
  					data <- list(CAUCHY_SD=CauchySd, SEASONALITY=SEASONALITY,
  							MIN_POW=MIN_POW, MAX_POW=MAX_POW, 
  							MIN_SIGMA=MIN_SIGMA, MIN_VAL=MIN_VAL, 
  							MIN_NU=MIN_NU,  MAX_NU=MAX_NU, 
  							POWX_ALPHA=POWX_ALPHA, POWX_BETA=POWX_BETA,
  							POW_TREND_ALPHA=POW_TREND_ALPHA, POW_TREND_BETA=POW_TREND_BETA,
								y=y, N=n, SKEW=SKEW) # to be passed on to Stan
  					
  					avgRHat=1e20; irep=1
  					for (irep in 1:MAX_NUM_OF_REPEATS) {
							initializations <- list();
							for (irr in 1:STAN_CLUSTER_SIZE) {
								initializations[[irr]]=list( 
										initSu=rnorm(SEASONALITY,1,0.1)
								)
							}
  						
  						samples1=
  								sampling(stanModel,   
  										data=data, 
											init=initializations,
  										pars=parameters,
  										iter=NUM_OF_ITER*2^(irep-1),
  										chains=STAN_CLUSTER_SIZE,
  										cores=STAN_CLUSTER_SIZE,
  										open_progress=F,
  										refresh = 1000,
                      control=list(adapt_delta = 0.9))
  						
  						ainfo=summary(samples1)
  						RHats=ainfo$summary[,10]
  						RHats=as.numeric(RHats[is.finite(RHats)])
  						currRHat=mean(RHats, na.rm=T)
  						if (!PRODUCTION_MODE)	print(samples1)  
  						Log.Info(paste("currRHat",currRHat))
  						if (currRHat<=MAX_RHAT_ALLOWED) {
  							samples=samples1
  							avgRHat=currRHat
  							break
  						} else {
  							if (currRHat<avgRHat) {#in the worst case this is at least once executed, because avgRHat is initialized high
  								samples=samples1
  								avgRHat=currRHat
  							} else {
  								Log.Warn ("worse...")
  								Log.Warn(paste("currRHat",currRHat))
  							}
  							Log.Warn(paste("trying to do better...",item))
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
  					
						b=extract(samples)$b
						if(!is.null(b)) {
							lastB=b[,lastDataIndex]
							lastBM=mean(lastB)
							paste("lastB",round(lastBM,2))
							
							bSm <- extract(samples)$bSm
							#summary(bSm)
							bSmM=mean(bSm)
							paste("bSm",round(bSmM,2))
							
							locTrendFract <- extract(samples)$locTrendFract
							summary(locTrendFract)
							locTrendFractM=mean(locTrendFract)
							paste("locTrendFract",round(locTrendFractM,2))
						}
  					
  					if (SEASONALITY>1) {
  					  sSm <- extract(samples)$sSm
						  #summary(sSm)
						  sSmM=mean(sSm)
						  paste("sSm",sSmM)
						
						  s=extract(samples)$s
						  sM=apply(s,2,mean)
						  #summary(sM)
						  cat("s:"); print(sM)
						}
						
  					sigma <- extract(samples)$sigma
  					#summary(sigma)
  					sigmaM=mean(sigma)
  					paste("sigma",sigmaM)
            
  					powx <- extract(samples)$powx
  					if (!is.null(powx)) {
  						powxM=mean(powx)
  						paste("powx",powxM)	
  					}
  					
  					offsetSigma=extract(samples)$offsetSigma
  					if (!is.null(offsetSigma)) {
  						summary(offsetSigma)
  						offsetSigmaM=mean(offsetSigma)
  						paste("offsetSigma",round(offsetSigmaM,2))	
  					}
  					
  					powTrend <- extract(samples)$powTrend
            if (!is.null(powTrend)) {
              #summary(powTrend)
              powTrendM=mean(powTrend)
              paste("powTrend",round(powTrendM,2))
              
              coefTrend <- extract(samples)$coefTrend
              #summary(coefTrend)
              coefTrendM=mean(coefTrend)
              paste("coefTrend",round(coefTrendM,2))
            }
  					
  					levSm <- extract(samples)$levSm
  					#summary(levSm)
  					levSmM=mean(levSm)
  					paste("levSm",levSmM)
  					
  					nu <- extract(samples)$nu
  					#summary(nu)
  					nuM=mean(nu)
  					paste("nu",nuM)
						
						innovSm<- extract(samples)$innovSm
						if (!is.null(innovSm)) {
							innovSmM=mean(innovSm)
							paste("innovSm",round(innovSmM,2))
							
							smoothedInnovSize<-extract(samples)$smoothedInnovSize
							lastSmoothedInnovSize=smoothedInnovSize[,lastDataIndex]
							lastSmoothedInnovSizeM=mean(lastSmoothedInnovSize)
							paste("lastSmoothedInnovSize",round(lastSmoothedInnovSizeM,2))
						}
            
            tskew<- extract(samples)$tskew
            if (!is.null(tskew)) {
              #summary(tskew)
              tskewM=mean(tskew)
              paste("tskew",tskewM)
            }

  					
  					Log.Info(paste(item,anchorDate))
  					lastDetails=paste0("nu=",round(nuM,2), ", lSm=",round(levSmM,2), " sigma=",round(sigmaM,2))
            if (!is.null(powTrend)) {
              lastDetails=paste0(lastDetails, ", coefT=",round(coefTrendM,2), ", powT=",round(powTrendM,2))
            }
  					if (!is.null(powx)) {
  						lastDetails=paste0(lastDetails, ", powx=",round(powxM,2))
  					}
						if (!is.null(offsetSigma)) {
							lastDetails=paste0(lastDetails, ", offsetS=",round(offsetSigmaM,2))
						}
						if(!is.null(b)) {
							lastDetails=paste0(lastDetails,
								", lastB=",round(lastBM,1), ", bSm=",round(bSmM,2), 
								", lTFract=", round(locTrendFractM,2))
						}
						if (SEASONALITY>1) {
							lastDetails=paste0(lastDetails,", sSm=",round(sSmM,2))
						}
						if (!is.null(innovSm)) {
							lastDetails=paste0(lastDetails,
								", innovSize=",round(lastSmoothedInnovSizeM,1), ", innovSm=",round(innovSmM,2))
						}
            if (!is.null(tskew)) {
              lastDetails=paste0(lastDetails, ", tskew=",round(tskewM,2))
            }
  					Log.Info(lastDetails)
  					
            
            #simulate->forecast
						bSmS=0; bS=0; locTrendFractS=0; tskewS=0; coefTrendS=0; powTrendS=0; #these initializations are important, do not remove. For manually skewed distribution tskewS=0 not SKEW (SKEW is only used during fitting)
  					t=1; irun=1
  					if (SEASONALITY>1) {
						  sS=rep(1,SEASONALITY+MAX_FORECAST_HORIZON)
						}
  					yf=matrix(lastLevelM,nrow=NUM_OF_TRIALS, ncol=MAX_FORECAST_HORIZON)
  					for (irun in 1:NUM_OF_TRIALS) {
  						indx=sample(mc_range,1)
  						prevLevel=lastLevel[indx]
  						levSmS=levSm[indx]
  						nuS=nu[indx]
  						sigmaS=sigma[indx]
              if (!is.null(powTrend)) {
                powTrendS=powTrend[indx]
                coefTrendS=coefTrend[indx]
              }
							if(!is.null(b)) {
								bSmS=bSm[indx]
								bS=lastB[indx]
								locTrendFractS=locTrendFract[indx]
							}
							if (!is.null(innovSm)) {
								innovSmS=innovSm[indx]
								innovSize=lastSmoothedInnovSize[indx]
								offsetsigmaS=offsetSigma[indx]
							} else if (!is.null(powx)) {
								powxS=powx[indx]
								offsetsigmaS=offsetSigma[indx]
							}
              if (!is.null(tskew)) {
                tskewS=tskew[indx]
              }
							
							#t=2
							if (SEASONALITY>1) {
  						  sS[1:SEASONALITY]=s[indx,(ncol(s)-SEASONALITY+1):ncol(s)]
                for (t in 1:MAX_FORECAST_HORIZON) {
                  seasonA=sS[t]
                  expVal=(prevLevel+ coefTrendS*abs(prevLevel)^powTrendS)*seasonA;
                  if (!is.null(powx)) {
                    omega=sigmaS*(abs(prevLevel))^powxS+offsetsigmaS
                  } else if (!is.null(innovSm)) {
                    omega=sigmaS*innovSize + offsetsigmaS
                  } else {
                    omega=sigmaS
                  }
                  error=rst(n=1, xi=0 ,	omega=omega, alpha=tskewS, nu=nuS)
                  yf[irun,t]=min(MAX_VAL,max(MIN_VAL,expVal+error))
                  currLevel=max(MIN_VAL,levSmS*yf[irun,t]/seasonA + (1-levSmS)*prevLevel) ;
                  
                  if (currLevel>MIN_VAL) {
                    #innovSize=innovSmS*abs(error)/seasonA+(1-innovSmS)*innovSize;
                    prevLevel=currLevel
                  } 
                  sS[t+SEASONALITY] <- sS[t];
                }	#through horizons
              # yf[irun,]
              } else { #nonseasonal
  							for (t in 1:MAX_FORECAST_HORIZON) {
  								expVal=prevLevel + coefTrendS*(abs(prevLevel))^powTrendS + locTrendFractS*bS
  								if (!is.null(powx)) {
                    omega=sigmaS*(abs(prevLevel))^powxS+offsetsigmaS
  								} else if (!is.null(innovSm)) {
  									omega=sigmaS*innovSize + offsetsigmaS
  								} else {
  									omega=sigmaS
  								}
  								error=rst(n=1, xi=0, omega=omega, alpha=tskewS, nu=nuS)
  								
  								yf[irun,t]=min(MAX_VAL,max(MIN_VAL, expVal+error))
  								currLevel=max(MIN_VAL,levSmS*yf[irun,t] + (1-levSmS)*prevLevel) ;
  								
  								if (currLevel>MIN_VAL) {
  									bS= bSmS*(currLevel-prevLevel)+(1-bSmS)*bS
  									#innovSize=innovSmS*abs(error)+(1-innovSmS)*innovSize;
  									prevLevel=currLevel
  								} 
  							} #through horizons
  						# yf[irun,]
  						}
  					} #through trials (simulations)
  					avgYfs=apply(yf,2,quantile,probs=QUANTS)	
  				} else {#not enough data, so use RandomWalk/ETS. This require relatively recent version of the forecast package. Here number of points >MIN_NUMBER_OF_DATA_POINTS_MAKE_FORECAST
  					mod=ets(y)
  					forec=forecast(mod, h=MAX_FORECAST_HORIZON, level=coverageLevels)
  					avgYfs=matrix(0,nrow=length(QUANTS), ncol=MAX_FORECAST_HORIZON)
  					for (ip in 1:length(positionsOfLowerTailsInForecastLower)) {
  						ipos=positionsOfLowerTailsInForecastLower[[ip]]
  						avgYfs[ip,]=pmax(MIN_VAL,pmin(MAX_VAL,forec$lower[,ipos]))	
  					}
  					avgYfs[ip+1,]=forec$mean
  					for (i in 1:length(highPercentiles)) {
  						avgYfs[ip+1+i,]=pmax(MIN_VAL,pmin(MAX_VAL,forec$upper[,i]))
  					}
            lastDetails=mod$method
            if (is.null(lastDetails)) {
              lastDetails=paste(mod$components, collapse='_')
            }
            yf=NULL
  				}
  				forecastVect=pmax(MIN_VAL,pmin(MAX_VAL,avgYfs[2,]))	
  				
  				versionNumber=paste(gsub('-','',as.character(anchorDate),fixed=T),Version_Number,sep='')
  				save_df=NULL; predictionHorizon=1; currMonth=-1
  				for (predictionHorizon in (1:MAX_FORECAST_HORIZON)) {
  					targetDate=anchorDate+predictionHorizon*FREQ
  					if (is.null(save_df)){
  						save_df=createForecastDataFrame(oneItem_df) #in reader
  					} else {
  						save_df=rbind(save_df, createForecastDataFrame(oneItem_df) )
  					}
  				}
  				sqlSave(fchannel, save_df, tablename=FORECAST_TABLE, append=T, rownames=F,verbose=F)  
  				
  				#plot
  				imageFileNam=paste(gsub('[*& ]','_',item, fixed=F),anchorDate, sep='_')
  				imageFileName=paste(imageFileNam,'.png',sep='')
  				filePath=paste(IMAGES_DIR, imageFileName,sep='/')
  				unlink(filePath)
  				png(file=filePath, bg="white", width = imageWidth, height = imageHeight, pointsize=16)
  				
          if (is.null(yf)) {
            pathsMatrix=NULL
          } else {
            pathsMatrix=yf[sample(nrow(yf),NUM_OF_TRIALS_TO_SHOW),]
          }
  				plotQ(Forecast_Name, item, pastDates, y, futureDates, forecastVect, 
  						knownFutureValues=yFuture,
  						quants=QUANTS[-indexOfMedian], quantileMatrix=avgYfs[-indexOfMedian,], 
  						pathsMatrix=pathsMatrix,
  						indxOfQuantToDisp=indexOfHigh, xlabel=lastDetails, ylabel=VAR, logAxis=LOG_FOR_DISPLAY)
  				
  				dev.off() ;
  				InsertGraph(imageFileName, htmlFilePath, htmlFilePath2)
        } #enough data
			} #through anchorDates of one item
			updateTrueValues(oneItem_df) #in reader
		} #!is.null(anchorDates)
	} #enough days
	gc()
} #through items
EndReport(htmlFilePath, htmlFilePath2)
odbcCloseAll()
Log.Info(predictionAlgorithm)
