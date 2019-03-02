# Test SGT or S2GT on M4 hourly data
# We will run in parallel, utilizinge all cores on the server - the more of them the faster the whole process will be.
# We will  go through all of the over 400 series, but in around 5% cases, where the seasonality does not appear to be (24,168)
# we switch to SGT. See at the bottom for a piece of code that helped to uncover these irregular cases.
#(Running SGT-only can be enforced by specifying USE_S2GT=FALSE)
#
# You can see the progress by opening (and refreshing) M4Hourly.html in S2GT_M4 subdirectory of your working directory.


USE_S2GT=FALSE

library(Rlgt)
# install.packages("devtools")
# devtools::install_github("carlanetto/M4comp2018")
library(M4comp2018)
library(doParallel)

USE_DB=TRUE  #if using the database. 
#One of the main benefits is that you can restart calculations (just run with the seame runName)
# 1) create the destination tables first (look at the end of the file for an example script)
# 2) create ODBC 64-bit system name "slawek" pointing to the destination database. 
if (USE_DB) {
	library(RODBC)
}
#set.seed(12)
ODBC_SOURCE_NAME="slawek"  
options(width=180)
if (.Platform$OS.type=="windows")  memory.limit(10000)
imageWidth=1000; imageHeight=400

variable='Hourly' 
LBack=0 #not used currently, but needed for db compatibity with nn-based systems

runName='mult_Avg_F'
#if running with the same runName, the system will try to finish, presumably unfinished, run.
level_method="seasAvg"  #c("HW", "seasAvg","HW_sAvg")


if (USE_S2GT) {
	SEASONALITY=24
	SEASONALITY2=168
} else {
	SEASONALITY=168
	SEASONALITY2=1
}

#names(M4[[1]])
hourly=Filter(function(l) l$period == variable, M4)
NUM_OF_CASES=length(hourly) #running over 400 cases would take quite a few days :-)
#NUM_OF_CASES=10
str(hourly[[1]])
H=hourly[[1]]$h

unexpectedSeasonalityList=list()
unexpectedSeasonalityList[["174"]]=23
unexpectedSeasonalityList[["202"]]=23
unexpectedSeasonalityList[["222"]]=21
unexpectedSeasonalityList[["223"]]=16
unexpectedSeasonalityList[["224"]]=23
unexpectedSeasonalityList[["230"]]=23
unexpectedSeasonalityList[["240"]]=10
unexpectedSeasonalityList[["244"]]=17
unexpectedSeasonalityList[["247"]]=23
unexpectedSeasonalityList[["248"]]=23
unexpectedSeasonalityList[["249"]]=23
unexpectedSeasonalityList[["254"]]=23
unexpectedSeasonalityList[["255"]]=23
unexpectedSeasonalityList[["267"]]=23
unexpectedSeasonalityList[["268"]]=23
unexpectedSeasonalityList[["271"]]=23
unexpectedSeasonalityList[["272"]]=23
unexpectedSeasonalityList[["289"]]=23
unexpectedSeasonalityList[["325"]]=23



OUTPUT_DIR="S2GT_M4"
fullOutputDir=file.path(getwd(),  paste0(OUTPUT_DIR,variable,'_',runName))
print(paste("The output will go to",fullOutputDir))
if (!file.exists(fullOutputDir)){
	dir.create(fullOutputDir)
}
htmlFilePath=file.path(fullOutputDir,"M4Hourly.html")
unlink(htmlFilePath)
imagesDir=file.path(fullOutputDir,'images')
if (!file.exists(imagesDir)){
	dir.create(imagesDir)
}
cat(paste0("<html><head><title>",OUTPUT_DIR, "</title></head> <body bgcolor=#D0D0D0>"), file = htmlFilePath, append = FALSE)


numOfCores=parallel:::detectCores()
CLUSTER_SIZE=as.integer(numOfCores/4)  #Each cluster will use 4 cores (as per default in rlgt.control())

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

#stopCluster(cl) # this is useful if rerunning the script in the same session. First time it will report error - it is OK
#Sys.sleep(2)
unlink("rclust.out")
cl = makeCluster(CLUSTER_SIZE, outfile="rclust.out")
registerDoParallel(cl)
sink(); graphics.off();
stanOutputPath='stanOutput.txt' #in working directory
unlink(stanOutputPath)


now=Sys.time()
if (USE_DB) { #check in meta table if this is the continuation run
	fchannel<- odbcConnect(ODBC_SOURCE_NAME)
	odbcGetInfo(fchannel)
	doneAlready_query=paste0("select dateTimeOfPrediction from M4StanModels",
 			" where run='",runName,  
			"' and LBack=",LBack,
			" and variable='",variable,
			"' ")
	doneAlready_df=sqlQuery(fchannel,doneAlready_query, stringsAsFactors =F)
	if (nrow(doneAlready_df)>0) {
		now=doneAlready_df$dateTimeOfPrediction
	} else {
		save_df=data.frame(run=runName, LBack=LBack, variable=variable, dateTimeOfPrediction=now, comments='')
		sqlSave(fchannel, save_df, tablename='M4StanModels', append=T, rownames=F,verbose=F)
	}
}

firstTime=TRUE; i=1
ret_df=foreach(i=1:NUM_OF_CASES, .combine=rbind, .inorder=FALSE, .packages=c("Rlgt","RODBC")) %dopar% { # the long loop :-)	
	if (firstTime==TRUE) { #needs to be done inside every new slave process
		fchannel<- odbcConnect(ODBC_SOURCE_NAME)
		odbcGetInfo(fchannel)
		sink(file=stanOutputPath, append=TRUE, split=TRUE)
		firstTime=FALSE
	}
	series=hourly[[i]]$st
	actuals = hourly[[i]]$xx   #need to be here, whethere we calculate or just read from db
	if (!USE_S2GT) {#convert freq from 24 to 168
		tsp_info=tsp(actuals)
		newStart=(tsp_info[1]*tsp_info[3])/SEASONALITY
		actuals=ts(actuals, start=newStart, freq=SEASONALITY)
		#plot(actuals)
	}

	proceed=TRUE
	if (USE_DB) {
		doneAlready_query=paste0("select * from M4StanModels m, M4Stan d",
				" where run='",runName,  
				"' and LBack=",LBack,
				" and variable='",variable,
				"' and series='",series,
				"' and m.dateTimeOfPrediction=d.dateTimeOfPrediction
        order by horizon")
		doneAlready_df=sqlQuery(fchannel,doneAlready_query, stringsAsFactors =F)
		if (nrow(doneAlready_df)>0) {
			proceed=FALSE
		}
	}
	
	if (proceed) {
		print(paste("starting",series))
		trainData = hourly[[i]]$x 
		if (!USE_S2GT) {#convert freq from 24 to 168
			tsp_info=tsp(trainData)
			newStart=(tsp_info[1]*tsp_info[3])/SEASONALITY
			trainData=ts(trainData, start=newStart, freq=SEASONALITY)
			#plot(trainData)
		}
	
		if (is.null(unexpectedSeasonalityList[[as.character(i)]])) {
			rstanmodel <- rlgt(trainData,seasonality2=SEASONALITY2, #but if SEASONALITY2==1, then we are using SGT
				control=rlgt.control(NUM_OF_ITER=10000),   
				#seasonality.type="generalized",
				level.method=level_method,
				verbose=TRUE)
		} else {#use SGT
			trainData = as.numeric(trainData) #to remove incorrect frequency stamp
			actuals = as.numeric(actuals)
			rstanmodel <- rlgt(trainData, seasonality=unexpectedSeasonalityList[[as.character(i)]], 
				control=rlgt.control(NUM_OF_ITER=10000),
				verbose=TRUE) 
		}
		
		if (SEASONALITY2==1) {
			startParToDisplay=3
		} else {
			startParToDisplay=4
		}			
						
		forec= forecast(rstanmodel, h = H, level=c(95,98))
		#str(forec$model$params, max.level=1)
		
		imageFileName=paste(series,'.png',sep='')
		relPath=file.path('images',imageFileName)
		filePath=file.path(imagesDir, imageFileName)
		unlink(filePath)
		png(file=filePath, bg="white", width = imageWidth, height = imageHeight, pointsize=16)
		
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
		
		dev.off() ;
		cat(paste0('<img src="',relPath,'" alt="',series,'" height="',imageHeight,'" width="',imageWidth,'">'), file = htmlFilePath, append = TRUE)
	
		sMAPE=mean(abs(forec$mean-actuals)/(forec$mean+actuals))*200
		q975Loss=quantileLoss(forec$upper[,1], actuals, 0.975)	
		q99Loss=quantileLoss(forec$upper[,2], actuals, 0.99)
		q025Loss=quantileLoss(forec$lower[,1], actuals, 0.025)
		cat(paste0("<p> ",Sys.time(), " ",series,  " sMAPE=",signif(sMAPE,3),"% </p>"), file = htmlFilePath, append = TRUE)
		params_txt=NULL;
		for (ipar in startParToDisplay:length(forec$model$params)) {
			paramName=names(forec$model$params)[ipar]
			if (!is.null(params_txt)) params_txt=paste0(params_txt,", ")
			params_txt=paste0(params_txt, paramName, "=",signif(median(forec$model$params[[paramName]]),2))
		}
		cat(paste0("<p> Medians of params: ",params_txt,"</p>"), file = htmlFilePath, append = TRUE)
		
		if (USE_DB) {
			save_df=data.frame(dateTimeOfPrediction=now, series=series, horizon=1:H,
					actual=actuals, predQ50=forec$mean, predQ2_5=forec$lower[,1], predQ97_5=forec$upper[,1], predQ99=forec$upper[,2])
			sqlSave(fchannel, save_df, tablename='M4Stan', append=T, rownames=F,verbose=F)
		}
		
		data.frame(series=series, sMAPE=sMAPE, 
			q025Loss=q025Loss, q975Loss=q975Loss, q99Loss=q99Loss,
			numOfCases99pExceeded=sum(actuals>forec$upper[,2]),
			numOfCases975pExceeded=sum(actuals>forec$upper[,1]),
			numOfCases025pExceeded=sum(actuals<forec$lower[,1]))
	} else {
		#doneAlready_df
	  if (sum(abs(actuals-doneAlready_df$actual))>1e-4*mean(actuals))  {
			stop(paste0("error, diffs between db and new actuals for ",series))
		}
		sMAPE=mean(abs(doneAlready_df$predQ50-actuals)/(doneAlready_df$predQ50+actuals))*200
		q975Loss=quantileLoss(doneAlready_df$predQ97_5, actuals, 0.975)	
		q99Loss=quantileLoss(doneAlready_df$predQ99, actuals, 0.99)
		q025Loss=quantileLoss(doneAlready_df$predQ2_5, actuals, 0.025)
		
		data.frame(series=series, sMAPE=sMAPE, 
				q025Loss=q025Loss, q975Loss=q975Loss, q99Loss=q99Loss,
				numOfCases99pExceeded=sum(actuals>doneAlready_df$predQ99),
				numOfCases975pExceeded=sum(actuals>doneAlready_df$predQ97_5),
				numOfCases025pExceeded=sum(actuals<doneAlready_df$predQ2_5))
	}
}	
sMAPE=mean(ret_df$sMAPE)
q975Loss=mean(ret_df$q975Loss)
q99Loss=mean(ret_df$q99Loss)
q025Loss=mean(ret_df$q025Loss)
exceed99=mean(ret_df$numOfCases99pExceeded)/H*100
exceed975=mean(ret_df$numOfCases975pExceeded)/H*100
exceed025=mean(ret_df$numOfCases025pExceeded)/H*100
print(paste0("SUMMARY: Num of cases:", nrow(ret_df), ", sMAPE:",signif(sMAPE,4),
	', % of time 99p exceeded:',signif(exceed99,4),
  ', % of time 97.5p exceeded:',signif(exceed975,4), ', % of time 2.5p exceeded:',signif(exceed025,4)
	,', q99Loss:',signif(q99Loss,4),
	', q2.5Loss:',signif(q025Loss,4),', q95Loss:',signif(q975Loss,4) ))

	
	
	
#The unexpectedSeasonalityIndices was created by running following code and visually inspecting, accepting/rejecting proposals made by findfrequency() function
#library(forecast)
#par(mfrow=c(2,1))
#for (indx in 1:length(hourly)) {
#	y=hourly[[indx]]$x
#	freq=findfrequency(y)
#	if (freq!=24) {
#		print(paste(indx,frequency(y),freq))
#		plot(hourly[[indx]]$x, main=indx)
#		for (ii in 1:length(y)) {
#			abline(v=start(hourly[[indx]]$x)+ii)
#    }
#		acf(hourly[[indx]]$x)
#		Sys.sleep(10)
#	}
#}

#CREATE TABLE [dbo].[M4StanModels](
#		[run] [varchar](164) NOT NULL,
#		[LBack] [tinyint] NOT NULL,
#		[variable] [varchar](20) NOT NULL,
#		[dateTimeOfPrediction] [datetime] NOT NULL,
#		[comments] [varchar](300) NULL,
#		CONSTRAINT [M4StanModels_pk] PRIMARY KEY CLUSTERED 
#				(
#				[run] ASC,
#				[LBack] ASC,
#				[variable] ASC,
#       [dateTimeOfPrediction] asc
#		)
#) 
#CREATE TABLE [dbo].[M4Stan](
#		[dateTimeOfPrediction] [datetime] NOT NULL,
#		[series] [varchar](50) NOT NULL,
#		[horizon] [tinyint] NOT NULL,
#		[actual] [real] NOT NULL,
#		[predQ50] [real] NOT NULL,
#		[predQ2_5] [real] NOT NULL,
#		[predQ97_5] [real] NOT NULL,
#		[predQ99] [real] NOT NULL,
#		CONSTRAINT [M4py_PK] PRIMARY KEY CLUSTERED 
#				(
#				[dateTimeOfPrediction] ASC,
#				[series] ASC,
#				[horizon] ASC
#		)
#)

