# Test SGT - M3 data
# We will use all cores on the server - the more of them the faster the whole process will be.
# If US_DB is true, the forecasts and actuals are saved in an database (SQL Server in this case, but with small changes any ODBC- compatible database, e.g. MySql would do).


library(Mcomp)
library(Rlgt)
library(doParallel)
library(ggplot2)

USE_DB=TRUE  #if using the database, 
# 1) create the destination tables first (look at the end of the file for an example script)
# 2) create ODBC 64-bit system name "slawek" pointing to the destination database. 
if (USE_DB) {
	library(RODBC)
}

#set.seed(12)
ODBC_SOURCE_NAME="slawek"  
options(width=180)
imageWidth=1000; imageHeight=400

runName='mult_HW_s' #if using a database, make it different for every run
variable='MONTHLY' #or "YEARLY", "QUARTERLY", "MONTHLY" or "OTHER". 
LBack=0

OUTPUT_DIR="GT_M3_"
fullOutputDir=file.path(getwd(),  paste0(OUTPUT_DIR,variable,'_',runName))
print(paste("The output will go to",fullOutputDir))
if (!file.exists(fullOutputDir)){
	dir.create(fullOutputDir)
}
htmlFilePath=file.path(fullOutputDir,paste0("M3_",variable,".html"))
unlink(htmlFilePath)
imagesDir=file.path(fullOutputDir,'images')
if (!file.exists(imagesDir)){
	dir.create(imagesDir)
}
cat(paste0("<html><head><title>",OUTPUT_DIR, "</title></head> <body bgcolor=#D0D0D0>"), file = htmlFilePath, append = FALSE)


numOfCores=parallel:::detectCores()
CLUSTER_SIZE=as.integer(numOfCores/4)  #Each cluster will use 4 cores (as per default in rlgt.control())

M3.data <- subset(M3,variable)
NUM_OF_CASES=length(M3.data)
#NUM_OF_CASES=20
H=length(M3.data[[1]]$xx)
seasonality=frequency(M3.data[[1]]$x)
startParamToDisplay=2
if (seasonality>1) {
	startParamToDisplay=3	
}
	

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
	#fchannel<- odbcConnect(ODBC_SOURCE_NAME)
  fchannel= odbcConnect(ODBC_SOURCE_NAME, uid='SA',pwd='your_pwd')
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
		sqlSave(fchannel, save_df, tablename='MStanModels', append=T, rownames=F,verbose=F)
	}
}

firstTime=TRUE; i=1
ret_df=foreach(i=1:NUM_OF_CASES, .combine=rbind, .inorder=FALSE, .packages=c("Rlgt","RODBC")) %dopar% { # the long loop :-)	
	if (firstTime==TRUE) { #needs to be done inside every new slave process
		#fchannel<- odbcConnect(ODBC_SOURCE_NAME)
	  fchannel= odbcConnect(ODBC_SOURCE_NAME, uid='SA',pwd='your_pwd')
		odbcGetInfo(fchannel)
		sink(file=stanOutputPath, append=TRUE, split=TRUE)
		firstTime=FALSE
	}
	series=M3.data[[i]]$sn
	actuals <- M3.data[[i]]$xx   #need to be here, whethere we calculate or just read from db
	
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
		trainData <- M3.data[[i]]$x
		rstanmodel <- rlgt(trainData, 
				#seasonality.type="generalized",  #<--------------------------------------------------------------------------------------------
				#level.method="seasAvg",  #c("HW", "seasAvg","HW_sAvg"),
			control=rlgt.control(NUM_OF_ITER=5000), #we do not need to specify seasonality, as it is extracted from M3.data[[i]]$x
			verbose=TRUE)
		forec= forecast(rstanmodel, h = H, level=c(90,95,98))
		#str(forec, max.level=1)
		#forec$lower  #5%  2.5%   1%
		#forec$upper #95% 97.5%   99%
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
		q05Loss=quantileLoss(forec$lower[,1], actuals, 0.05)
		q95Loss=quantileLoss(forec$upper[,1], actuals, 0.95)	
		q025Loss=quantileLoss(forec$lower[,2], actuals, 0.025)
		q975Loss=quantileLoss(forec$upper[,2], actuals, 0.975)
		q01Loss=quantileLoss(forec$lower[,3], actuals, 0.01)
		q99Loss=quantileLoss(forec$upper[,3], actuals, 0.99)
		cat(paste0("<p> ",Sys.time(), " ",series,  " sMAPE=",signif(sMAPE,3),"% </p>"), file = htmlFilePath, append = TRUE)
		params_txt=NULL;ipar=3
		for (ipar in startParamToDisplay:length(forec$model$params)) {
			paramName=names(forec$model$params)[ipar]
			if (!is.null(params_txt)) params_txt=paste0(params_txt,", ")
			params_txt=paste0(params_txt, paramName, ":",signif(median(forec$model$params[[paramName]]),2))
		}
		cat(paste0("<p> ",params_txt,"</p>"), file = htmlFilePath, append = TRUE)
		
		if (USE_DB) {
			save_df=data.frame(dateTimeOfPrediction=now, series=series, horizon=1:H,
					actual=actuals, predQ50=forec$mean, 
					predQ5=forec$lower[,1], predQ95=forec$upper[,1],
					predQ2_5=forec$lower[,2], predQ97_5=forec$upper[,2], 
					predQ1=forec$lower[,3], predQ99=forec$upper[,3])
			sqlSave(fchannel, save_df, tablename='M4Stan', append=T, rownames=F,verbose=F)
		}
		
		data.frame(series=series, sMAPE=sMAPE, 
		  q05Loss=q05Loss, q95Loss=q95Loss, 
			q025Loss=q025Loss, q975Loss=q975Loss, 
			q01Loss=q01Loss, q99Loss=q99Loss,
			numOfCases05pExceeded=sum(actuals<forec$lower[,1]),
			numOfCases95pExceeded=sum(actuals>forec$upper[,1]),
			numOfCases025pExceeded=sum(actuals<forec$lower[,2]),
			numOfCases975pExceeded=sum(actuals>forec$upper[,2]),
			numOfCases01pExceeded=sum(actuals<forec$lower[,3]),
			numOfCases99pExceeded=sum(actuals>forec$upper[,3]))
	} else {
		#doneAlready_df
	  if (sum(abs(actuals-doneAlready_df$actual))>1e-4*mean(actuals))  {
			stop(paste0("error, diffs between db and new actuals for ",series))
		}
		sMAPE=mean(abs(doneAlready_df$predQ50-actuals)/(doneAlready_df$predQ50+actuals))*200
		q05Loss=quantileLoss(doneAlready_df$predQ5, actuals, 0.05)
		q95Loss=quantileLoss(doneAlready_df$predQ95, actuals, 0.95)	
		q025Loss=quantileLoss(doneAlready_df$predQ2_5, actuals, 0.025)
		q975Loss=quantileLoss(doneAlready_df$predQ97_5, actuals, 0.975)	
		q01Loss=quantileLoss(doneAlready_df$predQ1, actuals, 0.01)
		q99Loss=quantileLoss(doneAlready_df$predQ99, actuals, 0.99)
		
		data.frame(series=series, sMAPE=sMAPE, 
		           q05Loss=q05Loss, q95Loss=q95Loss, 
		           q025Loss=q025Loss, q975Loss=q975Loss, 
		           q01Loss=q01Loss, q99Loss=q99Loss,
		           numOfCases05pExceeded=sum(actuals<doneAlready_df$predQ5),
		           numOfCases95pExceeded=sum(actuals>doneAlready_df$predQ95),
		           numOfCases025pExceeded=sum(actuals<doneAlready_df$predQ2_5),
		           numOfCases975pExceeded=sum(actuals>doneAlready_df$predQ97_5),
		           numOfCases01pExceeded=sum(actuals<doneAlready_df$predQ1),
		           numOfCases99pExceeded=sum(actuals>doneAlready_df$predQ99)
		          )
	} #done already
}	#through cases

sMAPE=mean(ret_df$sMAPE)
q05Loss=mean(ret_df$q05Loss)
q95Loss=mean(ret_df$q95Loss)
q025Loss=mean(ret_df$q025Loss)
q975Loss=mean(ret_df$q975Loss)
q01Loss=mean(ret_df$q01Loss)
q99Loss=mean(ret_df$q99Loss)

exceed05=mean(ret_df$numOfCases05pExceeded)/H*100
exceed95=mean(ret_df$numOfCases95pExceeded)/H*100
exceed025=mean(ret_df$numOfCases025pExceeded)/H*100
exceed975=mean(ret_df$numOfCases975pExceeded)/H*100
exceed01=mean(ret_df$numOfCases01pExceeded)/H*100
exceed99=mean(ret_df$numOfCases99pExceeded)/H*100

print(paste0("SUMMARY: Num of cases:", nrow(ret_df), ", sMAPE:",signif(sMAPE,4),
  ', % of time exceeded 1p:',signif(exceed01,4),
  ', 2.5p:',signif(exceed025,4),
  ', 5p:',signif(exceed05,4),
  ', 95p:',signif(exceed95,4), 
  ', 97.5p:',signif(exceed975,4), 
	', 99p:',signif(exceed99,4),
  ', qLoss 1p:',signif(q01Loss,4),
  ', 2.5p:',signif(q025Loss,4),
  ', 5p:',signif(q05Loss,4),
  ', 95p:',signif(q95Loss,4),
  ', 97.5p:',signif(q975Loss,4),
	', 99p:',signif(q99Loss,4)
   ))


#CREATE TABLE [dbo].[MStanModels](
#		[run] [varchar](164) NOT NULL,
#		[LBack] [tinyint] NOT NULL,
#		[variable] [varchar](20) NOT NULL,
#		[dateTimeOfPrediction] [datetime] NOT NULL,
#		[comments] [varchar](300) NULL,
#		CONSTRAINT [MStanModels_pk] PRIMARY KEY CLUSTERED 
#				(
#				[run] ASC,
#				[LBack] ASC,
#				[variable] ASC,
#       [dateTimeOfPrediction] asc
#		)
#) 
#CREATE TABLE [dbo].[MStan](
#[dateTimeOfPrediction] [datetime] NOT NULL,
#[series] [varchar](50) NOT NULL,
#[horizon] [tinyint] NOT NULL,
#[actual] [real] NOT NULL,
#[predQ50] [real] NOT NULL,
#[predQ5] [real] NOT NULL,
#[predQ95] [real] NOT NULL,
#[predQ2_5] [real] NOT NULL,
#[predQ97_5] [real] NOT NULL,
#[predQ1] [real] NOT NULL,
#[predQ99] [real] NOT NULL,
#CONSTRAINT [Mpy_PK] PRIMARY KEY CLUSTERED 
#(
#  [dateTimeOfPrediction] ASC,
#  [series] ASC,
#  [horizon] ASC
#))

