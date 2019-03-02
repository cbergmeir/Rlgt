# Test SGT or S2GT on M4 hourly data
# We will run in par alle, utilizinge all cores on the server - the more of them the faster the whole process will be.
# We will  go through all of the over 400 series, but in around 5% cases, where the seasonality does not appear to be (24,168)
# we switch to SGT. See at the bottom for a piece of code that helped to uncover these irregular cases.
# You can see the progress by opening (and refreshing) M4Hourly.html in S2GT_M4 subdirectory of your working directory.
# New: The code may run in mostly S2GT mode (as explained above)
# or purely SGT, using as the seasonality 168, just change the logical variable below

USE_S2GT=FALSE

library(Rlgt)
# install.packages("devtools")
# devtools::install_github("carlanetto/M4comp2018")
library(M4comp2018)
library(doParallel)

#set.seed(12)
options(width=180)
if (.Platform$OS.type=="windows")  memory.limit(10000)
imageWidth=1000; imageHeight=400

H=48
if (USE_S2GT) {
	SEASONALITY=24
	SEASONALITY2=168
} else {
	SEASONALITY=168
	SEASONALITY2=1
}

#names(M4[[1]])
hourly=Filter(function(l) l$period == "Hourly", M4)
NUM_OF_CASES=length(hourly) #running over 400 cases would take quite a few days :-)
#NUM_OF_CASES=10

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
fullOutputDir=file.path(tempdir(),OUTPUT_DIR)
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

i=1
ret_df=foreach(i=1:NUM_OF_CASES, .combine=rbind, .inorder=FALSE, .packages=c("Rlgt")) %dopar% { # the long loop :-)	
	if (i==1) { #needs to be done inside every new slave process
		sink(file=stanOutputPath, append=TRUE, split=TRUE)
	}
	series=hourly[[i]]$st
	print(paste("starting",series))
	
	trainData = hourly[[i]]$x 
	actuals = hourly[[i]]$xx  
	
	if (is.null(unexpectedSeasonalityList[[as.character(i)]])) {
		rstanmodel <- rlgt(trainData,seasonality2=SEASONALITY2, #but if SEASONALITY2==1, then we are using SGT
				control=rlgt.control(NUM_OF_ITER=10000),   
				#seasonality.type="generalized",
				level.method="HW_sAvg",  #c("HW", "seasAvg","HW_sAvg"),
				verbose=TRUE)
		startParToDisplay=4
	} else {
		trainData = as.numeric(trainData) #to remove incorrect frequency stamp
		actuals = as.numeric(actuals)
		rstanmodel <- rlgt(trainData, seasonality=unexpectedSeasonalityList[[as.character(i)]], 
				control=rlgt.control(NUM_OF_ITER=10000),
				verbose=TRUE) #use SGT
		startParToDisplay=3
	}
					
	forec= forecast(rstanmodel, h = H, level=c(90,98))
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
	q95Loss=quantileLoss(forec$upper[,1], actuals, 0.95)	
	q99Loss=quantileLoss(forec$upper[,2], actuals, 0.99)
	q5Loss=quantileLoss(forec$lower[,1], actuals, 0.05)
	cat(paste0("<p> ",Sys.time(), " ",series,  " sMAPE=",signif(sMAPE,3),"% </p>"), file = htmlFilePath, append = TRUE)
	params_txt=NULL;
	for (ipar in startParToDisplay:length(forec$model$params)) {
		paramName=names(forec$model$params)[ipar]
		if (!is.null(params_txt)) params_txt=paste0(params_txt,", ")
		params_txt=paste0(params_txt, paramName, "=",signif(median(forec$model$params[[paramName]]),2))
	}
	cat(paste0("<p> Medians of params: ",params_txt,"</p>"), file = htmlFilePath, append = TRUE)
	
	data.frame(series=series, sMAPE=sMAPE, 
		q5Loss=q5Loss, q95Loss=q95Loss, q99Loss=q99Loss,
		numOfCases95pExceeded=sum(actuals>forec$upper[,1]),
		numOfCases5pExceeded=sum(actuals<forec$lower[,1]))
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

