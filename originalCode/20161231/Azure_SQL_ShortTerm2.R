# regions, short-term, X structure, using SGT
# Nov 2016


SQL.SERVER.NAME='.'
DATABASE.NAME='slawek'
##################

R_DIR="C:/progs/r/short_term/Azure_SQL_ShortTerm_Forecast"
LIB_DIR="C:/progs/r/short_term/libs"

source(paste(R_DIR,'Azure_SQL_ShortTerm2_Reader.R',sep='/'))
source(paste(LIB_DIR,'prodDriverConsts.R',sep='/'))
algorithm='LGT'
source(paste(LIB_DIR,paste0(algorithm,'.R'),sep='/'))#this will take a few minutes perhaps
OUTPUT_DIR="c:/reports/SQLS_ShortTerm"
if(!dir.exists(OUTPUT_DIR)){
	dir.create(OUTPUT_DIR, recursive = TRUE)
}


FORECAST_TABLE='Azure_SQL_ShortTerm_Forecast'
reportTitle='SQL Azure in Regions per Type'
Forecast_Name = 'AzureShortTerm'
Service="Azure"
Location_Level=3
Location_LevelName="region"
Item_Level=2
Item_LevelName="Platform"
Item_Parent_Name=""
Version_Number="00"
Customer_Segment = ""


START_DATE=as.Date('2014-10-01')
VAR='DTU'
Forecast_Unit = VAR
MAX_SIZE_OF_LOOK_BACK_WINDOW=365*6/12;
ANCHOR_STEP=30 #daily data
FREQ=1
FORECAST_HORIZONS=(1:133)
MAX_FORECAST_HORIZON=max(FORECAST_HORIZONS)
PERCENTILES=c(5,50,55,60,65,70,75,80,85,90,95,99) #so even you can chage here, 
# and even if I make the reader dynamically writing them (not done yet),
# the forecast table needs to be in perfect sync with these. So do not touch lightly :-)
QUANTS=PERCENTILES/100 
#LOG_FOR_DISPLAY='y' #if you want logarithmic y axis 

POW_TREND_ALPHA=4
SKEW=0
MIN_NU=5
if (SKEW!=0) {
	predictionAlgorithm=paste(paste0(algorithm,SKEW),paste0('TA',POW_TREND_ALPHA),paste0('MinNu',MIN_NU),sep='_');	
} else {
	predictionAlgorithm=paste(algorithm,paste0('TA',POW_TREND_ALPHA),paste0('MinNu',MIN_NU),sep='_');
}
preProcessingAlgorithm="6m"

#################################
readerParams=list()
readerParams[['VAR']]=VAR
bigSum_df=reader(readerParams)
#doForecast(bigSum_df)
source(paste(LIB_DIR,'prodDriver.R',sep='/'))

