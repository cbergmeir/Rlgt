# constants to accompany prodDriver.R
# They are defaults that can be overwritten in the main program, e.g. Azure_StorageXIO_ShortTerm2.R
# Author: slsmyl
# Dec 2016
###############################################################################

library(rstan)

PRODUCTION_MODE=F #if true, no backtesting is done

#high-level prarllelism
numOfCores=parallel:::detectCores()
if (numOfCores==8) {
	CLUSTER_SIZE=3
} else if (numOfCores==16) {
	CLUSTER_SIZE=6
} else {
	CLUSTER_SIZE=as.integer(numOfCores/3)
} 
STAN_CLUSTER_SIZE=4 #so at peak there will be CLUSTER_SIZE*STAN_CLUSTER_SIZE rscript processes running, but not all the time (some chains stop earlier), so some oversubscription is OK
rstan_options(auto_write = TRUE)
options(mc.cores = STAN_CLUSTER_SIZE, width=180)

MIN_MEMORY_MB=3999
if (memory.limit()<MIN_MEMORY_MB) {
	memory.limit(MIN_MEMORY_MB)
}	
memory.limit()

CLUST_OUT_DIR=OUTPUT_DIR
CLUST_OUT_PATH=paste(CLUST_OUT_DIR, "rclust.txt", sep='/')
if(!dir.exists(CLUST_OUT_DIR)){
	dir.create(CLUST_OUT_DIR, recursive = TRUE)
}

Sys.setenv(TZ = "GMT");
Sys.timezone()
today=Sys.Date();
START_DATE=today-5*365; #default
now=Sys.time()

imageHeight=600
imageWidth=1300		

LOG_FOR_DISPLAY='' #e.g. put 'y' if you want logarithm on y on the image 
MIN_DAYS_TO_MAKE_FORECAST=60
MIN_NUMBER_OF_DATA_POINTS_MAKE_FORECAST=9 #really minimum, Random Walk or some basic ets() will be used 
SHOW_ALG_DETAILS=T
Version_Number="00"

###
NUM_OF_TRIALS=2000
NUM_OF_TRIALS_TO_SHOW=150
if (NUM_OF_TRIALS<NUM_OF_TRIALS_TO_SHOW)
  NUM_OF_TRIALS_TO_SHOW=NUM_OF_TRIALS

NUM_OF_ITER=5000 #sampling, you can try to reduce to say 2.5k to reduce calc time, but I am not usre if you get much
MAX_NUM_OF_REPEATS=3
MAX_RHAT_ALLOWED=1.005

MAX_NU=15 
MIN_NU=2
MIN_SIGMA=0.01
MIN_VAL=0.01
MAX_VAL=1e38 # roughly what SQL Server real can take
MIN_POW=-0.5; #you could change it to 0 if you do not like the negative trend pow
MAX_POW=1
CAUCHY_SD_DIV=200 # could be 50-300, does not matter too much

POW_TREND_ALPHA=1 #to make the forecast more curved, make it larger like 3 to 6
POW_TREND_BETA=1

POWX_ALPHA=1  
POWX_BETA=1  #if the powx fitted is too often too high (i.e.> 0.6) you can try to tame it down by increasing POWX_BETA to say 3 to 5, or more (but if say 10 or 20 it becomes too strong-armed) 

SKEW=0 #would be overwritten in the main program if a skewed model is used
FREQ=1 #for weekly data, set to 7
SEASONALITY=1 

STEP_BACK_WINDOW=0
