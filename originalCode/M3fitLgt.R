# Author: Slawek Smyl 
# Oct 2016
# An example fitting the (non-seasonal) LGT to one of the M3 yearly curves
###############################################################################
# install.packages("Mcomp")
# install.packages("rstan") #but first install RTools

require(Mcomp)
library(rstan)
library(sn)

CAUCHY_SD_DIV=200
STAN_CLUSTER_SIZE=4
rstan_options(auto_write = TRUE)
options(mc.cores = STAN_CLUSTER_SIZE, width=180)

QUANTS=c(0.05,0.5,0.95)
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
imageWidth=1200; imageHeight=600
NUM_OF_TRIALS=5000  
NUM_OF_TRIALS_TO_SHOW=100


#non-seasonal model
parameters = c("l", "b", "nu", "sigma", "levSm",  "bSm", 
		"powx", "coefTrend",  "powTrend", "offsetSigma", "locTrendFract")
modelX = '
	data {  
		real<lower=0> CAUCHY_SD;
		real MIN_POW;  real MAX_POW;
		real<lower=0> MIN_SIGMA;
		real<lower=1> MIN_NU; real<lower=1> MAX_NU;
		int<lower=1> N;
		vector<lower=0>[N] y;
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
		
		l[1] = y[1]; b[1] = 0;
		powTrend= (MAX_POW-MIN_POW)*powTrendBeta+MIN_POW;
		
		for (t in 2:N) {
			l[t]  = levSm*y[t] + (1-levSm)*l[t-1] ;
			b[t]  = bSm*(l[t]-l[t-1]) + (1-bSm)*b[t-1] ;
		}
	}
	model {
		sigma ~ cauchy(0,CAUCHY_SD) T[0,];
		offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
		coefTrend ~ cauchy(0,CAUCHY_SD);
		
		for (t in 2:N) {
		y[t] ~ student_t(nu, l[t-1]+coefTrend*fabs(l[t-1])^powTrend+locTrendFract*b[t-1], 
		  sigma*fabs(l[t-1])^powx+ offsetSigma);
	}
}
'  
stanModel = stan_model(model_code=modelX)
#str(stanModel)

fitLgt <- function(y,maxPredictionHorizon=6) {
  n=length(y)
  CauchySd=max(y)/CAUCHY_SD_DIV
  data = list(CAUCHY_SD=CauchySd, 
      y=y, N=n) # to be passed on to Stan  
  
  avgRHat=1e20; irep=1
  for (irep in 1:MAX_NUM_OF_REPEATS) {
    #initializations = list();
    samples1=
        sampling(#control=list(adapt_delta = 0.9, max_treedepth=11),
            stanModel,   
            data=data, 
            #init=initializations,
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
  
  levSm = extract(samples)$levSm
  #summary(levSm)
  levSmM=mean(levSm)
  paste("levSm",levSmM)
  
  nu = extract(samples)$nu
  #summary(nu)
  nuM=mean(nu)
  paste("nu",nuM)
  
  coefTrend = extract(samples)$coefTrend
  #summary(coefTrend)
  coefTrendM=mean(coefTrend)
  paste("coefTrend",coefTrendM)
  
  powTrend = extract(samples)$powTrend
  #summary(powTrend)
  powTrendM=mean(powTrend)
  paste("powTrend",powTrendM)
  
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
  
  b=extract(samples)$b
  bM=apply(b,2,mean)
  lastB=b[,lastDataIndex]
  lastBM=mean(lastB)
  paste("lastB",round(lastBM,2))
  
  bSm = extract(samples)$bSm
  #summary(bSm)
  bSmM=mean(bSm)
  paste("bSm",round(bSmM,2))
  
  locTrendFract = extract(samples)$locTrendFract
  #print(summary(locTrendFract))
  locTrendFractM=mean(locTrendFract)
  paste("locTrendFract",round(locTrendFractM,2))
  
  lastDetails=paste("coefT=",round(coefTrendM,2), ", powT=",round(powTrendM,2),
      ", sigma=",round(sigmaM,2),", powx=",round(powxM,2),", offsetS=",round(offsetSigmaM,2), 
      ", nu=",round(nuM,2), ", lSm=",round(levSmM,2), 
      ", lastB=",round(lastBM,1), ", bSm=",round(bSmM,2), 
      ", lTFract=", round(locTrendFractM,2), sep='')				
  print(lastDetails)
  
  
  t=1; irun=1
  yf=matrix(lastLevelM,nrow=NUM_OF_TRIALS, ncol=maxPredictionHorizon)
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
  } #through trials
  yf		
}


rangeOfYearlySeries=1:645
idr=sample(rangeOfYearlySeries,1) # choose series

series=M3[[idr]]$sn
x=as.numeric(M3[[idr]]$x)
#n=M3[[idr]]$n
y=x+rnorm(n,0,sd=abs(min(x))*0.0001) # I found that adding a bit of jitter is helping Stan in case of some flat series		
# summary(y)		
yy=as.numeric(M3[[idr]]$xx)
maxPredictionHorizon=M3[[idr]]$h
plot(c(y,yy), main=series,  col=c(rep(1,n),rep(2,maxPredictionHorizon)), ylab='',xlab='years') #just on screen for debugging

yf=fitLgt(y)
avgYfs=apply(yf,2,quantile,probs=QUANTS)

ymax=max(c(y,yy),max(avgYfs[3,]))
ymin=min(min(c(y,yy)),min(avgYfs[1,]))
plot(c(y,yy), main=series,  col=c(rep(1,n),rep(2,maxPredictionHorizon)), ylab='',xlab=lastDetails,ylim=c(ymin,ymax) )
for (irun in 1:NUM_OF_TRIALS_TO_SHOW) 
	lines((n+1):(n+maxPredictionHorizon),yf[irun,], col='gray')
lines((n+1):(n+maxPredictionHorizon),avgYfs[1,], col='pink',lwd=2)
lines((n+1):(n+maxPredictionHorizon),avgYfs[2,], col='blue',lwd=2)
lines((n+1):(n+maxPredictionHorizon),avgYfs[3,], col='pink',lwd=2)
lines(c(y,yy), main=series,  col=c(rep(1,n),rep(2,maxPredictionHorizon)),type='b')

	