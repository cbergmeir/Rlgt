# Non-Seasonal Local Global Trend (LGT) algorithm
# Dec 2016
# Slawek
#############################################

library(rstan)

predictionAlgorithm=paste('LGT');

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
		real<lower=0> POW_TREND_ALPHA; real<lower=0> POW_TREND_BETA; 
		real<lower=0> POWX_ALPHA; real<lower=0> POWX_BETA; 
	}
	parameters {
		real<lower=MIN_NU,upper=MAX_NU> nu; 
		real<lower=0> sigma;
		real <lower=0,upper=1>levSm;
		real <lower=0,upper=1> bSm;
		real <lower=0,upper=1> powx;
		real bInit;
		real <lower=0,upper=1> powTrendBeta;
		real coefTrend;
		real <lower=MIN_SIGMA> offsetSigma;
		real <lower=-0.25,upper=1> locTrendFract;
	} 
	transformed parameters {
		real <lower=MIN_POW,upper=MAX_POW>powTrend;
		vector[N] l; vector[N] b;
		
		l[1] = y[1]; b[1] = bInit;
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
		powTrendBeta ~ beta(POW_TREND_ALPHA, POW_TREND_BETA);
    powx ~ beta(POWX_ALPHA, POWX_BETA);
    bInit ~ normal(0,CAUCHY_SD);
		
		for (t in 2:N) {
			y[t] ~ student_t(nu, l[t-1]+coefTrend*fabs(l[t-1])^powTrend+locTrendFract*b[t-1], 
				sigma*fabs(l[t-1])^powx+ offsetSigma);
		}
	}
'    
stanModel = stan_model(model_code=modelX)
#str(stanModel)
