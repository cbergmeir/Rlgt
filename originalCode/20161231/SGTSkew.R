# Seasonal Global Trend (SGT) algorithm, skew
# Dec 2016
# Slawek
#############################################

library(rstan)

predictionAlgorithm=paste('SGTSkew');

#seasonal, standard SGT with function for error size
parameters = c("l", "s", "sSm","nu", "sigma", "levSm", 
		"powx", "coefTrend", "powTrend", "offsetSigma")
modelS = '
  functions { 
	  real skew_student_t_log(real y, real nu, real mu, real sigma, real skew) {
		  real z; real zc; 
		  if (sigma <= 0)
		    reject("Scale has to be positive.  Found sigma=", sigma);

      z= (y-mu)/sigma;
      zc= skew*z*sqrt((nu+1)/(nu+square(z)));
		  return -log(sigma) + student_t_lpdf(z | nu, 0, 1)
                       + student_t_lcdf(zc | nu+1, 0, 1);
		}
	}
	data {  
		int<lower=2> SEASONALITY;
		real<lower=0> CAUCHY_SD;
		real MIN_POW;  real MAX_POW;
		real<lower=0> MIN_SIGMA;
		real<lower=1> MIN_NU; real<lower=1> MAX_NU;
		int<lower=1> N;
		vector<lower=0>[N] y;
		real<lower=0> POW_TREND_ALPHA; real<lower=0> POW_TREND_BETA; 
		real<lower=0> POWX_ALPHA; real<lower=0> POWX_BETA; 
		real SKEW; 
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
		vector[SEASONALITY] initSu; //unnormalized
	} 
	transformed parameters {
		real <lower=MIN_POW,upper=MAX_POW>powTrend;
		vector[N] l;
		vector[SEASONALITY] inits;
		vector[N+SEASONALITY] s;
		real sumsu;
		
		sumsu = 0;
		for (i in 1:SEASONALITY) 
			sumsu = sumsu+ initSu[i];
		for (i in 1:SEASONALITY) 
			inits[i] = initSu[i]*SEASONALITY/sumsu;
		
		for (i in 1:SEASONALITY) {
			s[i] = inits[i];
		}
		s[SEASONALITY+1] = inits[1];
		
		l[1] = y[1]/s[1];
		powTrend= (MAX_POW-MIN_POW)*powTrendBeta+MIN_POW;
		
		for (t in 2:N) {
			l[t]  = levSm*y[t]/(s[t]) + (1-levSm)*l[t-1] ;  
			s[t+SEASONALITY] = sSm*y[t]/l[t]+(1-sSm)*s[t];
		}
	}
	model {
		real expVal;

		sigma ~ cauchy(0,CAUCHY_SD) T[0,];
		offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
		coefTrend ~ cauchy(0, CAUCHY_SD);
	  powTrendBeta ~ beta(POW_TREND_ALPHA, POW_TREND_BETA);
    powx ~ beta(POWX_ALPHA, POWX_BETA);
		
		for (t in 1:SEASONALITY) {
			initSu[t] ~ normal (1, 0.3) T[0.01,];
		}
		
		for (t in 2:N) {
		  expVal = (l[t-1]+ coefTrend*fabs(l[t-1])^powTrend)*s[t];
		  y[t] ~ skew_student_t(nu, expVal, sigma*fabs(expVal)^powx+ offsetSigma, SKEW);
		}
	}
' 
stanModel = stan_model(model_code=modelS)
#str(stanModelS)	