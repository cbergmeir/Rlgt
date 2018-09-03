// Dual Seasonal Global Trend (S2GT) algorithm, using smoothed Error size (S2GTe)

data {  
	int<lower=2> SEASONALITY;
	int<lower=2> SEASONALITY2;
	real<lower=0> CAUCHY_SD;
	real MIN_POW_TREND;  real MAX_POW_TREND;
	real<lower=0> MIN_SIGMA;
	real<lower=1> MIN_NU; real<lower=1> MAX_NU;
	int<lower=SEASONALITY+SEASONALITY2> N;
	vector<lower=0>[N] y;
	real<lower=0> POW_TREND_ALPHA; real<lower=0> POW_TREND_BETA;
	real<lower=0> POW_SEASON_ALPHA; real<lower=0> POW_SEASON_BETA;
	int<lower=0,upper=1> USE_REGRESSION;
	int<lower=0,upper=1> USE_GENERALIZED_SEASONALITY;
	int<lower=1> J;
	matrix[N, J] xreg;  
	vector<lower=0>[J] REG_CAUCHY_SD;
}
transformed data {
  real<lower=0> reg0CauchySd=mean(REG_CAUCHY_SD)*10;
}
parameters {
 	vector[J]  regCoef; real regOffset;
	real<lower=MIN_NU,upper=MAX_NU> nu; 
	real<lower=0> sigma;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1>sSm;
	real <lower=0,upper=1>s2Sm;
	real <lower=0,upper=1>powx;
	real <lower=0,upper=1> powTrendBeta;
	real coefTrend;
	real <lower=MIN_SIGMA> offsetSigma;
	real <lower=0,upper=1>innovSm;
	real <lower=0> innovSizeInit;
	vector[SEASONALITY] initSu;
	real <lower=0,upper=1> powSeason;
	vector[SEASONALITY2] initSu2;
	real <lower=0,upper=1> powSeason2;
} 
transformed parameters {
	real <lower=MIN_POW_TREND,upper=MAX_POW_TREND>powTrend;
	vector<lower=0>[N] l;
	vector[N+SEASONALITY] s;
	vector[N+SEASONALITY2] s2;
	real r; //regression component
	vector<lower=0>[N] expVal; 
	vector<lower=0>[N] smoothedInnovSize;
	real sumsu;
	
	r=0;
	if (USE_REGRESSION==1)
		r = xreg[1,:] * regCoef + regOffset;
		
	if (USE_GENERALIZED_SEASONALITY==1) {
		for (i in 1:SEASONALITY) 
    		s[i] = initSu[i];
    	for (i in 1:SEASONALITY2) 
    		s2[i] = initSu2[i];	
    	l[1] = y[1] - r;
	} else {
		sumsu = 0;
		for (i in 1:SEASONALITY) 
			sumsu = sumsu+ initSu[i];
		for (i in 1:SEASONALITY) 
			s[i] = initSu[i]*SEASONALITY/sumsu;	
		
		sumsu = 0;
		for (i in 1:SEASONALITY2) 
			sumsu = sumsu+ initSu2[i];
		for (i in 1:SEASONALITY2) 
			s2[i] = initSu2[i]*SEASONALITY2/sumsu;
			
		l[1] = (y[1]-r)/(s[1]*s2[1]);
	}
	s[SEASONALITY+1] = s[1];
	s2[SEASONALITY2+1] = s2[1];
	
	smoothedInnovSize[1]=innovSizeInit;
	
	powTrend= (MAX_POW_TREND-MIN_POW_TREND)*powTrendBeta+MIN_POW_TREND;
	expVal[1] = y[1];
	
	for (t in 2:N) {
		if (USE_REGRESSION==1)
			r = xreg[t,:] * regCoef + regOffset;
		if (USE_GENERALIZED_SEASONALITY==1) {
		    l[t]  = levSm*(y[t] - s[t]*l[t-1]^powSeason - s2[t]*l[t-1]^powSeason2 -r) + (1-levSm)*l[t-1] ;  //As usually, we skip global trend in the level update formula. Why? Becasue it works better :-)
    		s[t+SEASONALITY]= sSm*(y[t] - l[t-1] - coefTrend*l[t-1]^powTrend - s2[t]*l[t-1]^powSeason2 - r)/l[t-1]^powSeason + (1-sSm)*s[t]; 
    		s2[t+SEASONALITY2]= s2Sm*(y[t] - l[t-1] - coefTrend*l[t-1]^powTrend - s[t]*l[t-1]^powSeason - r)/l[t-1]^powSeason2 + (1-s2Sm)*s2[t]; 
    		expVal[t]=l[t-1]+ coefTrend*l[t-1]^powTrend + s[t]*l[t-1]^powSeason + s2[t]*l[t-1]^powSeason2 + r;
    		smoothedInnovSize[t]=innovSm*fabs(y[t]-expVal[t])+(1-innovSm)*smoothedInnovSize[t-1];
		} else {	
			l[t]  = levSm*(y[t]-r)/(s[t]*s2[t]) + (1-levSm)*l[t-1];
			s[t+SEASONALITY] = sSm*(y[t]-r)/(l[t]*s2[t])+(1-sSm)*s[t];
			s2[t+SEASONALITY2] = s2Sm*(y[t]-r)/(l[t]*s[t])+(1-s2Sm)*s2[t];
			expVal[t]=(l[t-1]+ coefTrend*l[t-1]^powTrend)*s[t]*s2[t] + r;
			smoothedInnovSize[t]=innovSm*fabs(y[t]-expVal[t])/s[t]+(1-innovSm)*smoothedInnovSize[t-1];
		}
	}
}
model {
	sigma ~ cauchy(0,CAUCHY_SD) T[0,];
	offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];	
	coefTrend ~ cauchy(0, CAUCHY_SD);
	innovSizeInit~ cauchy(y[1]/100,CAUCHY_SD) T[0,];
	powTrendBeta ~ beta(POW_TREND_ALPHA, POW_TREND_BETA);
	powSeason ~ beta(POW_SEASON_ALPHA, POW_SEASON_BETA);
	powSeason2 ~ beta(POW_SEASON_ALPHA, POW_SEASON_BETA);
	
	if (USE_GENERALIZED_SEASONALITY==1) {
		for (t in 1:SEASONALITY)
			initSu[t] ~ cauchy (0, y[t]*0.1);	
		for (t in 1:SEASONALITY2)
			initSu2[t] ~ cauchy (0, y[t]*0.1);		
	} else {
		for (t in 1:SEASONALITY) 
    		initSu[t] ~ cauchy (1, 0.3) T[0.01,];
		for (t in 1:SEASONALITY2)
			initSu2[t] ~ cauchy (1, 0.3) T[0.01,];
	}
	regCoef ~ cauchy(0, REG_CAUCHY_SD);
	regOffset ~ cauchy(0, reg0CauchySd);
	
	for (t in 2:N) {
	  y[t] ~ student_t(nu, expVal[t], sigma*smoothedInnovSize[t-1] + offsetSigma);
	}
}
