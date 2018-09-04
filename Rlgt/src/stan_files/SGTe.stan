// Seasonal Global Trend algorithm, using smoothed Error size (SGTe)

data {  
	int<lower=2> SEASONALITY;
	real<lower=0> CAUCHY_SD;
	real MIN_POW_TREND;  real MAX_POW_TREND;
	real<lower=0> MIN_SIGMA;
	real<lower=1> MIN_NU; real<lower=1> MAX_NU;
	int<lower=SEASONALITY> N;
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
	real <lower=0,upper=1> powTrendBeta;
	real coefTrend;
	real <lower=MIN_SIGMA> offsetSigma;
	real <lower=0,upper=1>innovSm;
	real <lower=0> innovSizeInit;
	vector[SEASONALITY] initSu;
	real <lower=0,upper=1> powSeason;
} 
transformed parameters {
	real <lower=MIN_POW_TREND,upper=MAX_POW_TREND>powTrend;
	vector<lower=0>[N] l;
	vector[N+SEASONALITY] s;
	vector[N] r; //regression component
	vector<lower=0>[N] expVal; 
	vector<lower=0>[N] smoothedInnovSize;
	real sumsu;
	
	if (USE_REGRESSION==1)
		r = xreg * regCoef + regOffset;
	else 
		r=rep_vector(0, N);	
		
	if (USE_GENERALIZED_SEASONALITY==1) {
		for (i in 1:SEASONALITY) 
    		s[i] = initSu[i];
    	l[1] = y[1] - r[1];
	} else {
		sumsu = 0;
		for (i in 1:SEASONALITY) 
			sumsu = sumsu+ initSu[i];
		for (i in 1:SEASONALITY) 
			s[i] = initSu[i]*SEASONALITY/sumsu;	
		l[1] = (y[1]-r[1])/s[1];
	}
	s[SEASONALITY+1] = s[1];
	
	smoothedInnovSize[1]=innovSizeInit;
		
	powTrend= (MAX_POW_TREND-MIN_POW_TREND)*powTrendBeta+MIN_POW_TREND;
	expVal[1] = y[1];

	for (t in 2:N) {
		if (USE_GENERALIZED_SEASONALITY==1) {
		    l[t]  = levSm*(y[t] - s[t]*l[t-1]^powSeason -r[t]) + (1-levSm)*l[t-1] ;  //As usually, we skip global trend in the level update formula. Why? Becasue it works better :-)
    		s[t+SEASONALITY]= sSm*(y[t] - l[t-1]- coefTrend*l[t-1]^powTrend -r[t])/l[t-1]^powSeason + (1-sSm)*s[t]; 
    		expVal[t]=l[t-1]+ coefTrend*l[t-1]^powTrend + s[t]*l[t-1]^powSeason + r[t];
    		smoothedInnovSize[t]=innovSm*fabs(y[t]-expVal[t])+(1-innovSm)*smoothedInnovSize[t-1];
		} else {	
			l[t]  = levSm*(y[t]-r[t])/(s[t]) + (1-levSm)*l[t-1];
			s[t+SEASONALITY] = sSm*(y[t]-r[t])/l[t]+(1-sSm)*s[t];
			expVal[t]=(l[t-1]+ coefTrend*l[t-1]^powTrend)*s[t] + r[t];
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
	
	if (USE_GENERALIZED_SEASONALITY==1) 
		for (t in 1:SEASONALITY)
			initSu[t] ~ cauchy (0, y[t]*0.1);	
	else
		for (t in 1:SEASONALITY) 
    		initSu[t] ~ cauchy (1, 0.3) T[0.01,];
    		
	regCoef ~ cauchy(0, REG_CAUCHY_SD);
	regOffset ~ cauchy(0, reg0CauchySd);
	
	
	for (t in 2:N) {
	  y[t] ~ student_t(nu, expVal[t], sigma*smoothedInnovSize[t-1] + offsetSigma);
	}
}
