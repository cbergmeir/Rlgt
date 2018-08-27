// Dual Seasonal Global Trend (S2GT) algorithm, v2. Proposed level equals avg of last MAX_SEASONALITY points 

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
}
transformed data {
  int<lower=2> MAX_SEASONALITY;

  MAX_SEASONALITY=SEASONALITY;
  if (MAX_SEASONALITY<SEASONALITY2)
    MAX_SEASONALITY=SEASONALITY2;
}
parameters {
	real<lower=MIN_NU,upper=MAX_NU> nu; 
	real<lower=0> sigma;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1>sSm;
	real <lower=0,upper=1>s2Sm;
	real <lower=0,upper=1>powx;
	real <lower=0,upper=1> powTrendBeta;
	real coefTrend;
	real <lower=MIN_SIGMA> offsetSigma;
	vector[SEASONALITY] initSu; //unnormalized
	vector[SEASONALITY2] initSu2; //unnormalized
} 
transformed parameters {
	real <lower=MIN_POW_TREND,upper=MAX_POW_TREND>powTrend;
	vector<lower=0>[N] l;
	vector<lower=0>[N+SEASONALITY] s;
	vector<lower=0>[N+SEASONALITY2] s2;
	real sumsu;
	real newLevelP;
	real movingSum;
	
	sumsu = 0;
	for (i in 1:SEASONALITY) 
		sumsu = sumsu+ initSu[i];
	for (i in 1:SEASONALITY) 
		s[i] = initSu[i]*SEASONALITY/sumsu;
	s[SEASONALITY+1] = s[1];
	
	sumsu = 0;
	for (i in 1:SEASONALITY2) 
		sumsu = sumsu+ initSu2[i];
	for (i in 1:SEASONALITY2) 
		s2[i] = initSu2[i]*SEASONALITY2/sumsu;
	s2[SEASONALITY2+1] = s2[1];
	
	//l[1] = y[1]/(s[1]*s2[1]);
	powTrend= (MAX_POW_TREND-MIN_POW_TREND)*powTrendBeta+MIN_POW_TREND;
	
	movingSum=y[1];
	for (t in 2:MAX_SEASONALITY) 
		movingSum=movingSum+y[t];
	newLevelP=movingSum/MAX_SEASONALITY;
	l[1] =newLevelP;
	
	for (t in 2:N) {
		if (t>MAX_SEASONALITY) {
			movingSum=movingSum+y[t]-y[t-MAX_SEASONALITY];
			newLevelP=movingSum/MAX_SEASONALITY;
			l[t]  = levSm*newLevelP + (1-levSm)*l[t-1] ;
		} else {
			l[t]=newLevelP; //same starting level
		} 
		s[t+SEASONALITY] = sSm*y[t]/(l[t]*s2[t])+(1-sSm)*s[t];
		s2[t+SEASONALITY2] = s2Sm*y[t]/(l[t]*s[t])+(1-s2Sm)*s2[t];
	}
}
model {
	real expVal;

	sigma ~ cauchy(0,CAUCHY_SD) T[0,];
	offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];	
	coefTrend ~ cauchy(0, CAUCHY_SD);
	powTrendBeta ~ beta(POW_TREND_ALPHA, POW_TREND_BETA);

	for (t in 1:SEASONALITY)
		initSu[t] ~ cauchy (1, 0.3) T[0.01,];
	for (t in 1:SEASONALITY2)
		initSu2[t] ~ cauchy (1, 0.3) T[0.01,];
	
	for (t in 2:N) {
	  expVal = (l[t-1]+ coefTrend*l[t-1]^powTrend)*s[t]*s2[t];
	  y[t] ~ student_t(nu, expVal, sigma*expVal^powx+ offsetSigma);
	}
}
