// Seasonal Holt-Winter Dampen algorithm

data {  
    int<lower=2> SEASONALITY;
	real<lower=0> CAUCHY_SD;
	real<lower=0> MIN_SIGMA;
	int<lower=1> N;
	vector<lower=0>[N] y; 
    real<lower=1> MIN_NU; real<lower=1> MAX_NU;
    real<lower=0> POW_SIGMA_ALPHA; real<lower=0> POW_SIGMA_BETA; 
}
parameters {
    real<lower=MIN_NU,upper=MAX_NU> nu;
	real<lower=0.8,upper=0.98> psi; 
	real<lower=0> sigma;
    real <lower=0,upper=1>powx;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1> bSm;
    real <lower=0,upper=1>sSm;
	real bInit;
    real <lower=MIN_SIGMA> offsetSigma;
    vector[SEASONALITY] initSu; //unnormalized
} 
transformed parameters {
	vector[N] l; vector[N] b; 
    vector[SEASONALITY] inits;
	vector[N+SEASONALITY] s;
	real sumsu;
	
	l[1] = y[1]; b[1] = bInit;
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
	
	for (t in 2:N) {
		l[t]  = levSm*y[t]/(s[t]) + (1-levSm)*(l[t-1]+psi*b[t-1]) ;
        b[t]  = bSm*(l[t]-l[t-1]) + (1-bSm)*psi*b[t-1] ;  
		s[t+SEASONALITY] = sSm*y[t]/(l[t])+(1-sSm)*s[t];
	}
}

model {
    real expVal;
	sigma ~ cauchy(0,CAUCHY_SD) T[0,];
    bInit ~ normal(0,CAUCHY_SD);
    powx ~ beta(POW_SIGMA_ALPHA, POW_SIGMA_BETA);
    offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
	for (t in 1:SEASONALITY) {
		initSu[t] ~ normal (1, 0.3) T[0.01,];
	}
	for (t in 2:N) {
        expVal = (l[t-1]+psi*b[t-1])*s[t];
		y[t] ~ student_t(nu, expVal, sigma*fabs(expVal)^powx+ offsetSigma);
	}
}

