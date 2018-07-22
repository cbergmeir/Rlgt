// Non-Seasonal Holt-Winter Dampen algorithm

data {  
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
    real <lower=0,upper=1> powx;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1> bSm;
	real bInit;
	real <lower=MIN_SIGMA> offsetSigma;
} 
transformed parameters {
	vector[N] l; vector[N] b; 
	l[1] = y[1]; b[1] = bInit;

	
	for (t in 2:N) {
		l[t]  = levSm*y[t] + (1-levSm)*(l[t-1]+psi*b[t-1]) ;
		b[t]  = bSm*(l[t]-l[t-1]) + (1-bSm)*psi*b[t-1] ;
	}
}
model {
	sigma ~ cauchy(0,CAUCHY_SD) T[0,];
    offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
    bInit ~ normal(0,CAUCHY_SD);
    powx ~ beta(POW_SIGMA_ALPHA, POW_SIGMA_BETA);
	
	for (t in 2:N) {
		y[t] ~ student_t(nu, l[t-1]+psi*b[t-1], sigma*fabs(l[t-1])^powx+ offsetSigma);
	}
}

