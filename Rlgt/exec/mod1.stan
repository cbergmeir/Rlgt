// Non-Seasonal Local Global Trend (LGT) algorithm

data {  
	real<lower=0> CAUCHY_SD;
	real<lower=0> MIN_SIGMA;
	int<lower=1> N;
	vector<lower=0>[N] y; 
}
parameters {
	real<lower=0.8,upper=0.98> psi; 
	real<lower=0> sigma;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1> bSm;
	real bInit;
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
    bInit ~ normal(0,CAUCHY_SD);
	
	for (t in 2:N) {
		y[t] ~ normal(l[t-1]+psi*b[t-1], sigma);
	}
}

