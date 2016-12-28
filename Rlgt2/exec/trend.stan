data {
	real<lower=0> CAUCHY_SD;
	real MIN_POW;  real MAX_POW;
	real<lower=0> MIN_SIGMA;
	//real<lower=1> MIN_NU; real<lower=1> MAX_NU;
	int<lower=1> N;
	vector<lower=0>[N] y;
}

parameters {
	real<lower=MIN_SIGMA> sigma;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1> bSm;
	real <lower=0,upper=1> powTrendBeta;
	real coefTrend;
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
	sigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
	coefTrend ~ cauchy(0,CAUCHY_SD);

	for (t in 2:N) {
		y[t] ~ normal(l[t-1]+coefTrend*fabs(l[t-1])^powTrend+locTrendFract*b[t-1], sigma);
	}
}
