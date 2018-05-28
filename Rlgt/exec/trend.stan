//Local and Global Trend, homoscedastic model with Normal distribution of error

data {
	real<lower=0> CAUCHY_SD;
	real MIN_POW_TREND;  real MAX_POW_TREND;
	real<lower=0> MIN_SIGMA;
	int<lower=1> N;
	vector<lower=0>[N] y;
}

parameters {
	real<lower=MIN_SIGMA> sigma;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1> bSm;
	real <lower=0,upper=1> powTrendBeta;
	real coefTrend;
	real <lower=-0.25,upper=1> locTrendFract;
}

transformed parameters {
	real <lower=MIN_POW_TREND,upper=MAX_POW_TREND>powTrend;
	vector[N] l; vector[N] b;

	l[1] = y[1]; b[1] = 0;
	powTrend= (MAX_POW_TREND-MIN_POW_TREND)*powTrendBeta+MIN_POW_TREND;

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
