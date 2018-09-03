// non-seasonal Local and Global Trend (LGT) algorithm

data {  
	real<lower=0> CAUCHY_SD;
	real MIN_POW_TREND;  real MAX_POW_TREND;
	real<lower=0> MIN_SIGMA;
	real<lower=1> MIN_NU; real<lower=1> MAX_NU;
	int<lower=1> N;
	vector<lower=0>[N] y;
	real<lower=0> POW_TREND_ALPHA; real<lower=0> POW_TREND_BETA;
	int<lower=0,upper=1> USE_REGRESSION;
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
	real <lower=0,upper=1> bSm;
	real <lower=0,upper=1> powx;
	real bInit;
	real <lower=0,upper=1> powTrendBeta;
	real coefTrend;
	real <lower=MIN_SIGMA> offsetSigma;
	real <lower=0,upper=1> locTrendFract;
} 
transformed parameters {
	real <lower=MIN_POW_TREND,upper=MAX_POW_TREND>powTrend;
	vector<lower=0>[N] l; 
	vector[N] b;
	real r; //regression component
	vector<lower=0>[N] expVal; 
	
	l[1] = y[1]; b[1] = bInit;
	powTrend= (MAX_POW_TREND-MIN_POW_TREND)*powTrendBeta+MIN_POW_TREND;
	expVal[1] = y[1];
	r=0;

	for (t in 2:N) {
		if (USE_REGRESSION==1)
			r = xreg[t,:] * regCoef + regOffset;
		expVal[t]=l[t-1]+coefTrend*l[t-1]^powTrend+locTrendFract*b[t-1]+r;
		l[t] = levSm*(y[t]-r) + (1-levSm)*l[t-1] ;  
		b[t] = bSm*(l[t]-l[t-1]) + (1-bSm)*b[t-1] ;
	}
}
model {
	sigma ~ cauchy(0,CAUCHY_SD) T[0,];
	offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
	coefTrend ~ cauchy(0,CAUCHY_SD);
	powTrendBeta ~ beta(POW_TREND_ALPHA, POW_TREND_BETA);
	bInit ~ cauchy(0,CAUCHY_SD);
	regCoef ~ cauchy(0, REG_CAUCHY_SD);
	regOffset ~ cauchy(0, reg0CauchySd);
	
	for (t in 2:N) {
		y[t] ~ student_t(nu, expVal[t], sigma * l[t-1] ^ powx + offsetSigma) ;
	}
}
