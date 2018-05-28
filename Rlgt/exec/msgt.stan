data {
	int<lower=1> NUM_OF_S2_DAYS;
	real<lower=0> CAUCHY_SD;
	real MIN_POW;  real MAX_POW;
	real<lower=0> MIN_SIGMA;
	real<lower=1> MIN_NU; real<lower=1> MAX_NU;
	int<lower=1> N;
	vector<lower=0>[N] y;
	int<lower=0,upper=1> S2_INDEX[N];
	real<lower=0> TREND_ALPHA; real<lower=0> TREND_BETA;
	real<lower=0> POWX_ALPHA;  real<lower=0> POWX_BETA;
	int<lower=2> S2_SEASONALITY; //must be multiple of 7
	vector[S2_SEASONALITY] INIT_S2U; //unnormalized
}

transformed data {
	int<lower=2> YEARLY_SEASONALITY;
	int<lower=1> NWS; //length of vector of yearly seasonalities
    real NWS_real; //quick hack to avoid R check warning

	YEARLY_SEASONALITY=52;

	NWS_real =N/7.0+1+YEARLY_SEASONALITY;
    while (NWS<NWS_real){ //just to convert real back to int
        NWS=NWS+1;
    }
}

parameters {
	real<lower=MIN_NU,upper=MAX_NU> nu;
	real<lower=0> sigma;
	real <lower=0,upper=1>powx;
	real <lower=MIN_SIGMA> offsetSigma;
	real<lower=0> sigma2;
	real <lower=0,upper=1>powx2;
	real <lower=MIN_SIGMA> offsetSigma2;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1>sSm;
	real <lower=0,upper=1>s2Sm;
	real <lower=0,upper=1>sySm;
	real <lower=0,upper=1> powTrendBeta;
	real coefTrend;
	vector[7] initSu; //unnormalized
	vector[S2_SEASONALITY] initS2u; //unnormalized
	vector[YEARLY_SEASONALITY] initSyu; //unnormalized
}

transformed parameters {
	real <lower=MIN_POW,upper=MAX_POW>powTrend;
	vector[N] l;

	vector[7] inits;
	vector[N+7+1] s;

	vector[S2_SEASONALITY] inits2;
	vector[NUM_OF_S2_DAYS+S2_SEASONALITY+1] s2;

	vector[YEARLY_SEASONALITY] initSy;
	vector[NWS] sy;

	real sumsu; real sumy;

	sumsu = 0;
	for (i in 1:7)
		sumsu = sumsu+ initSu[i];
	for (i in 1:7)
		inits[i] = initSu[i]*7/sumsu;
	for (i in 1:7) {
		s[i] = inits[i];
	}
	s[7+1] = inits[1];


	sumsu = 0;
	for (i in 1:S2_SEASONALITY)
		sumsu = sumsu+ initS2u[i];
	for (i in 1:S2_SEASONALITY)
		inits2[i] = initS2u[i]*S2_SEASONALITY/sumsu;
	for (i in 1:S2_SEASONALITY) {
		s2[i] = inits2[i]*INIT_S2U[i];
	}

	sumsu = 0;
	for (i in 1:YEARLY_SEASONALITY)
		sumsu = sumsu+ initSyu[i];
	for (i in 1:YEARLY_SEASONALITY)
		initSy[i] = initSyu[i]*YEARLY_SEASONALITY/sumsu;
	for (i in 1:YEARLY_SEASONALITY) {
		sy[i] = initSy[i];
	}


	if (S2_INDEX[1]==0) {
		l[1] = y[1]/(s[1]*sy[1]);
	} else {
		l[1] = y[1]/(s2[1]*sy[1]);
	}
	powTrend= (MAX_POW-MIN_POW)*powTrendBeta+MIN_POW;

	{
		int is2; int iy;
		sumy=0; is2=0; iy=1;
		for (t in 2:N) {
			if (S2_INDEX[t]==0) {
				l[t]  = levSm*y[t]/(s[t]*sy[iy]) + (1-levSm)*l[t-1] ;
				s[t+7] = sSm*y[t]/(l[t]*sy[iy])+(1-sSm)*s[t];
				sumy=sumy+y[t]/(l[t]*s[t]);
			} else {//this code will not work if the first day of the series (in-sample) is in inside the Thanksgiving week. This possibility needs be removed in the calling R code.
				is2=is2+1;
				l[t]  = levSm*y[t]/(s2[is2]*sy[iy]) + (1-levSm)*l[t-1] ;
				s2[is2+S2_SEASONALITY] = s2Sm*y[t]/(l[t]*sy[iy])+(1-s2Sm)*s2[is2];
				s[t+7] = s[t];
				sumy=sumy+y[t]/(l[t]*s2[is2]);
			}
			if ((t-1)%7==0) {
				sy[iy+YEARLY_SEASONALITY] = sySm*sumy/7+(1-sySm)*sy[iy];
				sumy=0;
				iy=iy+1;
			}
		}
	}
}

model {
	real expVal;real riy; int is2; int iy;

	powTrendBeta ~ beta(TREND_ALPHA, TREND_BETA);
	powx ~  beta(POWX_ALPHA, POWX_BETA);
	sigma ~ cauchy(0,CAUCHY_SD) T[0,];
	offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];
	sigma2 ~ cauchy(0,5*CAUCHY_SD) T[0,];
	offsetSigma2 ~ cauchy(MIN_SIGMA,5*CAUCHY_SD) T[MIN_SIGMA,];
	coefTrend ~ cauchy(0, CAUCHY_SD);

	for (t in 1:7)
		initSu[t] ~ normal (1, 0.2) T[0.05,];
	for (t in 1:S2_SEASONALITY)
		initS2u[t] ~ normal (1, 0.2) T[0.1,];//this is just a disturbance on top of externally-provided initialization
	for (t in 1:YEARLY_SEASONALITY)
		initSyu[t]~ normal (1, 0.3) T[0.05,];

	is2=0;
	for (t in 2:N) {
		riy=(t-1)/7.0+1;
        while (iy<riy){ //just to convert real back to int
            iy=iy+1;
        }
		if (S2_INDEX[t]==0) {
			expVal = (l[t-1]+ coefTrend*fabs(l[t-1])^powTrend)*s[t]*sy[iy];
			y[t] ~ student_t(nu, expVal, sigma*fabs(expVal)^powx+ offsetSigma);
		} else {
			is2=is2+1;
			expVal = (l[t-1]+ coefTrend*fabs(l[t-1])^powTrend)*s2[is2]*sy[iy];
			y[t] ~ student_t(nu, expVal, sigma2*fabs(expVal)^powx2+ offsetSigma2);
		}
	}
}

