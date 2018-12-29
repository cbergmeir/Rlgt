// Seasonal Global Trend (SGT) algorithm

data {  
	int<lower=2> SEASONALITY;
	real<lower=SEASONALITY,upper=SEASONALITY+1> SEASONALITY_F;  //possibly non-integer seasonality. Normally SEASONALITY_F=SEASONALITY
	real<lower=0> CAUCHY_SD;
	real MIN_POW_TREND;  real MAX_POW_TREND;
	real<lower=0> MIN_SIGMA;
	real<lower=1> MIN_NU; real<lower=1> MAX_NU;
	int<lower=SEASONALITY+1> N;
	vector<lower=0>[N] y;
	real<lower=0> POW_TREND_ALPHA; real<lower=0> POW_TREND_BETA;
	real<lower=0> POW_SEASON_ALPHA; real<lower=0> POW_SEASON_BETA;
	int<lower=0,upper=1> USE_REGRESSION;
	int<lower=0,upper=1> SEASONALITY_TYPE;  //0- multiplicative, 1- generalized
	int<lower=0,upper=1> USE_SMOOTHED_ERROR;
	int<lower=0> NUM_OF_SEASON_INIT_CYCLES;
	int<lower=0,upper=3> LEVEL_CALC_METHOD;  //0-classical(HW), 2-average, 3- avg of 0 and 2  (yes 1 is missing :-) )
	int<lower=1> J;
	matrix[N, J] xreg;  
	vector<lower=0>[J] REG_CAUCHY_SD;
}
transformed data {
  real <lower=0,upper=1> fractSeasonality;
	real<lower=0> reg0CauchySd=mean(REG_CAUCHY_SD)*10;

	if (SEASONALITY_F>SEASONALITY) {
		fractSeasonality=SEASONALITY_F-SEASONALITY;
		//print("Non-integer seasonality used.");
	}
	else
		fractSeasonality=0;
}
parameters {
 	vector[J]  regCoef; real regOffset;
	real<lower=MIN_NU,upper=MAX_NU> nu; 
	real<lower=0> sigma;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1>llevSm;   //used for LEVEL_CALC_METHOD==3
	real <lower=0,upper=1>sSm;
	real <lower=0,upper=1>powx;
	real <lower=0,upper=1> powTrendBeta;
	real coefTrend;
	real <lower=MIN_SIGMA> offsetSigma;
	real <lower=0,upper=1>innovSm;
	real <lower=0> innovSizeInit;
	vector[SEASONALITY] initS;
	real <lower=0,upper=1> powSeason;
} 
transformed parameters {
	real <lower=MIN_POW_TREND,upper=MAX_POW_TREND>powTrend;
	vector<lower=0>[N] l;
	vector<lower=0>[N] l0;
	vector[N+SEASONALITY+1] s;  //1 extra in case of non-integer seasonality
	vector[N] r; //regression component
	vector<lower=0>[N] expVal; 
	vector<lower=0>[N] smoothedInnovSize;
	real seasonalityP;
	real sumsu;
	real newLevelP;
	real movingSum;
	
	if (USE_REGRESSION)
		r = xreg * regCoef + regOffset;
	else 
		r=rep_vector(0, N);	
		
	if (SEASONALITY_TYPE==0){
		sumsu = 0;
		for (i in 1:SEASONALITY) 
			sumsu = sumsu+ exp(initS[i]); 
		for (i in 1:SEASONALITY) {
			s[i] = exp(initS[i])*SEASONALITY/sumsu;
			//print(i," ",s[i]);
		}	
		l[1] = (y[1]-r[1])/s[1];  
	} else  {  //generalized
		for (i in 1:SEASONALITY) 
			s[i] = initS[i];
		l[1] = y[1] - r[1];  
	} 

	s[N+SEASONALITY+1]=1;  //in case of integer seasonality, the last value is not filled and Stan does not like it
	s[SEASONALITY+1] = s[1];
	s[SEASONALITY+2] = s[2]; //needed in case of non-integer seasonality, otherwise s[SEASONALITY+2] will get overwritten
	
	if (USE_SMOOTHED_ERROR)
	  smoothedInnovSize[1]=innovSizeInit;
	else
	  smoothedInnovSize[1]=1;  //has to have some value, not NA
	
	powTrend= (MAX_POW_TREND-MIN_POW_TREND)*powTrendBeta+MIN_POW_TREND;
	expVal[1] = y[1];
	
	if (LEVEL_CALC_METHOD==3) 
		l0[1]=l[1];
	else 
		l0=rep_vector(0, N);	//not used
	
	if (LEVEL_CALC_METHOD>0) {
		movingSum=y[1]-r[1];
		for (t in 2:SEASONALITY) 
			movingSum=movingSum+y[t]-r[t];
	}

	for (t in 2:N) {	
		if (LEVEL_CALC_METHOD>0 && t>SEASONALITY) 
			movingSum=movingSum+(y[t]-r[t])-(y[t-SEASONALITY]-r[t-SEASONALITY]);
					
		//expVal and level	
		if (SEASONALITY_TYPE==0) {//HW
			expVal[t]=(l[t-1]+ coefTrend*l[t-1]^powTrend)*s[t] + r[t];   //expVal[t] can't use y[t] or  anything derived from it
			newLevelP=(y[t]-r[t])/s[t];
		} else { //if (SEASONALITY_TYPE==1) 
			expVal[t]=l[t-1]+ coefTrend*l[t-1]^powTrend + s[t]*l[t-1]^powSeason + r[t];  //expVal[t] can't use y[t] or  anything derived from it 
			newLevelP=y[t] - r[t] - s[t]*l[t-1]^powSeason ;
		}
		 	
		//level cont
		if (LEVEL_CALC_METHOD==0) 
			l[t]  = levSm*newLevelP + (1-levSm)*l[t-1];
		else if (LEVEL_CALC_METHOD==2) {
		  if (t<=SEASONALITY)
		  	l[t]  = levSm*newLevelP + (1-levSm)*l[t-1];
		  else  
				l[t]= levSm*movingSum/SEASONALITY+(1-levSm)*l[t-1];
		}	
		else  if (LEVEL_CALC_METHOD==3) {
			l0[t]  = levSm*newLevelP + (1-levSm)*l0[t-1];
			if (t<=SEASONALITY)
			  l[t]=l0[t];
			else
			  l[t]  = llevSm*l0[t]+(1-llevSm)*movingSum/SEASONALITY;
		}
		 
		//seasonality
		if (SEASONALITY_TYPE==0) {//HW
			seasonalityP = sSm*(y[t]-r[t])/l[t]+(1-sSm)*s[t];
		} else { //if (SEASONALITY_TYPE==1) {//generalized
			seasonalityP=sSm*(y[t] - l[t] -r[t])/l[t-1]^powSeason + (1-sSm)*s[t];
		}
		//seasonality cont			
		if (fractSeasonality>0) {
			s[t+SEASONALITY+1]=seasonalityP;  //with fractSeasonality weight
			s[t+SEASONALITY]=fractSeasonality*s[t+SEASONALITY]+(1-fractSeasonality)*seasonalityP;
		} else
			s[t+SEASONALITY]=seasonalityP;
			
		//size of error	 
		if (USE_SMOOTHED_ERROR)
			smoothedInnovSize[t]=innovSm*fabs(y[t]-expVal[t])+(1-innovSm)*smoothedInnovSize[t-1];
		else	
			smoothedInnovSize[t]=1;  //has to have some value, not NA 
	}
}
model {
	sigma ~ cauchy(0,CAUCHY_SD) T[0,];
	offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];	
	coefTrend ~ cauchy(0, CAUCHY_SD);
	powTrendBeta ~ beta(POW_TREND_ALPHA, POW_TREND_BETA);
	levSm ~ beta(1, 2);
	
	if(USE_SMOOTHED_ERROR)
		innovSizeInit~ cauchy(y[1]/100,CAUCHY_SD) T[0,];
		
	if (USE_REGRESSION) {
		regCoef ~ cauchy(0, REG_CAUCHY_SD);
		regOffset ~ cauchy(0, reg0CauchySd);
	}		
	
	if (SEASONALITY_TYPE==0) {//HW
		for (t in 1:SEASONALITY) 
    		initS[t] ~ cauchy (0, 4);  //exp(8)=3000
	} else { //if (SEASONALITY_TYPE==1) {//generalized
		powSeason ~ beta(POW_SEASON_ALPHA, POW_SEASON_BETA); 
		for (t in 1:SEASONALITY)
			initS[t] ~ cauchy (0, y[t]*0.3);
	} 
		
	for (t in 2:N) {
	  if (USE_SMOOTHED_ERROR==0)
	  	y[t] ~ student_t(nu, expVal[t], sigma*expVal[t]^powx+ offsetSigma);  // expVal[t]^powx , not l[t-1]^powx
	  else
	  	y[t] ~ student_t(nu, expVal[t], sigma*smoothedInnovSize[t-1] + offsetSigma);
	}
}

