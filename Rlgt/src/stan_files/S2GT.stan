// Dual Seasonal Global Trend (S2GT) algorithm. 

data {  
	int<lower=2> SEASONALITY;
	real<lower=SEASONALITY,upper=SEASONALITY+1> SEASONALITY_F;  //possibly non-integer seasonality. Normally SEASONALITY_F=SEASONALITY
	int<lower=SEASONALITY+1> SEASONALITY2;  //here checking that SEASONALITY<SEASONALITY2
	real<lower=SEASONALITY2,upper=SEASONALITY2+1> SEASONALITY2_F;  //possibly non-integer seasonality. Normally SEASONALITY2_F=SEASONALITY2
	real<lower=0> CAUCHY_SD;
	real MIN_POW_TREND;  real MAX_POW_TREND;
	real<lower=0> MIN_SIGMA;
	real<lower=1> MIN_NU; real<lower=1> MAX_NU;
	int<lower=SEASONALITY+SEASONALITY2> N;
	vector<lower=0>[N] y;
	real<lower=0> POW_TREND_ALPHA; real<lower=0> POW_TREND_BETA;
	real<lower=0> POW_SEASON_ALPHA; real<lower=0> POW_SEASON_BETA;
	int<lower=0,upper=1> USE_REGRESSION;
	int<lower=0,upper=1> USE_GENERALIZED_SEASONALITY;
	int<lower=0,upper=1> USE_SMOOTHED_ERROR;
	int<lower=0> NUM_OF_SEASON_INIT_CYCLES;
	int<lower=0,upper=1> LEVEL_CALC_METHOD;  //0-classical, 1-avg over largest SEASONALITY 
	int<lower=1> J;
	matrix[N, J] xreg;  
	vector<lower=0>[J] REG_CAUCHY_SD;
}
transformed data {
	real <lower=0,upper=1> fractSeasonality;
	real <lower=0,upper=1> fractSeasonality2;
	real<lower=0> reg0CauchySd=mean(REG_CAUCHY_SD)*10;
	vector[SEASONALITY] firstRatios;    
	vector[SEASONALITY2] firstRatios2;
	real sumy; int j;   
	
	for (i in 1:SEASONALITY)
		firstRatios[i]=1;	
					
	sumy = 0; j=1;
 	while(j<=NUM_OF_SEASON_INIT_CYCLES && j*SEASONALITY<=N)  {
    sumy=0; 
		for (i in 1:SEASONALITY) 
			sumy = sumy+ y[(j-1)*SEASONALITY+i];
		for (i in 1:SEASONALITY) 
		  if (j==1) 
		  	firstRatios[i] = y[(j-1)*SEASONALITY+i]*SEASONALITY/sumy;		//at this stage we do not have access to the regression
		  else
			firstRatios[i] = firstRatios[i]+y[(j-1)*SEASONALITY+i]*SEASONALITY/sumy;	   
		j=j+1;
	}
	if (j>1) {
		j=j-1;
		for (i in 1:SEASONALITY) {
			firstRatios[i]=firstRatios[i]/j;
			//print(firstRatios[i]) ;	
		}
	}	
	
	
	for (i in 1:SEASONALITY2)
		firstRatios2[i]=1;	
					
	sumy = 0; j=1;
 	while(j<=NUM_OF_SEASON_INIT_CYCLES && j*SEASONALITY2<=N)  {
    	sumy=0; 
		for (i in 1:SEASONALITY2) 
			sumy = sumy+ y[(j-1)*SEASONALITY2+i];
		for (i in 1:SEASONALITY2) 
		  if (j==1) 
		  	firstRatios2[i] = y[(j-1)*SEASONALITY2+i]*SEASONALITY2/sumy/firstRatios[(i-1)%SEASONALITY+1];		//at this stage we do not have access to the regression
		  else
			firstRatios2[i] = firstRatios2[i]+y[(j-1)*SEASONALITY2+i]*SEASONALITY2/sumy/firstRatios[(i-1)%SEASONALITY+1];
		j=j+1;
	}
	if (j>1) {
		j=j-1;
		for (i in 1:SEASONALITY2) {
			firstRatios2[i]=firstRatios2[i]/j;
			//print(i, " ",firstRatios2[i]);
		}
	}
		
	if (SEASONALITY_F>SEASONALITY) {
		fractSeasonality=SEASONALITY_F-SEASONALITY;
		print("Non-integer seasonality used.");
	}
	else
		fractSeasonality=0;
		
	if (SEASONALITY2_F>SEASONALITY2) {
		fractSeasonality2=SEASONALITY2_F-SEASONALITY2;
		print("Non-integer seasonality2 used.");
	}
	else
		fractSeasonality2=0;	
}
parameters {
 	vector[J]  regCoef; real regOffset;
	real<lower=MIN_NU,upper=MAX_NU> nu; 
	real<lower=0> sigma;
	real <lower=0,upper=1>levSm;
	real <lower=0,upper=1>sSm;
	real <lower=0,upper=1>s2Sm;
	real <lower=0,upper=1>powx;
	real <lower=0,upper=1> powTrendBeta;
	real coefTrend;
	real <lower=MIN_SIGMA> offsetSigma;
	real <lower=0,upper=1>innovSm;
	real <lower=0> innovSizeInit;
	vector[SEASONALITY] initSu;
	real <lower=0,upper=1> powSeason;
	vector[SEASONALITY2] initSu2;
	real <lower=0,upper=1> powSeason2;
} 
transformed parameters {
	real <lower=MIN_POW_TREND,upper=MAX_POW_TREND>powTrend;
	vector<lower=0>[N] l;
	vector[N+SEASONALITY+1] s;  //1 extra in case of non-integer seasonality
	vector[N+SEASONALITY2+1] s2;
	vector[N] r; //regression component
	vector<lower=0>[N] expVal; 
	vector<lower=0>[N] smoothedInnovSize;
	real seasonalityP; real seasonalityP2;
	real sumsu;
	real newLevelP;
	real movingSum=0;
	
	if (USE_REGRESSION)
		r = xreg * regCoef + regOffset;
	else 
		r=rep_vector(0, N);	
		
	if (USE_GENERALIZED_SEASONALITY) {
		for (i in 1:SEASONALITY) 
    		s[i] = initSu[i];
    	for (i in 1:SEASONALITY2) 
    		s2[i] = initSu2[i];	
    	l[1] = y[1] - r[1];
	} else {
		sumsu = 0;
		for (i in 1:SEASONALITY) 
			sumsu = sumsu+ initSu[i];
		for (i in 1:SEASONALITY) 
			s[i] = firstRatios[i]*initSu[i]*SEASONALITY/sumsu;	
		
		sumsu = 0;
		for (i in 1:SEASONALITY2) 
			sumsu = sumsu+ initSu2[i];
		for (i in 1:SEASONALITY2) 
			s2[i] = firstRatios2[i]*initSu2[i]*SEASONALITY2/sumsu;
			
		l[1] = (y[1]-r[1])/(s[1]*s2[1]);  //initialization for LEVEL_CALC_METHOD==0, it will get overwritten otherwise
	}
	s[N+SEASONALITY+1]=1;  //for integer seasonality the last value is not filled and Stan does not like it
	s[SEASONALITY+1] = s[1];
	s[SEASONALITY+2] = s[2]; //needed in case of non-integer seasonality, otherwise s[SEASONALITY+2] will get overwritten
	
	s2[N+SEASONALITY2+1]=1;  //for integer seasonality the last value is not filled and Stan does not like it
	s2[SEASONALITY2+1] = s2[1];
	s2[SEASONALITY2+2] = s2[2];
	
	if (USE_SMOOTHED_ERROR)
	  smoothedInnovSize[1]=innovSizeInit;
	else
	  smoothedInnovSize[1]=1;
	
	powTrend= (MAX_POW_TREND-MIN_POW_TREND)*powTrendBeta+MIN_POW_TREND;
	expVal[1] = y[1];
	
	if (LEVEL_CALC_METHOD>0) {
		movingSum=y[1]-r[1];
		for (t in 2:SEASONALITY2) 
			movingSum=movingSum+y[t]-r[t];
		newLevelP=movingSum/SEASONALITY2;
	
		for (t in 1:SEASONALITY2)
			l[t] = newLevelP;
	}

	
	for (t in 2:N) {
		if (LEVEL_CALC_METHOD>0 && t>SEASONALITY2) 
			movingSum=movingSum+(y[t]-r[t])-(y[t-SEASONALITY2]-r[t-SEASONALITY2]);
		
		if (USE_GENERALIZED_SEASONALITY) {
			if (LEVEL_CALC_METHOD==0)
				newLevelP=y[t] - s[t]*l[t-1]^powSeason - s2[t]*l[t-1]^powSeason2 -r[t];
			else 
				newLevelP=movingSum/SEASONALITY2;
			l[t]  = levSm*newLevelP + (1-levSm)*l[t-1];	
			
			seasonalityP=sSm*(y[t] - l[t] - s2[t]*l[t]^powSeason2 - r[t])/l[t]^powSeason + (1-sSm)*s[t]; 
    		seasonalityP2= s2Sm*(y[t] - l[t] - s[t]*l[t]^powSeason - r[t])/l[t]^powSeason2 + (1-s2Sm)*s2[t]; 
    		expVal[t]=l[t-1]+ coefTrend*l[t-1]^powTrend + s[t]*l[t-1]^powSeason + s2[t]*l[t-1]^powSeason2 + r[t];
		} else {	
			if (LEVEL_CALC_METHOD==0)
				newLevelP=(y[t]-r[t])/(s[t]*s2[t]);
			else 
				newLevelP=movingSum/SEASONALITY2;
			l[t]  = levSm*newLevelP + (1-levSm)*l[t-1];	
			
			seasonalityP = sSm*(y[t]-r[t])/(l[t]*s2[t])+(1-sSm)*s[t];
			seasonalityP2 = s2Sm*(y[t]-r[t])/(l[t]*s[t])+(1-s2Sm)*s2[t];
			expVal[t]=(l[t-1]+ coefTrend*l[t-1]^powTrend)*s[t]*s2[t] + r[t];
		}
		if (fractSeasonality>0) {
    		s[t+SEASONALITY+1]=seasonalityP;  //with fractSeasonality weight
    		s[t+SEASONALITY]=fractSeasonality*s[t+SEASONALITY]+(1-fractSeasonality)*seasonalityP;
    	} else
    		s[t+SEASONALITY]=seasonalityP;
    		
		if (fractSeasonality2>0) {
    		s2[t+SEASONALITY2+1]=seasonalityP2;  //with fractSeasonality weight
    		s2[t+SEASONALITY2]=fractSeasonality2*s2[t+SEASONALITY2]+(1-fractSeasonality2)*seasonalityP2;
    	} else
    		s2[t+SEASONALITY2]=seasonalityP2;
    	
		if (USE_SMOOTHED_ERROR)
			smoothedInnovSize[t]=innovSm*fabs(y[t]-expVal[t])+(1-innovSm)*smoothedInnovSize[t-1];
		else	
			smoothedInnovSize[t]=1;
		//print(innovSm, " ",t," ",smoothedInnovSize[t]);
	}
}
model {
	sigma ~ cauchy(0,CAUCHY_SD) T[0,];
	offsetSigma ~ cauchy(MIN_SIGMA,CAUCHY_SD) T[MIN_SIGMA,];	
	coefTrend ~ cauchy(0, CAUCHY_SD);
	powTrendBeta ~ beta(POW_TREND_ALPHA, POW_TREND_BETA);
	
	if(USE_SMOOTHED_ERROR)
		innovSizeInit~ cauchy(y[1]/100,CAUCHY_SD) T[0,];
		
	if (USE_REGRESSION) {
		regCoef ~ cauchy(0, REG_CAUCHY_SD);
		regOffset ~ cauchy(0, reg0CauchySd);
	}		
	if (USE_GENERALIZED_SEASONALITY) {
		powSeason ~ beta(POW_SEASON_ALPHA, POW_SEASON_BETA); 
		powSeason2 ~ beta(POW_SEASON_ALPHA, POW_SEASON_BETA);
		for (t in 1:SEASONALITY)
			initSu[t] ~ cauchy (0, y[t]*0.1);	
		for (t in 1:SEASONALITY2)
			initSu2[t] ~ cauchy (0, y[t]*0.1);		
	} else {
		for (t in 1:SEASONALITY) 
    		initSu[t] ~ cauchy (1, 0.3) T[0.01,];
		for (t in 1:SEASONALITY2)
			initSu2[t] ~ cauchy (1, 0.3) T[0.01,];
	}
	
	for (t in 2:N) {
	  if (USE_SMOOTHED_ERROR==0)
	  	y[t] ~ student_t(nu, expVal[t], sigma*expVal[t]^powx+ offsetSigma);
	  else
	  	y[t] ~ student_t(nu, expVal[t], sigma*smoothedInnovSize[t-1] + offsetSigma);
	}
}
