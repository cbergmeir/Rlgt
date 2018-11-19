#' @title Sets and initializes the control parameters
#' @description This function initializes and sets the control parameters, i.e. 
#' hyperparameter values which control the prior distribution of the \code{\link{rlgtfit}} model. 
#' The purpose of this function is mainly to provide a default value for each of the hyperparameters. 
#' The function also accepts a customised set of values of the parameters as provided in the input of this function. 
#' This function is used in conjunction with the \code{\link{rlgt}} function. 
#' @param	ADAPT_DELTA Target Metropolis acceptance rate. See Stan manual. Suggested range is between (0.85-0.97).
#' @param	MAX_TREE_DEPTH NUTS maximum tree depth. See Stan manual for more details. Suggested range is between (10-15), defaut is 12.
#' @param NUM_OF_CHAINS Number of MCMC chains. Suggested range is 3 to 4. Default is 4.
#' @param NUM_OF_CORES Number of cores used for calculations. It can be smaller than NUM_OF_CHAINS, 
#' but for best computational speed, it should be equal to NUM_OF_CHAINS. Default is 4. 
#' @param ADD_JITTER Whether to add a very small amount (sd=min(y)*0.0001) of jitter to the input series. 
#' It is sometimes useful in cases of series with some perfectly flat sections. Default is TRUE. 
#' @param CAUCHY_SD_DIV Cauchy distribution is used for some parameters with non-obvious range. The error size hyperparameter 
#' of this distribution is calculated by dividing the max value of the time series by this constant. 
#' Suggested range is between (100,300). Default 150.
#' @param NUM_OF_ITER Starting number of iterations for each chain. Suggested range is between (2000,10000). Default is 5000.
#' Generally, the longer the series, the smaller is the value to reach convergence. 
#' Some models e.g. those with "innov" error size method are more difficult to fit and require more iterations. 
#' @param MAX_NUM_OF_REPEATS Maximum number of the sampling procedure repeats if the fit is unsatisfactorily, i.e. avgRHat>MAX_RHAT_ALLOWED.
#' Each round will double the number of iterations which could potentially double the running time. 
#' Suggested range is between (2,4). Default is 2.  
#' @param MAX_RHAT_ALLOWED Maximum average value of Rhat's that suggests a good fit, i.e. the treshold 
#' below which the fit is considered as acceptable. Consult Stan's manual for more details on Rhat. 
#' Suggested range is between (1.005,1.02). Default is 1.006. 
#' @param NUM_OF_SEASON_INIT_CYCLES For seasonal models, number of seasonality periods used for establishing initial seasonality coefficients. Default is 3.
#' @param MIN_NU Minimum degrees of freedom of the Student's distribution that is used in most models. Suggested range(1.2, 5). Default 2.
#' @param MAX_NU Maximum degrees of freedom of the Student's distribution. Suggested range is between (15,30). Default 20. 
#' @param MIN_POW_TREND Minimum value of the global trend power coefficient. Suggested range is between (-1,0). Default -.5
#' @param MAX_POW_TREND Maximum value of the global trend power coefficient. It should be 1 to allow the model to approach exponential growth when needed.
#' Default is 1.
#' @param POW_TREND_ALPHA Alpha parameter of Beta prior distribution.  
#' To make the forecast more upward curved, so to nudge it towards larger values, make the parameter larger. Suggested range is between (1,6)
#' Default 1.
#' @param	POW_TREND_BETA Beta parameter of  Beta prior distribution for the global trend power coefficient. 1 by default, see also above.
#' @param	POW_SEASON_ALPHA Alpha parameter of Beta distribution that is the prior of the power coefficient in the formula of the generalized seasonality in gSGT model. 
#' 1 by default, increasing it (say, to 3 or 5) will push the seasonality towards multiplicative behavior.
#' @param	POW_SEASON_BETA Beta parameter of Beta distribution that is the prior of the power coefficient in the formula of the generalized seasonality in gSGT model. 
#' 1 by default.
#' @param MIN_SIGMA Minimum size of the fitted sigma, applied for numerical stability. Must be positive. 1e-10 by default.
#' @param MIN_VAL Minimum value that forecast can take. Must be positive. 1e-30 by default.
#' @param MAX_VAL Maximum value the forecast can take. 1e38 by default.

# @returnType list
#' @return list of control parameters
#' @export
rlgt.control <- function(
	ADAPT_DELTA=0.9, 
	MAX_TREE_DEPTH=12,
	NUM_OF_CHAINS=4,
	NUM_OF_CORES=4,
	ADD_JITTER=TRUE,
	CAUCHY_SD_DIV=150,
	NUM_OF_ITER=5000,
	MAX_NUM_OF_REPEATS=2,
	MAX_RHAT_ALLOWED=1.006,
	NUM_OF_SEASON_INIT_CYCLES=3,
	MIN_NU=2,
  MAX_NU=20,
  MIN_POW_TREND=-0.5,
  MAX_POW_TREND=1,
  POW_TREND_ALPHA=1,
  POW_TREND_BETA=1,
  POW_SEASON_ALPHA=1,  
  POW_SEASON_BETA=1,
	MIN_SIGMA=1e-10, 
	MIN_VAL=1e-30, 
	MAX_VAL=1e38  
) {
  list(
		NUM_OF_CHAINS=NUM_OF_CHAINS,
		NUM_OF_CORES=NUM_OF_CORES,
		ADD_JITTER=ADD_JITTER, 
		CAUCHY_SD_DIV=CAUCHY_SD_DIV,
    MAX_RHAT_ALLOWED=MAX_RHAT_ALLOWED,
		NUM_OF_SEASON_INIT_CYCLES=NUM_OF_SEASON_INIT_CYCLES,
    NUM_OF_ITER=NUM_OF_ITER,
    MAX_NUM_OF_REPEATS=MAX_NUM_OF_REPEATS,
    MIN_SIGMA=MIN_SIGMA,
    MIN_NU=MIN_NU,
    MAX_NU=MAX_NU,
    MIN_POW_TREND=MIN_POW_TREND,
    MAX_POW_TREND=MAX_POW_TREND,
    POW_TREND_ALPHA=POW_TREND_ALPHA,
    POW_TREND_BETA=POW_TREND_BETA,
    POW_SEASON_ALPHA=POW_SEASON_ALPHA,
    POW_SEASON_BETA=POW_SEASON_BETA,
    ADAPT_DELTA=ADAPT_DELTA,
    MAX_TREE_DEPTH=MAX_TREE_DEPTH,
		MIN_VAL=MIN_VAL,
		MAX_VAL=MAX_VAL
  )
}
