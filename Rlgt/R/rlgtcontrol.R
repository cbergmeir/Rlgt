#' @title Sets and initializes the control parameters
#' @description This function initializes and sets the control parameters, i.e. 
#' hyperparameter values which control the prior distribution of the Rlgt model. 
#' The purpose of this function is mainly to provide a default value for each of the hyperparameter. 
#' The function also accepts a customised set of values of the parameters as provided in the input of this function. 
#' This function is used in conjunction of used with the \code{\link{rlgt}} function. 
#' @param MAX_RHAT_ALLOWED Maximum average value of Rhat that suggests a good fit, i.e. the treshold 
#' below which the fit is considered as acceptable. Consult Stan's manual for more details. 
#' Suggested range is between (1.005,1.02), also see MAX_NUM_OF_REPEATS description below.
#' @param NUM_OF_ITER Number of iterations for each chain. Suggested range is between (1000,5000). 
#' Generally, the longer the series, the smaller is the value to reach convergence.
#' See also MAX_NUM_OF_REPEATS description below.
#' @param MAX_NUM_OF_REPEATS Maximum number of the sampling procedure repeats if the fit is unsatisfactorily, i.e. avgRHat>MAX_RHAT_ALLOWED.
#' Each round will double the number of iterationsm which could potentially double the running time. Suggested range is between (2,4)
#' @param CAUCHY_SD_DIV Cauchy distribution is used for parameters with non-obvious range. The error size hyperparameter 
#' of this distribution is calculated by dividing max value of the time series by this constant. Suggested range is between (100,300)
#' @param MIN_SIGMA Minimum size of the fitted sigma, applied for numerical stability. Must be positive. 
#' @param MIN_NU Minimum degrees of freedom of the Student's distribution, that is used in most models. Suggested range(1.2, 5)
#' @param MAX_NU Maximum degrees of freedom of the Student's distribution. Suggested range is between (15,30) 
#' @param MIN_POW_TREND Minimum value of the global trend power coefficient. Suggested range is between (-1,0) 
#' @param MAX_POW_TREND Maximum value of the global trend power coefficient. It should stay at 1 to allow the model to approach exponential growth when needed.
#' @param POW_TREND_ALPHA Alpha hyperparameter of Beta prior distribution.  
#' To make the forecast more curved, make it larger. Suggested range is between (1,6)
#' @param	POW_TREND_BETA Beta hyperparameter of  Beta prior distribution for the global trend power coefficient. 1 by default, see also above.
#' @param	POW_SEASON_ALPHA Alpha parameter of Beta distribution that is the prior of the power coefficient in the formula of the generalized seasonality in gSGT model. 
#' 1 by default, increasing it will push the seasonality towards multiplicative behavior. Seel also below.
#' @param	POW_SEASON_BETA Beta parameter of Beta distribution that is the prior of the power coefficient in the formula of the generalized seasonality in gSGT model. 1 by default.
#' To push the seasonality to resemble more multiplicative seasonality, encourage larger power coefficient by increasing POW_SEASON_ALPHA, say, to 3 or 5.
#' @param	ADAPT_DELTA Target Metropolis acceptance rate. See Stan manual. Suggested range is between (0.85-0.97). 
#' @param	MAX_TREE_DEPTH NUTS maximum tree depth. See Stan manual for more details. Suggested range is between (10-12).
#' @param	SEASONALITY Number of seasons pertaining to the main seasonality pattern in the data, e.g. 24 for hourly seasonality. Used only by seasonal models.
#' @param	SEASONALITY2 Number of seasons pertaining to the secondary seasonality pattern in the data, e.g. 168 for hourly seasonality. Used only by dual-seasonality models.
# @returnType list
#' @return list of control parameters
#' @export
rlgt.control <- function(
  MAX_RHAT_ALLOWED=1.005, 
  NUM_OF_ITER=2500,
  MAX_NUM_OF_REPEATS=3,
  CAUCHY_SD_DIV=200,
  MIN_SIGMA=0.001, 
  MIN_NU=2,
  MAX_NU=20,
  MIN_POW_TREND=-0.5,
  MAX_POW_TREND=1,
  POW_TREND_ALPHA=1,
  POW_TREND_BETA=1,
  POW_SEASON_ALPHA=1,  
  POW_SEASON_BETA=1, 
  ADAPT_DELTA=0.9, 
  MAX_TREE_DEPTH=11,
  SEASONALITY=1,
	SEASONALITY2=1
) {
  
  list(CAUCHY_SD_DIV=CAUCHY_SD_DIV,
    MAX_RHAT_ALLOWED=MAX_RHAT_ALLOWED,
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
    SEASONALITY=SEASONALITY,
		SEASONALITY2=SEASONALITY2
  )
}
