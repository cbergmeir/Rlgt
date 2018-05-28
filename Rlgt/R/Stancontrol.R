

#'
#' @title Sets and initializes the main parameters of the algorithm
#' @description ANother building block: This function initialize and set the parameters of the algorithm. This function is created mainly to provide a default value for each parameter.
#' It generates a list of parameters, to be used with the \code{\link{fit.lgt}} function. 
#' @param MAX_RHAT_ALLOWED Maximum average Rhat that suggests a good fit, see Stan's manual. Suggested range(1.005,1.02), see also MAX_NUM_OF_REPEATS description below.
#' @param NUM_OF_ITER Number of iterations for each chain. Suggested range(1000,5000). Generally, the longer the series, the smaller the vallue will do. 
#' See also MAX_NUM_OF_REPEATS description below.
#' @param MAX_NUM_OF_REPEATS Maximum number of the sampling procedure repeats if the fit is unsatisfactorily (avgRHat>MAX_RHAT_ALLOWED).
#' Each round doubles the number of iterations. Suggested range(2,4)
#' @param CAUCHY_SD_DIV For parameters with non-obvious range Cauchy distribution is used. The error size of this distribution 
#' is calculated by dividing max value of the time series by this constant. Suggested range(100,300)
#' @param MIN_SIGMA Minimum size of the fitted sigma, applied for numerical stability. Must bve positive. 
#' @param MIN_NU Minimum degrees of freedom of the Student's distribution, that is used in most models. Suggested range(1.2, 5)
#' @param MAX_NU Maximum degrees of freedom of the Student's distribution. Suggested range(15,30) 
#' @param MIN_POW_TREND Minimum value of power of trend coefficient. Suggested range(-1,0) 
#' @param MAX_POW_TREND Maximum value of power of trend coefficient. It should stay 1 to allow the model to approach exponential growth when needed.
#' @param POW_TREND_ALPHA Alpha parameter of Beta distribution that is the prior of the power coefficient in the formula of trend parameter.
#' To make the forecast more curved, make it larger. Suggested range(1,6)
#' @param	POW_TREND_BETA Beta parameter of Beta distribution that is the prior of the power of trend parameter. 1 by default, see also above.
#' @param	POW_SIGMA_ALPHA Alpha parameter of Beta distribution that is the prior of the power coefficient in the formula of the error size. 1 by default, see also below.
#' @param	POW_SIGMA_BETA Beta parameter of Beta distribution that is the prior of the power coefficient in the formula of the error size.
#' If the powSigma fitted is considered too often too high (i.e.> 0.7) you can attempt to tame it down by increasing POW_SIGMA_BETA.  Suggested range(1,4).
#' ADAPT_DELTA Target Metropolis acceptance rate. See Stan manual. Suggested range (0.8-0.97). 
#' MAX_TREE_DEPTH NUTS maximum tree depth. See Stan manual. Suggested range (10-12).
#' @param	SEASONALITY E.g. 12 for monthly seasonality. 1 for non-seasonal models
#' @param	SKEW Skew of error distribution used by manually-skewed models. 0 be default. 
#' @param MAX_TREE_DEPTH Description
#' @param ADAPT_DELTA Description
#' Setting it negative makes negative innovations having smaller impact on the fitting than the positive ones,
#' which would have the effect of making a model "more optimistic". Suggested range (-0.5, 0.5).
#' @returnType list
#' @return list of control parameters
#' @export
lgt.control <- function(
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
  POW_SIGMA_ALPHA=1,  
  POW_SIGMA_BETA=1, 
  ADAPT_DELTA=0.9, 
  MAX_TREE_DEPTH=11,
  SEASONALITY=1,
  SKEW=0 
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
    POW_SIGMA_ALPHA=POW_SIGMA_ALPHA,
    POW_SIGMA_BETA=POW_SIGMA_BETA,
    ADAPT_DELTA=ADAPT_DELTA,
    MAX_TREE_DEPTH=MAX_TREE_DEPTH,
    SEASONALITY=SEASONALITY,
    SKEW=SKEW
  )
}
