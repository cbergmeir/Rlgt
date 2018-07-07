#' @name Rlgt-package
#' @aliases Rlgt
#' @docType package
#' @keywords forecasting, exponential smoothing
# @exportPattern "^[[:alpha:]]+"
#' @useDynLib Rlgt, .registration=TRUE
#' @import methods
#' @import Rcpp
#' @examples
#' x <- 1
#' @docType package
#' @title Getting started with the Rlgt package
#' 
#' @description An implementation of innovative Bayesian ETS models named
#' LGT (for non-seasonal time series data) and SGT (for time series data).
#'  These models have been tested on M3-competition dataset in which they 
#'  outperform all of the models originally participating in the competition.
#' The following sections will briefly outline the mathematical constructions 
#' and the rationale for these models.
#' 
#' @section LGT(Local and Global Trend):
#' In terms of mathematical notation, the model can be fully represented as follow 
#' \eqn{y_{t+1} \sim Student (v,y_{t+1}, \sigma _{t+1})}  (eq.1.1)  
#' \eqn{  y_{t+1}=l_{t}+ \gamma l_{t}^{ \rho }+ \lambda b_{t}   }  (eq. 1.2)
#' \eqn{  l_{t}= \alpha y_{t}+ \left( 1- \alpha  \right)  \left(  l_{t-1} \right)    }  (eq. 1.3)
#' \eqn{  b_{t+1}= \beta  \left( l_{t+1}-l_{t} \right) + \left( 1- \beta  \right) b_{t}  }  (eq. 1.4)
#' \eqn{   \sigma _{t+1}= \sigma l_{t}^{ \tau}+ \varsigma    } (eq. 1.5)
#' and the notations are defined as follow:
#' 
#' \eqn{  y_{t}   } : value of the dependent variable of interest at time t
#' \eqn{  y_{t+1}   } : predicted value of y at time t+1 given information up to time t
#' \eqn{   \sigma _{t}   } : variance of the distribution at time t
#' \eqn{  l_{t}   } : level at time t
#' \eqn{  b_{t}   } : local trend at time t
#' 
#' Parameters of the model which need to be estimated
#' \eqn{  v   } : degrees of freedom of the t-distribution
#' \eqn{   \gamma    } : coefficient of the global trend
#' \eqn{   \rho    } : power coefficient of the global trend
#' \eqn{   \lambda    }: damping coefficient of the local trend
#' \eqn{   \alpha    } : smoothing parameter for the level
#' \eqn{   \beta    } : smoothing parameter for the local trend
#' \eqn{   \sigma    } : coefficient of the heteroscedastic standard deviation
#' \eqn{   \tau} : power coefficient of the heteroscedastic standard deviation
#' \eqn{   \varsigma} : base/ minimum value of the standard deviation

#' 
#' 
#' 
NULL

