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
#' \deqn{y_{t+1} \sim Student (v,y_{t+1}, \sigma _{t+1}) \quad (eq.1.1)}    
#' \deqn{  y_{t+1}=l_{t}+ \gamma l_{t}^{ \rho }+ \lambda b_{t}  \quad  (eq. 1.2)} 
#' \deqn{  l_{t}= \alpha y_{t}+ \left( 1- \alpha  \right)  \left(  l_{t-1} \right)  } 
#' \deqn{  b_{t+1}= \beta  \left( l_{t+1}-l_{t} \right) + \left( 1- \beta  \right) b_{t}  \quad  (eq. 1.4)}
#' \deqn{   \sigma _{t+1}= \sigma l_{t}^{ \tau}+ \varsigma   \quad  (eq. 1.5) }
#' and the notations are defined as follow:
#' 
#' \deqn{  y_{t} \quad: value of the dependent variable of interest at time t}
#' \deqn{  y_{t+1} \quad: predicted value of y at time t+1 given information up to time t}
#' \deqn{   \sigma _{t}\quad: variance of the distribution at time t}
#' \deqn{  l_{t}  \quad: level at time t}
#' \deqn{  b_{t}\quad: local trend at time t}
#' 
#' Parameters of the model which need to be estimated
#' \deqn{  v  \quad: degrees of freedom of the t-distribution}
#' \deqn{   \gamma   \quad : coefficient of the global trend}
#' \deqn{   \rho   \quad  : power coefficient of the global trend}
#' \deqn{   \lambda  \quad  : damping coefficient of the local trend}
#' \deqn{   \alpha   \quad  : smoothing parameter for the level}
#' \deqn{   \beta   \quad  : smoothing parameter for the local trend}
#' \deqn{   \sigma   \quad  : coefficient of the heteroscedastic standard deviation}
#' \deqn{   \tau \quad : power coefficient of the heteroscedastic standard deviation}
#' \deqn{   \varsigma \quad : base/ minimum value of the standard deviation}


NULL

