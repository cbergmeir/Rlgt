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
#' 
#' 
#' 
#' 
NULL

