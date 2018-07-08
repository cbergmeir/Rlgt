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
#' \deqn{  l_{t}= \alpha y_{t}+ \left( 1- \alpha  \right)  \left(  l_{t-1} \right) \quad (eq. 1.3) } 
#' \deqn{  b_{t+1}= \beta  \left( l_{t+1}-l_{t} \right) + \left( 1- \beta  \right) b_{t}  \quad  (eq. 1.4)}
#' \deqn{   \sigma _{t+1}= \sigma l_{t}^{ \tau}+ \varsigma   \quad  (eq. 1.5) }
#' and the notations are defined as follow:
#' \describe{
#' \item{\eqn{y_{t}}}{value of the dependent variable of interest at time t}
#' \item{\eqn{y_{t+1}}}{predicted value of y at time t+1 given information up to time t}
#' \item{\eqn{\sigma _{t}}}{variance of the distribution at time t}
#' \item{\eqn{l_{t}}}{level at time t}
#' \item{\eqn{b_{t}}}{local trend at time t}
#' }
#' Parameters of the model which need to be estimated
#' \describe{
#' \item{\eqn{v}}{degrees of freedom of the t-distribution}
#' \item{\eqn{\gamma}}{coefficient of the global trend}
#' \item{\eqn{\rho}}{power coefficient of the global trend}
#' \item{\eqn{\lambda}}{damping coefficient of the local trend}
#' \item{\eqn{\alpha}}{smoothing parameter for the level}
#' \item{\eqn{\beta}}{smoothing parameter for the local trend}
#' \item{\eqn{\sigma}}{coefficient of the heteroscedastic standard deviation}
#' \item{\eqn{\tau}}{power coefficient of the heteroscedastic standard deviation}
#' \item{\eqn{\varsigma}}{base/ minimum value of the standard deviation}
#' }


#' The rationale of each individual equation of the model is discussed below.
#' \itemize{
#' \item \strong{eq. 1.1. Student’s t error distribution}
#' The data value follows Student’s t-distribution around the 
#' expected value of the data with a time-varying standard deviation. 
#' The Student-t distribution can be seen as a generalisation of the normal 
#' distribution to allow for a fat-tailed error distribution
#' 
#' \item \strong{eq. 1.5. Heteroscedastic Error }
#' In addition to accounting for possible fat-tailed error distribution, 
#' the error function also allows the variance of the error to change as the level changes.  
#' This is achieved by allowing the scale (deviation) parameter of the assumed 
#' Student’s t-distribution to vary in proportion to the current level of the time series. 
#' This will account for common situations, where the magnitude of the error will increase 
#' as the value of the data points increases.
#' Moreover, this relationship does not necessarily be linear as the parameter  
#' \eqn{\tau}  controls the growth of the variance of the error term. 
#' In practice, the values taken by the parameter is often limited between 0 and 1. 
#' A value of 0 corresponds to constant variance, i.e. homoscedasticity, whereas a 
#' value of 1 approximates the behaviour of the multiplicative error ETS model. 
#' The value between 0 and 1 will describe an error function which grows in relation 
#' to the increase in data value, but at a slower pace than linear growth.

#' \item \strong{eq. 1.2: One-step prediction forecast}
#' There are three distinct terms that constitute this Bayesian ETS model: a level term, 
#' and a couple of different trends. The first term eqn{l_{t}} is the level term, 
#' while the second term  \eqn{l_{t}^{ \rho }}  refers to the global trend which increases 
#' with the level of the dependent variables  \eqn{l_{t}}  at a constant rate  
#' \eqn{ \gamma  } . Similar to the heteroscedasticity in eq. 1.5, the change in 
#' the value of the dependant variable does not necessarily grow linearly 
#' with respect to the level of the variable. This relationship can be tuned in 
#' by the use of exponential parameter  \eqn {  \rho  } . The interpretation of this 
#' global trend is also analogous to the time-varying error term, i.e. 
#' the value of  \eqn{ 0< \rho <1 }  indicates a global trend which grows 
#' faster than the additive models but slower than the multiplicative model.
#' The last term of the right-hand side of eq.1.2. refers to the usual local dampen trend in ETS model. 
#' However, there is an addition of dampen parameter  \eqn{  \lambda  } , constrained such that  
#' \eqn{ -1<  \lambda  <1 } , to reduce the strength of the local trend model.

#' \item \strong{eq. 1.3: Level adjustment equation}
#' This equation is defined according to the classical linear trend ETS model. 
#' The level at time t ( \eqn{ l_{t} } ) is calculated as a weighted average of 
#' the current observation  \eqn{ y_{t} }  and the previous level at lag 1 ( \eqn{ l_{t-1} } ) 
#' with smoothing parameter  \eqn{  \alpha  }.

#' \item \strong{eq. 1.4: Local trend adjustment equation}
#' Similarly, the evolution of the local trend  \eqn{ b_{t}  } 
#' is identically defined to the linear trend method.\  The local trend 
#' at time t  \eqn{ b_{t} }  is obtained as a weighted average of the 
#' difference in level terms  \eqn{  \left( l_{t}-l_{t-1} \right)  }  
#' and the trend at time t-1 ( \eqn{ b \_  \left( t-1 \right) )  } 
#'  with the smoothing parameter  \eqn{  \beta  } .
#' }

#' \strong{Parameters’ Prior Distributions in L/SGT Models}
#' The default prior distributions of the parameters are given below:

#' \describe{
#' \item{\eqn{ \sigma, \gamma, \varsigma }} {Cauchy distribution with 0 location value and the scale parameter equals to 1/200 of the maximum value of y}
#' \item{\eqn { b}}{Normally distributed with a mean of 0 and standard deviation of 1/200 of the maximum value of y.}
#' \item{\eqn{ \phi}} {Uniformly distributed between -1 and 1}
#' \item{\eqn{\alpha, \beta, \zeta}}{Uniform between 0 and 1}
#' \item{\eqn{ s_{t}}}{ (i.i.d) normal with a mean of 1 and standard deviation of 0.3 before being normalised}
#' \item{\eqn{\rho}}{Uniform between -0.5 and 1.0}
#' \item{\eqn{v}}{Uniformly distributed between 2 and 20}
#' \item{\eqn{\tau}}{Beta distribution with shape parameters  \eqn{\alpha =1}  and \eqn{\beta =1 }}
#' }

NULL

