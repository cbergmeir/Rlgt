#' @name Rlgt-package
#' @aliases Rlgt
#' @docType package
#' @keywords forecasting exponential smoothing
# @exportPattern "^[[:alpha:]]+"
#' @useDynLib Rlgt, .registration=TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#' @importFrom stats frequency
#' @importFrom stats runif
#' 
#' @docType package
#' @title Getting started with the Rlgt package
#' @description An implementation of Bayesian ETS models named
#' LGT (for non-seasonal time series data) and SGT (for seasonal time series data).
#'  These models have been tested on the M3 competition dataset in which they 
#'  outperform all of the models originally participating in the competition.
#' 
#' @section Getting started:
#' 
#' The best way to get started with the package is to have a look at the vignettes and the various demos that ship with the package. 
#' There is a vignette with examples of how to use the various methods included in the package, and a vignette that discusses some of the
#' theoretical background.
#'
#' As to the demos, you can find their source code in the "demo" subfolder in the package sources (available on CRAN). 
#' There are some basic demos and other more advanced ones that run on subsets of the M3 dataset and run potentially for hours.
#'
#' The package contains models for seasonal and non-seasonal data, allows for external regressors, and different 
#' error distributions. In the following, we briefly also present some of the theoretical background of the methods.
#'   
#' @section LGT (Local and Global Trend):
#' The LGT model is constructed based on Holt’s linear trend method. 
#' The model is designed to allow for a more general term of error by allowing 
#' for heteroscedasticity and an addition of constant "global" trend in the model.
#' 
#' \subsection{Model Equations}{
#' 
#' In terms of mathematical notation, the model can be fully represented as follows:
#' 
#' \deqn{y_{t+1} \sim Student (\nu,y_{t+1}, \sigma _{t+1}) \quad (eq.1.1)}{y{t+1} ~ Student (\nu, yf{t+1}, \sigma{t+1})  ....(eq.1.1)}   
#' \deqn{\widehat{y}_{t+1}=l_{t}+ \gamma l_{t}^{ \rho }+ \lambda b_{t}  \quad  (eq. 1.2)}{yf{t+1} = l{t}+ \gamma*l{t}^\rho + \lambda*b{t}   ....(eq. 1.2)}
#' \deqn{l_{t+1}= \alpha y_{t+1}+ \left( 1- \alpha  \right)  \left(  l_{t} \right) \quad (eq. 1.3)}{l{t+1} = \alpha*y{t+1} + (1-\alpha)*l{t}    ....(eq. 1.3)} 
#' \deqn{b_{t+1}= \beta  \left( l_{t+1}-l_{t} \right) + \left( 1- \beta  \right) b_{t}  \quad  (eq. 1.4)}{b{t+1} = \beta*(l{t+1}-l{t}) + (1-\beta)*b{t}     ....(eq. 1.4)}
#' \deqn{ \widehat{\sigma}_{t+1}= \sigma l_{t}^{ \tau}+ \xi   \quad  (eq. 1.5) }{\sigmaf{t+1} = \sigma*l_{t}^(\tau) + \xi      ....(eq. 1.5)}
#' }
#' 
#' \subsection{Notations}{
#' 
#' \describe{
#' \item{\eqn{y_{t}}{y{t}}}{value of the dependent variable of interest at time t}
#' \item{\eqn{\widehat{y}_{t+1}}{yf{t+1}}}{forecasted value of y at time t+1 given information up to time t}
#' \item{\eqn{\widehat{\sigma}_{t+1}}{\sigmaf{t+1}}}{forecasted deviation at time t+1 given information up to time t}
#' \item{\eqn{l_{t}}{l{t}}}{level at time t}
#' \item{\eqn{b_{t}}{b{t}}}{local trend at time t}
#' }
#' }
#' 
#' 
#' \subsection{Parameters}{
#' 
#' \describe{
#' \item{\eqn{\nu}}{degrees of freedom of the t-distribution}
#' \item{\eqn{\gamma}}{coefficient of the global trend}
#' \item{\eqn{\rho}}{power coefficient of the global trend}
#' \item{\eqn{\lambda}}{damping coefficient of the local trend}
#' \item{\eqn{\alpha}}{smoothing parameter for the level term}
#' \item{\eqn{\beta}}{smoothing parameter for the local trend term}
#' \item{\eqn{\sigma}}{coefficient of the heteroscedastic standard deviation}
#' \item{\eqn{\tau}}{power coefficient of the heteroscedastic standard deviation}
#' \item{\eqn{\xi}}{base or minimum value of the standard deviation}
#' }
#' }
#' 
#' 
#' @section SGT (Seasonal, Global Trend):
#' 
#' The SGT model was designed as a seasonal counterpart to the LGT model. 
#' Similar to LGT, this model is devised to allow for a global trend term and heteroscedastic error.
#' 
#' \subsection{Model Equations}{
#' 
#' \deqn{ y_{t+1} \sim Student \left( \nu,\widehat{y}_{t+1}, \sigma _{t+1} \right)  \quad (eq. 2.1) }{y{t+1} ~ Student(\nu, yf{t+1}, \sigma{t+1})   ....(eq. 2.1)} 
#' \deqn{ \widehat{y}_{t+1}= \left( l_{t}+ \gamma l_{t}^{ \rho } \right)  s_{t+1} \quad (eq. 2.2)}{yf{t+1} = (l{t}+ \gamma*l{t}^\rho) * s{t+1}    ....(eq. 2.2)} 
#' \deqn{ l_{t+1}= \alpha  \frac{y_{t+1}}{s_{t+1}}+ \left( 1- \alpha  \right)  \left( l_{t} \right) \quad (eq. 2.3)}{l_{t+1} = \alpha*y{t+1}/s{t+1} + (1-\alpha)*l{t}   ....(eq. 2.3)}  
#' \deqn{ s_{t+m+1}= \zeta  \frac{y_{t+1}}{l_{t+1}}+ \left( 1- \zeta  \right) s_{t+1}  \quad (eq. 2.4)}{s{t+m+1} = \zeta*y{t+1}/l{t+1} + (1-\zeta)*s{t+1}   ....(eq. 2.4)}
#' \deqn{ \widehat{\sigma}_{t+1}= \sigma \widehat{y}_{t+1}^{ \tau}+ \xi \quad (eq. 2.5)}{\sigmaf{t+1} = \sigma*yf{t+1}^\tau + \xi    ....(eq. 2.5)}
#' }
#' 
#' \subsection{Additional Notations}{
#' 
#' \describe{
#' \item{\eqn{s_{t}}{s{t}}}{seasonality factor at time t}
#' \item{\eqn{ m }}{number of seasons in the data (e.g. 12 for monthly, 4 for quarterly)}
#' }
#' }
#' 
#' \subsection{Additional Parameters}{
#' 
#' \describe{
#'  \item{\eqn{\zeta}}{smoothing parameter for the seasonality terms}
#' }
#' 
#' }
#' 
#' @section S2GT (Double Seasonal, Global Trend):
#' 
#' S2GT is designed as an extension to SGT for time series data which exhibit two seasonality patterns. 
#' 
#' \subsection{Model Equations}{
#' 
#' \deqn{ y_{t+1} \sim Student \left( \nu,\widehat{y}_{t+1}, \sigma _{t+1} \right) \quad (eq. 3.1)}{y{t+1} ~ Student (\nu, yf{t+1}, \sigma{t+1})    (eq. 3.1)} 
#' \deqn{ \widehat{y}_{t+1}=\left( l_{t}+ \gamma l_{t}^{ \rho } \right) s_{t+1}w_{t+1} \quad (eq. 3.2)}{yf{t+1} = (l{t}+\gamma*l{t}^\rho) * s{t+1} * w{t+1}    (eq. 3.2)} 
#' \deqn{ l_{t}= \alpha  \frac{y_{t}}{s_{t}w_{t}}+ \left( 1- \alpha  \right) \left( l_{t-1} \right)  \quad (eq. 3.3)}{l{t}= \alpha*y{t}/(s{t}*w{t}) + (1- \alpha)*(l{t-1})    (eq. 3.3)} 
#' \deqn{ s_{t+m}= \zeta  \frac{y_{t}}{l_{t}w_{t}}+ \left( 1- \zeta  \right) s_{t} \quad (eq. 3.4)}{s{t+m} = \zeta*y{t}/(l{t}*w{t}) + (1-\zeta)*s{t}    (eq. 3.4)} 
#' \deqn{ w_{t+d}= \delta  \frac{y_{t}}{l_{t}s_{t}}+ \left( 1- \delta  \right) w_{t} \quad (eq. 3.5)}{w{t+d} = \delta*y{t}/(l{t}*s{t}) + (1-\delta)*w{t}    (eq. 3.5)} 
#' \deqn{ \widehat{\sigma} _{t+1}= \sigma y_{t+1}^{ \tau}+ \xi  \quad (eq. 3.6)}{\sigmaf{t+1} = \sigma*y{t+1}^\tau + \xi    (eq. 3.6)} 
#' }
#' 
#' \subsection{Additional Notations}{
#' \describe{
#' \item{\eqn{w_{t}}{w{t}}}{second seasonality factor prevailing at time t}
#' \item{\eqn{d}}{number of (second) seasons in a complete period (e.g. 12 for monthly, 4 for quarterly)}
#' }
#' }
#' 
#' \subsection{Additional Parameters}{
#' \describe{
#' \item{\eqn{\delta}}{smoothing parameters for the second seasonality factors}
#' }
#' }
#' @section Prior Distributions for LGT and SGT models:
#' 
#' The default prior distributions of the parameters are given below:
#' 
#' \describe{
#' \item{\eqn{ \sigma,\gamma,\xi }{\sigma, \gamma,\xi }}{Cauchy distribution with 0 location value and the scale parameter equals to 1/200 of the maximum value of y}
#' \item{\eqn{b}}{Normally distributed with a mean of 0 and standard deviation of 1/200 of the maximum value of y.}
#' \item{\eqn{\phi}}{Uniformly distributed between -1 and 1}
#' \item{\eqn{\alpha, \beta, \zeta, \delta}}{Uniform between 0 and 1}
#' \item{\eqn{ s_{t}, w_{t}}{ s{t}, w{t}} for \eqn{t<= m,d}{ (i.i.d) normal with a mean of 1 and standard deviation of 0.3 before being normalised}
#' \item{\eqn{\rho}}{Uniform between -0.5 and 1.0}
#' \item{\eqn{\nu}}{Uniformly distributed between 2 and 20}
#' \item{\eqn{\tau}}{Beta distribution with shape parameters  \eqn{\alpha =1}  and \eqn{\beta =1 }}
#' }
#' 
#' Note that some of these prior distributions can be adjusted by the users in the \code{rlgt} function
#' 
# @references 
# Stan Development Team (2017). RStan: the R interface to Stan. R package version 2.16.2. http://mc-stan.org
# 
NULL

