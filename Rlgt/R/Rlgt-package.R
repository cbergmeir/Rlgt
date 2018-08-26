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
#' 
#' @docType package
#' @title Getting started with the Rlgt package
#' @description An implementation of innovative Bayesian ETS models named
#' LGT (for non-seasonal time series data) and SGT (for time series data).
#'  These models have been tested on M3-competition dataset in which they 
#'  outperform all of the models originally participating in the competition.
#' The following sections will briefly outline the mathematical constructions 
#' and the rationale for these models.
#' 
#' @section LGT(Local and Global Trend):
#' \subsection{Model Equations}{
#' 
#' In terms of mathematical notation, the model can be fully represented as follow:
#' 
#' \deqn{y_{t+1} \sim Student (v,y_{t+1}, \sigma _{t+1}) \quad (eq.1.1)}{y{t+1} ~ Student (v, y{t+1}, \sigma{t+1})     (eq.1.1)}   
#' \deqn{y_{t+1}=l_{t}+ \gamma l_{t}^{ \rho }+ \lambda b_{t}  \quad  (eq. 1.2)}{y{t+1} = l{t}+ \gamma *l{t}^\rho + \lambda*b{t}    (eq. 1.2)}
#' \deqn{l_{t}= \alpha y_{t}+ \left( 1- \alpha  \right)  \left(  l_{t-1} \right) \quad (eq. 1.3)}{l{t} = \alpha *y{t} + (1-\alpha)*l{t-1}    (eq. 1.3)} 
#' \deqn{b_{t+1}= \beta  \left( l_{t+1}-l_{t} \right) + \left( 1- \beta  \right) b_{t}  \quad  (eq. 1.4)}{b{t+1} = \beta*(l{t+1}-l{t}) + (1-\beta)*b{t}     (eq. 1.4)}
#' \deqn{\sigma _{t+1}= \sigma l_{t}^{ \tau}+ \varsigma   \quad  (eq. 1.5) }{\sigma{t+1} = \sigma*l_{t}^(\tau) + \varsigma      (eq. 1.5)}
#' }
#' 
#' \subsection{Notations}{
#' 
#' \describe{
#' \item{\eqn{y_{t}}{y{t}}}{value of the dependent variable of interest at time t}
#' \item{\eqn{y_{t+1}}{y{t+1}}}{predicted value of y at time t+1 given information up to time t}
#' \item{\eqn{\sigma_{t}}{\sigma{t}}}{variance of the distribution at time t}
#' \item{\eqn{l_{t}}{l{t}}}{level at time t}
#' \item{\eqn{b_{t}}{b{t}}}{local trend at time t}
#' }
#' }
#' 
#' 
#' \subsection{Parameters}{
#' 
#' \describe{
#' \item{\eqn{v}}{degrees of freedom of the t-distribution}
#' \item{\eqn{\gamma}}{coefficient of the global trend}
#' \item{\eqn{\rho}}{power coefficient of the global trend}
#' \item{\eqn{\lambda}}{damping coefficient of the local trend}
#' \item{\eqn{\alpha}}{smoothing parameter for the level term}
#' \item{\eqn{\beta}}{smoothing parameter for the local trend term}
#' \item{\eqn{\sigma}}{coefficient of the heteroscedastic standard deviation}
#' \item{\eqn{\tau}}{power coefficient of the heteroscedastic standard deviation}
#' \item{\eqn{\varsigma}}{base or minimum value of the standard deviation}
#' }
#' }
#' 
#' 
#' @section SGT (Seasonal, Global Trend):
#' 
#' SGT model was designed as a seasonal counterpart to the LGT model. 
#' Similar to LGT, this model is devised to allow for global trend term and heteroscedastic error.
#' 
#' \subsection{Model Equations}{
#' 
#' \deqn{ y_{t+1} \sim Student \left( v,y_{t+1}, \sigma _{t+1} \right)  \quad (eq. 2.1) }{y{t+1} ~ Student(v, y{t+1}, \sigma{t+1})     (eq. 2.1)} 
#' \deqn{ y_{t+1}= \left( l_{t}+ \gamma l_{t}^{ \rho } \right)  s_{t+1-m} \quad (eq. 2.2)}{y_{t+1} = (l{t}+ \gamma*l{t}^\rho) * s{t+1-m}     (eq. 2.2)} 
#' \deqn{ l_{t}= \alpha  \frac{y_{t}}{s_{t}}+ \left( 1- \alpha  \right)  \left( l_{t-1} \right) \quad (eq. 2.3)}{l_{t} = \alpha*y{t}/s{t} + (1-\alpha)*(l{t-1})    (eq. 2.3)}  
#' \deqn{ s_{t+m}= \zeta  \frac{y_{t}}{l_{t}}+ \left( 1- \zeta  \right) s_{t}  \quad (eq. 2.4)}{s{t+m} = \zeta*y_{t}/l_{t} + (1-\zeta)*s_{t}     (eq. 2.4)}
#' \deqn{ \sigma _{t+1}= \sigma y_{t+1}^{ \tau}+ \varsigma \quad (eq. 2.5)}{\sigma{t+1} = \sigma*y{t+1}^{\tau} + \varsigma     (eq. 2.5)}
#' }
#' 
#' \subsection{Additional Notations}{
#' 
#' \describe{
#' \item{\eqn{s_{t}}}{seasonality factor at time t}
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
#' \subsection{Rationale for Mathematical Equations}{
#' 
#' Most of the components of the model is similar to the non-seasonal model above. 
#' A couple of modifications have been applied to this model in comparison to LGT:
#' 
#' \enumerate{
#' \item Removal of local dampen trend
#' \item Addition of multiplicative seasonality term
#' }
#' 
#' The mathematical equations have also been adjusted accordingly:
#' 
#' \itemize{
#' \item \strong{Eq. 2.2.One-step Ahead Prediction Equation}
#' 
#' Comparing eq. 2.2 with eq. 1.2, it is apparent that the local trend term  
#' \eqn{\lambda b_{t}}  has been removed from the SGT model. 
#' Consequently, the equation governing the evolution of local trend over time (eq. 2.4) is 
#' no longer needed. Based on empirical evidence, the role of the local trend was 
#' found to be insignificant for seasonal data in M3-competition dataset.
#' 
#' Moreover, there is also an addition of multiplicative seasonality term modelled 
#' after Holt-Winters seasonality term. In this case, the data value is decomposed 
#' into level value and seasonal value (i.e.  { y_{t}=l_{t} s_{t} } ). Thus, it follows 
#' that the one-step prediction value of  { y_{t+1} }  (eq. 2.2) is obtained by 
#' multiplying the predicted level term  \eqn{l_{t}+ \gamma l_{t}^{ \rho }} with the 
#' seasonal term  \eqn{s_{t+1}}. 
#'
#' \item \strong{Eq. 2.4.Seasonal Factors Adjustment Equation}
#' 
#' The evolution of the seasonal component in eq. 2.4. is also based on standard 
#' ETS model, i.e. a weighted average of predicted seasonal level  
#' \eqn{\left( \frac{y_{t}}{l_{t}} \right)}  and previous seasonal 
#' value from the past observations \eqn{s_{t}} .
#' 
#' \item \strong{Eq. 2.3.Level Adjustment Equation}
#' 
#' The evolution of level term in eq. 2.3 is defined analogously to the level 
#' term in LGT model. The level term of time t  \eqn{\left( l_{t} \right)}  is 
#' calculated as a weighted average of predicted current level value at time t  
#' \eqn{ \frac{y_{t}}{s_{t}} }  and the previous estimate of the level at time t 
#' given information up to t-1, i.e.  \eqn{ l_{t-1}+ \gamma l_{t-1}^{ \rho }.}  
#' 
#'  \item \strong{Eq. 2.5.Heteroscedastic Error}
#'  
#' There is also a straightforward modification to the error term in SGT eq. 2.5 
#' due to the introduction of seasonal factors. Compared to eq. 1.5, the error term
#'  in SGT is allowed to vary in proportion to \eqn{  y  } instead of  \eqn{ l_{t} } . 
#'  The reason for that is because the error term is likely to become larger d
#'  uring the peak seasons (i.e. seasons with high seasonality coefficient). 
#'  Hence, it is intuitive to link the dispersion of the error to the predicted 
#'  value of the dependent variable. Note that in the case of LGT, the predicted 
#'  value of y is also a function of the previous level  \eqn{ l_{t-1} }  and 
#'  thus, the error term can be defined in terms of  \eqn{ l_{t-1} }  instead.
#' }
#' }
#' 
#' @section S2GT (Double Seasonal, Global Trend):
#' 
#' S2GT is designed as an extension to SGT in time-series data which exhibit two seasonality patterns. 
#' The additional second seasonality factor generalises the model to capture a number of seasonalities 
#' that exist in the data. A classic example of this type of data is the hourly electricity consumption data 
#' which are affected by the time of the day as well as the day in the week. 
#' The mathematical modelling is based on the classical Double Seasonal Exponential Smoothing 
#' model by Taylor (2003).
#' 
#' \subsection{Model Equations}{
#' 
#' \deqn{y_{t+1} \sim Student \left( v,y_{t+1}, \sigma _{t+1} \right) \quad (eq. 3.1)} 
#' \deqn{ y_{t+1}=\left( l_{t}+ \gamma l_{t}^{ \rho } \right) s_{t+1}w_{t+1} \quad (eq. 3.2)} 
#' \deqn{ l_{t}= \alpha  \frac{y_{t}}{s_{t}w_{t}}+ \left( 1- \alpha  \right) \left( l_{t-1} \right)  \quad (eq. 3.3)} 
#' \deqn{ s_{t+m}= \zeta  \frac{y_{t}}{l_{t}w_{t}}+ \left( 1- \zeta  \right) s_{t} \quad (eq. 3.4)} 
#' \deqn{ w_{t+d}= \delta  \frac{y_{t}}{l_{t}s_{t}}+ \left( 1- \delta  \right) w_{t} \quad (eq. 3.5)} 
#' \deqn{  \sigma _{t+1}= \sigma y_{t+1}^{ \tau}+ \varsigma  \quad (eq. 3.6)} 
#' 
#' }
#' \subsection{Additional Notations}{
#' \describe{
#' \item{\eqn{ w_{t} }}{second seasonality factor prevailing at time t}
#' \item{\eqn{ d }}{number of (second) seasons in a complete period (e.g. 12 for monthly, 4 for quarterly}
#' }
#' }
#' 
#' \subsection{Additional Parameters}{
#' \describe{
#' \item{\eqn{\delta}}{smoothing parameters for the second seasonality factors}
#' }
#' }
#' 
#' \subsection{Rationale for Mathematical Equations}{
#' There is a minor difference in this model compared to the SGT model:
#' 
#' \enumerate{
#' \item Addition of second seasonality factors
#' }
#' 
#' The mathematical equations have been adjusted as follow:
#' 
#' \itemize{
#' 
#' \item \strong{Eq. 3.2. One-step Ahead Prediction Equation}
#'
#' Comparing eq. 2.2 with eq. 2.2, there is an addition of the second seasonality factor in multiplicative form. 
#' Therefore, both seasonality factors will affect the proportionality between the expected 
#' value of the observed data and the level. 
#' 
#' \item \strong{Equation 3.3. Level Adjustment Equation}
#' 
#' The evolution of level term in eq. 3.3 is largely identical to eq.2.3. 
#' The only difference is that since the observed value is now seen as a product 
#' of the level terms with both seasonality factors, the predicted current level 
#' value at time t is now calculated as   \eqn{ \frac{y_{t}}{l_{t}s_{t}} }. 
#' The weighted average in the right-hand side of the equation is then adjusted accordingly.
#' 
#' \item \strong{Eq. 3.4. First Seasonal Factors Adjustment Equation}
#' 
#' The evolution of the seasonal component in eq. 3.4. is also similar to eq. 2.4. 
#' Again, the predicted seasonal level based on current observation has been changed to  
#' \eqn{ \frac{y_{t}}{l_{t}w_{t}} }  to undo the effect of the other seasonality 
#' factor on the observed data value. 
#' 
#' \item \strong{Eq. 3.5. Second Seasonal Factors Adjustment Equation}
#' 
#' The adjustment equation for the second seasonality term is identical to eq. 3.4. 
#' with the slight change in the parameters' notations and values. 
#' }
#' }
#' 
#' @section Prior Distributions for LGT and SGT models:
#' 
#' The default prior distributions of the parameters are given below:
#' 
#' \describe{
#' \item{\eqn{ \sigma,\gamma,\varsigma }}{Cauchy distribution with 0 location value and the scale parameter equals to 1/200 of the maximum value of y}
#' \item{\eqn{b}}{Normally distributed with a mean of 0 and standard deviation of 1/200 of the maximum value of y.}
#' \item{\eqn{\phi}}{Uniformly distributed between -1 and 1}
#' \item{\eqn{\alpha,\beta,\zeta, \delta}}{Uniform between 0 and 1}
#' \item{\eqn{ s_{t}, w_{t}}}{ (i.i.d) normal with a mean of 1 and standard deviation of 0.3 before being normalised}
#' \item{\eqn{\rho}}{Uniform between -0.5 and 1.0}
#' \item{\eqn{v}}{Uniformly distributed between 2 and 20}
#' \item{\eqn{\tau}}{Beta distribution with shape parameters  \eqn{\alpha =1}  and \eqn{\beta =1 }}
#' }
#' 
#' Note that some of these prior distributions can be adjusted by the users in the \code{rlgt} function
#' 
#' @references 
#' Stan Development Team (2017). RStan: the R interface to Stan. R package version 2.16.2. http://mc-stan.org
#' 
NULL

