# Bayesian exponential smoothing -- R/C++ implementation (no Stan needed)
#
# @title Fit a Bayesian exponential smoothing model
# @param y.full A vector containing the time series to smooth
# @param n.samples Number of posterior samples to generate.
# @param burnin Number of burn-in samples.
# @param m The frequency of seasonality. For seasonal model, m (seasonality) > 1. A non-seasonal model will be fitted by default.
# @section Details:
# Draws a series of samples from the posterior distribution of a Bayesian exponential smoothing model.
# 
# @return An object containing the results of the sampling process, plus some additional information.
# 
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats binomial glm predict rgamma rmultinom rt
# 
# @examples 
# \dontrun{
## Build data and test
# library(Mcomp)
# M3.data <- subset(M3,"yearly")
# 
# train.data = list()
# future.data = list()
# for (i in 1:645)
# {
#   train.data[[i]] = as.numeric(M3.data[[i]]$x)
#   future.data[[i]] = as.numeric(M3.data[[i]]$xx)  
# }
# 
# ## Test -- change below to test more series
# w.series = 1:20
# # w.series = 1:645        # uncomment to test all series
# 
# # use 10,000 posterior samples; change n.samples to 20,000 to test that as well if you want
# s = system.time({rv=blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4)})
# 
# s                         # overall timing info
# s[[3]] / length(w.series) # per series time
# 
# mean(rv$sMAPE)            # performance in terms of mean sMAPE
# mean(rv$InCI)/6           # coverage of prediction intervals -- should be close to 95%
# }
# 
# @export
blgt <- function(y.full, burnin = 1e4, n.samples = 1e4, m = 1, homoscedastic = F)
{
  # nu proposal
  # nu.prop = c(0.1,0.2,0.4,0.6,0.8,1,1.15,1.35,1.6,1.95, 2.4, 3, 4, 5.6, 8.84, 18.63, 1e3)
  # nu.prop = c(0.47,0.53,0.6,0.68,0.77,0.875,1,1.15,1.35,1.6,1.95, 2.4, 3, 4, 5.6, 8.84, 18.63, 1e3)
  nu.prop = c(1.6,1.95, 2.4, 3, 4, 5.6, 8.84, 18.63, 1e3)
  
  
  # Process data
  max.y  = max(y.full)
  n.full = length(y.full)
  y      = y.full[2:length(y.full)]
  y.diff = y.full[2:length(y.full)] - y[1:(length(y.full)-1)] # TODO: Tell Daniel about probable error here
  n      = length(y)
  
  # add seasonal flag
  if (m > 1) {
    seasonal = T
    # print("Fitting a seasonal model...")
  } else {
    # print("Fitting a non-seasonal model now...")
    # print("If you want to fit a seasonal model, set Seasonality to be at least 2.")
    seasonal = F
  }
  
  # Initialise variables
  b1        = y.diff[1]
  
  log.s         = rep(0,m)
  s.ix          = ((0:(n)) %% m) + 1
  s.ix          = s.ix[2:(n+1)]
  lambda2.log.s = rep(1, m-1)
  nu.log.s      = rep(1, m-1)
  tau2.log.s    = 10
  xi.log.s      = 1  
  
  alpha     = 0.7
  zeta      = 0.7
  beta      = 0.7
  rho       = 0.5
  
  if (seasonal) {
    rv.smooth     = blgt.sexpsmooth(y.full, alpha, beta, zeta, y.full[1], b1, log.s)
  } else {
    rv.smooth = blgt.expsmooth(y.full, alpha, beta, y.full[1], 0)
  }

  l         = rv.smooth$l[1:n]
  b         = rv.smooth$b[1:n]
  # blgt.expsmooth will return s = 1 for non-seasonal version
  s         = rv.smooth$s[2:(n+1)]

  lambda2.w1 = 1
  lambda2.w2 = 1
  lambda2.b1 = 1
  
  phi        = 1
  tau        = 0
  chi2       = 0.5
  chi2.l2.scale = 1 #sd(y.full)
  chi2.lambda2 = 1
  sigma2     = 0
  xi2        = 1
  nu         = 5
  omega2     = matrix(1, n, 1)
  ar         = 0
  e1         = 0

  l1         = y[1]
  
  e          = matrix(0, n, 1)
  
  # Hyperparameters
  max.y.diff = max(abs(y.diff))

  # Initialise weights with least-squares solution
  X     = matrix(0, n, 2)
  X[,1] = l^rho
  X[,2] = b
  w     = solve(t(X) %*% X, t(X) %*% (y-l))
  mu    = s*(l + X %*% w)

  w.s   = 1
  
  # Sample storage
  rv = list()
  rv$n.samples = n.samples
  rv$sigma2    = matrix(0, n.samples, 1)
  rv$xi2       = matrix(0, n.samples, 1)
  rv$phi       = matrix(0, n.samples, 1)
  rv$chi2      = matrix(0, n.samples, 1)
  rv$chi2.lambda2 = matrix(0, n.samples, 1)
  rv$w         = matrix(0, n.samples, 2)
  rv$alpha     = matrix(0, n.samples, 1)
  rv$beta      = matrix(0, n.samples, 1)
  rv$zeta      = matrix(0, n.samples, 1)
  rv$rho       = matrix(0, n.samples, 1)
  rv$tau       = matrix(0, n.samples, 1)
  rv$nu        = matrix(0, n.samples, 1)
  rv$l1        = matrix(0, n.samples, 1)
  rv$b1        = matrix(0, n.samples, 1)
  rv$lt        = matrix(0, n.samples, 1)
  rv$bt        = matrix(0, n.samples, 1)
  rv$et        = matrix(0, n.samples, 1)
  rv$log.s     = matrix(0, n.samples, m)
  rv$y.on.l    = matrix(0, n.samples, m)
  rv$L         = matrix(0, n.samples, 1)
  rv$log.s1    = matrix(0, n.samples, m)
  rv$w.s       = matrix(0, n.samples, 1)
  rv$l2.log.s  = matrix(0, n.samples, m-1)
  rv$t2.log.s  = matrix(0, n.samples, 1)
  rv$s.ix      = s.ix
  rv$m         = m
  rv$L         = matrix(0, n.samples, 1)
  rv$y         = y.full
  
  mu.hat       = matrix(0, n, 1)
  l.hat        = matrix(0, n, 1)
  b.hat        = matrix(0, n, 1)
  
  if (seasonal) {
    theta      = c(mgrad.logit(alpha,c(0,1)), mgrad.logit(zeta,c(0,1)));  
  } else {
    theta      = c(mgrad.logit(alpha,c(0,1)), mgrad.logit(beta,c(0,1)));  
  }
  
  # Sampling setup
  #theta = ...
  mh.tune.theta = mgrad.initialise(2, bayes.exp.L.mh, bayes.exp.h.mh, 75, 1e-20, 1e20, n.samples, burnin, F)
  if (seasonal)
    mh.tune.log.s = mgrad.initialise(m-1, bayes.exp.L.s.mh, bayes.exp.h.s.mh, 75, 1e-20, 1e20, n.samples, burnin, F)
  
  #nu.prop  = seq(from = 4, to = 30, length.out = 20)
  #nu.prop = c(4, 5.6, 8.84, 18.63, 1e3)
  rho.prop = seq(from = -0.5, to = 1, length.out = 30)
  tau.prop = seq(from=0, to=1, length.out=30)
  #tau.prop = -1
  #tau.prop = exp(seq(from=-3,to=0,length.out=30))
  phi.prop = seq(from=0, to=1, length.out=30)
  
  iter     = 0
  k        = 0
  
  # Hyperparameters
  sigma2.B     = 0
  sigma2.lv    = 1
  chi2.scale   = max.y/100
  
  w1.scale     = (max.y/100)
  w2.scale     = 1
  w1.delta     = 1
  w2.delta     = 1

  b1.scale     = max.y/100
  b1.delta     = 1
  
  sample.tau = !homoscedastic
  sample.phi = !homoscedastic

  # Precompute random variables
  rng.gam = matrix(0, burnin+n.samples, 6)
  rng.gam[,1] = rgamma(burnin + n.samples, shape = (w1.delta+1)/2, rate = 1)
  rng.gam[,2] = rgamma(burnin + n.samples, shape = (w2.delta+1)/2, rate = 1)
  rng.gam[,3] = rgamma(burnin + n.samples, shape = (b1.delta+1)/2, rate = 1)
  rng.gam[,4] = rgamma(burnin + n.samples, shape = (n)/2, rate = 1)
  rng.gam[,5] = rgamma(burnin + n.samples, shape = 1, rate = 1)
  rng.gam[,6] = rgamma(burnin + n.samples, shape = (n+1)/2, rate = 1)
  
  # Main sampling loop
  while (k < n.samples)
  {
    ## Noise model
    e[1]   = y[1] - mu[1] + ar*e1
    e[2:n] = (y[2:n] - mu[2:n]) + ar*(y[1:n-1]-mu[1:n-1])
    l2.tau = l^(2*tau)

    # Sample chi2
    v2   = phi^2 + (1-phi)^2*l2.tau
    V    = sum(e^2/v2/omega2)/2 + sigma2.B/2
    # improper prior for chi2
    chi2 = V / rng.gam[iter+1,4]
    #chi2 = V / stats::rgamma(1, shape=n/2, scale=1)
    
    # # chi2 with half-cauchy prior
    # chi2 = (V + 1/chi2.lambda2) / rng.gam[iter+1,6]
    # scale = chi2.l2.scale^2 + 1/chi2
    # chi2.lambda2 = scale / stats::rexp(1)

    # If using half-Cauchy ...    
    #V = sum(e^2/v2/omega2)/2 + 1/sigma2.lv
    #chi2 = V / rng.gam[iter+1,4]
    
    #V = 1/chi2 + 1/chi2.scale^2
    #sigma2.lv = V / rng.gam[iter+1,5]
    
    # Sample the l.v.s
    #shape = (nu + 1)/2
    #scale = (e^2/(chi2*v2) + nu)/2
    #omega2 = as.matrix(1 / stats::rgamma(n, shape=shape, scale=1/scale), n, 1)
    omega2 = as.matrix( ((e^2/(chi2*v2) + nu)/2) / stats::rgamma(n, shape=(nu+1)/2, scale=1), n, 1)
    
    # Convert to xi2/sigma2
    xi2    = chi2*phi^2
    sigma2 = chi2*(1-phi)^2

    ## Sample nu
    logOmega2 = sum(log(omega2))
    OneOverOmega2 = sum(1/omega2)
    
    L.prop = -n*(nu.prop/2)*log(nu.prop/2) + (nu.prop/2+1)*logOmega2 + (nu.prop/2)*OneOverOmega2 + n*lgamma(nu.prop/2)
    #L.prop = L.prop - (2/2)*log( (nu.prop+1)/(nu.prop+3) ); 
    #L.prop = L.prop - (1/2)*log( psigamma(nu.prop/2,1) 
    #                - psigamma((nu.prop+1)/2, 1) - 2*(nu.prop+5)/nu.prop/(nu.prop+1)/(nu.prop+3) ); 
    nu = bayes.exp.grid.sample(nu.prop, L.prop)

    ## Linear combination model
    #l2.tau = l^(2*tau)
    v2 = omega2*(xi2 + l2.tau*sigma2)
    l.rho = l^rho
    
    # Sample w(2)
    e  = y - s*(l + w[2]*b)
    w[1] = bayes.exp.sample.w(e, s*l.rho, v2, lambda2.w1*w1.scale^2)

    # Sample w(2)
    e = y - s*(l + l.rho*w[1])
    w[2] = bayes.exp.sample.w(e, s*b, v2, lambda2.w2*w2.scale^2, range = c(-100,1))
    if (seasonal) {
      w[2] = 0
    }

    # Form mu
    mu   = w.s*s*(l + l.rho*w[1] + b*w[2])
    
    # Cauchy for w(1) and w(2)
    lambda2.w1 = (w[1]^2/2/w1.scale^2 + w1.delta/2) / rng.gam[iter+1,1]
    lambda2.w2 = (w[2]^2/2/w2.scale^2 + w2.delta/2) / rng.gam[iter+1,2]

    ## Sample b1
    if (seasonal) {
      b1 = 0
    } else {
      db.db1  = (1-beta)^(0:(n-1))
      dmu.db1 = w[2]*db.db1
      g       = sum(-(y-mu)*dmu.db1/v2)
      H       = sum(dmu.db1^2/v2)
      Q2      = lambda2.b1*b1.scale^2
      
      mu.b1   = Q2*(H*b1 - g)/(H*Q2 + 1)
      s2.b1   = 1/(H + 1/Q2)
      
      # Sample b1 & lambda2.b1
      b1 = rnorm(1, mu.b1, sd=sqrt(s2.b1))
      lambda2.b1 = ((b1^2/b1.scale^2 + b1.delta)/2) / rng.gam[iter+1,3]
    }

    ## Sample alpha/beta
    approx.c  = c(1,1) * 12
    approx.mu = 0
    rv.sample = mgrad.Sample(theta, approx.c, y.full, mh.tune.theta, w, sigma2, xi2, omega2, nu, tau, b1, rho, ar, e1, log.s, seasonal)
    
    if (rv.sample$accept)    
    {

      alpha = mgrad.ilogit(rv.sample$theta[1]+approx.mu,c(0,1))$x
      if (seasonal) {
        zeta = mgrad.ilogit(rv.sample$theta[2]+approx.mu,c(0,1))$x
      } else {
        beta  = mgrad.ilogit(rv.sample$theta[2]+approx.mu,c(0,1))$x 
      }
    }
    mh.tune.theta = rv.sample$tune
    theta   = rv.sample$theta
    
    if (seasonal)
    {
      approx.c  = lambda2.log.s*tau2.log.s
      #approx.c = rep(2.5,m-1)
      #rv.sample = mgrad.Sample(log.s[1:(m-1)], approx.c, y.full, mh.tune.log.s, alpha, beta, w, sigma2, xi2, omega2, nu, tau, b1, rho, ar, e1, c(1,s.ix), w.s)
      rv.sample = mgrad.Sample(log.s[1:(m-1)], approx.c, y.full, mh.tune.log.s, alpha, beta, zeta, w, sigma2, xi2, omega2, nu, tau, b1, rho, ar, e1, seasonal)

      mh.tune.log.s = rv.sample$tune
      log.s         = c(rv.sample$theta, -sum(rv.sample$theta))
      #s             = exp(log.s)[s.ix]


      # Sample the horseshoe stuff
      rv.hs = blgt.horseshoe.lambda2(log.s[1:(m-1)], tau2.log.s, nu.log.s)
      
      lambda2.log.s = rv.hs$lambda2
      nu.log.s      = rv.hs$nu
      
      rv.hs = blgt.horseshoe.tau2(log.s[1:(m-1)], lambda2.log.s, xi.log.s)
      tau2.log.s    = rv.hs$tau2
      xi.log.s      = rv.hs$xi
    }
    
    # Update
    if (seasonal) {
      rv.smooth = blgt.sexpsmooth(y.full, alpha, 0, zeta, l1, 0, log.s)
    } else {
      rv.smooth = blgt.expsmooth(y.full, alpha, beta, l1, b1)
    }
    lt = rv.smooth$l[n.full]
    bt = rv.smooth$b[n.full]
    l  = rv.smooth$l[1:n]
    b  = rv.smooth$b[1:n]
    s  = rv.smooth$s[2:n.full]
    
    # Sample tau
    log.l = log(l)
    e = y - s*(l + l.rho*w[1] + b*w[2])
    if (sample.tau)
    {
      P = matrix(0, length(tau.prop), 2)
      P[,1] = phi
      P[,2] = tau.prop
      tau = rexpsmooth.grid.sample.tau.phi(P,chi2,e,log.l,omega2,nu)[2]
    }
    
    # Sample phi
    if (sample.phi)
    {
      P = matrix(0, length(phi.prop), 2)
      P[,1] = phi.prop
      P[,2] = tau
      phi = rexpsmooth.grid.sample.tau.phi(P,chi2,e,log.l,omega2,nu)[1]
    }
    
    ## Sample rho
    log.l   = log(l)
    rho = rexpsmooth.grid.sample.rho(rho.prop, y - s*(l+w[2]*b), xi2+sigma2*exp(log.l*2*tau), log.l, w[1], nu, s)$theta

    ## Store sample
    iter = iter+1
    if (iter > burnin)
    {
      k = k+1
      
      # Store
      rv$sigma2[k] = sigma2
      rv$xi2[k]    = xi2
      rv$w[k,]     = w
      rv$alpha[k]  = alpha
      rv$beta[k]   = beta
      rv$zeta[k]   = zeta
      rv$rho[k]    = rho
      rv$tau[k]    = tau
      rv$nu[k]     = nu
      rv$phi[k]    = phi
      rv$chi2[k]   = chi2
      rv$b1[k]     = b1
      rv$lt[k]     = lt
      rv$bt[k]     = bt
      rv$et[k]     = e[n]
      rv$log.s1[k,] = log.s
      rv$log.s[k,]  = rv.smooth$log_s[(n.full-m+1):n.full]
      if (seasonal)
        rv$y.on.l[k,] = log( y[(n-m+1):n] / c(l[(n-m+2):n],lt) )

      rv$l2.log.s[k,] = lambda2.log.s
      rv$t2.log.s[k]  = tau2.log.s
      
      mu.hat = mu.hat + mu
    }
  }
  
  rv$mu.hat = mu.hat / n.samples
  
  rv
}

########################################################################
# Sample a regression coefficient
bayes.exp.sample.w <- function(y, x, v2, lambda2, range = NULL)
{
  Q2 = lambda2
  XY = sum(y*x/v2)
  X2 = sum(x^2/v2)
  
  mu.w = (Q2*XY)/(Q2*X2+1);
  s2.w = 1/(X2 + 1/Q2);
  
  if (!is.null(range))
  {
    w = rtruncnorm(1, a=range[1], b=range[2], mean=mu.w, sd=sqrt(s2.w))
  }
  else
  {
    w = rnorm(1, mu.w, sqrt(s2.w)) 
  }
  
  w
}


########################################################################
# Sample horseshoe L.V.s
blgt.horseshoe.lambda2 <- function(b, tau2, nu)
{
  p = length(b)
  
  # Sample lambda2
  scale   = 1/nu + b^2 / 2 / tau2
  lambda2 = scale / stats::rexp(p)
  
  scale = 1 + 1/lambda2
  nu    = scale / stats::rexp(p)  
  
  return(list(lambda2=lambda2,nu=nu))
}

########################################################################
# Sample tau2
blgt.horseshoe.tau2 <- function(b, lambda2, xi)
{
  p = length(b)

  shape = (p+1)/2
  scale = 1/xi + sum(b^2 / lambda2) / 2
  tau2  = 1 / stats::rgamma(1, shape=shape, scale=1/scale)
  
  # Sample xi
  scale = 1 + 1/tau2
  xi    = scale / stats::rexp(1)
  
  return(list(tau2=tau2,xi=xi))
}


########################################################################
# Simple grid sample
bayes.exp.grid.sample <- function(theta.prop, L.prop)
{
  p.prop = exp(-L.prop + min(L.prop))
  theta.prop[which(rmultinom(1, 1, p.prop)==1)]
}

########################################################################
# Likelihood for alpha & beta
# @export
bayes.exp.L.mh <- function(theta, y.full, L.stats, aux.stats, w, sigma2, xi2, omega2, nu, tau, b1, rho, ar, e1, log.s, seasonal)
{
  n.full = length(y.full)
  n      = n.full-1
  
  # Extract parameters
  alpha = mgrad.ilogit(theta[1],c(0,1))$x
  if (seasonal) {
    zeta = mgrad.ilogit(theta[2],c(0,1))$x
  } else {
    beta  = mgrad.ilogit(theta[2],c(0,1))$x 
  }
  
  l1 = y.full[1]
  y  = y.full[2:n.full]
  if (!seasonal) {
    rv = blgt.expsmooth(y.full, alpha, beta, l1, b1)
  } else {
    rv = blgt.sexpsmooth(y.full, alpha, 0, zeta, l1, 0, log.s)

  }
  
  
  # Compute the likelihood of the new model
  l = rv$l[1:n]
  b = rv$b[1:n]
  s = rv$s[2:n.full]
  mu = s*(l + l^rho*w[1] + b*w[2])
  e = y - mu
  
  v2 = xi2 + sigma2*l^(2*tau)
  #E = [e(1) + e1*ar(1); e(2:end) + e(1:end-1)*ar(1)];

  E = matrix(0,length(l),1)
  E[1]   = e[1] + e1*ar
  E[2:n] = e[2:n] + e[1:n-1]*ar
  
  q = e^2/nu/v2 + 1
  L = (nu+1)/2*sum(log(q)) + (1/2)*sum(log(v2))
  
  ## Gradients
  # Gradient of neg-log-likelihood
  #de2.dbeta  = -2*rv$db_dbeta[1:n]*w[2]*e
  #de2.dalpha = 2*(-rv$db_dalpha[1:n]*w[2] - l^(rho-1)*rv$dl_dalpha[1:n]*rho*w[1] - rv$dl_dalpha[1:n])*e
  
  #A = (nu+1)/2
  #G = nu*xi2
  #dL.dbeta  = A*de2.dbeta / (G*(e^2/G+1))
  #dL.dalpha = A*de2.dalpha / (G*(e^2/G+1))
  
  #g = c(0,0)
  #g[1] = sum(dL.dalpha)
  #g[2] = sum(dL.dbeta)
  
  #g[1] = g[1] * exp(theta[1]) / (exp(theta[1]) + 1)^2
  #g[2] = g[2] * exp(theta[2]) / (exp(theta[2]) + 1)^2
  
  # Gradients of neg-log-prior
  #a = 1;
  #b = 1/2;
  
  #sum((b+a)*log(exp(theta[1:2])+1) - a*theta[1:2])
  #g = g + (a+b)*exp(theta[1:2]) / (exp(theta[1:2])+1) - a/theta[1:2]
  
  # Gradients of approx. prior
  #g = g - theta[1:2]/approx.c
  
  #L = sum(e^2/v2)/2 + (1/2)*sum(log(v2))
  
  # Gradient for alpha
  
  list(L = L, g = c(0,0), L.stats = rv, aux.stats = NULL)
}

########################################################################
# Prior for alpha & beta
# @export
bayes.exp.h.mh <- function(theta, c)
{
  # Beta prior on (alpha,beta)
  a = 1;
  b = 1/2;
  
  sum((b+a)*log(exp(theta[1:2])+1) - a*theta[1:2])
}
########################################################################
# Likelihood for seasonal adjustments
bayes.exp.L.s.mh <- function(theta, y.full, L.stats, aux.stats, alpha, beta, zeta, w, sigma2, xi2, omega2, nu, tau, b1, rho, ar, e1, seasonal)
{
  n.full = length(y.full)
  n      = n.full-1
  m      = length(theta)
 
  #t = theta[1:m]
  #t = t - mean(t)
  #t = exp(t)
  t = exp(c(theta[1:m], -sum(theta[1:m])))
  t = pmin(t,1e6)
  t = pmax(t,1/1e6)
  
  # Seasonal adjustments
  #s = exp(theta[1:m])[s.ix]
  l1 = y.full[1]
  y  = y.full[2:n.full]
  rv = blgt.sexpsmooth(y.full, alpha, 0, zeta, l1, 0, log(t))
  
  s = rv$s[2:n.full]
  
  # Compute the likelihood of the new model
  l = rv$l[1:n]
  b = rv$b[1:n]
  mu = s*(l + l^rho*w[1] + b*w[2])
  e = y - mu
  
  v2 = xi2 + sigma2*l^(2*tau)

  E = matrix(0,length(l),1)
  E[1]   = e[1] + e1*ar
  E[2:n] = e[2:n] + e[1:n-1]*ar
  
  q = e^2/nu/v2 + 1
  L = (nu+1)/2*sum(log(q)) + (1/2)*sum(log(v2))

  # Compute gradients
  dL.ds = matrix(0, n, m)
  mu = rv$l + w[1]*rv$l^rho
  mu.s = mu[1:(n.full-1)] * s
  l.rho.2 = rho*w[1]*rv$l^(rho-1)
  E2 = (y-mu.s)^2
  for (j in 1:m)
  {
    dmu = rv$dl_ds[,j] + l.rho.2*rv$dl_ds[,j]
    dmu.s = s*dmu[1:(n.full-1)] + mu[1:(n.full-1)]*rv$dlogs_ds[2:n.full,j]*s
    dv2.ds = 2*tau*sigma2*l^(2*tau-1)*rv$dl_ds[1:(n.full-1),j]
    dE2 = -2*(y - mu.s)*dmu.s
    # dL.ds[,j] = (nu+1)*dE2 / (2*nu*xi2*(E2/nu/xi2+1))
    dL.ds[,j] = -nu*dv2.ds/2/v2 + (nu+1)*(nu*dv2.ds + dE2)/2/(nu*v2+E2)
  }    

  g = colSums(dL.ds)
  # v2 can be small when phi = 0, tau = -1
  g[which(is.nan(g))] = 0

  #
  list(L = L, g = g, L.stats = rv, aux.stats = NULL)
}

########################################################################
# Prior for seasonal adjusters
bayes.exp.h.s.mh <- function(theta, c)
{
  # Normal prior on log-seasonal adjusters
  m = length(theta)
  #theta = theta-mean(theta)
  sum(theta[1:m]^2/c)/2
  
  #sum(-theta + log(exp(2*theta)+1))
  
  # Cauchy prior
  # sum(log(1+theta^2))
  
  #a = 1;
  #b = 1/2;
  
  #sum((b+a)*log(exp(theta[1:2])+1) - a*theta[1:2])
}


########################################################################
# Forecast
# @export
blgt.forecast <- function(rv, h, ns = 1e6)
{
  # Sample
  n = length(rv$y)
  n.samples = rv$n.samples
  m = rv$m
  I = sample(n.samples, ns, replace=T)
  
  yp = matrix(rv$y[n], ns, 1)
  yf = matrix(0, ns, h)
  
  # Initialise
  bS = rv$bt[I]
  prev.level = rv$lt[I]
  error = rv$et[I]

  seasonal = F
  if (m > 1) {
    seasonal = T
    log.s = matrix(0, ns, m + h)
    log.s[,1:m] = rv$log.s[I,]
  
    y.on.l = matrix(0, ns, m + h)
    y.on.l[,1:m] = rv$y.on.l[I,]
  }
  # Forecast
  for (i in 1:h)
  {
    # Update s
    if (seasonal) {
      log.s[,i+m] = rv$zeta[I]*y.on.l[,i] + (1-rv$zeta[I])*log.s[,i]
      # bound log.s within (-2,2)
      log.s[,i+m] = 4 * 1 / (1 + exp(-log.s[,i+m])) - 2
      # or take the original value
      # idx <- which(log.s[,i+m] < -0.5 | log.s[,i+m] > 0.5)
      # log.s[idx, i+m] <- log.s[idx, i]
    }

    l = prev.level
    if (seasonal) {
      exp.val = exp(log.s[,i+m])*(l + rv$w[I,1]*l^rv$rho[I] + rv$w[I,2]*bS)
    } else {
      exp.val = l + rv$w[I,1]*l^rv$rho[I] + rv$w[I,2]*bS
    }
    
    scale = sqrt(rv$xi2[I] + rv$sigma2[I]*l^(2*rv$tau[I]))
    
    e = rt(ns, rv$nu[I]) * scale

    yf[,i] = pmax(pmin(exp.val + e, 1e38), 1e-3)    
    
    if (seasonal) {
      cur.level = pmax(1e-3, rv$alpha[I]*yf[,i]/exp(log.s[,i+m]) + (1-rv$alpha[I])*prev.level)
    } else {
      cur.level = pmax(1e-3, rv$alpha[I]*yf[,i] + (1-rv$alpha[I])*prev.level)
    }
    idx <- which(yf[,i] <= 1e-3)
    cur.level[idx] <- prev.level[idx]
    
    if (seasonal) {
      y.on.l[,i+m] = log(yf[,i]/cur.level)
      y.on.l[idx,i+m] <- y.on.l[idx,i]
    }
    
    bS = rv$beta[I]*(cur.level-prev.level) + (1-rv$beta[I])*bS
    prev.level = cur.level
  }
  
  # Return forecasts
  list(yf.med  = apply(yf, 2, function(x){quantile(x,0.5)}),
       yf.CI05 = apply(yf, 2, function(x){quantile(x,0.025)}),
       yf.CI95 = apply(yf, 2, function(x){quantile(x,0.975)}),
       yf      = yf)
}

########################################################################
# Compute the sMAPE
# @export
blgt.sMAPE <- function(yp, yt)
{
  mean( abs(yp - yt)/(yp + yt) )*200
}

blgt.MASE <- function(yp, yt, train, m) {
  mae <- mean( abs(yp - yt) )
  
  n <- length(train)
  sNaive <- mean( abs( train[1:(n-m)] - train[(m+1):n] ))
  
  mae / sNaive
}

########################################################################
#' @title Rlgt LSGT Gibbs run in parallel
#' @description  Fit a list of series and produce forecast, then calculate the accuracy.
#' @param train A list of training series.
#' @param future A list of corresponding future values of the series.
#' @param n.samples The number of samples to sample from the posterior (the default is 2e4).
#' @param burnin The number of burn in samples (the default is 1e4).
#' @param parallel Whether run in parallel or not (Boolean value only, default \code{TRUE}). 
#' @param m The seasonality period, with default \code{m=1}, i.e., no seasonality specified.
#' @param homoscedastic Run with homoscedastic or heteroscedastic version of the Gibbs sampler version. By default it is set to \code{FALSE}, i.e., run a heteroscedastic model.
#' @return returns a forecast object compatible with the forecast package in R
#' @examples
#' \dontrun{demo(exampleScript)}
#' \dontrun{
#' ## Build data and test
#' library(Mcomp)
#' M3.data <- subset(M3,"yearly")
#' 
#' train.data = list()
#' future.data = list()
#' for (i in 1:645)
#' {
#'   train.data[[i]] = as.numeric(M3.data[[i]]$x)
#'   future.data[[i]] = as.numeric(M3.data[[i]]$xx)  
#' }
#' 
#' ## Test -- change below to test more series
#' w.series = 1:20
#' # w.series = 1:645        # uncomment to test all series
#' 
#' # use 10,000 posterior samples; change n.samples to 20,000 to test that as well if you want
#' s = system.time({
#'   rv=blgt.multi.forecast(train.data[w.series], future.data[w.series], n.samples=1e4)
#' })
#' 
#' s                         # overall timing info
#' s[[3]] / length(w.series) # per series time
#' 
#' mean(rv$sMAPE)            # performance in terms of mean sMAPE
#' mean(rv$InCI)/6           # coverage of prediction intervals -- should be close to 95%
#' }
#' @export
blgt.multi.forecast <- function(train, future, n.samples = 2e4, burnin = 1e4, parallel = T, m = 1, homoscedastic = F)
{
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop(
      "Package \"doParallel\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop(
      "Package \"foreach\" must be installed to use this function.",
      call. = FALSE
    )
  }
  #
  n.cores = Inf
  n.series = length(train)

  # If parallel
  if (parallel)
  {
    # Find out how many cores are available
    n.available.cores = parallel::detectCores()[1] - 1
    if (is.na(n.cores))
    {
      n.cores = n.available.cores
    }
    n.cores = min(n.available.cores, n.cores)

    # If more than one core available
    if (n.cores > 1)
    {
      #n.samples = ceiling(n.samples/n.cores)
      ix = list()
      q = ceiling(n.series/n.cores)

      j = 1
      for (i in 1:n.cores)
      {
        if (i < n.cores)
        {
          ix[[i]] = (j:(j+q-1))
        }
        else
        {
          if (j <= n.series)
            ix[[i]] = (j:n.series)
          else
            n.cores = n.cores - 1
        }
        j = j +q;
      }
    }
    else
    {
      warning('Only 1 additional core available -- no parallelisation will occur')
    }
  }
  else
  {
    n.cores = 1
  }

  # Main sampling loop
  if (n.cores > 1)
  {
    # Register a cluster
    cores   = parallel::detectCores()
    cluster = parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cluster)

    # Sample each chain
    rv.p = foreach::`%dopar%` (foreach::foreach(i=1:n.cores, .packages = "Rlgt"),
      {
        rv = vector("list", length(ix[[i]]))
        rv$sMAPE = rep(0, length(ix[[i]]))
        rv$InCI  = rep(0, length(ix[[i]]))
        for (j in 1:length(ix[[i]]))
        {
          k = ix[[i]][j]
          rv[[j]] = list()
          rv[[j]]$model    = blgt(train[[k]], burnin = burnin, n.samples = n.samples, m = m, homoscedastic = homoscedastic)
          rv[[j]]$forecast = blgt.forecast(rv[[j]]$model,length(future[[k]]),1e5)
          rv[[j]]$sMAPE    = blgt.sMAPE(rv[[j]]$forecast$yf.med, future[[k]])
          rv[[j]]$MASE     = blgt.MASE(rv[[j]]$forecast$yf.med, future[[k]], train[[k]], m)
          rv[[j]]$InCI     = sum( (future[[k]] > rv[[j]]$forecast$yf.CI05) & (future[[k]] < rv[[j]]$forecast$yf.CI95) )
          rv[[j]]$model    = NULL
        }
        rv
      })

    # Done -- close the cluster
    parallel::stopCluster(cluster)
  }

  ## Single core
  else
  {
    rv = list()
    rv$model = vector("list", n.series)
    rv$forecast = vector("list", n.series)
    rv$sMAPE = rep(0, n.series)
    rv$MASE = rep(0, n.series)
    rv$InCI  = rep(0, n.series)

    for (j in 1:n.series)
    {
      rv$model[[j]]    = blgt(train[[j]], burnin = burnin, n.samples = n.samples, m = m, homoscedastic = homoscedastic)
      rv$forecast[[j]] = blgt.forecast(rv$model[[j]],length(future[[j]]),1e5)
      rv$sMAPE[j]      = blgt.sMAPE(rv$forecast[[j]]$yf.med, future[[j]])
      rv$MASE[j]       = blgt.MASE(rv$forecast[[j]]$yf.med, future[[j]], train[[j]], m)
      rv$InCI[j]       = sum( (future[[j]] > rv$forecast[[j]]$yf.CI05) & (future[[j]] < rv$forecast[[j]]$yf.CI95) )
      rv$model[[j]]    = NULL
    }
  }

  # If multiple cores, repack into one list
  if (n.cores > 1)
  {
    rv = list()
    # rv$model = vector("list", n.series)
    rv$forecast = vector("list", n.series)
    rv$sMAPE = rep(0, n.series)
    rv$MASE = rep(0, n.series)
    rv$InCI  = rep(0, n.series)
    for (i in 1:n.cores)
    {
      for (j in 1:length(ix[[i]]))
      {
        k = ix[[i]][j]
        # rv$model[[k]] = rv.p[[i]][[j]]$model
        rv$forecast[[k]] = rv.p[[i]][[j]]$forecast
        rv$sMAPE[k] = rv.p[[i]][[j]]$sMAPE
        rv$MASE[k] = rv.p[[i]][[j]]$MASE
        rv$InCI[k] = rv.p[[i]][[j]]$InCI
      }
    }
  }

  rv
}


#function tune = mgrad_Initialise(p, L, h, window, delta_min, delta_max, nsamples, burnin, display)
# #' @export
mgrad.initialise <- function(p, L, h, window, delta.min, delta.max, n.samples, burnin, display)
{
  # Create a tuning object containing relevant data
  tune = list()
  
  # Likelihood and prior functions
  tune$p         = p
  
  tune$L         = L
  tune$L.stats   = NULL
  tune$aux.stats = NULL
  
  tune$h         = h
  
  # Step-size setup
  tune$M         = 0
  tune$window    = window
  
  tune$delta.max = delta.max;
  tune$delta.min = delta.min;
  
  # Start in phase 1
  tune$iter             = 0
  tune$n.samples        = n.samples
  tune$burnin           = burnin
  tune$W.phase          = 1
  tune$W.burnin         = 0
  tune$n.burnin.windows = floor(burnin/window)
  tune$m.window         = rep(0, tune$n.burnin.windows)
  tune$n.window         = rep(0, tune$n.burnin.windows)
  tune$delta.window     = rep(0, tune$n.burnin.windows)
  tune$display          = display
  
  tune$phase.cnt       = c(0,0)
  
  tune$M               = 0
  tune$D               = 0
  
  tune$b.tune          = NULL
  
  # Start at the maximum delta (phase 1)
  tune$delta           = delta.max
  
  tune
}


# #' @export
mgrad.tune <- function(tune)
{
  #MH_TUNE update tuning parameters for a Metropolis-Hastings sampler
  #   tune = mh_Tune(...) updates the adaptive step-size for a 
  #   Metropolis-Hastings sampler
  #
  #   The input arguments are:
  #       tune   - [1 x 1] a Metropolis-Hastings tuning structure
  #
  #   Return values:
  #       tune   - [1 x 1] updated Metropolis-Hastings tuning structure
  #
  #   (c) Copyright Enes Makalic and Daniel F. Schmidt, 2020
  
  DELTA.MAX   = exp(80)
  DELTA.MIN   = exp(-40)
  NUM.PHASE.1 = 15
  
  # Perform a tuning step, if necessary
  if ( (tune$iter %% tune$window) == 0)
  {
    # Store the measured acceptance probability
    tune$W.burnin                    = tune$W.burnin + 1;
    tune$delta.window[tune$W.burnin] = log(tune$delta);
    tune$m.window[tune$W.burnin]     = tune$M+1
    tune$n.window[tune$W.burnin]     = tune$D+2
    
    # If in phase 1, we are exploring the space uniformly
    if (tune$W.phase == 1)
    {
      tune$phase.cnt[1] = tune$phase.cnt[1]+1;
      #d = linspace(log(DELTA_MIN),log(DELTA_MAX),NUM_PHASE_1);
      d = seq(log(DELTA.MIN), log(DELTA.MAX), length.out = NUM.PHASE.1)
      
      tune$delta = exp(d[tune$phase.cnt[1]]);
      
      if (tune$delta < tune$delta.min)
      {
        tune$delta.min = tune$delta
      }
      if (tune$delta > tune$delta.max)
      {
        tune$delta.max = tune$delta
      }
      
      # If we have exhausted 
      if (tune$phase.cnt[1] == NUM.PHASE.1)
      {
        tune$W.phase = 2
      }
    }
    # Else in phase 2, we are probing randomly guided by model
    else {
      tune$phase.cnt[2] = tune$phase.cnt[2]+1
      
      # Fit a logistic regression to the response and generate new random probe point
      #tune$b.tune = glmfit(tune.delta_window(1:tune.W_burnin), [tune.m_window(1:tune.W_burnin), tune.n_window(1:tune.W_burnin)], 'binomial');
      YY          = matrix(c(tune$m.window[1:tune$W.burnin], tune$n.window[1:tune$W.burnin]), tune$W.burnin, 2)
      tune$b.tune = suppressWarnings(glm(y ~ x, data=data.frame(y=YY[,1]/YY[,2],x=tune$delta.window[1:tune$W.burnin]), family=binomial, weights=YY[,2]))
      
      probe.p     = runif(1)*0.7 + 0.15
      tune$delta  = exp( -(log(1/probe.p-1) + tune$b.tune$coefficients[1])/tune$b.tune$coefficients[2] )
      tune$delta  = min(tune$delta, DELTA.MAX);
      tune$delta  = max(tune$delta, DELTA.MIN);
      
      if (tune$delta > tune$delta.max)
      {
        tune$delta.max = tune$delta
      }
      else if (tune$delta < tune$delta.min)
      {
        tune$delta.min = tune$delta
      }
      
      if (tune$delta == tune$delta.max || tune$delta == tune$delta.min)
      {
        tune$delta = exp(runif(1)*(log(tune$delta.max) - log(tune$delta.min)) + log(tune$delta.min))
      }
    }
    
    #
    tune$M = 0
    tune$D = 0
  }
  
  # If we have reached last sample of burn-in, select a suitable delta
  if (tune$iter == tune$burnin)
  {
    # If the algorithm has not grown and shrunk delta, give an error
    if (tune$phase.cnt[2] < 100)
    {
      #error('Metropolis-Hastings sampler has not explored the step-size space sufficiently; please increase the number of burnin samples');
    }
    
    #tune$b.tune = glmfit(tune.delta_window(1:tune.W_burnin), [tune.m_window(1:tune.W_burnin), tune.n_window(1:tune.W_burnin)], 'binomial');
    YY          = matrix(c(tune$m.window[1:tune$W.burnin], tune$n.window[1:tune$W.burnin]), tune$W.burnin, 2)
    tune$b.tune = suppressWarnings(glm(y ~ x, data=data.frame(y=YY[,1]/YY[,2],x=tune$delta.window[1:tune$W.burnin]), family=binomial, weights=YY[,2]))
    
    # Select the final delta to use
    tune$delta  = exp( -(log(1/0.55-1) + tune$b.tune$coefficients[1])/tune$b.tune$coefficients[2] )
    #if (tune.delta == 0 || isinf(tune.delta))
    #    tune.delta = log(tune.delta_max)/2;
    #end
    
    if (tune$display)
    {
      df = data.frame(log.delta=tune$delta.window[1:tune$W.burnin], p=YY[,1]/YY[,2])
      df.2 = data.frame(y=predict(tune$b.tune,newdata=data.frame(x=seq(log(tune$delta.min)-5,log(tune$delta.max)+5,length.out=100)),type="response"), x=seq(log(tune$delta.min)-5,log(tune$delta.max)+5,length.out=100))
      #print(ggplot(df, aes(x=log.delta,y=p)) + geom_point() + geom_line(data=df.2,aes(x=x,y=y,color="red")) + labs(title=paste("mGrad Burnin Tuning: Final delta = ", sprintf("%.3g",tune$delta), sep=""), x="log(delta)", y="Estimated Probability of Acceptance") + theme(legend.position = "none"))
      
      # Produce a diagonostic plot
      #figure(1);
      #clf;
      #plot(tune.delta_window, tune.m_window./tune.n_window, '.');
      #hold on;
      #de = linspace(log(tune.delta_min)-5, log(tune.delta_max)+5);
      #plot(de, 1./(1+exp(-(de*tune.b_tune(2) + tune.b_tune(1)))), 'LineWidth', 1.5);
      #grid on;
      #xlabel('$\log(\delta)$','Interpreter','Latex');
      #ylabel('Estimated Acceptance Probabi lity','Interpreter','Latex');
      #xlim([log(tune.delta_min), log(tune.delta_max)]);
      #
      #title(sprintf('mGrad Burnin Tuning: Final $\\delta=%.3g$', tune.delta),'Interpreter','Latex');
    }
  }
  
  # Done
  tune    
}

# #' @export
mgrad.Sample <- function(theta, c, data, tune, ...)
{
  # Generate the proposal
  L = tune$L(theta, data, tune$L.stats, tune$aux.stats, ...)
  tune$L.stats = L$L.stats
  tune$aux.stats = L$aux.stats
  #[L, g, tune.Lstats, tune.auxstats] = tune.L(theta, data, tune.Lstats, tune.auxstats, varargin{:});
  delta.p = rep(1,tune$p) * tune$delta
  
  # Proposal
  prop.v  = c*delta.p*(4*c + delta.p)/(2*c + delta.p)^2
  prop.mu = c*(delta.p*-L$g + 2*theta)/(2*c + delta.p)
  
  # theta_new = normrnd( prop_mu, sqrt(prop_v) );
  #theta.new = randn(tune.p,1).*sqrt(prop_v) + prop_mu;
  theta.new = rnorm(tune$p) * sqrt(prop.v) + prop.mu
  
  #L.new = tune.L(theta.new, data, NA, tune$aux.stats, varargin{:});
  L.new = tune$L(theta.new, data, NULL, tune$aux.stats, ...)
  
  delta.new = rep(1,tune$p) * tune$delta
  prop.mu.new = c*(delta.new * -L.new$g + 2*theta.new)/(2*c + delta.new)
  
  # Accept/reject?
  # Log-priors
  L.h     = tune$h(theta, c)
  L.h.new = tune$h(theta.new, c)
  
  # Log-proposals
  L.prop = sum( (theta.new - prop.mu)^2/2/prop.v );
  L.prop.new = sum( (theta - prop.mu.new)^2/2/prop.v );  
  
  accept = F
  tune$D = tune$D + 1
  
  d = L.new$L + L.h.new + L.prop.new - L$L - L.h - L.prop
  if (is.nan(d))
  {
    # Check which element is 'nan' and report an error
    if (is.nan(L.new$L))
    {
      stop("Likelihood at new proposal is NaN")
    }
    if (is.nan(L.prop.new))
    {
      stop("L.prop.new is NaN -- likely that gradients are NaN")
    }
  }
  
  if (runif(1) < exp( - (L.new$L + L.h.new + L.prop.new - L$L - L.h - L.prop) ))
  {
    theta = theta.new
    
    tune$L.stats = L.new$L.stats
    tune$aux.stats = L.new$aux.stats
    
    #fprintf('Diff=%g; L_old = %g; L_new=%g; L_prop=%g; L_prop_new=%g\n', exp(-(L_new + L_h_new + L_prop_new - L - L_h - L_prop)), L, L_new, L_prop, L_prop_new);
    
    accept = T
    tune$M = tune$M + 1
  }
  
  # Are we tuning?
  tune$iter = tune$iter + 1
  if (tune$iter <= tune$burnin)
  {
    tune = mgrad.tune(tune)
  }
  
  # Return values
  rv = list()
  rv$theta  = theta
  rv$tune   = tune
  rv$accept = accept
  if (accept)
  {
    rv$L = L.new$L
  }
  else
  {
    rv$L = L$L
  }
  
  rv
}

# #' @export
mgrad.logit <- function(x, bnds)
{
  # Map to (0,1)
  x = (x-bnds[1]) / (bnds[2]-bnds[1])
  
  # Map to (-inf, +inf)
  log(x/(1-x))
}

# #' @export
mgrad.ilogit <- function(x, bnds)
{
  # Map to (0,1)
  y = exp(x)/(exp(x)+1)
  y = (y+1e-16)/(1+2e-16)
  
  # Map to [bnds(1), bnds(2)]
  rv = list()
  rv$x = y * (bnds[2]-bnds[1]) + bnds[1]
  rv$g = (bnds[2]-bnds[1]) * exp(x)/(exp(x)+1)^2
  
  rv
}

# @export
blgt.expsmooth <- function(y, alpha, beta, l1, b1)
{
  rcpp_expsmooth(y,alpha,beta,l1,b1)
}

blgt.sexpsmooth <- function(y, alpha, beta, zeta, l1, b1, log.s)
{
  rcpp_sexpsmooth(y,alpha,beta,zeta,l1,b1,log.s)
}

# #' @export
rexpsmooth.grid.sample.tau.phi <- function(theta.prop, chi2, e, log.l, omega2, nu)
{
  rcpp_GridSampleTauPhi(theta.prop, runif(1), chi2, e, log.l, omega2, nu)
}

# #' @export
rexpsmooth.grid.sample.phi <- function(phi.prop, chi2, tau, e, logl, nu)
{
  rcpp_GridSamplePhi(phi.prop, runif(1), chi2, tau, e, logl, nu)
}

# #' @export
rexpsmooth.grid.sample.rho <- function(rho.prop, ytilde, chi2, log.l, w1, nu, s, rho.scale = NULL)
{
  if (is.null(rho.scale))
  {
    rho.scale = rep(1, length(rho.prop))
  }
  rcpp_GridSampleRho(rho.prop, runif(1), ytilde, chi2, log.l, w1, nu, s, rho.scale)
}

# #' @export
rexpsmooth.grid.sample.rho.gaussian.mix <- function(rho.prop, ytilde, v2, log.l, w1)
{
  rcpp_GridSampleRhoGaussianMix(rho.prop, runif(1), ytilde, v2, log.l, w1)
}
