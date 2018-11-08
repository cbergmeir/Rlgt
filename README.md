# README #

This package is based on code from Slawek Smyl to implement LGT, a local- and global trend exponential smoothing forecasting method using Rstan for model fitting.

### Installation ###

Rlgt: R CMD build runs the "cleanup" script that runs roxygen and creates the source files. Then, R CMD INSTALL can be used to install the .tar.gz package.

RlgtLik: Implementation of LGT into code extracted from the "forecast" package.

### TODOs ###

1. We need a paragraph to say what is so good about the package/model? Main selling points? What can we do that nobody else can?

An implementation of various Bayesian Exponential Smoothing (ETS) models which
    have been found to outperform all of the original models in the M3 competition.
    These models include LGT (Local-Global Trend), SGT (Seasonal Global Trend), and
    their variations. The Bayesian model fitting is based on the RStan package.

Local-global trend, multiple seasonalities, external variables, high accuracy, prediction intervals, non-normal errors (fat tailed), heteroscedasticity, very flexible...

trend between multiplicative and additive (between linear and non-linear)
seasonality is also generalized between multiplicative and additive.

2. In the vignette, we need one/some of the examples from the demo and go through them step by step (?) Which examples are suitable?

LGT_REG Runs LGT forecast with and without a regression component.
lynx demo with the lynx dataset.
SGT_REG Runs SGT forecast with and without a regression component.


LGT_M3 Runs through a subset of M3 yearly data, showing two methods of passing the non-seasonal series data.
LGT&SGT_M3 Using M3 series, it demos several possibilities or sub-versions of seasonal and non-seasonal models.
LGT&SGT    Implementation of LGT and SGT models on the R buit-in datasets.
S2GT_M4 It demos several ways of passing data to dual seasonality models, and several dual seasonality models, using M4 hourly data set.
S2GT_M4Hourly_parallel Parallel implementation of S2GT (and occasionally SGT) models on hourly time-series of M4-competition dataset.
SGT_M3 It uses quarterly subset of M3 data to demo ways of passing seasonal data.
SGT_M3parallel Parallel implementation of SGT model on monthly time-series data from M3-competition dataset.
SGT_M4Weekly_parallel Parallel implementation of SGT model with non-integer seasonality on weekly time-series data from M4-competition dataset.

3. Checking the demos:

running the lynx demo, I get:

SAMPLING FOR MODEL 'SGT' NOW (CHAIN 1).
Chain 1: Initialization between (-2, 2) failed after 100 attempts.
[1] "Error in sampler$call_sampler(args_list[[i]]) : Initialization failed."
error occurred during calling the sampler; sampling not done

--> if I reduce the "seasonality" parameter it works

S2GT_M4 demo gives me the same error, doesn't run.

--> Seems to be a problem with newest Stan version

4. Rlgt-package: I added some pointers to the demos but should be expanded with content discussed above.

5. DESCRIPTION FILE: revise description with "main selling points"

6. need to describe umcsent.example dataset in data.R


