# README #

This package is based on code from Slawek Smyl to implement LGT, a local- and global trend exponential smoothing forecasting method using Rstan for model fitting.

### Installation ###

Rlgt: R CMD build runs the "cleanup" script that runs roxygen and creates the source files. Then, R CMD INSTALL can be used to install the .tar.gz package.

RlgtLik: Implementation of LGT into code extracted from the "forecast" package.

### TODOs ###

In Rlgt:
implement seasonal models into Rlgt2 change init mechanism. Implement S3 methods according to Stan developer guide.
