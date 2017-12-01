pkgname <- "Rlgt2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Rlgt2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Rlgt2-package")
### * Rlgt2-package

flush(stderr()); flush(stdout())

### Name: Rlgt2-package
### Title: Getting started with the Rlgt2 package
### Aliases: Rlgt2-package Rlgt2
### Keywords: exponential forecasting, smoothing

### ** Examples

x <- 1



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
