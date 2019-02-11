#' This is a method of the \code{link{rlgtfit}} class to produce posterior intervals
#' 
#' @title rlgtfit posterior interval
#' @param object an object of class rlgtfit
#' @param prob percentile level to be generated (multiple values can be accepted as a vector)
#' @param type currently only central is available
#' @param ... currently not in use
#' @return confidence interval
#' @S3method posterior_interval rlgtfit
#' @method posterior_interval rlgtfit
#' @importFrom rstantools posterior_interval 
#' @examples 
# \dontrun{
#' # The following is a toy example that runs within a few seconds. To get good 
#' # fitting results the number of iterations should be set to at least 2000, and 
#' # 4 chains should be used (the default). To speed up computation the number of 
#' # cores should also be adjusted (default is 4).
#' 
#' rlgt_model <- rlgt(lynx, 
#'        control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=50, NUM_OF_CHAINS = 1, 
#'                             NUM_OF_CORES = 1), verbose=TRUE)
#'
#' # print the model details
#' posterior_interval(rlgt_model)
#}
#' @export

posterior_interval.rlgtfit <- function(object,
    prob = 0.9,
    type = "central",
    ...) {
    if (!identical(type, "central"))
      stop("Currently the only option for 'type' is 'central'.",
        call. = FALSE)
    mat <- as.matrix(object[["samples"]])
    
    posterior_interval(mat, prob = prob)
  }
