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
#'\dontrun{
#' rlgt_model <- rlgt(lynx,
#'      control=rlgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), verbose=TRUE)
#'
#' # print the model details
#' posterior_interval(rlgt_model)
#'}
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
