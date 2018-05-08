

#' This is a method of lgt object to produce posterior interval
#' 
#' @title lgt posterior interval
#' @param object an object of class lgt
#' @param prob percentile level to be generated (multiple values can be accepted as a vector)
#' @param type currently only central is available
#' @param ... currently not in use
#' @return confidence interval
#' @S3method posterior_interval lgt
#' @method posterior_interval lgt
#' @importFrom rstantools posterior_interval 
#' @author wibowo
#' @examples 
#'\dontrun{
#' lgt_model <- fit.lgt(lynx, model="LGT", nCores=4, nChains=4,
#' control=lgt.control(MAX_NUM_OF_REPEATS=1, NUM_OF_ITER=2000), 
#' verbose=TRUE)

#' # print the model details
#' posterior_interval(lgt_model)
#'}
#' @export

posterior_interval.lgt <- function(object,
    prob = 0.9,
    type = "central",
    ...) {
    if (!identical(type, "central"))
      stop("Currently the only option for 'type' is 'central'.",
        call. = FALSE)
    mat <- as.matrix(object[["samples"]])
    
    posterior_interval(mat, prob = prob)
  }
