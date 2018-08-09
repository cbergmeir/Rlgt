#' Print out some characteristics of a \code{\link{rlgtfit}} model.
#'  
#' @title Generic print function for lgt models
#' @param x an rlgtfit object
#' @param ... additional function parameters (currently not used)
#' @aliases summary.rlgt
#' @importFrom stats median
#' @export
#' @S3method print rlgtfit
#' @method print rlgtfit
#' @rdname lgt
print.rlgtfit <- function(x, ...) {
  if(!inherits(x, "rlgtfit")) stop("not a legitimate result of an Rlgt model")
  
  print(lapply(x$params,median))
  
  invisible(x)
}
