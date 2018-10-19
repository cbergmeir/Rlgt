#' Print out some characteristics of an \code{\link{rlgtfit}} model.
#'  
#' @title Generic print function for rlgtfit models
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
  
  if (x$use.regression) {
    temp <- lapply(x$params,median)
    temp[['regCoef']] <- colMeans(x$params[['regCoef']])
    print(temp)
  } else {
    print(lapply(x$params,median))
  }
  
  invisible(x)
}
