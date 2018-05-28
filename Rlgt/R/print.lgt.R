
#' Print out some characteristics of a \code{\link{lgt}} model.
#'  
#' @title Generic print function for lgt models
#' @param x the \code{\link{lgt}} model
#' @param ... additional function parameters (currently not used)
#' @aliases summary.lgt
#' @export
# @S3method print lgt
# @method print lgt
# @rdname lgt
print.lgt <- function(x, ...) {
  if(!inherits(x, "lgt")) stop("not a legitimate lgt result")
  
  print(x$paramMeans)
  
  invisible(x)
}
