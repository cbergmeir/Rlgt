
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
  
  
  #    lastDetails=paste("coefT=",round(coefTrendM,2), ", powT=",round(powTrendM,2),
  #      ", sigma=",round(sigmaM,2),", powx=",round(powxM,2),", offsetS=",round(offsetSigmaM,2), 
  #      ", nu=",round(nuM,2), ", lSm=",round(levSmM,2), 
  #      ", lastB=",round(lastBM,1), ", bSm=",round(bSmM,2), 
  #      ", lTFract=", round(locTrendFractM,2), sep='')				
  print(x$paramMeans)
  
  #  cat(sprintf("NumTotalEvalEA: %d\n", x$numEvalEA))
  #  cat(sprintf("NumTotalEvalLS: %d\n", x$numEvalLS))  
  #  
  #  ratio_effort <- x$numEvalEA/(x$numEvalEA+x$numEvalLS)
  #  cat(sprintf("RatioEffort EA/LS: [%.0f/%.0f]\n", 100*ratio_effort, 100*(1-ratio_effort)))
  #  
  #  ratio_alg <- x$improvementEA/(x$improvementEA+x$improvementLS)
  #  cat(sprintf("RatioImprovement EA/LS: [%.0f/%.0f]\n", 100*ratio_alg, 100*(1-ratio_alg)))
  #  #cat(sprintf(("Restarts: %d\n", restarts))
  #  
  #  if((x$numTotalEA != 0) && (x$numTotalLS != 0)) {
  #    cat(sprintf("PercentageNumImprovement[EA]: %d%%\n", round((x$numImprovementEA*100)/x$numTotalEA)))
  #    cat(sprintf("PercentageNumImprovement[LS]: %d%%\n", round((x$numImprovementLS*100)/x$numTotalLS)))    
  #  }
  #  
  #  cat(sprintf("Time[EA]: %.2f\n", x$timeMsEA))
  #  cat(sprintf("Time[LS]: %.2f\n", x$timeMsLS))
  #  cat(sprintf("Time[MA]: %.2f\n", x$timeMsMA))
  #  cat(sprintf("RatioTime[EA/MA]: %.2f\n", 100*x$timeMsEA/x$timeMsMA))
  #  cat(sprintf("RatioTime[LS/MA]: %.2f\n", 100*x$timeMsLS/x$timeMsMA))
  #  
  #  
  #  cat("Fitness:\n",sep="")
  #  print(x$fitness)
  #  cat("Solution:\n",sep="") 
  #  print(x$sol)
  
  invisible(x)
}

