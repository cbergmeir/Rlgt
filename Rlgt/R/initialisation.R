
#' Initialize a model from the Rlgt family
#' 
#' A building block:to validate the model type and generate the corresponding list of parameters for the model.
#' 
#' @param model.type type of the forecasting model selected, a character object
#' @return an RLGT skeleton model
#' 
#' @importFrom rstan stan_model
#' @export
initModel <- function(model.type = NULL){
  
  if(is.null(model.type)) {
    print("No model type was provided, generating an LGT model.")
    model.type <- "LGT"
  }
  
  model <- list()
  
  if(model.type=="LGT")  {
    #Non-Seasonal Local Global Trend model
    model[["parameters"]] <- c("l", "b", "nu", "sigma", "levSm",  "bSm", 
      "powx", "coefTrend",  "powTrend", "offsetSigma", "locTrendFract")
    
    model[["model"]] <- stanmodels$lgt
    class(model) <- c("RlgtStanModelLGT")
  }	
	else if(model.type=="SGT") {
		#Seasonal Global Trend model
		model[["parameters"]] <- c("l", "s", "sSm","nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		
		model[["model"]] <- stanmodels$sgt
		class(model) <- c("RlgtStanModelSGT")
	} 
	else if(model.type=="gSGT") {
		#generalized Seasonality Global Trend model
		model[["parameters"]] <- c("l", "s", "sSm","nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma", "powSeason")
		
		model[["model"]] <- stanmodels$gSGT
		class(model) <- c("RlgtStanModelgSGT")
	} 
	else if(model.type=="S2GT")  {
		#Non-Seasonal Local Global Trend model
		model[["parameters"]] <- c("l", "s", "sSm", "s2", "s2Sm", "nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		model[["model"]] <- stanmodels$S2GT
		class(model) <- c("RlgtStanModelS2GT")
	}

  
  class(model) <- c("RlgtStanModel", class(model))
  model  
  
} 


#' @title rlgtfit class
#' @description A building block: a constructor function for the "rlgtfit" class. 
#' This class will be used as an output to the \code{\link{rlgt}} function.
#' @param y time-series data for training (provided as a vector or a ts object).
#' @param model.type the type of rlgt model
#' @param rlgtmodel an rlgt model.
#' @param params list of parameters of the model (to be fitted).
#' @param control list of control parameters, i.e. hyperparameter values 
#' for the model's prior distribution. See \code{\link{rlgt.control}}
#' @param samples stanfit object representing the MCMC samples
#' @return an rlgtfit instance

rlgtfit <- function(y, model.type, rlgtmodel, params, control, samples) {
	# we can add our own integrity checks
	value <- list(x = y, model.type = model.type,
	              has.regression = has.regression,
	              model = lgtmodel, params = params, 
	              control = control, samples = samples)
	
	# class can be set using class() or attr() function
	attr(value, "class") <- "rlgtfit"
	value
}




