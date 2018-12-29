
#' Initialize a model from the Rlgt family
#' 
#' This is an internal function that usually won't be called by users directly. It validates the model type and generates the corresponding list of parameters for the model.
#' 
#' @param model.type type of the forecasting model selected, a character object
#' @param use.regression binary parameter indicating whether additional regressors will be used for forecasting in multivariate settings.
#' @param seasonalityMethodId Seasonality method Id (0- HW, 1- generalized).
#' @param levelMethodId Level method Id.
#' @param useSmoothingMethodForError if the non-standard function for error size should be used, one based on smoothed innovations or surprises 
#' @return an Rlgt skeleton model
#' 
#' @importFrom rstan stan_model
#' @export
initModel <- function(model.type=NULL, use.regression=FALSE, 
		seasonalityMethodId=0, levelMethodId=0,   
		useSmoothingMethodForError=FALSE) {
  
  if(is.null(model.type)) {
    print("No model type was provided, generating an LGT model.")
    model.type <- "LGT"
  }
  
  model <- list()
  
	#model[["parameters"]] are those parameters we are interested in extracting. So e.g. if the model does not use regression, although the Stan code lists regression coefs, we are not listing them here
  if (model.type=="LGT")  {
    #Non-Seasonal Local Global Trend model
		if (useSmoothingMethodForError) {
			model[["parameters"]] <- c("l", "b", "smoothedInnovSize", "innovSm",
					"coefTrend",  "powTrend", "locTrendFract", "sigma", "offsetSigma",
					"levSm", "bSm", "nu")		
		} else {
			model[["parameters"]] <- c("l", "b",
					"coefTrend",  "powTrend", "locTrendFract", "sigma", "offsetSigma",
					"levSm", "bSm", "nu", "powx")	
		}	
		model[["model"]] <- stanmodels$LGT
		class(model) <- c("RlgtStanModelLGT")
  }	
	else if(model.type=="SGT") {
		#Seasonal Global Trend model
		if (useSmoothingMethodForError) {
			model[["parameters"]] <- c("l", "s", "smoothedInnovSize", "innovSm",
					"coefTrend", "powTrend", "sigma", "offsetSigma",
					"levSm", "sSm", "nu")
		} else {
			model[["parameters"]] <- c("l", "s", 
					"coefTrend", "powTrend", "sigma", "offsetSigma",
					"levSm", "sSm", "nu", "powx")
		}
		model[["model"]] <- stanmodels$SGT
		class(model) <- c("RlgtStanModelSGT")
	}  
	else if(model.type=="S2GT")  { #dual seasonal
		if (useSmoothingMethodForError) {
			model[["parameters"]] <- c("l", "s", "s2", "smoothedInnovSize", "innovSm",
					"coefTrend", "powTrend", "sigma", "offsetSigma", 
					"levSm", "sSm", "s2Sm", "nu")			
		} else {
			model[["parameters"]] <- c("l", "s", "s2", 
					"coefTrend", "powTrend", "sigma", "offsetSigma",
					"levSm", "sSm", "s2Sm", "nu", "powx")
		}
		if (seasonalityMethodId==1) { #powSeason added later
			model[["parameters"]] <- c(model[["parameters"]], "powSeason2")
		}
		
		model[["model"]] <- stanmodels$S2GT
		class(model) <- c("RlgtStanModelS2GT")
	}

  if (use.regression) {
		model[["parameters"]] <- c(model[["parameters"]], "regCoef", "regOffset", "r")     
  }
	if (seasonalityMethodId==1) {
		model[["parameters"]] <- c(model[["parameters"]],"powSeason")
	}  
	if (levelMethodId==3) {
		model[["parameters"]] <- c("l0", model[["parameters"]],"llevSm")
	}
  
  class(model) <- c("RlgtStanModel", class(model))
  model  
  
} 


#' @title rlgtfit class
#' @description A constructor function for objects of class \code{rlgtfit}, the main class of the package. Objects of this class 
#' are output from the \code{\link{rlgt}} function. This constructor will usually not be called by users directly.
#' @param y time series data for training (provided as a vector or a ts object).
#' @param model.type the type of rlgt model, one of: "LGT", "SGT", "S2GT"
#' @param use.regression whether the data has any additional variables to be used with forecasting, i.e. multivariate time-series.
#' @param seasonalityMethodId Seasonality method Id (0- HW, 1- generalized).
#' @param levelMethodId Level method Id.
#' @param useSmoothingMethodForError if the non-standard function for error size should be used, one based on smoothed innovations or surprises 
#' @param seasonality This specification of seasonality will be overridden by frequency of y, if y is of ts or msts class. 
#' 1 by default, i.e. no seasonality.
#' @param seasonality2 Second seasonality. If larger than 1, a dual seasonality model will be used. 
#' This specification of seasonality will be overridden by the second seasonality of y, if y is of msts class. 
#' 1 by default, i.e. no seasonality or single seasonality.
#' @param rlgtmodel an rlgt model.
#' @param params list of parameters of the model (to be fitted).
#' @param control list of control parameters, i.e. hyperparameter values 
#' for the model's prior distribution. See \code{\link{rlgt.control}}
#' @param samples stanfit object representing the MCMC samples
#' @return an rlgtfit instance

rlgtfit <- function(y, model.type, use.regression,
		seasonalityMethodId, levelMethodId,  
		useSmoothingMethodForError=FALSE,
		seasonality, seasonality2,
    rlgtmodel, params, control, samples) {
	# we can add our own integrity checks
	value <- list(x = y, model.type = model.type,
	              use.regression = use.regression,
								seasonalityMethodId=seasonalityMethodId, levelMethodId=levelMethodId,
								useSmoothingMethodForError=useSmoothingMethodForError,
								seasonality=seasonality, seasonality2=seasonality2,
	              model = rlgtmodel, params = params, 
	              control = control, samples = samples)
	
	# class can be set using class() or attr() function
	attr(value, "class") <- "rlgtfit"
	value
}
