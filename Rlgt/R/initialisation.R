
#' Initialize a model from the Rlgt family
#' 
#' A building block:to validate the model type and generate the corresponding list of parameters for the model.
#' 
#' @param model.type type of the forecasting model selected, a character object
#' @param use.regression binary parameter indicating whether additional regressors will be used for forecasting in multivariate settings.
#' @param useGeneralizedSeasonality If generalized seasonality is to be used. Default FALSE.
#' @return an RLGT skeleton model
#' 
#' @importFrom rstan stan_model
#' @export
initModel <- function(model.type=NULL, use.regression=FALSE, useGeneralizedSeasonality=FALSE) {
  
  if(is.null(model.type)) {
    print("No model type was provided, generating an LGT model.")
    model.type <- "LGT"
  }
  
  model <- list()
  
	#model[["parameters"]] are those parameters we are interested in extracting. So e.g. if the model does not use regression, although the Stan code lists regression coefs, we are not listing them here
  if(model.type=="LGT")  {
    #Non-Seasonal Local Global Trend model
    model[["parameters"]] <- c("l", "b",
			"coefTrend",  "powTrend",  "locTrendFract",
			"nu",  "levSm",  "bSm",
      "powx", "sigma", "offsetSigma"
			)
		model[["model"]] <- stanmodels$LGT
		class(model) <- c("RlgtStanModelLGT")
  }
	else if (model.type=="LGTe") {
		#Non-Seasonal Local Global Trend model with smoothed error size
		model[["parameters"]] <- c("l", "b", "smoothedInnovSize", 
				"coefTrend",  "powTrend", "locTrendFract",
				"nu",  "levSm",  "bSm",
				"innovSm", "sigma", "offsetSigma"
				)
		model[["model"]] <- stanmodels$LGTe
		class(model) <- c("RlgtStanModelLGTe")
	} 	
	else if(model.type=="SGT") {
		#Seasonal Global Trend model
		model[["parameters"]] <- c("l", "s", "sSm","nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		model[["model"]] <- stanmodels$SGT
		class(model) <- c("RlgtStanModelSGT")
	}  
	else if (model.type=="SGTe") {
		#Seasonal Global Trend model with smoothed error size 
		model[["parameters"]] <- c("l", "s", "smoothedInnovSize", 
				"coefTrend", "powTrend", "sigma", "offsetSigma",
				"levSm", "sSm", "innovSm", "nu")
		model[["model"]] <- stanmodels$SGTe
		class(model) <- c("RlgtStanModelSGTe")
	}
	else if(model.type=="S2GT")  {
		#Non-Seasonal Local Global Trend model
		model[["parameters"]] <- c("l", "s", "s2", "sSm", "s2Sm", "nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		if (useGeneralizedSeasonality) { #powSeason added later
			model[["parameters"]]<-c(model[["parameters"]],"powSeason2")
		}
		model[["model"]] <- stanmodels$S2GT
		class(model) <- c("RlgtStanModelS2GT")
	}
	else if(model.type=="S2GTe")  {
		#Non-Seasonal Local Global Trend model
		model[["parameters"]] <- c("l", "s", "s2", "sSm", "s2Sm", "nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		if (useGeneralizedSeasonality) { #powSeason added later
			model[["parameters"]]<-c(model[["parameters"]],"powSeason2")
		}
		model[["model"]] <- stanmodels$S2GTe
		class(model) <- c("RlgtStanModelS2GTe")
	}
	
  if (use.regression) {
		model[["parameters"]]=c(model[["parameters"]],"regCoef","regOffset","r")     
  }
	if (useGeneralizedSeasonality) {
		model[["parameters"]]<-c(model[["parameters"]],"powSeason")
	}
  
  class(model) <- c("RlgtStanModel", class(model))
  model  
  
} 


#' @title rlgtfit class
#' @description A building block: a constructor function for the "rlgtfit" class. 
#' This class will be used as an output to the \code{\link{rlgt}} function.
#' @param y time-series data for training (provided as a vector or a ts object).
#' @param model.type the type of rlgt model
#' @param use.regression whether the data has any additional variables to be used with forecasting, i.e. multivariate time-series.
#' @param useGeneralizedSeasonality If generalized seasonality is to be used.
#' @param levelMethodId Used with dual seasonality models
#' @param seasonality 
#' @param seasonality2 
#' @param rlgtmodel an rlgt model.
#' @param params list of parameters of the model (to be fitted).
#' @param control list of control parameters, i.e. hyperparameter values 
#' for the model's prior distribution. See \code{\link{rlgt.control}}
#' @param samples stanfit object representing the MCMC samples
#' @return an rlgtfit instance

rlgtfit <- function(y, model.type, use.regression,
		useGeneralizedSeasonality,  levelMethodId=levelMethodId, 
		seasonality, seasonality2,
    rlgtmodel, params, control, samples) {
	# we can add our own integrity checks
	value <- list(x = y, model.type = model.type,
	              use.regression = use.regression,
								useGeneralizedSeasonality=useGeneralizedSeasonality,
								levelMethodId=levelMethodId,  
								seasonality=seasonality, seasonality2=seasonality2,
	              model = rlgtmodel, params = params, 
	              control = control, samples = samples)
	
	# class can be set using class() or attr() function
	attr(value, "class") <- "rlgtfit"
	value
}
