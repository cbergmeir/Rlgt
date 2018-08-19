
#' Initialize a model from the Rlgt family
#' 
#' A building block:to validate the model type and generate the corresponding list of parameters for the model.
#' 
#' @param model.type type of the forecasting model selected, a character object
#' @param use.regression binary parameter indicating whether additional regressors will be used for forecasting in multivariate settings.
#' @return an RLGT skeleton model
#' 
#' @importFrom rstan stan_model
#' @export
initModel <- function(model.type=NULL, use.regression=FALSE){
  
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
	else if (model.type=="LGTe") {
		#Non-Seasonal Local Global Trend model with smoothed error size
		model[["parameters"]] <- c("l", "b", "smoothedInnovSize", 
				"coefTrend",  "powTrend", "sigma", "offsetSigma",
				"bSm", "locTrendFract", "levSm", "innovSm", "nu")
		
		model[["model"]] <- stanmodels$LGTe
		class(model) <- c("RlgtStanModelLGTe")
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
		model[["parameters"]] <- c("l", "s", "sSm", "s2", "s2Sm", "nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		model[["model"]] <- stanmodels$S2GT
		class(model) <- c("RlgtStanModelS2GT")
	}

  if (use.regression) {
    # append regression parameters
    model[["parameters"]] <- c(model[["parameters"]], "regCoef")
    # only support LGT and SGT
    if (model.type == "LGT")  {
      #Non-Seasonal Local Global Trend model
      model[["model"]] <- stanmodels$lgt_reg
      class(model) <- c("RlgtStanModelLGT_REG")
    } else if(model.type == "SGT") {
      #Seasonal Global Trend model
      model[["model"]] <- stanmodels$sgt_reg
      class(model) <- c("RlgtStanModelSGT_REG")
    } 
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
#' @param rlgtmodel an rlgt model.
#' @param params list of parameters of the model (to be fitted).
#' @param control list of control parameters, i.e. hyperparameter values 
#' for the model's prior distribution. See \code{\link{rlgt.control}}
#' @param samples stanfit object representing the MCMC samples
#' @return an rlgtfit instance

rlgtfit <- function(y, model.type, use.regression,
                    rlgtmodel, params, control, samples) {
	# we can add our own integrity checks
	value <- list(x = y, model.type = model.type,
	              use.regression = use.regression,
	              model = rlgtmodel, params = params, 
	              control = control, samples = samples)
	
	# class can be set using class() or attr() function
	attr(value, "class") <- "rlgtfit"
	value
}




