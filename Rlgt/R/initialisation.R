
#' Initialize a model from Global Trend (GT) family
#' 
#' This function is only used as a building block for the main functions
#' 
#' @title Initialize a model
#' @param model.type type of the forecasting model selected, a character object
#' @return SkeletonModel
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
	else if(model.type=="S2GT")  {
		#Non-Seasonal Local Global Trend model
		model[["parameters"]] <- c("l", "s", "sSm", "s2", "s2Sm", "nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		model[["model"]] <- stanmodels$S2GT
		class(model) <- c("RlgtStanModelS2GT")
	}
	else if (model.type=="SGTe") {
		#Seasonal Global Trend model with smoothed error size 
		model[["parameters"]] <- c("l", "s", "smoothedInnovSize", 
				"coefTrend", "powTrend", "sigma", "offsetSigma",
				"levSm", "sSm", "innovSm", "nu")
		model[["model"]] <- stanmodels$SGTe
		class(model) <- c("RlgtStanModelSGTe")
	}
  else if (model.type=="Dampen") {
    #dampen ETS fitted with Bayesian method
    model[["parameters"]] <- c("l", "b","sigma", "bSm", "levSm", "psi")
    
    model[["model"]] <- stanmodels$Dampen
    class(model) <- c("RlgtStanModelDampen")
  }
  else if (model.type=="SDampen") {
    #dampen seasonal ETS fitted with Bayesian method
    model[["parameters"]] <- c("l", "b", "s", "sSm", "sigma", "bSm", "levSm", "psi")
    
    model[["model"]] <- stanmodels$SDampen
    class(model) <- c("RlgtStanModelDampen")
  }
  else if (model.type=="TDampen") {
    #dampen ETS fitted with Bayesian method
    model[["parameters"]] <- c("l", "b","sigma", "bSm", "levSm", "psi","nu", "powx", "offsetSigma")
    
    model[["model"]] <- stanmodels$TDampen
    class(model) <- c("RlgtStanModelDampen")
  }
  else if (model.type=="TSDampen") {
    #dampen seasonal ETS fitted with Bayesian method
    model[["parameters"]] <- c("l", "b", "s", "sSm", "sigma", "bSm", "levSm", "psi", "nu", "powx", "offsetSigma")
    
    model[["model"]] <- stanmodels$TSDampen
    class(model) <- c("RlgtStanModelDampen")
  }	
  else if (model.type=="Trend") {
    #trend-only, Gaussian error, homoscedastic version of LGT
    model[["parameters"]] <- c("l", "b", "sigma", "levSm",  "bSm", "coefTrend", "powTrend", "locTrendFract")
    
    model[["model"]] <- stanmodels$trend
    class(model) <- c("RlgtStanModelTrend")
  } else if(model.type=="LGT_Reg")  {
    #Non-Seasonal Local Global Trend model
    model[["parameters"]] <- c("l", "b", "nu", "sigma", "levSm",  "bSm", 
                               "powx", "coefTrend",  "powTrend", "offsetSigma", 
                               "locTrendFract", "xreg_coef")
    model[["model"]] <- stanmodels$lgt_reg
    class(model) <- c("RlgtStanModelLGT")
  }	else if(model.type=="SGT_Reg") {
    #Seasonal Global Trend model
    model[["parameters"]] <- c("l", "s", "sSm","nu", "sigma", "levSm", 
                               "powx", "coefTrend", "powTrend", "offsetSigma",
                               "xreg_coef")
    model[["model"]] <- stanmodels$sgt_reg
    class(model) <- c("RlgtStanModelSGT")
  } 
  
  class(model) <- c("RlgtStanModel", class(model))
  model  
  
} 




#' @title lgt class
#' @description Another building block: a constructor function for the "lgt" class
#' @param y the time series data
#' @param lgtmodel type of lgtmodel selected
#' @param params list of parameters of the models (to be fitted)
#' @param control list of control parameters (chosen by the user)
#' @param samples stanfit object representing the MCMC samples
#' @return lgt instance

lgt <- function(y, model.type, has.regression, 
                lgtmodel, params, control, samples) {
	# we can add our own integrity checks
	value <- list(x = y, model.type = model.type,
	              has.regression = has.regression,
	              model = lgtmodel, params = params, 
	              control = control, samples = samples)
	
	# class can be set using class() or attr() function
	attr(value, "class") <- "lgt"
	value
}


dampen <- function(y, model.type = model.type, 
                   lgtmodel, params, paramMean, 
                   seasonality, samples) {
	# we can add our own integrity checks
	value <- list(x = y, model.type = model.type, 
	              model = lgtmodel, params = params, 
	              paramMeans = paramMean, 
	              SEASONALITY = seasonality, 
	              samples = samples)
	
	# class can be set using class() or attr() function
	attr(value, "class") <- "dampen"
	value
}

