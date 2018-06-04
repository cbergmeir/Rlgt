
#' Initialize a model from Global Trend (GT) family
#' 
#' This function is only used as a building block for the main functions
#' 
#' @title Initialize a model
#' @param modelType type of the forecasting model selected, a character object
#' @return SkeletonModel
#' 
#' @importFrom rstan stan_model
#' @export
initModel <- function(modelType = NULL){
  
  if(is.null(modelType)) {
    print("No model type was provided, generating an LGT model.")
    modelType <- "LGT"
  }
  
  model <- list()
  
  if(modelType=="LGT")  {
    #Non-Seasonal Local Global Trend model
    model[["parameters"]] <- c("l", "b", "nu", "sigma", "levSm",  "bSm", 
      "powx", "coefTrend",  "powTrend", "offsetSigma", "locTrendFract")
    
    model[["model"]] <- stanmodels$lgt
    class(model) <- c("RlgtStanModelLGT")
  }	
	else if (modelType=="LGTe") {
		#Non-Seasonal Local Global Trend model with smoothed error size
		model[["parameters"]] <- c("l", "b", "smoothedInnovSize", 
				"coefTrend",  "powTrend", "sigma", "offsetSigma",
				"bSm", "locTrendFract", "levSm", "innovSm", "nu")
		
		model[["model"]] <- stanmodels$LGTe
		class(model) <- c("RlgtStanModelLGTe")
	} 
	else if(modelType=="SGT") {
		#Seasonal Global Trend model
		model[["parameters"]] <- c("l", "s", "sSm","nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		
		model[["model"]] <- stanmodels$sgt
		class(model) <- c("RlgtStanModelSGT")
	} 
	else if(modelType=="S2GT")  {
		#Non-Seasonal Local Global Trend model
		model[["parameters"]] <- c("l", "s", "sSm", "s2", "s2Sm", "nu", "sigma", "levSm", 
				"powx", "coefTrend", "powTrend", "offsetSigma")
		model[["model"]] <- stanmodels$S2GT
		class(model) <- c("RlgtStanModelS2GT")
	}
	else if (modelType=="SGTe") {
		#Seasonal Global Trend model with smoothed error size 
		model[["parameters"]] <- c("l", "s", "smoothedInnovSize", 
				"coefTrend", "powTrend", "sigma", "offsetSigma",
				"levSm", "sSm", "innovSm", "nu")
		model[["model"]] <- stanmodels$SGTe
		class(model) <- c("RlgtStanModelSGTe")
	}
  else if (modelType=="Dampen") {
    #dampen ETS fitted with Bayesian method
    model[["parameters"]] <- c("l", "b","sigma", "bSm", "levSm", "psi")
    
    model[["model"]] <- stanmodels$Dampen
    class(model) <- c("RlgtStanModelDampen")
  }
  else if (modelType=="SDampen") {
    #dampen seasonal ETS fitted with Bayesian method
    model[["parameters"]] <- c("l", "b", "s", "sSm", "sigma", "bSm", "levSm", "psi")
    
    model[["model"]] <- stanmodels$SDampen
    class(model) <- c("RlgtStanModelDampen")
  }
  else if (modelType=="TDampen") {
    #dampen ETS fitted with Bayesian method
    model[["parameters"]] <- c("l", "b","sigma", "bSm", "levSm", "psi","nu", "powx", "offsetSigma")
    
    model[["model"]] <- stanmodels$TDampen
    class(model) <- c("RlgtStanModelDampen")
  }
  else if (modelType=="TSDampen") {
    #dampen seasonal ETS fitted with Bayesian method
    model[["parameters"]] <- c("l", "b", "s", "sSm", "sigma", "bSm", "levSm", "psi", "nu", "powx", "offsetSigma")
    
    model[["model"]] <- stanmodels$TSDampen
    class(model) <- c("RlgtStanModelDampen")
  }	
  else if (modelType=="Trend") {
    #trend-only, Gaussian error, homoscedastic version of LGT
    model[["parameters"]] <- c("l", "b", "sigma", "levSm",  "bSm", "coefTrend", "powTrend", "locTrendFract")
    
    model[["model"]] <- stanmodels$trend
    class(model) <- c("RlgtStanModelTrend")
  } 
  
  class(model) <- c("RlgtStanModel", class(model))
  model  
  
} 

