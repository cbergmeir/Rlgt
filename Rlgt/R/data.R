#' Weekly Initial Claims of US Unemployment Benefits & Google Trends Queries
#'
#' A dataset containing the weekly initial claims for US unemployment benefits against 
#' a few related Google trend queries from Jan 2010 - June 2018. 
#' This aims to mimick the dataset from Scott and Varian (2014).
#'
#' @usage data("iclaims.example")
#' @name iclaims.example
#' @format A data frame with 443 rows and 5 variables with log-transformation
#' \describe{
#'   \item{week}{date of records starting by Mondays with US calendar format}
#'   \item{claims}{weekly initial claims of unemployment benefits in thousands}
#'   \item{trend.unemploy}{normalized trend queries retreived from gtrendsR API}
#'   \item{trend.filling}{normalized trend queries retreived from gtrendsR API}
#'   \item{trend.job}{normalized trend queries retreived from gtrendsR API}
#' }
#' @docType data
#' @keywords datasets
#' @references 
#' U.S. Employment and Training Administration, Initial Claims [ICNSA], retrieved from FRED, Federal Reserve Bank of St. Louis; 
#' \url{https://fred.stlouisfed.org/series/ICNSA}, October 27, 2018.
#' @references
#' Trend queries from Google search engine. 
#' \url{https://trends.google.com/trends/?geo=US}
#' @references 
#' An interface for retrieving and displaying the information returned online by Google Trends is provided. Trends (number of hits) over the time as well as geographic representation of the results can be displayed.
#' \url{https://CRAN.R-project.org/package=gtrendsR} 
#' @references 
#' Scott, S. L. and Varian, H. R. (2014). Predicting the Present with Bayesian Structural Time Series.
#' International Journal of Mathematical Modeling and Optimization 5 4–23.
#' \url{http://people.ischool.berkeley.edu/~hal/Papers/2013/pred-present-with-bsts.pdf} 
#' 
NULL

#' University of Michigan Monthly Survey of Consumer Sentiment  & Google Trends Queries
#'
#' A dataset containing monthly University of Michigan survey of Consumer Sentiment
#' along a few related google trend queries Jan from 2014 - June 2018. 
#' This aims to mimick the dataset from Scott and Varian (2014).
#'
#' @usage data("umcsent.example")
#' @name umcsent.example
#' @format A data frame with 174 rows and 8 variables with log-transformation
#' \describe{
#'   \item{date}{first date of each month in US calendar format}
#'   \item{consumer.sent}{monthly initial claims of University of Michigan: Consumer Sentiment}
#'   \item{search.engine}{normalized trend queries retreived from gtrendsR API}
#'   \item{financial.planning}{normalized trend queries retreived from gtrendsR API}
#'   \item{bus.news}{normalized trend queries retreived from gtrendsR API}
#'   \item{investing}{normalized trend queries retreived from gtrendsR API}
#'   \item{energy.utilities}{normalized trend queries retreived from gtrendsR API}
#' }
#' @docType data
#' @keywords datasets
#' @references 
#' University of Michigan, University of Michigan: Consumer Sentiment [UMCSENT], retrieved from FRED, Federal Reserve Bank of St. Louis; 
#' \url{https://fred.stlouisfed.org/series/UMCSENT}, November 17, 2018.
#' @references
#' Trends queries from google search engine. 
#' \url{https://trends.google.com/trends/?geo=US}
#' @references 
#' An interface for retrieving and displaying the information returned online by Google Trends is provided. Trends (number of hits) over the time as well as geographic representation of the results can be displayed.
#' \url{https://CRAN.R-project.org/package=gtrendsR} 
#' @references 
#' Scott, S. L. and Varian, H. R. (2012). Bayesian Variable Selection for Nowcasting Economic Time Series.
#' \url{https://www.aeaweb.org/conference/2013/retrieve.php?pdfid=447} 
#' 
NULL
