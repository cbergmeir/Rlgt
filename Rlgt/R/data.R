#' Weekly Initial Claims of US Unemployment Benefits & Google Trends Queries
#'
#' A dataset containing the weekly initial claims for US unemployment benefits aginst 
#' a few related google trend queries
#' since 2010
#'
#' @usage data("iclaims.example")
#' @name iclaims.example
#' @format A data frame with 443 rows and 5 variables:
#' \describe{
#'   \item{week}{date of records starting by Mondays with US calendar format}
#'   \item{claims}{weekly initial claims of unemployment benefits in thousands with log-transformation}
#'   \item{trend.unemploy}{normalized trend queries retreived from gtrendsR API with log-transformation}
#'   \item{trend.filling}{normalized trend queries retreived from gtrendsR API with log-transformation}
#'   \item{trend.job}{normalized trend queries retreived from gtrendsR API with log-transformation}
#' }
#' @source \url{https://fred.stlouisfed.org/series/ICNSA}
#' @source \url{https://trends.google.com/trends/?geo=US}
#' @source \url{https://CRAN.R-project.org/package=gtrendsR}
#' @docType data
#' @keywords datasets
#' @references 
#' U.S. Employment and Training Administration, Initial Claims [ICNSA], retrieved from FRED, Federal Reserve Bank of St. Louis; \url{https://fred.stlouisfed.org/series/ICNSA}, October 27, 2018.
#' @references
#' Trends queries from google search engine. 
#' {https://trends.google.com/trends/?geo=US}
#' @references 
#' An interface for retrieving and displaying the information returned online by Google Trends is provided. Trends (number of hits) over the time as well as geographic representation of the results can be displayed.
#' \url{https://CRAN.R-project.org/package=gtrendsR} 
#' 
NULL
