library(Rlgt)
library(gtrendsR)
library(dplyr)
library(lubridate)
iclaim <- read.csv('./data-raw/ICNSA.csv', stringsAsFactors = FALSE)
iclaim$DATE <- as.Date(iclaim$DATE)

x <- gtrends(keyword = "unemployment", geo = "US")$interest_over_time

search.unemploy1 <- gtrends(keyword = "unemployment", geo = "US", 
                           time = "2010-01-01 2014-12-31")$interest_over_time
search.unemploy2 <- gtrends(keyword = "unemployment", geo = "US", 
                            time = "2015-01-01 2018-06-30")$interest_over_time
search.unemploy <- rbind(search.unemploy1, search.unemploy2)
search.filling1 <- gtrends(keyword = "filling", geo = "US", 
                          time = "2010-01-01 2014-12-31")$interest_over_time
search.filling2 <- gtrends(keyword = "filling", geo = "US", 
                          time = "2015-01-01 2018-06-30")$interest_over_time
search.filling <- rbind(search.filling1, search.filling2)
search.job1 <- gtrends(keyword = "job", geo = "US", 
                      time = "2010-01-01 2014-12-31")$interest_over_time
search.job2 <- gtrends(keyword = "job", geo = "US", 
                      time = "2015-01-01 2018-06-30")$interest_over_time
search.job <- rbind(search.job1, search.job2)

search.unemploy$date <- as.Date(search.unemploy$date)
search.filling$date  <- as.Date(search.filling$date)
search.job$date      <- as.Date(search.job$date)

search.unemploy.clean <- search.unemploy %>%
  dplyr::select(date, hits) %>%
  rename(trend.unemploy = hits)
search.filling.clean <- search.filling %>%
  dplyr::select(date, hits) %>%
  rename(trend.filling = hits)
search.job.clean <- search.job %>%
  dplyr::select(date, hits) %>%
  rename(trend.job = hits)

iclaim_example <- iclaim %>%
  rename(claims = ICNSA , week = DATE) %>%
  mutate(week = floor_date(week, unit = "weeks", week_start = 7) + 7) %>%
  inner_join(search.unemploy.clean, by = c('week' = 'date')) %>%
  inner_join(search.filling.clean, by = c('week' = 'date')) %>%
  inner_join(search.job.clean, by = c('week' = 'date')) 

#set.seed(12)
options(width=180)

# predict initial unemployment claims of US based on google search data
curr_series <- "iclaimsNSA"
data.train  <- log(iclaim_example$claims)
sizeTestSet <- length(data.train)

mod <- list()
forecasts <- list()
#--------------------------------
#Fit LGT model
# debug(rlgt)
mod[["LGT"]] <- rlgt(data.train, model.type ="LGT", nCores=4, nChains=4,
                     control=lgt.control(MAX_NUM_OF_REPEATS=10, 
                                         NUM_OF_ITER=2000), 
                     verbose=TRUE)
# print the model details
print(mod[["LGT"]])

# print the interval for all vars
posterior_interval(mod[["LGT"]])

forecasts[["LGT"]] <- forecast(mod[["LGT"]], h = sizeTestSet/2, 
                               level=c(80, 95, 98))
plot(forecasts[["LGT"]], main=paste(curr_series,'by LGT'))

# Use AirPassanger data as an example for a seasonal dataset
seasonal_data <- AirPassengers
curr_series <- "AirPassengers"
data.train <- seasonal_data 
sizeTestSet <- frequency(AirPassengers)
#--------------------------------
#Fit SGT model
mod[["SGT"]] <- rlgt(data.train, model.type="SGT", nCores=4, nChains=4,
                     control=lgt.control(MAX_NUM_OF_REPEATS=3, NUM_OF_ITER=1000), 
                     verbose=TRUE)
# print the model details
print(mod[["SGT"]])

# print the interval for all vars
posterior_interval(mod[["SGT"]])

forecasts[["SGT"]] <- forecast(mod[["SGT"]], h = sizeTestSet, level=c(80, 95, 98))
plot(forecasts[["SGT"]],main=paste(curr_series,'by SGT'))



