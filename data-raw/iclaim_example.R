library(Rlgt)
library(gtrendsR)
library(dplyr)
library(lubridate)
iclaim <- read.csv('../data-raw/ICNSA.csv', stringsAsFactors = FALSE)
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

iclaims.example.prestransform <- iclaim %>%
  rename(claims = ICNSA, week = DATE) %>%
  mutate(week = floor_date(week, unit = "weeks", week_start = 7) + 7,
         claims = claims / 1000) %>%
  inner_join(search.unemploy.clean, by = c('week' = 'date')) %>%
  inner_join(search.filling.clean, by = c('week' = 'date')) %>%
  inner_join(search.job.clean, by = c('week' = 'date')) 

iclaims.example <- iclaims.example.prestransform %>%
  mutate_at(.vars = c("consumer.sent", "search.engine",
                      "financial.planning", "bus.news",
                      "investing", "energy.utilities"),
            .funs = log)

save(iclaims.example, file = './Rlgt/data/iclaims.example.rda')