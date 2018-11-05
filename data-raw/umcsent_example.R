library(gtrendsR)
library(dplyr)
library(lubridate)
umscent <- read.csv('../data-raw/UMCSENT.csv', stringsAsFactors = FALSE)
umscent$DATE <- as.Date(umscent$DATE)

# search engine
search.engine <- gtrends(keyword = NA, geo = "US",
                         category = 485, time = "all")$interest_over_time
# financial planning
financial.planning <- gtrends(keyword = NA, geo = "US",
                              category = 903, time = "all")$interest_over_time
# business news
bus.news <- gtrends(keyword = NA, geo = "US",
                    category = 784, time = "all")$interest_over_time
# investing
investing <- gtrends(keyword = NA, geo = "US",
                     category = 107, time = "all")$interest_over_time
# energy utilities
energy.utilities <- gtrends(keyword = NA, geo = "US",
                            category = 233, time = "all")$interest_over_time

search.engine$date         <- as.Date(search.engine$date)
financial.planning$date    <- as.Date(financial.planning$date)
bus.news$date              <- as.Date(bus.news$date)
investing$date             <- as.Date(investing$date)
energy.utilities$date      <- as.Date(energy.utilities$date)

search.engine <- search.engine %>%
  dplyr::select(date, hits) %>%
  rename(search.engine = hits)
financial.planning <- financial.planning %>%
  dplyr::select(date, hits) %>%
  rename(financial.planning = hits)
bus.news <- bus.news %>%
  dplyr::select(date, hits) %>%
  rename(bus.news = hits)
investing <- investing %>%
  dplyr::select(date, hits) %>%
  rename(investing = hits)
energy.utilities <- energy.utilities %>%
  dplyr::select(date, hits) %>%
  rename(energy.utilities = hits)

umcsent.example.prestransform <- umscent %>%
  rename(consumer.sent = UMCSENT, date = DATE) %>%
  inner_join(search.engine, by = c('date')) %>%
  inner_join(financial.planning, by = c('date')) %>%
  inner_join(bus.news, by = c('date')) %>%
  inner_join(investing, by = c('date')) %>%
  inner_join(energy.utilities, by = c('date'))

curr_series <- "consumer.sent"
umcsent.example <- umcsent.example.prestransform %>%
  mutate_at(.vars = c("consumer.sent", "search.engine",
                      "financial.planning", "bus.news",
                      "investing", "energy.utilities"),
            .funs = log)

save(umcsent.example, file = './data/umcsent.example.rda')
