#' ---
#' title: "East Africa rainfall data analysis"
#' author: "Joris Wiethase"
#' date: "Nov 28, 2019"
#' ---
  
# Start with clean environment
rm(list = ls(all=TRUE))  

# Load packages
suppressWarnings(suppressMessages(suppressPackageStartupMessages({
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(lubridate)
library(scales)
library(zoo)
})))

# Set global ggplot theme
theme_set(ggthemes::theme_few(base_size = 12))

# Function for outlier detection
is_outlier <- function(x) {
  return(x < quantile(x, 0.25, na.rm = TRUE) - 1.5 * IQR(x, na.rm = TRUE) |
           x > quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE))
}

#'## Data import and preparation
# Import CHIRPS data set. Downloaded from http://ClimateSERV.servirglobal.net/ 
# as average rainfall per day[mm] data set 
df <- read.csv("raw_data/CHIRPS_av_day_small/b57196d3-2259-4709-bf40-e173270c4a4b.csv", skip = 1) 

# Edit column names
colnames(df) <- c("date", "prec_day") 

# Extract time values 
df_chirps <- df %>% 
  mutate(date = mdy(date),
         Month = month(date),
         month_ab = month.abb[Month],
         Week = week(date),
         Day = yday(date),
         Year = year(date),
         pentad = ceiling( (Day - leap_year(Year)*(Day > 59)) / 5 ),
         month_ab = factor(month_ab, levels = c( "Aug", "Sep", "Oct", "Nov",
                                                 "Dec", "Jan", "Feb", "Mar", 
                                                 "Apr", "May", "Jun", "Jul")),
         decade = NA)

# Define decades
df_chirps$decade[df_chirps$Year > 1980 & df_chirps$Year <= 1990] <- "1981-1990"
df_chirps$decade[df_chirps$Year > 1990 & df_chirps$Year <= 2000] <- "1991-2000"
df_chirps$decade[df_chirps$Year > 2000 & df_chirps$Year <= 2010] <- "2001-2010"
df_chirps$decade[df_chirps$Year > 2010 & df_chirps$Year <= 2020] <- "2011-2018"

#'## Data exploration

#' Re-arrange months, average monthly rainfall
df_monthly <- df_chirps %>% 
  group_by(Year, month_ab, decade) %>% 
  summarize(rain_month = sum(prec_day)) %>% 
  ungroup() %>% 
  mutate(month_ab = factor(month_ab, levels = c( "Aug", "Sep", "Oct", "Nov", 
                                                 "Dec", "Jan", "Feb", "Mar", 
                                                 "Apr", "May", "Jun", "Jul"))) %>% 
  group_by(month_ab, decade) %>% 
  summarize(median_month = median(rain_month))

#' Plot monthly average rainfall for every decade
ggplot(df_monthly) +
  geom_bar(aes(x = month_ab, y = median_month), stat = "identity") +
  facet_wrap(facets = "decade") + 
  theme_bw() 

#' Summarize rainfall for every month per year
df_chirps_m <- df_chirps %>% 
  group_by(Year, Month, decade, month_ab) %>% 
  summarize(prec_month = round(sum(prec_day), digits = 1)) 

#' Plot total monthly rainfall per year
ggplot(df_chirps_m) +
  geom_bar(aes(x = month_ab, y = prec_month), stat = "identity") +
  facet_wrap(facets = "Year") + 
  theme_bw()  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#'## Determine onset and cessation dates of rainy season for our study area, 
#'## following Liebmann et al. 2013, Dunning et al. 2016

# From Dunning et al. 2016, onset is ~Dec 30, and cessation is ~April 15, which 
# seems to agree with CHIRPS rainfall data. Calculate season dates using climatological 
# accumulated rainfall anomalies:
  
#' Using decades: Climatological mean rainfall for each day of the calendar year, 
#' averaged over one decade
df_chirps_av <- df_chirps %>% 
  group_by(decade, Day) %>% 
  summarize(av_day = mean(prec_day))

rain_seas_av <- df_chirps_av %>% 
  group_by(decade) %>% 
  mutate(an_mean_day = mean(av_day),    # Long-term climatological daily median rainfall 
                                        # (or overall daily average per decade)
         rain_ann = av_day-an_mean_day, # The total precipitation per day, minus the decadal
                                        # average for the same day = daily median rainfall anomaly
         an_accum = cumsum(rain_ann))

#' Plot the daily median rainfall vs rain anamoly
ggplot(rain_seas_av) +
  geom_line(aes(x = Day, y = av_day), col = "red") +
  geom_line(aes(x = Day, y = rain_ann), col = "blue") +
  facet_wrap(facets = "decade") + 
  ylab("Mean rainfall (mm)") 

#' Plot anomalous accumulation 
ggplot(rain_seas_av) +
  geom_line(aes(x = Day, y = an_accum, col = decade)) +
  ylab("Cumulative daily mean rainfall anomaly (mm)") 


#' Climatological wet season seems span the year break. Better to shift to strat at day 200

# Desired day-of-year sequence
seq <- c(200:366, 1:199)

#' Re-arrange the day of year
re_df_chirps_av <- df_chirps_av %>% 
  group_by(decade) %>% 
  arrange(factor(Day, levels = seq))

#' Calculate anomalous accumulation again
re_rain_seas_av <- re_df_chirps_av %>% 
  group_by(decade) %>% 
  mutate(an_mean_day = mean(av_day),         
         rain_ann = av_day-an_mean_day,      
         an_accum = cumsum(rain_ann),
         Day = factor(Day, levels = Day),
         onset = min(an_accum),
         cessation = max(an_accum),
         day_on = Day[an_accum == onset],
         day_cess = Day[an_accum == cessation])

#' Plot anomalous accumulation. Highlight the onset and cessation dates
ggplot(re_rain_seas_av, aes(x = Day, y = an_accum, group = decade, col = decade)) +
  geom_point(stat='summary', fun.y=sum, pch = ".") +
  stat_summary(fun.y=sum, geom="line") +
  ylab("Cumulative daily median rainfall anomaly (mm)")  +
  scale_x_discrete(breaks = seq[c(F, F, F, F, F, T)]) +
  scale_color_manual(values = c("red", "green", "blue", "orange")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#' The plot suggests that there might be two rainy seasons (starting around day 7-30).
#' Calculate biannual seasons following Dunning et al. 2016.

#' Smooth line using 21 day running median (S(d))

#' Dunning et al. 2016 use 30 day running mean. We choose a 15 day runing mean,
#' since the seasons seem to be more subtly separated 

smooth_df <- re_rain_seas_av %>% 
  group_by(decade) %>% 
  arrange(Day) %>% 
  mutate(roll_mean = rollmean(an_accum, 15, na.pad = T))

#' Identify days where S(d) lower/higher than the four preceding/following day
k <- 5
smooth_df <- smooth_df %>% 
  group_by(decade) %>% 
  mutate(Season = rollapplyr(roll_mean, k, function(x) all(x[k] > x[-k]), fill = NA)) %>% 
  mutate(Season = ifelse(Season == TRUE, "wet", "dry"))

#' The two longest seasons are considered to be the two seasons of interest.
#' Calculate length of seasons
smooth_df <- smooth_df %>% 
  # Group by decade and different seasons
  group_by(decade, run = data.table::rleid(Season)) %>% 
  # Get length of each season
  mutate(seas_length = n()) %>% 
  # Filter very short seasons (smaller than 5)
  filter(seas_length > 3) %>% 
  ungroup() %>% 
  # Get length of each season again, with subsetted data
  group_by(decade) %>% 
  mutate(run = data.table::rleid(Season)) %>% 
  group_by(decade, run) %>% 
  mutate(seas_length = n())
  
#' Get the onset and cessation dates for every rainy season
seas_dates <- smooth_df %>% 
  # Subset rainy seasons, of at least a week length
  filter(Season == "wet") %>% 
  group_by(decade, run) %>% 
  summarize(long_day_on = as.numeric(as.character(dplyr::first(Day))),
            long_day_cess = as.numeric(as.character(dplyr::last(Day))), 
            seas = seas_length[1]) %>% 
  group_by(decade) %>% 
  # Only keep long rain season
  filter(seas == max(seas)) %>% 
  mutate(date_on = as.Date(long_day_on, origin = "2018-01-01"),   # Some example year
         month_on = month.abb[month(date_on)],
         date_cess = as.Date(long_day_cess, origin = "2018-01-01"),
         month_cess = month.abb[month(date_cess)])

#' Plot the accumulation curve
p1 <- ggplot(smooth_df, aes(x = Day, y = roll_mean, group = decade, col = Season)) +
  geom_point(stat='summary', fun.y=sum, pch = ".") +
  stat_summary(fun.y=sum, geom="line", size = 1.5) +
  facet_wrap(facets = "decade")+
  ylab("Cumulative daily mean rainfall anomaly (mm)")  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(breaks = seq[c(F, F, F, F, F, F, F, F, F, F, F, F, F, T)]); p1



#'### SUMMARY season dates:
print(as.data.frame(seas_dates))

#' Cessation of rainy season is vairly unifom, and averages to 124
mean(seas_dates$day_cess)

#' Onset is more variable, and averages to 24 where long rainy season is distinguishable
mean(seas_dates$day_on[seas_dates$decade != "2001-2010"])

#'## Determine experimental values for rainfall and dry days
  
#' Rainy season days for our study area are ~24-124. Subset the main rainfall data frame accordingly
df_rain_long <- df_chirps %>% 
  group_by(Year) %>% 
  dplyr::filter(Day %in% c(24:124)) %>% 
  #group_by(Year) %>% 
  # # Remove some weekly rainfall outliers before later calculating the average rainfall 
  # # (Up for debate)
  #filter(is_outlier(prec_day) == FALSE) %>%      
  ungroup() %>% 
  group_by(Year, Week) %>% 
  # Divide total rainfall by number of weeks in a given wet season in a given year
  summarize(prec_week = sum(prec_day)) 

# Look at the distribution of weekly rainfall for every year

ggplot(df_rain_long[df_rain_long$Year %in% c(1991, 1992, 1993, 1994),]) +
  geom_density(aes(prec_week)) +
  geom_vline(data = plyr::ddply(df_rain_long[df_rain_long$Year %in% c(1991, 1992, 1993, 1994),], "Year", summarize, x = median(prec_week)), aes(xintercept=x, col = "red"), lty = 2) +
  geom_vline(data = plyr::ddply(df_rain_long[df_rain_long$Year %in% c(1991, 1992, 1993, 1994),], "Year", summarize, x = mean(prec_week)), aes(xintercept=x, col = "blue"), lty = 2) +
  facet_wrap(facets = "Year", scales = "free")  +
  scale_colour_manual(name = '', 
                      values =c('red'='red','blue'='blue'), labels = c('Mean','Median'), guide = "legend")



# As expected, the rainfall data is quite skewed. Rather than average, we should 
# use the median here to extract a single weekly rainfall value. However, using the 
# median turns rainfall values to zero for some, heavily skewed years (i.e. 1992).
# For now, go with the mean.
df_rain_long <- df_rain_long %>% 
  ungroup() %>% 
  group_by(Year) %>% 
  summarize(prec_week_long = mean(prec_week))

ggplot(df_rain_long) +
  geom_density(aes(prec_week_long)) 

#' Density plot of rainfall, with 5% and 95% quantile, as well as mean and median highlighted
jpeg("figures/plot2.jpg", width = 500, height = 500)
p2 <- ggplot() +
  geom_density(data = df_rain_long, aes(prec_week_long)) +
  geom_vline(xintercept = quantile(df_rain_long$prec_week_long, 0.05), 
             col = "red", alpha = .5, lty = 2, lwd = 1) +
  geom_vline(xintercept = quantile(df_rain_long$prec_week_long, 0.95), 
             col = "red", alpha = .5, lty = 2, lwd = 1) +
  geom_vline(xintercept = median(df_rain_long$prec_week_long), 
             col = "blue", alpha = .5, lty = 2, lwd = 1) +
  xlab("Weekly wet season rainfall per year (mm)"); p2
dev.off()
# Average and median are both quite similar
  

jpeg("figures/plot_comb.jpg", width = 1000, height = 500, quality = 100)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()


  
#' Calculate 5% and 95% quantile as well as median for rainfall
quantile(df_rain_long$prec_week_long, 0.05)
median(df_rain_long$prec_week_long)
quantile(df_rain_long$prec_week_long, 0.95)

#' Get length of consecutive dry days, as well as overall number of dry days, for every year
chirps_dry <- df_chirps %>% 
  group_by(Year) %>% 
  dplyr::filter(Day %in% c(24:124)) %>% 
  arrange(date) %>% 
  group_by(Year) %>% 
  do(dry_duration = rle(.$prec_day)$lengths[which(rle(.$prec_day)$values == 0)]) %>% 
  unnest() %>% 
  group_by(Year) %>% 
  mutate(total_dry = sum(dry_duration)
         ,outlier_dur = is_outlier(dry_duration)
         ) %>%
  dplyr::filter(outlier_dur == "FALSE")

#' Density plot for number of consecutive dry days, with 5% and 95% quantile highlighted
ggplot() +
  geom_density(data = chirps_dry, aes(dry_duration)) +
  geom_vline(xintercept = quantile(chirps_dry$dry_duration, 0.05), 
             col = "red", alpha = .5, lty = 2) +
  geom_vline(xintercept = quantile(chirps_dry$dry_duration, 0.95), 
             col = "red", alpha = .5, lty = 2) +
  xlab("Number of consecutive dry days during wet season") 

#' Calculate 5% and 95% quantile as well as median for number of consecutive dry days
quantile(chirps_dry$dry_duration, 0.05)
median(chirps_dry$dry_duration)
quantile(chirps_dry$dry_duration, 0.95)


#'#### SUMMARY: Our experimental values are:
median(chirps_dry$dry_duration)
quantile(chirps_dry$dry_duration, 0.95)

#' 39% decrease, 49% increase compared to control
1-quantile(df_rain_long$prec_week_long, 0.05)/median(df_rain_long$prec_week_long)
1-quantile(df_rain_long$prec_week_long, 0.95)/median(df_rain_long$prec_week_long)


#' What if we estimate weekly rainfall using proportions of long rain season rain?

#' Calculate the proportional rainfall per year
seas_rain <- df_chirps %>% 
  group_by(Year) %>% 
  summarize(wet_rain = sum(prec_day[Day %in% c(24:124)])/sum(prec_day))

median(seas_rain$wet_rain)
# 60% of rain falls in the long rain season

#' Calculate the weekly rainfall, using the proportional rainfall 
prop_rain_long <- df_chirps %>% 
  group_by(Year) %>% 
  summarize(year_rain = sum(prec_day),
            long_rain = 0.6*year_rain,
            nweek_longrain = length(unique(Week[Day %in% c(24:124)])),
            long_rain_week = long_rain/nweek_longrain)

ggplot(prop_rain_long) +
  geom_density(aes(long_rain_week)) +
  geom_vline(xintercept = quantile(prop_rain_long$long_rain_week, 0.05), 
             col = "red", alpha = .5, lty = 2) +
  geom_vline(xintercept = quantile(prop_rain_long$long_rain_week, 0.95), 
             col = "red", alpha = .5, lty = 2) +
  geom_vline(xintercept = mean(prop_rain_long$long_rain_week), 
             col = "blue", alpha = .5, lty = 2) +
  geom_vline(xintercept = median(prop_rain_long$long_rain_week), 
             col = "orange", alpha = .5, lty = 2) +
  xlab("Weekly wet season rainfall per year (mm)") 


#' 29% decrease, 42% increase compared to control
1-quantile(prop_rain_long$long_rain_week, 0.05)/median(prop_rain_long$long_rain_week)
median(prop_rain_long$long_rain_week)
1-quantile(prop_rain_long$long_rain_week, 0.95)/median(prop_rain_long$long_rain_week)

