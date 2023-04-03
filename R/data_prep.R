library(tidyverse)
df <- read.csv(file = 'data/solar_analysis_data.csv', header = TRUE, stringsAsFactors = FALSE)

#select the data columns we need

# GAP codes are NA outside of protected areas. set these to 5
df$GAP_Sts[is.na(df$GAP_Sts)] <- 5

# Get the state from the county census tract FIPS
df$fips <- as.character(df$FIPS)

# filter records to those in DC, DE, MD, NY, PA, VA
temp <- filter(df, state %in% c('Delaware', 'District of Columbia', 'Maryland', 'New York', 'Pennsylvania', 'Virginia'))%>%
  # filter (presumably) random polygons that got drawn in water
  filter(!is.na(FIPS))%>%
  filter(!is.na(impervious16))%>%
  filter(!is.na(cultivated16))%>%
  mutate(solar = as.numeric(year < 2022),
         statef = as.factor(state),
         statei = as.numeric(statef),
         ttd = ifelse(solar, year - 2016, NA),
         parcel_area = min(Shape_Area_x, Shape_Area_y, na.rm = TRUE),
         index_parcel = min(index_parcel_x, index_parcel_y, na.rm = TRUE),
         county = min(county_x, county_y, na.rm = TRUE),
         d = as.numeric(solar == 0)
  )%>%
  group_by(state)%>%
  mutate(tmax = max(ttd, na.rm = TRUE),
         l = min(ttd, na.rm = TRUE))%>%
  ungroup()%>%
  # as a final step, after we've calculated ttd and l with per-state first years as reference
  # set any panels with >10% impervious surface as panels that existed in 2016
  mutate(
    year = ifelse(impervious16 >= 10 & solar, 2016, year)
  )%>%
  select(year,
         ttd,
         tmax,
         solar,
         area,
         ssurgo,
         cultivated,
         cultivated16,
         impervious,
         impervious16,
         open,
         open16,
         tree_cover,
         tree_cover16,
         slope,
         lat,
         parcel_area,
         housing,
         income,
         pdensity,
         GAP_Sts,
         road_dist,
         line_dist,
         statei,
         statef,
         l,
         d)

# fill in missing income and ssurgo with means from solar and non solar
temp$ssurgo[is.na(temp$ssurgo) & temp$solar == 1] <- mean(temp$ssurgo[temp$solar == 1], na.rm = TRUE)
temp$ssurgo[is.na(temp$ssurgo) & temp$solar == 0] <- mean(temp$ssurgo[temp$solar == 0], na.rm = TRUE)
temp$income[is.na(temp$income) & temp$solar == 1] <- mean(temp$income[temp$solar == 1], na.rm = TRUE)
temp$income[is.na(temp$income) & temp$solar == 0] <- mean(temp$income[temp$solar == 0], na.rm = TRUE)

saveRDS(as.data.frame(temp), file = 'data/cleaned_solar_analysis_data.rds')


nlcdClasses <- c(
  'Water' = 11,
  'Dev., Open' = 21,
  'Dev., Low' = 22,
  'Dev., Med' = 23,
  'Dev., High' = 24,
  'Barren' = 31,
  'Forest, Decid.' = 41,
  'Forest, Conif.' = 42,
  'Fores, Mix' = 43,
  'Shrub' = 52,
  'Grass' = 71,
  'Pasture' = 81,
  'Crops' = 82,
  'Wetland, Herb' = 90,
  'Wetland, Wood' = 95
)

# Create Data for interactive web docs (markdown, shiny, etc.)
timeDF <- read.csv(file = 'data/data_raw/timeline.csv', header = TRUE)%>%
  pivot_longer(c(X2017, X2018, X2019, X2020, X2021),
               names_to = 'year',
               names_prefix = 'X',
               values_to = 'area')%>%
  mutate(std_area = area/Area)

saveRDS(timeDF, file = 'www/data/timelinedat.rds')

solarLCDF <- read.csv(file = 'data/data_raw/solar_lc_hist.csv', header = TRUE)%>%
  filter(state != 'DC')%>%
  mutate(area = 'Solar arrays',
         norm = count/sum(count, na.rm = TRUE))%>%
  group_by(state)%>%
  mutate(stnorm = count/sum(count, na.rm = TRUE))%>%
  ungroup()

# studyAreaLcDF <- read.csv(file = 'data/data_raw/lc_hist_studyarea.csv', header= TRUE)%>%
studyAreaLCDF <- read.csv(file = 'data/data_raw/lc_hist_padus.csv', header = TRUE)
colnames(studyAreaLCDF) <- c('state', 'landcover', 'gap', 'count', 'gapless')
studyAreaLCDF <- filter(studyAreaLCDF, state!='DC')%>%
  mutate(area = 'CPK study area',
         norm = count/sum(count, na.rm = TRUE),
         gapless_norm = gapless/sum(gapless, na.rm = TRUE))%>%
  group_by(state)%>%
  mutate(stnorm = count/sum(count, na.rm = TRUE),
         gapless_stnorm = gapless/sum(gapless), na.rm = TRUE)%>%
  ungroup()

state_lcDF <- bind_rows(solarLCDF, studyAreaLCDF)
state_lcDF$class <- vapply(state_lcDF$landcover,
                           function(x){
                             if(is.na(x)){
                               ''
                             }
                             else{
                               names(nlcdClasses)[grepl(x, nlcdClasses)]
                             }
                             },
                           FUN.VALUE = character(1))

saveRDS(state_lcDF, file = 'www/data/stateLandcover.rds')

lcDF <- group_by(state_lcDF, landcover, area)%>%
  summarize(
    count =sum(count, na.rm = TRUE),
    norm = sum(norm, na.rm = TRUE),
    gapless = sum(gapless, na.rm = TRUE),
    gapless_norm = sum(gapless_norm, na.rm = TRUE),
    class = first(class))

saveRDS(lcDF, file = 'www/data/landcover.rds')
