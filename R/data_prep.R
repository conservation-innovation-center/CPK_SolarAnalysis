library(tidyverse)
df <- read.csv(file = 'data/solar_analysis_data.csv', header = TRUE, stringsAsFactors = FALSE)

#select the data columns we need

# GAP codes are NA outside of protected areas. set these to 5
df$GAP_Sts[is.na(df$GAP_Sts)] <- 5

# Get the state from the county census tract FIPS
df$fips <- as.character(df$FIPS)

# filter records to those in DC, DE, MD, NY, PA, VA
cleanDF <- filter(df, state %in% c('Delaware', 'District of Columbia', 'Maryland', 'New York', 'Pennsylvania', 'Virginia'))%>%
  # filter (presumably) random polygons that got drawn in water
  filter(!is.na(FIPS))%>%
  filter(!is.na(impervious16))%>%
  filter(!is.na(cultivated16))%>%
  mutate(solar = as.numeric(year < 2022),
         statef = as.factor(state),
         statei = as.numeric(statef),
         ttd = ifelse(solar, year - 2016, NA),
         parel_area = min(Shape_Area_x, Shape_Area_y, na.rm = TRUE),
         index_parcel = min(index_parcel_x, index_parcel_y, na.rm = TRUE),
         county = min(county_x, county_y, na.rm = TRUE),
         d = as.numeric(solar == 0)
  )%>%
  group_by(state)%>%
  mutate(tmax = 2022 - min(year))%>%
  ungroup()%>%
  select(ttd,
         tmax,
         solar,
         area,
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
         parel_area,
         housing,
         income,
         POPULATION,
         GAP_Sts,
         road_dist,
         line_dist,
         statei,
         statef,
         d)

saveRDS(as.data.frame(cleanDF), file = 'data/cleaned_solar_analysis_data.rds')
