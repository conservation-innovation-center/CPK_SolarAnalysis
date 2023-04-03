# IMPORT LIBRARIES
library(plotly)
library(tidyverse)

# LOAD NECESSARY DATA
lcDF <- readRDS(file = 'www/data/landcover.rds')
stateLCDF <- readRDS(file = 'www/data/stateLandcover.rds')
timeDF <- readRDS(file = 'www/data/timelinedat.rds')
cleanDF <- readRDS(file = 'data/cleaned_solar_analysis_data.rds')
# transform predictor variables based on histograms
cleanDF$slope <- sqrt(cleanDF$slope)
cleanDF$line_dist <- log(cleanDF$line_dist+0.1)
cleanDF$road_dist <- log(cleanDF$road_dist+0.1)
cleanDF$GAP_Sts <- log(6-cleanDF$GAP_Sts)
cleanDF$pdensity <- log(cleanDF$pdensity + 0.1)

# DEFINE STANDARD AXIS FORMATTING
tickfont <- list(
  color = 'black',
  size = 16,
  family = 'Helvetica'
)

ax <- list(
  tickfont = tickfont,
  titlefont = list(
    color = 'black',
    size = 18,
    family = 'Helvetica'
  ),
  showgrid = FALSE,
  showline = FALSE,
  zeroline = TRUE
)



# FIGURE 1. TIMELINE
plot_ly(data = timeDF, type = 'scatter', mode = 'lines+markers')%>%
  add_lines(
    x = ~ as.integer(year),
    y = ~ area/100000,
    symbol = ~State,
    # colors = pal,
    line = list(
      color = 'black'
    ),
    marker = list(
      color = 'black',
      size = 12,
      symbol = rep(c('circle', 'x', 'diamond', 'triangle-up', 'square'), each = 5)
    )
  )%>%
  layout(
    xaxis = list(
      title = 'Year'
    )%>%append(ax),
    yaxis = list(
      title = 'Area (km<sup>2</sup>)'
    )%>%append(ax),
    legend = list(x = 0.1, y = 1, font = tickfont)
  )

plot_ly(data = timeDF, type = 'scatter', mode = 'lines+markers')%>%
  add_lines(
    x = ~ as.integer(year),
    y = ~ std_area,
    symbol = ~State,
    # colors = pal,
    line = list(
      color = 'black'
    ),
    marker = list(
      color = 'black',
      size = 12,
      symbol = rep(c('circle', 'x', 'diamond', 'triangle-up', 'square'), each = 5)
    )
  )%>%
  layout(
    xaxis = list(
      title = 'Year'
    )%>%append(ax),
    yaxis = list(
      title = 'Proportion of state area',
      showexponent = "all",
      exponentformat = "e"
    )%>%append(ax),
    legend = list(x = 0.1, y = 1, font = tickfont)
  )

# FIGURE 2 BAR PLOTS
# a) Whole Area
plot_ly(type = 'bar')%>%
  add_trace(
    data = lcDF[lcDF$area == 'CPK study area',],
    name = 'Study area',
    x = ~class,
    y = ~round(norm,4)*100,
    marker = list(
      color = 'white',
      line = list(
        color = 'black',
        width = 1
      )
    )
    )%>%
  add_trace(
    data = lcDF[lcDF$area == 'Solar arrays',],
    name = 'Solar arrays',
    x = ~class,
    y = ~round(norm,4)*100,
    marker = list(
      color = 'grey'
    )
  )%>%
  layout(
    title = list(
      text = 'a) Chesapeake study area',
      font = list(
        family = 'serif',
        color = 'black',
        size = 14
        ),
      x = 0.5,
      y = 0.95
    ),
    xaxis = list(
      title = 'NLCD class',
      tickangle = 45
    )%>%append(ax),
    yaxis = list(
      title = 'Percent landcover',
      range = c(0, 60)
    )%>%append(ax),
    legend = list(
      x = 0.7,
      y = 0.9,
      font = list(
        family = 'serif',
        size = 14,
        color = 'black'
      )
    ),
    margin = list(r = 0)
  )


states <- c('DE', 'MD', 'PA', 'NY', 'VA')
state_graphs <- lapply(states, function(st){
  p <- plot_ly(type = 'bar')%>%
    add_trace(
      data = stateLCDF[stateLCDF$area == 'CPK study area' & stateLCDF$state == st,],
      x = ~class,
      y = ~stnorm*100,
      name = 'Solar arrays',
      marker = list(
        color = 'white',
        line = list(
          color = 'black',
          width = 1
        )
      )
    )%>%
    add_trace(
      data = stateLCDF[stateLCDF$area == 'Solar arrays' & stateLCDF$state == st,],
      x = ~class,
      y = ~stnorm*100,
      name = 'CPK study area',
      marker = list(
        color = 'grey'
      )
    )%>%
    layout(
      title = list(
        text = paste(st)
      ),
      barmode = 'group',
      yaxis = list(
        range = c(0,60),
        title = 'Percent landcover'
      )%>%append(ax),
      xaxis = list(
        title = 'NLCD class'
      )%>%append(ax),
      showlegend = FALSE
    )
  return(p)
})

## VIOLIN PLOTS
outBinWeib <- readRDS(file = 'data/output/binWeib/fulllog_truncated_censored_both.rds')
gvsBinWeib <- readRDS(file = 'data/output/binweib/fulllog_truncated_censored_both_GVS.rds')
full <- c('impervious16',
          'open16',
          'tree_cover16',
          'cultivated16',
          'ssurgo',
          'slope', #square root(x)
          'GAP_Sts', #log(6-x)
          'line_dist', #square root(x)
          'road_dist', #square root(x)
          'pdensity',
          'income',
          'lat')

names(full) <- full

violinDF <- map_dfc(full, function(x){
  var <- paste('beta', x, sep = '_')
  outBinWeib$sims.list[[var]][13000:15000]
  })%>%
  pivot_longer(cols = everything(),
               names_to = 'variable',
               values_to = 'estimate')

violinDF <- lapply(full, function(x){
  var <- paste('beta', x, sep = '_')
  ind <- paste('w', x, sep = "_")
  vals <- gvsBinWeib$sims.list[[var]][gvsBinWeib$sims.list[[ind]] == 1]
  data.frame(variable = rep(var, length(vals)), estimate = vals, mn = rep(mean(vals, na.rm = TRUE), length(vals)))
})%>%
    bind_rows()%>%
    arrange(desc(mn))

plot_ly(
  data = violinDF,
  y = ~ estimate,
  x = ~ variable,
  type = 'violin',
  orientation = 'v',
  meanline = list(
    visible = TRUE
  ),
  line = list(
    color = 'black'
  ),
  points = FALSE,
  fillcolor = 'white'
)%>%
  layout(
    showlegend = FALSE,
    yaxis = list(
      title = 'Coefficient estimate'
    )%>%append(ax),
    xaxis = list(
      title = 'Covariate',
      tickmode = 'array',
      tickvals = c('beta_ssurgo', 'beta_cultivated16', 'beta_road_dist', 'beta_pdensity', 'beta_income', 'beta_open16', 'beta_tree_cover16', 'beta_GAP_Sts', 'beta_lat', 'beta_slope'),
      ticktext = c('Farm Score', 'Cultivated', 'Road Distance', 'Population', 'Income', 'Open', 'Tree Cover', 'GAP Status','Latitude', 'Slope')
    )%>%append(ax)
  )

alphas <- c('DE' = 'alpha[1]', 'DC' = 'alpha[2]', 'MD' = 'alpha[3]', 'NY' = 'alpha[4]', 'PA' = 'alpha[5]', 'VA' = 'alpha[6]')

alphaDF <- as.data.frame(outBinWeib$sims.list$alpha[13000:15000,])
colnames(alphaDF) <- c('DE', 'DC', 'MD', 'NY', 'PA', 'VA')
alphaDF <- pivot_longer(data = alphaDF,
             cols = everything(),
             names_to = 'variable',
             values_to = 'estimate')

plot_ly(
  data = filter(alphaDF),
  x = ~ estimate,
  y = ~ variable,
  type = 'violin',
  orientation = 'h',
  meanline = list(
    visible = TRUE
  ),
  # box = list(
  #   visible = TRUE
  # ),
  line = list(
    color = 'black'
  ),
  points = FALSE,
  fillcolor = 'white'
)%>%
  layout(
    showlegend = FALSE,
    xaxis = list(
      title = 'Coefficient estimate'
    )%>%append(ax),
    yaxis = list(
      title = 'Covariate'
    )%>%append(ax)
  )

salphas <- c('DE' = 'alpha[1]', 'DC' = 'alpha[2]', 'MD' = 'alpha[3]', 'NY' = 'alpha[4]', 'PA' = 'alpha[5]', 'VA' = 'alpha[6]')

shapeDF <- as.data.frame(gvsBinWeib$sims.list$shape[13000:15000,])
colnames(shapeDF) <- c('DE', 'DC', 'MD', 'NY', 'PA', 'VA')
shapeDF <- pivot_longer(data = shapeDF,
                        cols = everything(),
                        names_to = 'variable',
                        values_to = 'estimate')

plot_ly(
  data = filter(shapeDF, variable != 'DC'),
  y = ~ estimate,
  x = ~ variable,
  type = 'violin',
  orientation = 'v',
  meanline = list(
    visible = TRUE
  ),
  # box = list(
  #   visible = TRUE
  # ),
  line = list(
    color = 'black'
  ),
  points = FALSE,
  fillcolor = 'white'
)%>%
  add_lines(y = c(1,1), x = c('DE', 'VA'), line = list(dash = 'dash'))%>%
  layout(
    showlegend = FALSE,
    yaxis = list(
      title = 'Shape parameter estimate',
      range = c(0,4.5)
    )%>%append(ax),
    xaxis = list(
      title = 'State'
    )%>%append(ax)
  )

# RELATIONSHIPS

cleanDF <- readRDS(file = 'data/cleaned_solar_analysis_data.rds')
head(cleanDF)

# transform predictor variables based on histograms
cleanDF$slope <- sqrt(cleanDF$slope)
cleanDF$line_dist <- log(cleanDF$line_dist+0.1)
cleanDF$road_dist <- log(cleanDF$road_dist+0.1)
cleanDF$GAP_Sts <- log(6-cleanDF$GAP_Sts)
cleanDF$pdensity <- log(cleanDF$pdensity + 0.1)

var <- 'slope'
rg <- range(cleanDF[,var], na.rm = TRUE)
mn <- gvsProbs[[var]][2]
sd <- gvsProbs[[var]][3]

expected <- curve(exp(log(gvsBinWeib$mean$alpha_mu)+ mn*sqrt(x)), from = rg[1], to = rg[2])
upper <- curve(exp(log(gvsBinWeib$mean$alpha_mu+gvsBinWeib$sd$alpha_mu)+ (mn+sd)*sqrt(x)), from = rg[1], to = rg[2], add = TRUE, lty = 2)
lower <- curve(exp(log(gvsBinWeib$mean$alpha_mu-gvsBinWeib$sd$alpha_mu)+ (mn-sd)*sqrt(x)), from = rg[1], to = rg[2], add = TRUE, lty = 2)

plot_ly(type = 'scatter', mode = 'lines')%>%
  add_lines(
    x =~ expected$x,
    y = ~expected$y,
    line = list(
      color = 'black',
      width = 2
    ),
    name = 'Mean'
  )%>%
  add_lines(
    x =~ upper$x,
    y = ~upper$y,
    line = list(
      color = 'black',
      width = 2,
      dash = 'dash'
    ),
    name = '95% CI'
  )%>%
  add_lines(
    x =~ lower$x,
    y = ~lower$y,
    line = list(
      color = 'black',
      width = 2,
      dash = 'dash'
    ),
    name = '95% CI',
    showlegend = FALSE
  )%>%
  layout(
    xaxis = list(
      title = 'GAP Status'
    )%>%append(ax),
    yaxis = list(
      title = 'Lambdai'
    )%>%append(ax)
  )

## EFFECT HEATMAP
states <- c('DE', 'MD', 'PA', 'NY', 'VA')
classes <- sort(unique(state_lcDF$class), decreasing = TRUE)
### With GAP Included
state_chi_posthoc <- map(states, function(st){
  solarFreq <- filter(stateLCDF, area == 'Solar arrays', state == st, !is.na(landcover))%>%
    select(class, count)
  areaFreq <- filter(stateLCDF, area == 'CPK study area', state == st, !is.na(landcover))%>%
    select(class, stnorm)
  test_dat <- left_join(areaFreq, solarFreq, by = 'class')%>%
    arrange(desc(class))
  test_dat$count[is.na(test_dat$count)] <- 0
  lcChi <- chisq.test(test_dat$count, p = test_dat$stnorm)
  # adjustedResid <- lcChi$residuals/sd(lcChi$residuals)
  
  lcSim <- rmultinom(1000, sum(test_dat$count), p = test_dat$stnorm)
  diff <- test_dat$count - lcSim
  effect <- apply(diff, 1, mean)*(30^2)/100000
  p <- apply(diff < 0, 1, sum)/1000
  return(list(x2 = lcChi$statistic[[1]], df = lcChi$parameter[[1]], effects = effect, p = p))
})

effect_matrix <- matrix(data = NA, nrow = 5, ncol = 15)
for(i in 1:length(states)){
  print(length(state_chi_posthoc[[i]]$effects))
  effect_matrix[i,] <- state_chi_posthoc[[i]]$effects
}

solarFreq <- filter(lcDF, area == 'Solar arrays', !is.na(landcover))%>%
  arrange(desc(class))
areaFreq <- filter(lcDF, area == 'CPK study area', !is.na(landcover))%>%
  arrange(desc(class))%>%
  select(norm)
lcChi <- chisq.test(solarFreq$count, p = areaFreq$norm)
# adjustedResid <- lcChi$residuals/sd(lcChi$residuals)

lcSim <- rmultinom(1000, sum(solarFreq$count), p = areaFreq$norm)
diff <- solarFreq$count - lcSim
effect <- apply(diff, 1, mean)*(30^2)/100000
p <- apply(diff < 0, 1, sum)/1000

mat <- rbind(effect_matrix, effect[1:15])
xlab <- append(states, 'CBW')
plot_ly(z = t(mat),
        type = 'heatmap',
        x = xlab,
        y = classes[1:15],
        zmin = -20, zmax = 20,
        colors = colorRamp(c('black', 'white')),
        colorbar = list(
          title = list(
            text = 'Effect size',
            font = tickfont),
          tickfont = tickfont
          # orientation = "h",
          # x = 0.5,
          # xanchor = 'center',
          # y = 1,
          # yanchor = 'bottom'
        ))%>%
  layout(
    xaxis = list(
      title = 'State'
    )%>%append(ax),
    yaxis = list(
      title = 'NLCD class',
      # categoryorder = 'array',
      categoryarray = c(
        'Water',
        'Wetland, Wood',
        'Wetland, Herb',
        'Forest, Decid.',
        'Fores, Mix',
        'Shrub',
        'Grass',
        'Forest, Conif.',
        'Barren',
        'Pasture',
        'Crops',
        'Dev., Open',
        'Dev., Low',
        'Dev., Med',
        'Dev., High')
    )%>%append(ax)
  )%>%
  add_annotations(
    x = ~rep(xlab, each = 15),
    y = classes[1:15],
    text = ~as.character(t(round(mat,2))),
    font = list(
      color = 'white',
      family = 'serif',
      size = 12
    ),
    xref = 'x',
    yref = 'y',
    showarrow = FALSE)

### With GAP 1&2 Excluded
state_chi_posthoc <- map(states, function(st){
  solarFreq <- filter(state_lcDF, area == 'Solar arrays', state == st, !is.na(landcover))%>%
    select(class, count)
  areaFreq <- filter(state_lcDF, area == 'CPK study area', state == st, !is.na(landcover))%>%
    select(class, gapless_stnorm)
  test_dat <- left_join(areaFreq, solarFreq, by = 'class')%>%
    arrange(desc(class))
  test_dat$count[is.na(test_dat$count)] <- 0
  lcChi <- chisq.test(test_dat$count, p = test_dat$gapless_stnorm)
  # adjustedResid <- lcChi$residuals/sd(lcChi$residuals)
  
  lcSim <- rmultinom(1000, sum(test_dat$count), p = test_dat$gapless_stnorm)
  diff <- test_dat$count - lcSim
  effect <- apply(diff, 1, mean)*(30^2)/100000
  p <- apply(diff < 0, 1, sum)/1000
  return(list(x2 = lcChi$statistic[[1]], df = lcChi$parameter[[1]], effects = effect, p = p))
})

effect_matrix <- matrix(data = NA, nrow = 5, ncol = 15)
for(i in 1:length(states)){
  print(length(state_chi_posthoc[[i]]$effects))
  effect_matrix[i,] <- state_chi_posthoc[[i]]$effects
}

solarFreq <- filter(lcDF, area == 'Solar arrays', !is.na(landcover))%>%
  arrange(desc(class))
areaFreq <- filter(lcDF, area == 'CPK study area', !is.na(landcover))%>%
  arrange(desc(class))

lcChi <- chisq.test(solarFreq$count, p = areaFreq$gapless_norm)
# adjustedResid <- lcChi$residuals/sd(lcChi$residuals)

lcSim <- rmultinom(1000, sum(solarFreq$count), p = areaFreq$gapless_norm)
diff <- solarFreq$count - lcSim
effect <- apply(diff, 1, mean)*(30^2)/100000
p <- apply(diff < 0, 1, sum)/1000

mat <- rbind(effect_matrix, effect[1:15])
xlab <- append(states, 'CBW')
plot_ly(z = t(mat),
        type = 'heatmap',
        x = xlab,
        y = classes[1:15],
        zmin = -20, zmax = 20,
        colors = colorRamp(c('black', 'white')),
        colorbar = list(
          title = list(
            text = 'Effect size (km<sup>2</sup>)',
            font = tickfont),
          tickfont = tickfont
          # orientation = "h",
          # x = 0.5,
          # xanchor = 'center',
          # y = 1,
          # yanchor = 'bottom'
        ))%>%
  layout(
    xaxis = list(
      title = 'State'
    )%>%append(ax),
    yaxis = list(
      title = 'NLCD class',
      # categoryorder = 'array',
      categoryarray = c(
        'Water',
        'Wetland, Wood',
        'Wetland, Herb',
        'Forest, Decid.',
        'Fores, Mix',
        'Shrub',
        'Forest, Conif.',
        'Grass',
        'Pasture',
        'Crops',
        'Barren',
        'Dev., Open',
        'Dev., Low',
        'Dev., Med',
        'Dev., High')
    )%>%append(ax)
  )%>%
  add_annotations(
    x = ~rep(xlab, each = 15),
    y = classes[1:15],
    text = ~as.character(t(round(mat,2))),
    font = list(
      color = 'white',
      family = 'serif',
      size = 12
    ),
    xref = 'x',
    yref = 'y',
    showarrow = FALSE)

natural <- c(
  'Water',
  'Wetland, Wood',
  'Wetland, Herb',
  'Forest, Decid.',
  'Fores, Mix',
  'Shrub')

cultivated <- c(
  'Forest, Conif.',
  'Grass',
  'Pasture',
  'Crops'
)

modified <- c(
  'Barren',
  'Dev., Open',
  'Dev., Low',
  'Dev., Med',
  'Dev., High'
)

plot_ly(type = 'bar')%>%
  # add_trace(data = filter(state_lcDF, area == 'CPK study area'), x = ~state, y = ~gapless_stnorm, color = ~class)%>%
  # add_trace(data = filter(state_lcDF, area == 'Solar arrays')%>%
  #             mutate(category = ifelse(class %in% natural, 'Natural', ifelse(class %in% cultivated, 'Cultivated', ifelse(class%in%modified, 'Modified', 'Other')))),
  #           x = ~state, y = ~stnorm, color = ~category)%>%
  add_trace(data = filter(state_lcDF, area == 'CPK study area')%>%
              mutate(category = ifelse(class %in% natural, 'Natural', ifelse(class %in% cultivated, 'Cultivated', ifelse(class%in%modified, 'Modified', 'Other')))),
            x = ~state, y = ~stnorm, color = ~category)%>%
  layout(barmode = 'stack')

# BIODIVERSITY
biod <- read.csv(file = 'data/data_raw/biodiversity.csv', header = TRUE)
head(biod)
biod_long <- mutate(biod, all = birds + amphibians + reptile + mammals)%>%
  pivot_longer(cols = c(all, birds, amphibians, reptile, mammals), names_to = 'taxon', values_to = 'richness')%>%
  mutate(richness = round(ifelse(is.na(richness),0,richness)), solar = ifelse(year == 2030, 0, 1))%>%
  filter(state %in% c('Delaware', 'Maryland', 'New York', 'Pennsylvania', 'Virginia'))
head(biod_long)
plot_ly(data = filter(biod_long, taxon == 'all'), type = 'histogram', x = ~richness, color = ~state, histnorm = 'percent')
plot_ly(data = biod_long, type = 'box', x = ~state, y = ~log(richness+0.01), color = ~as.factor(solar))%>%
  layout(boxmode = 'group')

biod_state <- read.csv(file = 'data/data_raw/stateBiodiversity.csv', header = TRUE)%>%
  pivot_longer(cols = paste('X', 0:134, sep = ""), names_to = 'richness', values_to = 'freq')%>%
  mutate(freq = ifelse(is.na(freq), 0, freq), richness = as.integer(str_remove(richness, 'X')))%>%
  group_by(state, richness)%>%
  summarize(total_freq = sum(freq), .groups = 'drop_last')%>%
  mutate(total_p = total_freq/sum(total_freq))

head(biod_state)
plot_ly(biod_state, type ='bar', x = ~richness, y = ~total_p, color = ~state) 
