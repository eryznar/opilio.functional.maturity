library(tidyverse)
library(odbc)
library(DBI)
library(getPass)
library(keyring)
library(lifecycle)
library(data.table)
library(crabpack)
library(INLA)
library(sdmTMB)
library(glmmTMB)
library(broom)
library(sf)
library(gstat)
library(rnaturalearth)
library(raster)
library(concaveman)
library(png)
library(mgcv)
library(dplyr)
library(mgcv)
library(ggplot2)
library("rnaturalearth")
library(patchwork)
library(gratia)
library(MuMIn)
library(DHARMa)
library(mgcViz)


# Read in survey data
specEBS <- readRDS("./Data/snow_survey_specimenEBS.rda")
matEBS <- readRDS("./Data/snow_survey_maturityEBS.rda")

# Specify directory
dir <- "Y:/KOD_Research/Ryznar/Crab functional maturity"

# Read in prediction grid, convert to UTM
pred.grid <- readRDS(paste0(dir, "/Data/EBS_opilio_grid_5km_No_Land.rds")) %>%
              st_as_sf(., coords = c("Lon", "Lat"), crs = "+proj=longlat +datum=WGS84") %>%
              st_transform(., crs = "+proj=utm +zone=2") %>%
              cbind(st_coordinates(.)) %>%
              as.data.frame(.) %>%
              dplyr::select(Area_km2, X, Y) %>%
              #replicate_df(., "year", years) %>%
              mutate(X = X/1000, Y = Y/1000) %>%
              rename(LONGITUDE = X, LATITUDE = Y)

# Specify cutline coefficients
minima <- read.csv("./Output/opilio_cutline_minima.csv") %>%
  mutate(BETA0 = coef(lm(minima ~ midpoint))[1],
         BETA1 = coef(lm(minima ~ midpoint))[2])

BETA0 <- unique(minima$BETA0)
BETA1 <- unique(minima$BETA1)

# Functions
diagnose <- function(model){
  model.name <- deparse(substitute(model))
  
  ss <- summary(model) # model summary
  
  gam.check(model) # make sure smooth terms are appropriate
  
  plot(simulateResiduals(model)) # Checks uniformity, dispersion, outliers via DHARMa

    # plot facetted smooths
  sm.dat <- smooth_estimates(model) %>%
    pivot_longer(., cols = 6:9, names_to = "resp", values_to = "value")
  
  ggplot(sm.dat, aes(x = value, y = .estimate)) +
    geom_ribbon(sm.dat, mapping = aes(ymin = .estimate + 2 * .se, ymax = .estimate - 2 * .se), fill = "cadetblue", alpha = 0.25)+
    geom_line(color = "cadetblue", linewidth = 1.25) +
    facet_wrap(~ .smooth, scales = "free_x") +   # facet by smooth term name
    theme_bw()+
    ggtitle(model.name)+
    ylab("Partial effect")+
    xlab("Value")+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14)) -> sm.plot
  
  return(list(ss, sm.plot))
}

