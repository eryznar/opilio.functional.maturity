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

