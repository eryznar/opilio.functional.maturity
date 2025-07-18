library(tidyverse)
library(odbc)
library(DBI)
library(getPass)
library(keyring)
library(lifecycle)
library(data.table)
library(crabpack)
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

# Set years
years <- c(1989:2007, 2009:2013, 2015, 2017:2019, 2021:2024)

# Read in survey data
specEBS <- readRDS("./Data/snow_survey_specimenEBS.rda")
matEBS <- readRDS("./Data/snow_survey_maturityEBS.rda")

# Specify directory
dir <- "Y:/KOD_Research/Ryznar/Crab functional maturity"

data_dir <- "Y:/KOD_Survey/EBS Shelf/Data_Processing/Data/" # for survey data

# Specify cutline coefficients
minima <- read.csv("./Output/opilio_cutline_minima.csv") %>%
  mutate(BETA0 = coef(lm(minima ~ midpoint))[1],
         BETA1 = coef(lm(minima ~ midpoint))[2])

BETA0 <- unique(minima$BETA0)
BETA1 <- unique(minima$BETA1)

# Arithmetic minima
minima_arith <- read.csv("./Output/opilio_cutline_minima_arithmeticCW.csv") %>%
  mutate(BETA0 = coef(lm(minima ~ midpoint))[1],
         BETA1 = coef(lm(minima ~ midpoint))[2])

BETA0_arith <- unique(minima_arith$BETA0)
BETA1_arith <- unique(minima_arith$BETA1)

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


# Data
# minima data
minima <- read.csv("./Output/opilio_cutline_minima.csv")
minima_arith <- read.csv("./Output/opilio_cutline_minima_arithmeticCW.csv")

# Chela data compiled by Shannon
sh_chela <- read.csv(paste0(data_dir, "specimen_chela.csv")) %>% # already != HT 17, only shell 2, no special projects
  filter(SPECIES == "SNOW", HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE,
         YEAR %in% years) %>% # filter for males, sh2, only chela msrd, not HT17
  mutate(ratio = SIZE/CHELA_HEIGHT,
         LN_CH = log(CHELA_HEIGHT),
         LN_CW = log(SIZE),
         CW = SIZE) %>%
  filter(ratio > 2 & ratio < 35) %>% # filter extreme measurements
  dplyr::select(!c(ratio)) %>%
  mutate(cutoff = BETA0 + BETA1*LN_CW, # apply cutline model
         MATURE = case_when((LN_CH > cutoff) ~ 1,
                            TRUE ~ 0))

# crabpack chela
cp_chela <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
          #dplyr::select(colnames(chela_10.13)) %>%
          filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE,
                 YEAR %in% years) %>% # filter for males, sh2, only chela msrd, not HT17
          mutate(ratio = SIZE/CHELA_HEIGHT,
                  LN_CH = log(CHELA_HEIGHT),
                  LN_CW = log(SIZE),
                  CW = SIZE) %>%
          filter(ratio > 2 & ratio < 35) %>% # filter extreme measurements
          dplyr::select(!c(ratio)) %>%
          mutate(cutoff = BETA0 + BETA1*LN_CW, # apply cutline model
                 MATURE = case_when((LN_CH > cutoff) ~ 1,
                                    TRUE ~ 0))
