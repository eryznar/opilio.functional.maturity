library(tidyverse)
library(odbc)
library(DBI)
library(getPass)
library(keyring)
library(lifecycle)
library(data.table)
library(crabpack)


# Read in survey data
specEBS <- readRDS("./Data/snow_survey_specimenEBS.rda")
matEBS <- readRDS("./Data/snow_survey_maturityEBS.rda")