source("./Scripts/load_libs_params.R")




# Connect to Oracle via AFSC database (NOAA AFSC users only)
#channel <- crabpack::get_connected(db = "AFSC") #username: CRABBASE, pass: Chionoecetes_24
channel <- "API"

# Pull specimen data
species <- "SNOW"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             years = c(1975:2024),
                                             channel = channel)

saveRDS(specimen_data, "./Data/snow_survey_specimenEBS.rda")

# Separate chela data
chela <- specimen_data$specimen %>%
            filter(YEAR > 2005) # chela height only collected >=2006


# Get maturity
male_maturity_data <- crabpack::get_male_maturity(species = "SNOW",
                                                  region = "EBS",
                                                  channel = channel)

saveRDS(male_maturity_data, "./Data/snow_survey_maturityEBS.rda")

mat.dat <- male_maturity_data$male_mat_ratio
params <- male_maturity_data$model_parameters


# TANNER

# Pull specimen data
species <- "TANNER"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             years = c(1975:2024),
                                             channel = channel)

saveRDS(specimen_data, "./Data/tanner_survey_specimenEBS.rda")

# Separate chela data
chela <- specimen_data$specimen %>%
  filter(YEAR > 2005) # chela height only collected >=2006


# Get maturity
male_maturity_data <- crabpack::get_male_maturity(species = "TANNER",
                                                  region = "EBS",
                                                  channel = channel)

saveRDS(male_maturity_data, "./Data/tanner_survey_maturityEBS.rda")

mat.dat <- male_maturity_data$male_mat_ratio
params <- male_maturity_data$model_parameters
