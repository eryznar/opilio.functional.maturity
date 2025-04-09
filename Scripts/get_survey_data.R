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


# Get small and large mature male abundance/biomass
spec <- readRDS("./Data/snow_survey_specimenEBS.rda")
abund_55.65 <- calc_bioabund(crab_data = spec, 
                             sex = "male", 
                             region = "EBS", 
                             district = "ALL",
                             species = "SNOW", 
                             size_min = 55, 
                             size_max = 65, 
                             shell_condition = "new_hardshell") %>%
                mutate(bin = "Small (55-65mm)")

abund_95.105 <- calc_bioabund(crab_data = spec, 
                               sex = "male", 
                               region = "EBS", 
                               district = "ALL",
                               species = "SNOW", 
                               size_min = 95, 
                               size_max = 105, 
                               shell_condition = "new_hardshell") %>%
                mutate(bin = "Large (95-105mm)")

write.csv(rbind(abund_55.65, abund_95.105), "./Data/surveyabund_snowmales55.105.csv")


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
t.params <- male_maturity_data$model_parameters
