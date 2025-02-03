source("./Scripts/load_libs_params.R")




# Connect to Oracle via AFSC database (NOAA AFSC users only)
channel <- crabpack::get_connected(db = "AFSC") #username: CRABBASE, pass: Chionoecetes_24

# Pull specimen data
species <- "SNOW"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             years = c(1975:2024),
                                             channel = channel)

# Separate chela data
chela <- specimen_data$specimen %>%
            filter(YEAR > 2005) # chela height only collected >=2006


# Get maturity
male_maturity_data <- crabpack::get_male_maturity(species = "SNOW",
                                                  region = "EBS",
                                                  channel = channel)

mat.dat <- male_maturity_data$male_mat_ratio
params <- male_maturity_data$model_parameters
