# PURPOSE: to fit sdmTMB models for male snow crab

# Author: Emily Ryznar

# LOAD LIBS/PARAMS -----------------------------------------------------------------
source("./Scripts/load_libs_params.R")

# LOAD RESPONSE DATA ---------------------------------------------------------------
snow <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
            filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2) %>%
            mutate(CUTOFF = BETA0 + BETA1*(log(SIZE_1MM)),
                   MATURE = case_when((log(CHELA_HEIGHT) > CUTOFF) ~ 1,
                                      TRUE ~ 0)) 
  #==maturity coef # -2.8434, 0.2404
  survDAT$cutoff<- -2.8434+0.2404*survDAT$SIZE_1MM
  survDAT$mature<-survDAT$CHELA_HEIGHT>survDAT$cutoff
  
            st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
            st_transform(., crs = "+proj=utm +zone=2") %>%
            cbind(st_coordinates(.)) %>%
            as.data.frame(.) %>%
            mutate(X = X/1000, Y = Y/1000) %>%
            rename(LONGITUDE = X, LATITUDE = Y) %>%
            replace_na(list(PROP_MATURE = 0))
  

# LOAD FUNCTION TO FIT SDMTMB MODELS -----------------------------------------------
fit_models <- function(data, sex, years, dist, knots){
  
  # Filter data by params
  data %>%
    filter(YEAR %in% years) %>%
    mutate(year_fac = as.factor(YEAR)) -> data2
  
  # Make mesh
  mesh2 <- make_mesh(data2, c("LONGITUDE","LATITUDE"), n_knots = knots, type = "kmeans")
  
  if(dist == "DG"){
    # Fit prop mature model
    mod <- sdmTMB(PROP_MATURE ~ 1, # aggregate
                    spatial = "on",
                    #spatiotemporal = "iid",
                    mesh = mesh2,
                    family = delta_gamma(type = "poisson-link"),
                    #time = "YEAR",
                    anisotropy = TRUE,
                    data = data2)
    
   
  } else if(dist == "TW"){
    # Fit prop mature model
    mod <- sdmTMB(PROP_MATURE ~ 1, #the 0 is there so there is a factor predictor for each time slice
                    spatial = "on",
                    #spatiotemporal = "iid",
                    mesh = mesh2,
                    family = tweedie(link = "log"),
                    #time = "YEAR",
                    anisotropy = TRUE,
                    data = data2)
    
  } else{
    # Fit prop mature model
    mod <- sdmTMB(PROP_MATURE ~ 1, #the 0 is there so there is a factor predictor for each time slice
                    spatial = "on",
                    #spatiotemporal = "iid",
                    mesh = mesh2,
                    family = delta_lognormal(),
                    #time = "YEAR",
                    anisotropy = TRUE,
                    data = data2)
    
    
  }
  
  
  saveRDS(mod, paste0("./Output/Models/", sex, "_", dist, "_", knots, "_maturity.rda"))
   
  return(list(modTMB = mod, mesh = mesh2))
  
}

# RUN FUNCTION ---------------------------------------------------------------------
years <- c(1975:2019, 2021:2024)

out <- fit_models(snow, "male", years, "TW", 50)

# PREDICT --------------------------------------------------------------------------
pred <- predict(out$modTMB, newdata = pred.grid %>% 
                  dplyr::select(!year) %>%
                  distinct(), return_tmb_object = T)

dat <- pred$data %>%
  mutate(value = plogis(est1) * exp(est2))

ggplot(dat) +
  geom_tile(aes(y = LATITUDE, x = LONGITUDE, fill = value)) +
  scale_fill_viridis_c(name = "Proportion mature")+
  labs(y = "Latitude",
       x = "Longitude") +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal")
