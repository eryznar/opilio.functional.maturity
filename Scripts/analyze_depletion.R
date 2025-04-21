# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
# Directed fishery data
df.dat <- read.csv("./Data/opilio_directedfishery_catch.csv") %>%
  dplyr::select(Year, Retained_kt, Discarded_males_kt) %>%
  mutate(directedfish_biomass = Retained_kt + Discarded_males_kt) %>%
  dplyr::select(Year, directedfish_biomass)

# # Shannon chela data (just for 2010 and 2013 bc these are not in crabpack spec tables)
# chela_10.13 <- read.csv("./Data/snow_chela_UPDATED.csv") %>% # already != HT 17, only shell 2, no special projects
#                 filter(YEAR %in% c(2010, 2013)) %>%
#                 dplyr::select(!DATASET) 
                

# Survey specimen data
  readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
            #dplyr::select(colnames(chela_10.13)) %>%
            filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE) %>% # filter for males, sh2, only chela msrd, not HT17
            #rbind(., chela_10.13) %>% # bind with 2010, 2013 chela data from Shannon (chela weight tables)
            mutate(RATIO = SIZE/CHELA_HEIGHT) %>%
            filter(RATIO > 2 & RATIO < 35) %>% # filter extreme measurements
            dplyr::select(!c(RATIO)) %>%
            mutate(CUTOFF = BETA0 + BETA1*(log(SIZE_1MM)), # apply cutline model
                   MATURE = case_when((log(CHELA_HEIGHT) > CUTOFF) ~ 1,
                                      TRUE ~ 0),
                   CPUE = SAMPLING_FACTOR/AREA_SWEPT,
                   CPUE_KG = (CALCULATED_WEIGHT_1MM * SAMPLING_FACTOR)/1000) %>%
            dplyr::select(YEAR, SIZE_1MM, MATURE, SAMPLING_FACTOR, CALCULATED_WEIGHT_1MM, CPUE, CPUE_KG) %>%
            rename(Year = YEAR) %>%
            { . ->> SAM.dat; . } %>%
            mutate(bin = case_when((SIZE_1MM>=55 & SIZE_1MM<=65) ~ "Small (55-65mm)", # apply small and large bins for modeling below
                                   (SIZE_1MM>=95 & SIZE_1MM<=105) ~ "Large (95-105mm)")) %>%
            filter(is.na(bin) == FALSE) %>%
            group_by(Year, bin) %>%
            reframe(msrd_crab = n(),
                   msrd_crab_kg =sum(CALCULATED_WEIGHT_1MM)/1000, # multiply by # msrd crab?? No
                   mature_msrd_crab = sum(MATURE==1),
                   mature_msrd_crab_kg = sum(CALCULATED_WEIGHT_1MM[MATURE==1])/1000,
                   prop_mature = mature_msrd_crab/msrd_crab,
                   prop_mature_kg = mature_msrd_crab_kg/msrd_crab_kg) %>% # summing the sampling factor of only crab that are mature in each bin
            { . ->> bin.dat; . }
          
# Join binned survey specimen data with abundance from survey for same size bins
abund.dat <- read.csv("./Data/surveyabund_snowmales55.105.csv") %>%
                dplyr::select(YEAR, ABUNDANCE, BIOMASS_MT, bin) %>%
                rename(Year = YEAR, abundance = ABUNDANCE, biomass = BIOMASS_MT) %>%
                pivot_wider(., values_from = c("abundance", "biomass"), names_from = "bin") %>%
                rename(abundance_small = `abundance_Small (55-65mm)`, abundance_large = `abundance_Large (95-105mm)`,
                       biomass_small = `biomass_Small (55-65mm)`, biomass_large = `biomass_Large (95-105mm)`)

abund.bin.dat <- bin.dat %>%
                    dplyr::select(Year, bin, prop_mature, prop_mature_kg) %>%
                    pivot_wider(., names_from = "bin", values_from = c(prop_mature, prop_mature_kg)) %>%
                    rename(propmature_small = `prop_mature_Small (55-65mm)`, propmature_large = `prop_mature_Large (95-105mm)`,
                           propmature_small_kg = `prop_mature_kg_Small (55-65mm)`, propmature_large_kg = `prop_mature_kg_Large (95-105mm)`) %>%
                    right_join(., abund.dat) %>%
                    mutate(mature_abundance_small = (abundance_small * propmature_small)/1e6,
                           mature_abundance_large = (abundance_large * propmature_large)/1e6,
                           mature_biomass_small = (biomass_small * propmature_small)/1000,
                           mature_biomass_large = (biomass_large * propmature_small)/1000,
                           abundance_small = abundance_small/1e6,
                           abundance_large = abundance_large/1e6,
                           biomass_small = biomass_small/1000,
                           biomass_large= biomass_large/1000)
                

# Size at 50% maturity timeseries from maturity model
params <- read.csv("./Output/maturity_model_params.csv") %>%
  dplyr::select(YEAR, B_EST) %>%
  rename(Year = YEAR, size_at_mat = B_EST)

# Load March-April ice data
ice <- read.csv("./Data/ERA5ice_1972.2024.csv") %>%
          filter(name == "Mar-Apr ice cover") %>%
          dplyr::select(year, value) %>%
          rename(Year = year, MarApr_ice = value) %>%
          mutate(Year = Year - 0) # ice lag


# Bind all dataframes into df for modeling
model.dat <- right_join(df.dat, abund.bin.dat) %>%
            right_join(., ice) %>%
            #filter(is.na(size_at_mat) == FALSE) %>%
            right_join(., params) %>%
            mutate(mature_abundance_small = log(mature_abundance_small + 1),
                   mature_abundance_large = log(mature_abundance_large + 1),
                   mature_biomass_small = log(mature_biomass_small + 1),
                   mature_biomass_large = log(mature_biomass_large + 1),
                   size_at_mat = log(size_at_mat),
                   lag1_mature_abundance_small = lag(mature_abundance_small, n = 1),
                   lag1_mature_abundance_large = lag(mature_abundance_large, n = 1),
                   lag1_mature_biomass_small = lag(mature_biomass_small, n = 1),
                   lag1_mature_biomass_large = lag(mature_biomass_large, n = 1),
                   lag3_MarApr_ice = lag(MarApr_ice, n = 3))

write.csv(model.dat, "./Data/snow_male_GAM_modeldat.csv")

# FIT MODELS ---------------------------------------------------------------------------------------------
# proportion mature in small bin (55-65) ----
  # abundance covariates ----
  small.int.nolag.aa <- gam(propmature_small ~ s(directedfish_biomass, k = 4) + 
                           s(mature_abundance_small, mature_abundance_large, k = 4) +
                           s(MarApr_ice, k = 4),
                           family = betar(link = "logit"), 
                           method =  "REML", 
                           data = model.dat)
  
  
  small.main.nolag.aa <- gam(propmature_small ~ s(directedfish_biomass, k = 4) + 
                          s(mature_abundance_small, k = 4) +
                          s(mature_abundance_large, k = 4) +
                          s(MarApr_ice, k = 4),
                          family = betar(link = "logit"), 
                          method =  "REML", 
                          data = model.dat)
  
  small.int.lagabund.aa <- gam(propmature_small ~ s(directedfish_biomass, k = 4) + 
                            s(lag1_mature_abundance_small, lag1_mature_abundance_large, k = 4) +
                            s(MarApr_ice, k = 4),
                            family = betar(link = "logit"), 
                            method =  "REML", 
                            data = model.dat)
  
  small.main.lagabund.aa <- gam(propmature_small ~ s(directedfish_biomass, k = 4) + 
                            s(lag1_mature_abundance_small, k = 4) +
                            s(lag1_mature_abundance_large, k = 4) +
                            s(MarApr_ice, k = 4),
                            family = betar(link = "logit"), 
                            method =  "REML", 
                            data = model.dat)
  
  small.int.icelag.aa <- gam(propmature_small ~ s(directedfish_biomass, k = 4) + 
                            s(mature_abundance_small, mature_abundance_large, k = 4) +
                            s(lag3_MarApr_ice, k = 4),
                            family = betar(link = "logit"), 
                            method =  "REML", 
                            data = model.dat)
  
  small.main.icelag.aa <- gam(propmature_small ~ s(directedfish_biomass, k = 4) + 
                             s(mature_abundance_small, k = 4) +
                             s(mature_abundance_large, k = 4) +
                             s(lag3_MarApr_ice, k = 4),
                             family = betar(link = "logit"), 
                             method =  "REML", 
                             data = model.dat)
  
  small.int.iceabundlag.aa <- gam(propmature_small ~ s(directedfish_biomass, k = 4) + 
                            s(lag1_mature_abundance_small, lag1_mature_abundance_large, k = 4) +
                            s(lag3_MarApr_ice, k = 4),
                            family = betar(link = "logit"), 
                          data = model.dat)
  
  small.main.iceabundlag.aa <- gam(propmature_small ~ s(directedfish_biomass, k = 4) + 
                             s(lag1_mature_abundance_small, k = 4) +
                             s(lag1_mature_abundance_large, k = 4) +
                             s(lag3_MarApr_ice, k = 4),
                             family = betar(link = "logit"), 
                           method =  "REML", 
                           data = model.dat)
 
  # biomass covariates ----
  small.int.nolag.bb <- gam(propmature_small_kg ~ s(directedfish_biomass, k = 4) + 
                              s(mature_biomass_small, mature_biomass_large, k = 4) +
                              s(MarApr_ice, k = 4),
                            family = betar(link = "logit"), 
                            method =  "REML", 
                            data = model.dat)
  
  
  small.main.nolag.bb <- gam(propmature_small_kg ~ s(directedfish_biomass, k = 4) + 
                               s(mature_biomass_small, k = 4) +
                               s(mature_biomass_large, k = 4) +
                               s(MarApr_ice, k = 4),
                             family = betar(link = "logit"), 
                             method =  "REML", 
                             data = model.dat)
  
  small.int.lagabund.bb <- gam(propmature_small_kg~ s(directedfish_biomass, k = 4) + 
                                 s(lag1_mature_biomass_small, lag1_mature_biomass_large, k = 4) +
                                 s(MarApr_ice, k = 4),
                               family = betar(link = "logit"), 
                               method =  "REML", 
                               data = model.dat)
  
  small.main.lagabund.bb <- gam(propmature_small_kg ~ s(directedfish_biomass, k = 4) + 
                                  s(lag1_mature_biomass_small, k = 4) +
                                  s(lag1_mature_biomass_large, k = 4) +
                                  s(MarApr_ice, k = 4),
                                family = betar(link = "logit"), 
                                method =  "REML", 
                                data = model.dat)
  
  small.int.icelag.bb <- gam(propmature_small_kg ~ s(directedfish_biomass, k = 4) + 
                               s(mature_biomass_small, mature_biomass_large, k = 4) +
                               s(lag3_MarApr_ice, k = 4),
                             family = betar(link = "logit"), 
                             method =  "REML", 
                             data = model.dat)
  
  small.main.icelag.bb <- gam(propmature_small_kg~ s(directedfish_biomass, k = 4) + 
                                s(mature_biomass_small, k = 4) +
                                s(mature_biomass_large, k = 4) +
                                s(lag3_MarApr_ice, k = 4),
                              family = betar(link = "logit"), 
                              method =  "REML", 
                              data = model.dat)
  
  small.int.iceabundlag.bb <- gam(propmature_small_kg ~ s(directedfish_biomass, k = 4) + 
                                    s(lag1_mature_biomass_small, lag1_mature_biomass_large, k = 4) +
                                    s(lag3_MarApr_ice, k = 4),
                                  family = betar(link = "logit"), 
                                  data = model.dat)
  
  small.main.iceabundlag.bb <- gam(propmature_small_kg ~ s(directedfish_biomass, k = 4) + 
                                     s(lag1_mature_biomass_small, k = 4) +
                                     s(lag1_mature_biomass_large, k = 4) +
                                     s(lag3_MarApr_ice, k = 4),
                                   family = betar(link = "logit"), 
                                   method =  "REML", 
                                   data = model.dat)
  
  # model comparison ----
  AICc(small.int.nolag.aa, small.main.nolag.aa, small.int.lagabund.aa, small.main.lagabund.aa, 
       small.int.icelag.aa, small.main.icelag.aa, small.int.iceabundlag.aa, small.main.iceabundlag.aa,
       small.int.nolag.bb, small.main.nolag.bb, small.int.lagabund.bb, small.main.lagabund.bb, 
       small.int.icelag.bb, small.main.icelag.bb, small.int.iceabundlag.bb, small.main.iceabundlag.bb) %>%
    mutate(BEST = case_when((AICc == min(AICc)) ~ "Y",
                            TRUE ~ "N"))
  
  
  # diagnostics ----
  diagnose(small.int.nolag.bb)
  
  
  diagnose(small.main.nolag.bb)
  
  
  diagnose(small.int.nolag.aa)

  
  diagnose(small.main.nolag.aa)
 
  
# proportion mature in large bin (95-105mm) ----
  # abundance covariates ----
  large.int.nolag.aa<- gam(propmature_large ~ s(directedfish_biomass, k = 4) + 
                           s(mature_abundance_small, mature_abundance_large, k = 4) +
                           s(MarApr_ice, k = 4),
                         family = betar(link = "logit"), 
                         method =  "REML", 
                         data = model.dat)
  
  large.main.nolag.aa <- gam(propmature_large ~ s(directedfish_biomass, k = 4) + 
                            s(mature_abundance_small, k = 4) +
                            s(mature_abundance_large, k = 4) +
                            s(MarApr_ice, k = 4),
                          family = betar(link = "logit"), 
                          method =  "REML", 
                          data = model.dat)
  
  large.int.lagabund.aa <- gam(propmature_large ~ s(directedfish_biomass, k = 4) + 
                              s(lag1_mature_abundance_small, lag1_mature_abundance_large, k = 4) +
                              s(MarApr_ice, k = 4),
                            family = betar(link = "logit"), 
                            method =  "REML", 
                            data = model.dat)
  
  large.main.lagabund.aa <- gam(propmature_large ~ s(directedfish_biomass, k = 4) + 
                               s(lag1_mature_abundance_small, k = 4) +
                               s(lag1_mature_abundance_large, k = 4) +
                               s(MarApr_ice, k = 4),
                             family = betar(link = "logit"), 
                             method =  "REML", 
                             data = model.dat)
  
  large.int.icelag.aa <- gam(propmature_large ~ s(directedfish_biomass, k = 4) + 
                            s(mature_abundance_small, mature_abundance_large, k = 4) +
                            s(lag3_MarApr_ice, k = 4),
                          family = betar(link = "logit"), 
                          method =  "REML", 
                          data = model.dat)
  
  large.main.icelag.aa <- gam(propmature_large ~ s(directedfish_biomass, k = 4) + 
                             s(mature_abundance_small, k = 4) +
                             s(mature_abundance_large, k = 4) +
                             s(lag3_MarApr_ice, k = 4),
                           family = betar(link = "logit"), 
                           method =  "REML", 
                           data = model.dat)
  
  large.int.iceabundlag.aa <- gam(propmature_large ~ s(directedfish_biomass, k = 4) + 
                                 s(lag1_mature_abundance_small, lag1_mature_abundance_large, k = 4) +
                                 s(lag3_MarApr_ice, k = 4),
                               family = betar(link = "logit"), 
                               method =  "REML", 
                               data = model.dat)
  
  large.main.iceabundlag.aa <- gam(propmature_large ~ s(directedfish_biomass, k = 4) + 
                                  s(lag1_mature_abundance_small, k = 4) +
                                  s(lag1_mature_abundance_large, k = 4) +
                                  s(lag3_MarApr_ice, k = 4),
                                family = betar(link = "logit"), 
                                method =  "REML", 
                                data = model.dat)
  
  # biomass covariates ----
  large.int.nolag.bb<- gam(propmature_large_kg ~ s(directedfish_biomass, k = 4) + 
                             s(mature_biomass_small, mature_biomass_large, k = 4) +
                             s(MarApr_ice, k = 4),
                           family = betar(link = "logit"), 
                           method =  "REML", 
                           data = model.dat)
  
  large.main.nolag.bb <- gam(propmature_large_kg~ s(directedfish_biomass, k = 4) + 
                               s(mature_biomass_small, k = 4) +
                               s(mature_biomass_large, k = 4) +
                               s(MarApr_ice, k = 4),
                             family = betar(link = "logit"), 
                             method =  "REML", 
                             data = model.dat)
  
  large.int.lagabund.bb <- gam(propmature_large_kg ~ s(directedfish_biomass, k = 4) + 
                                 s(lag1_mature_biomass_small, lag1_mature_biomass_large, k = 4) +
                                 s(MarApr_ice, k = 4),
                               family = betar(link = "logit"), 
                               method =  "REML", 
                               data = model.dat)
  
  large.main.lagabund.bb <- gam(propmature_large_kg ~ s(directedfish_biomass, k = 4) + 
                                  s(lag1_mature_biomass_small, k = 4) +
                                  s(lag1_mature_biomass_large, k = 4) +
                                  s(MarApr_ice, k = 4),
                                family = betar(link = "logit"), 
                                method =  "REML", 
                                data = model.dat)
  
  large.int.icelag.bb <- gam(propmature_large_kg ~ s(directedfish_biomass, k = 4) + 
                               s(mature_biomass_small, mature_biomass_large, k = 4) +
                               s(lag3_MarApr_ice, k = 4),
                             family = betar(link = "logit"), 
                             method =  "REML", 
                             data = model.dat)
  
  large.main.icelag.bb <- gam(propmature_large_kg ~ s(directedfish_biomass, k = 4) + 
                                s(mature_biomass_small, k = 4) +
                                s(mature_biomass_large, k = 4) +
                                s(lag3_MarApr_ice, k = 4),
                              family = betar(link = "logit"), 
                              method =  "REML", 
                              data = model.dat)
  
  large.int.iceabundlag.bb <- gam(propmature_large_kg ~ s(directedfish_biomass, k = 4) + 
                                    s(lag1_mature_biomass_small, lag1_mature_biomass_large, k = 4) +
                                    s(lag3_MarApr_ice, k = 4),
                                  family = betar(link = "logit"), 
                                  method =  "REML", 
                                  data = model.dat)
  
  large.main.iceabundlag.bb <- gam(propmature_large_kg ~ s(directedfish_biomass, k = 4) + 
                                     s(lag1_mature_biomass_small, k = 4) +
                                     s(lag1_mature_biomass_large, k = 4) +
                                     s(lag3_MarApr_ice, k = 4),
                                   family = betar(link = "logit"), 
                                   method =  "REML", 
                                   data = model.dat)
  
  # model comparison ----
  AICc(large.int.nolag.aa, large.main.nolag.aa, large.int.lagabund.aa, large.main.lagabund.aa, 
       large.int.icelag.aa, large.main.icelag.aa, large.int.iceabundlag.aa, large.main.iceabundlag.aa,
       large.int.nolag.bb, large.main.nolag.bb, large.int.lagabund.bb, large.main.lagabund.bb, 
       large.int.icelag.bb, large.main.icelag.bb, large.int.iceabundlag.bb, large.main.iceabundlag.bb) %>%
    mutate(BEST = case_when((AICc == min(AICc)) ~ "Y",
                            TRUE ~ "N"))
  
  
  # diagnostics ----
  diagnose(large.int.nolag.bb)
  diagnose(large.main.nolag.bb)
 
  diagnose(large.int.nolag.aa)
  diagnose(large.main.nolag.aa)
  
  
# size at maturity ----
  # abundance covariates ----
  SAM.int.nolag.aa<- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                          s(mature_abundance_small, mature_abundance_large, k = 4) +
                          s(MarApr_ice, k = 4),
                          method =  "REML", 
                          data = model.dat)
  
  SAM.main.nolag.aa <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                               s(mature_abundance_small, k = 4) +
                               s(mature_abundance_large, k = 4) +
                               s(MarApr_ice, k = 4),
                             method =  "REML", 
                             data = model.dat)
  
  SAM.int.lagabund.aa <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                                 s(lag1_mature_abundance_small, lag1_mature_abundance_large, k = 4) +
                                 s(MarApr_ice, k = 4),
                               method =  "REML", 
                               data = model.dat)
  
  SAM.main.lagabund.aa <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                                  s(lag1_mature_abundance_small, k = 4) +
                                  s(lag1_mature_abundance_large, k = 4) +
                                  s(MarApr_ice, k = 4),
                                method =  "REML", 
                                data = model.dat)
  
  SAM.int.icelag.aa <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                               s(mature_abundance_small, mature_abundance_large, k = 4) +
                               s(lag3_MarApr_ice, k = 4),
                             method =  "REML", 
                             data = model.dat)
 
  
  SAM.main.icelag.aa <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                                s(mature_abundance_small, k = 4) +
                                s(mature_abundance_large, k = 4) +
                                s(lag3_MarApr_ice, k = 4),
                              method =  "REML", 
                              data = model.dat)
  
  SAM.int.iceabundlag.aa <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                                    s(lag1_mature_abundance_small, lag1_mature_abundance_large, k = 4) +
                                    s(lag3_MarApr_ice, k = 4),
                                  method =  "REML", 
                                  data = model.dat)
  
  SAM.main.iceabundlag.aa <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                                     s(lag1_mature_abundance_small, k = 4) +
                                     s(lag1_mature_abundance_large, k = 4) +
                                     s(lag3_MarApr_ice, k = 4),
                                   method =  "REML", 
                                   data = model.dat)
  
  # biomass covariates ----
  SAM.int.nolag.bb<- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                           s(mature_biomass_small, mature_biomass_large, k = 4) +
                           s(MarApr_ice, k = 4),
                         method =  "REML", 
                         data = model.dat)
  
  SAM.main.nolag.bb <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                             s(mature_biomass_small, k = 4) +
                             s(mature_biomass_large, k = 4) +
                             s(MarApr_ice, k = 4),
                           method =  "REML", 
                           data = model.dat)
  
  SAM.int.lagabund.bb <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                               s(lag1_mature_biomass_small, lag1_mature_biomass_large, k = 4) +
                               s(MarApr_ice, k = 4),
                             method =  "REML", 
                             data = model.dat)
  
  SAM.main.lagabund.bb <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                                s(lag1_mature_biomass_small, k = 4) +
                                s(lag1_mature_biomass_large, k = 4) +
                                s(MarApr_ice, k = 4),
                              method =  "REML", 
                              data = model.dat)
  
  SAM.int.icelag.bb <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                             s(mature_biomass_small, mature_biomass_large, k = 4) +
                             s(lag3_MarApr_ice, k = 4),
                           method =  "REML", 
                           data = model.dat)
  
  SAM.main.icelag.bb <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                              s(mature_biomass_small, k = 4) +
                              s(mature_biomass_large, k = 4) +
                              s(lag3_MarApr_ice, k = 4),
                            method =  "REML", 
                            data = model.dat)
  
  SAM.int.iceabundlag.bb <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                                  s(lag1_mature_biomass_small, lag1_mature_biomass_large, k = 4) +
                                  s(lag3_MarApr_ice, k = 4),
                                method =  "REML", 
                                data = model.dat)
  
  SAM.main.iceabundlag.bb <- gam(size_at_mat ~ s(directedfish_biomass, k = 4) + 
                                   s(lag1_mature_biomass_small, k = 4) +
                                   s(lag1_mature_biomass_large, k = 4) +
                                   s(lag3_MarApr_ice, k = 4),
                                 method =  "REML", 
                                 data = model.dat)
  
 # model comparison ----
  AICc(SAM.int.nolag.aa, SAM.main.nolag.aa, SAM.int.lagabund.aa, SAM.main.lagabund.aa, 
       SAM.int.icelag.aa, SAM.main.icelag.aa, SAM.int.iceabundlag.aa, SAM.main.iceabundlag.aa,
       SAM.int.nolag.bb, SAM.main.nolag.bb, SAM.int.lagabund.bb, SAM.main.lagabund.bb, 
       SAM.int.icelag.bb, SAM.main.icelag.bb, SAM.int.iceabundlag.bb, SAM.main.iceabundlag.bb) %>%
    mutate(BEST = case_when((AICc == min(AICc)) ~ "Y",
                            TRUE ~ "N"))

 # diagnostics ----
summary(SAM.main.icelag.aa)
gam.check(SAM.main.icelag.aa)
appraise(SAM.main.icelag.aa)
draw(SAM.main.icelag.aa)

summary(SAM.int.icelag.aa)
gam.check(SAM.int.icelag.aa)
appraise(SAM.int.icelag.aa)
draw(SAM.int.icelag.aa)

summary(SAM.main.icelag.bb)
gam.check(SAM.main.icelag.bb)
appraise(SAM.main.icelag.bb)
draw(SAM.main.icelag.bb)

summary(SAM.int.icelag.bb)
gam.check(SAM.int.icelag.bb)
appraise(SAM.int.icelag.bb)
draw(SAM.int.icelag.bb)

