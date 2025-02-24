# PURPOSE: to use morphometric maturity cutlines to define immature vs. mature for snow crab shell 2 males

# Author: Emily Ryznar, Jon Richar, Shannon Hennessey

# Load libs/params
source("./Scripts/load_libs_params.R")

# Read in maturity output from crabpack
ratio <- readRDS("./Data/snow_survey_maturityEBS.rda")$male_mat_ratio
pars <- readRDS("./Data/snow_survey_maturityEBS.rda")$model_parameters

# Read in minima data and run model to create cutlines
minima <- read.csv("./Output/opilio_cutline_minima.csv") %>%
  mutate(BETA0 = coef(lm(minima ~ midpoint))[1],
         BETA1 = coef(lm(minima ~ midpoint))[2])

# Read and filter survey data
haul <- readRDS("./Data/snow_survey_specimenEBS.rda")$haul %>%
  filter(#YEAR %in% c(2010, 2017, 2018, 2019) # years Jon uses for stock assessment estimates?,
          HAUL_TYPE !=17)

haul %>%
  dplyr::select(HAULJOIN, YEAR, MID_LATITUDE, MID_LONGITUDE) %>%
  rename(LATITUDE = MID_LATITUDE, LONGITUDE = MID_LONGITUDE)-> hauljoin

chela <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(SEX == 1,
         SHELL_CONDITION == 2, 
         HAUL_TYPE !=17, 
         #YEAR %in% c(2010, 2017, 2018, 2019) # years Jon uses for stock assessment estimates?
         is.na(CHELA_HEIGHT) == FALSE) %>%
  dplyr::select(HAULJOIN, SPECIES, REGION, DISTRICT, YEAR, SIZE, CHELA_HEIGHT, SAMPLING_FACTOR)

# Read in Jon data for 2010 and 2013
jon.dat <- read.csv("./Data/opilio_chela_height_TS.csv") %>%
  mutate(YEAR = as.numeric(substr(CRUISE, 1, 4)),
         SPECIES = "SNOW",
         REGION = "EBS",
         DISTRICT = "ALL") %>%
  dplyr::filter(YEAR %in% c(2010, 2013), 
                #HAUL_TYPE !=17, need to make sure this was filtered by Jon
                SEX == 1, 
                SHELL_CONDITION == 2) %>%
  dplyr::select(HAULJOIN, SPECIES, REGION, DISTRICT, YEAR, WIDTH, CHELA_HEIGHT, SAMPLING_FACTOR) %>%
  rename(SIZE = WIDTH)

# Join survey data and Jon's data from 2010 and 2013
chela <- rbind(jon.dat, chela)

# Pool and transform data using natural log
chela %>%
  mutate(LN_CW = log(SIZE),
         LN_CH = log(CHELA_HEIGHT)) -> chela2

# Subset data into size intervals at ln(CW) of 10mm and add in cutline params
chela2 %>%
  #filter(LN_CW >= 3.9 & LN_CW <= 4.6) %>% # these were filtered out by Jon
  mutate(BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),
         BIN2 = BIN,
         BETA0 = unique(minima$BETA0), # add in cutline params
         BETA1 = unique(minima$BETA1)) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         MIDPOINT = (UPPER + LOWER)/2) -> bin.dat

# Calculate proportion mature by size bin
prop.dat <- bin.dat %>%
            mutate(MATURITY = case_when(LN_CH >= (BETA0 + BETA1*LN_CW) ~ "MATURE",
                                        LN_CH < (BETA0 + BETA1*LN_CW) ~ "IMMATURE"),
                   WIDTH_1MM = floor(SIZE)) %>%
            unite("AREA", c(SPECIES, REGION, DISTRICT), sep = "_") %>%
            group_by(AREA, YEAR, HAULJOIN, BIN, MATURITY) %>%
            # Sum n_crabs per size bin
            reframe(COUNT = sum(SAMPLING_FACTOR))  %>%
            #dplyr::select(AREA, YEAR, HAULJOIN, BIN, MATURITY, COUNT) %>%
            # Expand/add 0-count size bins
            right_join(., right_join(hauljoin, expand_grid(BIN = unique(cut_width(c(0:200), width = 10, center = 5, closed = "left")),
                                             AREA = c("SNOW_EBS_ALL"), 
                                             MATURITY = c("MATURE", "IMMATURE"),
                                             YEAR = c(1989:2019, 2021:2024)), # also joining to unique hauljoins by year
                                  relationship = "many-to-many")) %>%
            replace_na(list(COUNT = 0)) %>%
            pivot_wider(names_from = MATURITY, values_from = COUNT) %>%
            mutate(TOTAL_CRAB= IMMATURE + MATURE,
                   PROP_MATURE = MATURE/TOTAL_CRAB) %>%
            separate(BIN, sep = ",", into = c("LOWER", "UPPER")) %>%
            mutate(LOWER = as.numeric(sub('.', '', LOWER)),
                   UPPER = as.numeric(gsub('.$', '', UPPER)),
                   MIDPOINT = (UPPER + LOWER)/2,
                   PROP_MATURE = case_when((LOWER <= 40) ~ 0, # set all crab in 40-50 bin or smaller to immature; previous code assigns 0 and 1 to crab if no crab were caught
                                           (LOWER >= 130) ~ 1, # set all crab in 130-140 bin or larger to mature
                                           TRUE ~ PROP_MATURE),
                   SIZE_BIN = paste0(LOWER, "-", UPPER)) %>%
            separate(AREA, sep = "_", into = c("SPECIES", "REGION", "DISTRICT")) %>%
            # Remove specific years: EBS snow - 2008, 2012, 2014, 2016, 2020; tanner - 1989, 2013, 2015, 2020, tannerE - 2011 too
            # filter(!((SPECIES == "TANNER" & YEAR %in% c(1989, 2013, 2015, 2020)) | (DISTRICT == "E166" & YEAR == 2011) |
            #            (SPECIES == "SNOW" & REGION == "EBS" & YEAR %in% c(2008, 2012, 2014, 2016, 2020)))) %>%
            # Format final output file
            rename(NUM_IMMATURE = IMMATURE,
                   NUM_MATURE = MATURE) %>%
            #replace_na(list(PROP_MATURE = 0)) %>% # SHOULD NOT REPLACE NAs with ZERO IF THERE WERE NO CRAB CAUGHT
            dplyr::select(SPECIES, REGION, DISTRICT, YEAR, SIZE_BIN, MIDPOINT, LATITUDE, LONGITUDE, NUM_IMMATURE, NUM_MATURE, TOTAL_CRAB, PROP_MATURE)
          
write.csv(prop.dat, "./Output/opilio_propmat_latlon.csv")

# Aggregate data by year
agg.prop.dat <- prop.dat %>%
                group_by(SPECIES, REGION, DISTRICT, YEAR, SIZE_BIN, MIDPOINT) %>%
                reframe(NUM_IMMATURE = sum(NUM_IMMATURE),
                        NUM_MATURE = sum(NUM_MATURE)) %>%
                mutate(TOTAL_CRAB = NUM_IMMATURE + NUM_MATURE,
                       PROP_MATURE = NUM_MATURE/TOTAL_CRAB,
                       LOWER = as.numeric(sub("-.*", "", SIZE_BIN)),
                       UPPER =  as.numeric(sub(".*-", "", SIZE_BIN)),
                       PROP_MATURE = case_when((LOWER <= 40 & TOTAL_CRAB !=0) ~ 0, # set all crab in 40-50 bin or smaller to immature; previous code assigns 0 and 1 to crab if no crab were caught
                                               (LOWER >= 130 & TOTAL_CRAB !=0) ~ 1, # set all crab in 130-140 bin or larger to mature
                                               TRUE ~ PROP_MATURE)) %>%
                #replace_na(list(PROP_MATURE = 0)) %>% # SHOULD NOT REPLACE NAs with ZERO IF THERE WERE NO CRAB CAUGHT
                dplyr::select(!c(LOWER, UPPER)) %>%
                group_by(YEAR) %>%
                arrange(., MIDPOINT, .by_group = TRUE) %>%
                ungroup()

write.csv(agg.prop.dat, "./Output/opilio_propmat_agg.csv")


## Fit size at maturity model --------------------------------------------------
# Jon filters out 2008, 2012, 2014, 2016 due to non-convergence?
# 2020 is omitted bc we didn't have a survey
# 2010 and 2013 are problematic for me I don't have any data (only via special projects?)

yrs <- unique(agg.prop.dat$YEAR)
yrs <- c(1991:2007,2009:2011, 2013, 2015, 2017:2019, 2021:2024)
params <- data.frame()
preds <- data.frame()

for(ii in 1:length(yrs)){
  print(paste("Fitting year", yrs[ii]))
  
  #filter by year
  agg.prop.dat %>%
    filter(YEAR == yrs[ii]) -> mod.dat
  
  #fit nls model
  mod <- nls(PROP_MATURE ~ (1/(1 + exp(-a*(as.numeric(MIDPOINT) - b)))),
              data = mod.dat,
              start = list(a = 0.10, b = 60.0),
              na.action = na.omit, 
              nls.control(maxiter = 5000))
  
  #pull out params
  A_EST <- summary(mod)$coefficients[1,1]
  A_SE <- summary(mod)$coefficients[1,2]
  B_EST <- summary(mod)$coefficients[2,1]
  B_SE <- summary(mod)$coefficients[2,2]
  
  #create summary df
  out <- data.frame(SPECIES = unique(mod.dat$SPECIES),
                    REGION = unique(mod.dat$REGION),
                    DISTRICT = unique(mod.dat$DISTRICT),
                    YEAR = yrs[ii],
                    A_EST = A_EST,
                    A_SE = A_SE,
                    B_EST = B_EST,
                    B_SE = B_SE)
  
  params <- rbind(params, out)
  
  #create prediction df
  out2 <- data.frame(SPECIES = unique(mod.dat$SPECIES),
                    REGION = unique(mod.dat$REGION),
                    DISTRICT = unique(mod.dat$DISTRICT),
                    YEAR = yrs[ii],
                    MIDPOINT = mod.dat$MIDPOINT,
                    PROP_MATURE = predict(mod, mod.dat))
  
  preds <- rbind(preds, out2)
  
}

write.csv(params, "./Output/maturity_model_params.csv")

##Plot maturity ogives
ggplot(preds, aes(MIDPOINT, PROP_MATURE, group = YEAR))+
  geom_line()+
  theme_bw()+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue")


