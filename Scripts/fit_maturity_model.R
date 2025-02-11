# PURPOSE: to use morphometric maturity cutlines to define immature vs. mature for snow crab shell 2 males

# Author: Emily Ryznar, Jon Richar, Shannon Hennessey

# Load libs/params
source("./Scripts/load_libs_params.R")

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
         is.na(CHELA_HEIGHT) == FALSE)


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
# create expansion grid
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
                   PROP_MATURE = case_when((LOWER <= 40) ~ 0, # set all crab in 40-50 bin or smaller to immature
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
            replace_na(list(PROP_MATURE = 0)) %>%
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
                       PROP_MATURE = case_when((LOWER <= 40) ~ 0, # set all crab in 40-50 bin or smaller to immature
                                               (LOWER >= 130) ~ 1, # set all crab in 130-140 bin or larger to mature
                                               TRUE ~ PROP_MATURE)) %>%
                replace_na(list(PROP_MATURE = 0)) %>%
                dplyr::select(!c(LOWER, UPPER))

write.csv(agg.prop.dat, "./Output/opilio_propmat_agg.csv")

## Fit size at maturity model --------------------------------------------------
latlon params <- agg.prop.dat %>%
                  group_by(YEAR) %>%
                  do(model = nls(PROP_MATURE ~ (1/(1 + exp(-a*(as.numeric(MIDPOINT) - b)))),
                                 data = .,
                                 start = list(a = 0.10, b = 60.0),
                                 na.action = na.omit)) %>% 
                  ungroup() %>%
                  mutate(nls = lapply(model, broom::tidy)) %>%
                  unnest(nls) %>%
                  pivot_wider(names_from = term, values_from = c(estimate, std.error, statistic, p.value)) %>%
                  rename(A_EST = estimate_a, 
                         A_SE = std.error_a,
                         B_EST = estimate_b,
                         B_SE = std.error_b) %>%
                  select(SPECIES, REGION, DISTRICT, YEAR, A_EST, A_SE, B_EST, B_SE) %>%
                  arrange(SPECIES, REGION, DISTRICT, YEAR)

  