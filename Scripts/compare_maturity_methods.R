#PURPOSE: 
# To compare different methods for calculating morphometric maturity, currently just for snow crab. 

# Author: Emily Ryznar

# NOTES:
# Make sure you're using updated chela database compiled by shannon, though can't use it for everything bc
# it doesn't have sampling factor...

# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")

# Set params ----
current_year <- 2024
years <- c(1989:2007, 2009, 2011, 2015, 2017:current_year)
years <- c(1989:current_year)


# Load data ----
# minima data
minima <- read.csv("./Output/opilio_cutline_minima.csv")

# Chela data compiled by Shannon
sh_chela <- read.csv("./Data/snow_chela_UPDATED.csv") %>% # already != HT 17, only shell 2, no special projects
          filter(YEAR %in% years) %>%
          filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE,
                 YEAR %in% years) %>% # filter for males, sh2, only chela msrd, not HT17
          mutate(ratio = SIZE/CHELA_HEIGHT,
                 LN_CH = log(CHELA_HEIGHT),
                 LN_CW = log(SIZE),
                 CW = SIZE) %>%
          filter(ratio > 2 & ratio < 35) %>% # filter extreme measurements
          dplyr::select(!c(ratio)) %>%
          mutate(cutoff = BETA0 + BETA1*CW, # apply cutline model
                 MATURE = case_when((LN_CH > cutoff) ~ 1,
                                    TRUE ~ 0))


# Crabpack specimen
cp_chela <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
          #dplyr::select(colnames(chela_10.13)) %>%
          filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE,
                 YEAR %in% years) %>% # filter for males, sh2, only chela msrd, not HT17
          mutate(ratio = SIZE/CHELA_HEIGHT,
                 LN_CH = log(CHELA_HEIGHT),
                 LN_CW = log(SIZE)) %>%
          filter(ratio > 2 & ratio < 35) %>% # filter extreme measurements
          dplyr::select(!c(ratio)) %>%
          mutate(cutoff = BETA0 + BETA1*LN_CW, # apply cutline model
                 MATURE = case_when((LN_CH > cutoff) ~ 1,
                                    TRUE ~ 0))

# Plot 1: cutline with 2024 highlighted ----
ggplot()+
  geom_point(sh_chela, mapping = aes(LN_CW, LN_CH, color = as.factor(MATURE)), alpha = 0.5)+
  theme_bw()+
  facet_wrap(~YEAR)+
  ylab("ln(chela height)")+
  xlab("ln(carapace width)")+
  geom_abline(slope = BETA1, intercept = BETA0, linewidth= 1, linetype = "dashed")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave("./Figures/year_facetted_cutline.png", width = 11, height = 8.5)

# Plot 2: legacy weighting of 0s and 1s by sampling factor per 10mm bin and 5mm bin ----
# Subset data into size intervals at ln(CW) and add in cutline params
# 10mm
cp_chela %>%
  mutate(BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),# 10mm
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         MIDPOINT = (UPPER + LOWER)/2,
         SIZE_1MM = floor(SIZE)) %>%
  group_by(YEAR, BIN, MIDPOINT, MATURE) %>%
  # Sum n_crabs per size bin
  reframe(COUNT = sum(SAMPLING_FACTOR)) %>%
  mutate(bin.size = "10mm") ->  bin.10mm

# 5mm
cp_chela %>%
  mutate(BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left"),# 5mm
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         MIDPOINT = (UPPER + LOWER)/2,
         SIZE_1MM = floor(SIZE)) %>%
  group_by(YEAR, BIN, MIDPOINT, MATURE) %>%
  # Sum n_crabs per size bin
  reframe(COUNT = sum(SAMPLING_FACTOR)) %>%
  mutate(bin.size = "5mm") ->  bin.5mm


rbind(bin.10mm, bin.5mm) -> plot.dat


  ggplot(plot.dat %>% filter(bin.size == "10mm"), aes(MIDPOINT, MATURE, size = COUNT)) +
    geom_point(alpha = 0.5)+
    facet_wrap(~YEAR)+
    ggtitle("Weighting method = sampling factor, 10mm bins")+
    scale_size_continuous(name = "Sampling factor")+
    theme_bw()+
     xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/binarymature_SFweight_10mm.png", height = 8.5, width = 11, units = "in")
  
  ggplot(plot.dat %>% filter(bin.size == "5mm"), aes(MIDPOINT, MATURE, size = COUNT)) +
    geom_point(alpha = 0.5, color = "cadetblue")+
    facet_wrap(~YEAR)+
    ggtitle("Weighting method = sampling factor, 5mm bins")+
    scale_size_continuous(name = "Sampling factor")+
    theme_bw()+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/binarymature_SFweight_5mm.png", height = 8.5, width = 11, units = "in")
  

# Plot 3: Proposed weighting of 0s and 1s by CPUE by 1mm bin ----
  # Get survey cpue by 1mm bin for ALL shell 2 males, not just chela msrd
  snow_cpue <- calc_cpue(crab_data = readRDS("./Data/snow_survey_specimenEBS.rda"),
                         species = "SNOW",
                         years = years,
                         sex = "male",
                         shell_condition = "new_hardshell",
                         bin_1mm = TRUE) %>%
               dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, CPUE)
   
  
  sh_chela %>%
    mutate(SIZE_1MM = floor(SIZE),
           MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                                TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
    right_join(., snow_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin 
     mutate(BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),# 10mm
            BIN2 = BIN) %>%
       separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
       mutate(LOWER = as.numeric(sub('.', '', LOWER)),
              UPPER = as.numeric(gsub('.$', '', UPPER)),
              MIDPOINT = (UPPER + LOWER)/2) %>%
       group_by(YEAR, BIN, MIDPOINT, MATURE, MAT_TEXT) %>%
       # Sum n_crabs per size bin
       reframe(CPUE = sum(CPUE)) %>%
       mutate(bin.size = "10mm") ->  bin.10mm
  
  
  sh_chela %>%
    mutate(SIZE_1MM = floor(SIZE),
           MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                                TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
    right_join(., snow_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin 
    mutate(BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left"),# 5mm
           BIN2 = BIN) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           MIDPOINT = (UPPER + LOWER)/2) %>%
    group_by(YEAR, BIN, MIDPOINT, MATURE, MAT_TEXT) %>%
    # Sum n_crabs per size bin
    reframe(CPUE = sum(CPUE)) %>%
    mutate(bin.size = "5mm") ->  bin.5mm
  
  
  
  
  rbind(bin.10mm, bin.5mm) %>%
    na.omit()-> plot.dat
  
  
  plot.dat %>% filter(!YEAR %in% c(2008, 2010, 2013, 2014), bin.size == "10mm") %>% mutate(CPUE = CPUE/1000) -> pd
  ggplot(pd, aes(MIDPOINT, MATURE, size = CPUE)) +
    geom_point(alpha = 0.5)+
    facet_wrap(~YEAR)+
    ggtitle("Weighting method = CPUE, 10mm bins")+
    scale_size_continuous(name = "CPUE\n(thousands)", breaks = seq(500, 2000, by = 500))+
    theme_bw()+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/binarymature_CPUEweight_10mm.png", height = 8.5, width = 11, units = "in")
  
  
  plot.dat %>% filter(!YEAR %in% c(2008, 2010, 2013, 2014), bin.size == "5mm") %>% mutate(CPUE = CPUE/1000) -> pd
  ggplot(pd, aes(MIDPOINT, MATURE, size = CPUE)) +
    geom_point(alpha = 0.5,  color = "cadetblue")+
    facet_wrap(~YEAR)+
    ggtitle("Weighting method = CPUE, 5mm bins")+
    scale_size_continuous(name = "CPUE\n(thousands)", breaks = seq(500, 2000, by = 500))+
    theme_bw()+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/binarymature_CPUEweight_5mm.png", height = 8.5, width = 11, units = "in")
  

  
  
# Plot 4: Proportion mature ogives using 5 and 10mm bins ----
  #10mm
  cp_chela %>%
    filter(YEAR %in% years) %>%
    mutate(BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),# 10mm
           BIN2 = BIN) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           MIDPOINT = (UPPER + LOWER)/2,
           SIZE_1MM = floor(SIZE)) %>%
    group_by(YEAR, MIDPOINT) %>%
    # Sum n_crabs per size bin
    reframe(TOTAL = sum(SAMPLING_FACTOR),
            NUM_MATURE = sum(SAMPLING_FACTOR[MATURE == 1]),
            PROP_MATURE = NUM_MATURE/TOTAL) %>%
    mutate(bin.size = "10mm", 
           PROP_MATURE = case_when((MIDPOINT <=35) ~ 0,
                                   (MIDPOINT >=135) ~ 1,
                                   TRUE ~ PROP_MATURE)) ->  bin.10mm.sf
  
  cp_chela %>%
    mutate(SIZE_1MM = floor(SIZE),
           MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                                TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
    right_join(., snow_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin 
    mutate(BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),# 10mm
           BIN2 = BIN) %>%
    filter(YEAR %in% years) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           MIDPOINT = (UPPER + LOWER)/2) %>%
    group_by(YEAR, MIDPOINT) %>%
    # Sum n_crabs per size bin
    reframe(TOTAL = sum(CPUE),
            NUM_MATURE = sum(CPUE[MATURE == 1]),
            PROP_MATURE = NUM_MATURE/TOTAL) %>%
    mutate(bin.size = "10mm", 
           PROP_MATURE = case_when((MIDPOINT <=35) ~ 0,
                                   (MIDPOINT >=135) ~ 1,
                                   TRUE ~ PROP_MATURE)) ->  bin.10mm.cpue
  
  
 plot.dat <- rbind(bin.10mm.sf %>% mutate(weight = "Sampling factor"), bin.10mm.cpue %>% mutate(weight = "CPUE"))

  ggplot(plot.dat, aes(MIDPOINT, PROP_MATURE, color = weight))+
    geom_point()+
    geom_line()+
    ggtitle("10mm bins")+
    facet_wrap(~YEAR)+
    theme_bw()+
    ylab("Proportion mature")+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/propmature_ogives_10mm.png", height = 8.5, width = 11, units = "in")
  
  
  #5mm
  cp_chela %>%
    filter(YEAR %in% years) %>%
    mutate(BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left"),# 10mm
           BIN2 = BIN) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           MIDPOINT = (UPPER + LOWER)/2,
           SIZE_1MM = floor(SIZE)) %>%
    group_by(YEAR, MIDPOINT) %>%
    # Sum n_crabs per size bin
    reframe(TOTAL = sum(SAMPLING_FACTOR),
            NUM_MATURE = sum(SAMPLING_FACTOR[MATURE == 1]),
            PROP_MATURE = NUM_MATURE/TOTAL) %>%
    mutate(bin.size = "10mm", 
           PROP_MATURE = case_when((MIDPOINT <=35) ~ 0,
                                   (MIDPOINT >=135) ~ 1,
                                   TRUE ~ PROP_MATURE)) ->  bin.5mm.sf
  
  cp_chela %>%
    mutate(SIZE_1MM = floor(SIZE),
           MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                                TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
    right_join(., snow_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin 
    mutate(BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left"),# 10mm
           BIN2 = BIN) %>%
    filter(YEAR %in% years) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           MIDPOINT = (UPPER + LOWER)/2) %>%
    group_by(YEAR, MIDPOINT) %>%
    # Sum n_crabs per size bin
    reframe(TOTAL = sum(CPUE),
            NUM_MATURE = sum(CPUE[MATURE == 1]),
            PROP_MATURE = NUM_MATURE/TOTAL) %>%
    mutate(bin.size = "10mm", 
           PROP_MATURE = case_when((MIDPOINT <=35) ~ 0,
                                   (MIDPOINT >=135) ~ 1,
                                   TRUE ~ PROP_MATURE)) ->  bin.5mm.cpue
  
  
  plot.dat <- rbind(bin.5mm.sf %>% mutate(weight = "Sampling factor"), bin.5mm.cpue %>% mutate(weight = "CPUE"))
  
  ggplot(plot.dat, aes(MIDPOINT, PROP_MATURE, color = weight))+
    geom_point()+
    geom_line()+
    ggtitle("5mm bins")+
    facet_wrap(~YEAR)+
    theme_bw()+
    ylab("Proportion mature")+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/propmature_ogives_5mm.png", height = 8.5, width = 11, units = "in")
  
  
  plot.dat <- rbind(bin.5mm.cpue %>% mutate(bin.size = "5mm"), bin.10mm.cpue)
  
  ggplot(plot.dat, aes(MIDPOINT, PROP_MATURE, color = bin.size))+
    geom_point()+
    geom_line()+
    ggtitle("Proportion mature at size, CPUE weighting")+
    facet_wrap(~YEAR)+
    theme_bw()+
    ylab("Proportion mature")+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/propmature_ogives_cpuebins.png", height = 8.5, width = 11, units = "in")
  

# Plot 5: Proportion mature ogives unweighted  ----
  #10mm
  sh_chela %>%
    filter(YEAR %in% years) %>%
    mutate(BIN = cut_width(SIZE, width = 10, center = 5, closed = "left"),# 10mm
           BIN2 = BIN) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           MIDPOINT = (UPPER + LOWER)/2,
           SIZE_1MM = floor(SIZE)) %>%
    group_by(YEAR, MIDPOINT) %>%
    # Sum n_crabs per size bin
    reframe(TOTAL = n(),
            NUM_MATURE = sum(MATURE == 1),
            PROP_MATURE = NUM_MATURE/TOTAL) %>%
    mutate(bin.size = "10mm", 
           PROP_MATURE = case_when((MIDPOINT <=35) ~ 0,
                                   (MIDPOINT >=135) ~ 1,
                                   TRUE ~ PROP_MATURE)) ->  bin.10mm.unweighted
  
  #5mm
  sh_chela %>%
    filter(YEAR %in% years) %>%
    mutate(BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left"),# 10mm
           BIN2 = BIN) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           MIDPOINT = (UPPER + LOWER)/2,
           SIZE_1MM = floor(SIZE)) %>%
    group_by(YEAR, MIDPOINT) %>%
    # Sum n_crabs per size bin
    reframe(TOTAL = n(),
            NUM_MATURE = sum(MATURE == 1),
            PROP_MATURE = NUM_MATURE/TOTAL) %>%
    mutate(bin.size = "10mm", 
           PROP_MATURE = case_when((MIDPOINT <=35) ~ 0,
                                   (MIDPOINT >=135) ~ 1,
                                   TRUE ~ PROP_MATURE)) ->  bin.5mm.unweighted
  
  
  #1mm
  sh_chela %>%
    filter(YEAR %in% years) %>%
    mutate(BIN = cut_width(SIZE, width = 1, center = 0.5, closed = "left"),# 10mm
           BIN2 = BIN) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           MIDPOINT = (UPPER + LOWER)/2,
           SIZE_1MM = floor(SIZE)) %>%
    group_by(YEAR, MIDPOINT) %>%
    # Sum n_crabs per size bin
    reframe(TOTAL = n(),
            NUM_MATURE = sum(MATURE == 1),
            PROP_MATURE = NUM_MATURE/TOTAL) %>%
    mutate(bin.size = "10mm", 
           PROP_MATURE = case_when((MIDPOINT <=35) ~ 0,
                                   (MIDPOINT >=135) ~ 1,
                                   TRUE ~ PROP_MATURE)) ->  bin.1mm.unweighted
 
   
  ggplot(bin.1mm.unweighted, aes(MIDPOINT, PROP_MATURE))+
    geom_point()+
    geom_line()+
    ggtitle("10mm bins")+
    facet_wrap(~YEAR)+
    theme_bw()+
    ylab("Proportion mature")+
    xlab("Bin width midpoint (mm)")
  
  ggsave("./Figures/propmature_ogives_10mm.png", height = 8.5, width = 11, units = "in")
  
  rbind(bin.10mm.unweighted %>% mutate(BIN_SIZE = "10mm"), 
        bin.5mm.unweighted %>% mutate(BIN_SIZE = "5mm"), 
        bin.1mm.unweighted %>% mutate(BIN_SIZE = "1mm")) %>%
    dplyr::select(!bin.size) -> unweighted.propmature
  
write.csv(unweighted.propmature, "./Output/snow_unweighted_propmature.csv")
# Plot 5: Sample size by year ----
  cp_chela %>%
    filter(YEAR %in% years)  %>%
    group_by(YEAR) %>%
    reframe(N = n()) -> plot.dat
  
  
 ggplot(plot.dat, aes(YEAR, N))+
   geom_bar(stat = "identity")+
   theme_bw()+
   scale_x_continuous(breaks = seq(1989, 2024, by = 2))+
   ylab("# chela measurements") +
   xlab("Year")
 
 ggsave("./Figures/chela_yearlysamplesize.png", height = 5, width = 7, units = "in")
 
  