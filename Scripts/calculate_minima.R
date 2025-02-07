# CH and CW measured to the enarest 1mm from 1990-2007 and 2009, and to the nearest 0.1mm in 2008 and 2010-2019

# Kernel density estimator to define distributions in CH to distinguish crab that 
# have and have not undergone terminal molt using density function

# 1) CH and CW measurements were obtained from males of all shell conditions from 2008, 2010, 2012, 2014, 2016, 2017
chela <- readRDS("./Data/tanner_survey_specimenEBS.rda")$specimen %>%
  filter(SEX == 1, YEAR %in% c(2008, 2010, 2012, 2014, 2016, 2017))

chela %>%
  filter(is.na(CHELA_HEIGHT) == FALSE) %>%
  group_by(YEAR) %>%
  reframe(N = n()) -> tt

# 2) Data were pooled and linearized using natural log
chela %>%
  mutate(LN_CW = log(SIZE),
         LN_CH = log(CHELA_HEIGHT)) -> chela2

# 3) Data were then subset into 24 ln(CW) size intervals at ln(CW) of 0.025
chela2 %>%
  filter(is.na(LN_CH) == FALSE, LN_CW >= 4.3 & LN_CW <= 4.925) %>% # these were filtered out by Jon
    mutate(BIN = cut_width(LN_CW, width = 0.025, center = 0.0125, closed = "left", dig.lab = 4),
           BIN2 = BIN) %>%
      separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
      mutate(LOWER = as.numeric(sub('.', '', LOWER)),
             UPPER = as.numeric(gsub('.$', '', UPPER)),
             MIDPOINT = (UPPER + LOWER)/2) -> bin.dat

# 4) Then sequentially applied KDE to the ln-chela height data for each interval, and the minima of the resulting
#    density distributions was used to define maturity classes in a given interval

bins <- unique(bin.dat$BIN)
bins <- "[4.525,4.55)"

min.dat <- data.frame()

for(ii in 1:length(bins)){
  # filter data by bin of interest
  bin.dat %>%
    filter(BIN == bins[1])-> opt.dat
  
  # calculate density for chela heights within that bin
  d <- density(opt.dat$LN_CH)
  
  # Find local maxima (modes)
  max <- which(diff(sign(diff(d$y))) == -2) + 1
  
  # Identify the two dominant modes to calculate minimum outside distribution tails
  ints <- round(sort(d$x[max])[1:2], 1)
  ints[2] <- ints[1] + 0.2
  
  # find minimum
  min <- optimize(approxfun(d$x, d$y), interval = ints)$minimum
  
  # create dataframe
  out <- data.frame(lwr = ints[1], 
                    upr = ints[2],
                    minima = min,
                    midpoint = opt.dat$MIDPOINT,
                    bin = bins[ii]) %>%
    distinct()
  
  min.dat <- rbind(min.dat, out)
  
}

min.dat <- min.dat[order(min.dat$midpoint),]

# 5) Two dominant modes within each extracted distribution were used to set boundaries for the region within which
#    the search algorithm would seek minima to prevent IDing minima on distribution tails
# 6) The cutline delineating the two maturity classes was then estimated as the best-fit linear regression of
#    the minima in the ln(CH) distributions (the division between immature/mature crabs within a CW interval) against
#    the midpoints of those CW intervals