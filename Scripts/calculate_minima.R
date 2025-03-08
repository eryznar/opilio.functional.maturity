# PURPOSE: to calculate minima between dominant distribution modes of opilio chela height. Minima will serve as 
# the distribution cutline to define morphometric maturity. 

# Author: Emily Ryznar, Jon Richar, Shannon Hennessey

# Load libs/params
source("./Scripts/load_libs_params.R")

# Load Jon chela data to get 2010 and 2013
jon.dat <- read.csv("./Data/opilio_chela_height_TS.csv") %>%
          mutate(YEAR = substr(CRUISE, 1, 4)) %>%
          dplyr::filter(YEAR %in% c(2010, 2013), 
                        #HAUL_TYPE !=17, need to make sure this was filtered by Jon
                        SEX == 1, 
                        SHELL_CONDITION == 2) %>%
  dplyr::select(HAULJOIN, YEAR, WIDTH, CHELA_HEIGHT, SAMPLING_FACTOR) %>%
  rename(SIZE = WIDTH)

# Obtain CH and CW measurements from shell 2 males from survey data
chela <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
          filter(SEX == 1,
                 SHELL_CONDITION == 2, 
                 HAUL_TYPE !=17, 
                 #YEAR %in% c(2010, 2017, 2018, 2019) # years Jon uses for stock assessment estimates?
                 is.na(CHELA_HEIGHT) == FALSE) %>%
  dplyr::select(HAULJOIN, YEAR, SIZE, CHELA_HEIGHT, SAMPLING_FACTOR)

# Join survey data and Jon's data from 2010 and 2013
chela <- rbind(jon.dat, chela)

# Pool and transform data using natural log
chela %>%
  mutate(LN_CW = log(SIZE),
         LN_CH = log(CHELA_HEIGHT)) -> chela2

# Subset data into size intervals at ln(CW) of 0.025
chela2 %>%
  filter(LN_CW >= 3.9 & LN_CW <= 4.6) %>% # these were filtered out by Jon
    mutate(BIN = cut_width(LN_CW, width = 0.025, center = 0.0125, closed = "left", dig.lab = 4),
           BIN2 = BIN) %>%
      separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
      mutate(LOWER = as.numeric(sub('.', '', LOWER)),
             UPPER = as.numeric(gsub('.$', '', UPPER)),
             MIDPOINT = (UPPER + LOWER)/2) -> bin.dat

# Sequentially apply KDE to the ln-chela height data for each interval, and identify minima of the resulting
# density distributions to define maturity classes in a given interval

bins <- unique(bin.dat$BIN)

min.dat <- data.frame()

for(ii in 1:length(bins)){
  # filter data by bin of interest
  bin.dat %>%
    filter(BIN == bins[ii])-> opt.dat
  
  # calculate density for chela heights within that bin
  d <- density(opt.dat$LN_CH, kernel = "gaussian")
  
  # Identify local minima and maxima:
  local_min_max <- function(x) {
    signs <- sign(diff(x))
    list(
      minima = which(diff(signs) > 0) + 1,
      maxima = which(diff(signs) < 0) + 1
    )
  }
  
  extrema <- local_min_max(d$y) # both local minima and maxima
 
  maxima <- d$x[extrema$maxima] # isolate maxima index locations
  d.max <- d$y[extrema$maxima] # identify density at maxima index locations
  
  cbind(maxima, d.max) %>%
    as.data.frame() %>%
    slice_max(d.max, n = 2) -> max.df # pull out maxima locations with top two greatest densities
  

  # Identify the two dominant modes based on maxima to calculate minimum outside distribution tails
  ints <- round(max.df$maxima, 1)
  ints <- ints[order(ints)]
  ints[2] <- ints[1]+0.2
 
  # find minimum in between dominant modes
  min <- optimize(approxfun(d$x, d$y), interval = ints)$minimum
  
  # create dataframe for output
  out <- data.frame(lwr = ints[1], 
                    upr = ints[2],
                    minima = min,
                    midpoint = opt.dat$MIDPOINT,
                    bin = bins[ii]) %>%
    distinct()
  
  min.dat <- rbind(min.dat, out)
  
}

write.csv(min.dat, "./Output/opilio_cutline_minima.csv")

# Plot to make sure minima were calculated correctly
right_join(bin.dat, min.dat %>% rename(BIN = bin)) -> plot.dat

ggplot()+
  geom_density(plot.dat, mapping = aes(LN_CH), linewidth = 1)+
  geom_vline(plot.dat, mapping = aes(xintercept = minima), color = "blue", linetype = "dashed", linewidth = 1)+
  facet_wrap(~BIN, scales = "free_x")+
  theme_bw()+
  ylab("Density")+
  xlab("ln(chela height)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> dens.plot

ggsave("./Figures/density_plots.png", width = 11, height = 8.5)

ggplot()+
  geom_histogram(plot.dat, mapping = aes(x = LN_CH), bins = 30, fill = "white", color = "black")+
  #geom_density(plot.dat, mapping = aes(LN_CH), color = "red", linetype = "dashed", linewidth = 1)+
  geom_vline(plot.dat, mapping = aes(xintercept = minima), color = "blue", linewidth = 1)+
  facet_wrap(~BIN, scales = "free_x")+
  theme_bw()+
  ylab("Density")+
  xlab("ln(chela height)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> hist.plot

# Fit linear model to evaluate ln(CH) minima (imm/mat division) and ln(CW) bin midpoint
mod <- lm(minima ~ midpoint, min.dat)

prd <- predict(mod, data = min.dat, se.fit=TRUE, interval="confidence", level=0.95)

cbind(min.dat %>% dplyr::select(!c(lwr, upr)), prd$fit) -> min.dat2

labs <- data.frame(x = 4, y = 3, lab = paste0("R-squared = ", round(summary(mod)$r.squared, 2), "\np < 0.001"))
                                              

ggplot()+
  geom_ribbon(min.dat2, mapping = aes(x = midpoint, ymin = lwr, ymax = upr), alpha = 0.4, fill = "cadetblue")+
  geom_line(min.dat2, mapping = aes(midpoint, fit), linewidth = 1, color = "cadetblue")+
  geom_point(min.dat2, mapping = aes(midpoint, minima))+
  #geom_density(plot.dat, mapping = aes(LN_CH), color = "red", linetype = "dashed", linewidth = 1)+
  theme_bw()+
  annotate("text", x = 4, y = 3, label = paste0("R-squared = ", round(summary(mod)$r.squared, 2), "\np < 0.001"))+
  ylab("Cutline (ln(chela height))")+
  xlab("Bin (ln(carapace width))")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> cut.plot

ggsave("./Figures/cutline_lm_plot.png", width = 6, height = 5)
