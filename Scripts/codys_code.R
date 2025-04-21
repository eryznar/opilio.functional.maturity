#survDAT<-read.csv("C:/Users/cody.szuwalski/Work/snow_2023_9/data/survey/EBSCrab_Haul/EBSCrab_Haul.csv",header=T,skip=5)
library(dplyr)
library(mgcv)
library(ggplot2)
library("rnaturalearth")
library(patchwork)

survDAT <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(#YEAR %in% c(2010, 2017, 2018, 2019) # years Jon uses for stock assessment estimates?,
    HAUL_TYPE !=17)


#==maturity coef # -2.8434, 0.2404 #cody's
survDAT$cutoff<- BETA0+BETA1*log(survDAT$SIZE_1MM)
survDAT$mature<-log(survDAT$CHELA_HEIGHT)>survDAT$cutoff

mod_dat<-survDAT %>%
  dplyr::select(YEAR, LATITUDE, LONGITUDE, SEX, SIZE, SAMPLING_FACTOR, mature)
mod_dat<-mod_dat[complete.cases(mod_dat),]
colnames(mod_dat)<-c("year","lat","lon","sex","size","sample","mature")

mod<-bam(data=filter(mod_dat,sex==1&size>40),mature~as.factor(year)+s(size)+s(lon,lat),family=binomial(link="logit"),weights=sample)
summary(mod)

year_effect<-data.frame(effect=c(0,summary(mod)$p.coeff),
                        year=unique(mod_dat$year))

yr_eff<-ggplot(year_effect)+
  geom_point(aes(x=year,y=effect),size=3)+theme_bw()+
  geom_abline(aes(slope=0,intercept=0),lty=2)

newdat<-data.frame(year=1989,
                   size=seq(40,140,5),
                   lon=weighted.mean(x=mod_dat$lon,w=mod_dat$sample),
                   lat=weighted.mean(x=mod_dat$lat,w=mod_dat$sample))

newdat$preds<-predict(mod,newdata=newdat,type='response',se.fit=TRUE)$fit
size_effect<-ggplot(newdat)+
  geom_line(aes(x=size,y=preds),lwd=1.2)+
  theme_bw()+expand_limits(y=0)+ylab("p(terminal_molt)")

library(maps)  # for map boundaries
# Set typical values for the other covariates:
typical_size <- 75
typical_year <- 1989

# Define grid boundaries using the range of your data for lat and long.
lat_range <- range(mod_dat$lat, na.rm = TRUE)
lon_range <- range(mod_dat$lon, na.rm = TRUE)

# Create a grid of lat and lon values
grid_data <- expand.grid(
  lat = seq(lat_range[1], lat_range[2], length.out = 200),
  lon = seq(lon_range[1], lon_range[2], length.out = 200)
)
# Add the typical values for the other covariates.
grid_data$size <- typical_size
grid_data$year <- typical_year

# Get predicted probability of molting on the response scale
grid_data$pred <- predict(mod, newdata = grid_data, type = "response")

world <- ne_countries(scale = "medium", returnclass = "sf")

p <- ggplot() +
  geom_raster(data = grid_data, aes(x = lon, y = lat, fill = pred)) +
  scale_fill_gradient(low = "blue", high = "red") +
  geom_sf(data=world) +
  coord_sf(xlim = c(-180,-160), ylim = c(50,65), expand = FALSE) +
  theme_bw()

yr_eff | size_effect | p
