# tutorial from https://github.com/coreysparks/DEM7263/blob/main/lectures/Lecture_5_GLMs_for_Spatial_Data.Rmd

library(haven)

url <- "/Users/agnesguo/Downloads/AHRF_2020-2021_SAS/AHRF2021.sas7bdat"

arf2020 <- read_sas(url)

library(spdep)
library(MASS)
library(spatialreg)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tigris)
library(sf)

# For our analysis, we will use the number of low birth weight infants and the total number of births as our numerator and denominator.
# 3-Year total birth (2017-2019): f1254617
# 3-Year total birth (2008-2010): f1254608
# 3-Year total low birth weight birth (2017-2019): f1255317
# 3-Year total birth (2008-2010): f1255308
# Percent Persons Age 0-17 in Pov (2019): f1332219
# the rural urban continuum code from USDA: f0002013
# the primary healthcare shortage area code from 2019: f0978719
# the per capital number of OB-GYN's in the county: f1168419/f1198419

# predictors for the analysis: child poverty rate from 2010, the rural urban continuum code from USDA, the primary healthcare shortage area code from 2010, 
# and the per capital number of OB-GYN's in the county. 
# Then we filter to have non-missing cases, which reduces the number of counties from 3,230 to 2,262.

pop_dat <-arf2020%>%
  mutate(cofips=f00004, 
         coname=f00010,
         state = f00011,
         births1416=f1254617,
         births0608=f1254608,
         lowbw1416=f1255317,
         lowbw0608=f1255308,
         childpov10= f1332210,
         rucc= as.factor(f0002013),
         hpsa10= as.factor(f0978719),
         obgyn10_pc= 1000*(f1168419/ f0453010),
         single_hh = 100* (f1160310/ f0453010),
         blk_pc =  f0453810,
         his_pc =  f0454210,
         unemploy =  f0679517 ,
         labor = f1451011,
         work_in = f1460611,
         lowedu =  1397515 ) %>%
  dplyr:: select(births1416, lowbw1416,births0608, lowbw0608,state, cofips, coname, childpov10, rucc, hpsa10, obgyn10_pc, single_hh, blk_pc, his_pc, unemploy, 
                 labor, work_in, lowedu)%>%
  filter(complete.cases(.))%>%
  as.data.frame()

filter(state %in% c("09", "17", "18", "19", "23", "25", "26", "27", "33", "34", "36", "39", "42", "44", "50", "55")) %>%
head(pop_dat)
summary(pop_dat)

options(tigris_class="sf")
usco<-counties(cb=T, year= 2019)
usco$cofips <-usco$GEOID
# focus on North US: Connecticut (09), Illinois(17), Indiana(18), Iowa(19), Maine(23), Massachusetts(25), Michigan(26), Minnesota(27),
# New Hampshire(33), New Jersey(34), New York(36), Ohio(39), Pennsylvania(42), Rhode Island(44), Vermont(50), Wisconsin(55)
usco_sub <- usco %>% filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78"))
ggplot(usco_sub) + geom_sf()

sts<-states(cb = T, year=2019)
sts<-st_boundary(sts)%>%
  filter(STATEFP %in% c("09", "17", "18", "19", "23", "25", "26", "27", "33", "34", "36", "39", "42", "44", "50", "55"))
ggplot(sts) + geom_sf()

pop_dat <- geo_join(usco_sub, pop_dat, by_sp = "cofips", by_df = "cofips", how = "left")

# exploratory analysis of spatial clustering
pop_dat$low_percent <- 100 *(pop_dat$lowbw1416/pop_dat$births1416)

# construct a k = 4 nearest neighbor list for low birth weight rate
knn <- knearneigh(coordinates(as_Spatial(pop_dat)), k = 4, longlat = T)
knn <- knn2nb(knn)
plot(st_geometry(pop_dat), main = "North US, K-Nearest Neighbors")
plot(knn, coordinates(as_Spatial(pop_dat)), col = "red", pch = ".", lwd = 0.5, add = T)

# low_percent map
ggplot(pop_dat) + 
  geom_sf(aes(fill = childpov10)) +
  scale_fill_viridis_c() 
  

0-pop_dat %>%
  filter(STATEFP %in% c("09", "17", "18", "19", "23", "25", "26", "27", "33", "34", "36", "39", "42", "44", "50", "55")) %>%
  mutate(lbrate=lowbw1416/births1416)%>%
  mutate(lb_group = cut(lbrate, breaks=quantile(lbrate, p=seq(0,1,length.out = 6), na.rm=T ), include.lowest=T ))%>%
  ggplot()+
  geom_sf(aes(fill= lb_group, color=NA))+
  scale_color_brewer(palette = "Blues")+
  scale_fill_brewer(palette = "Blues",na.value = "grey50")+
  geom_sf(data=sts, color="black")+
  coord_sf(crs = 2163) +
  ggtitle(label = "Proportion of births that were low birthweight, 2017-2019")


# Poisson regression
sub_dat <- filter(pop_dat, is.na(lowbw1416)==F)
fit_pois<- glm(lowbw1416 ~ offset(log(births1416+.00001)) + hpsa10, 
               family=poisson, 
               data=sub_dat)
summary(fit_pois)

exp(coef(fit_pois))


