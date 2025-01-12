# child poverty analysis
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

pop_dat <-arf2020%>%
  mutate(cofips=f00004, 
         coname=f00010,
         state = f00011,
         childpov10=  f1332219 ,
         rucc= as.factor(f0002013),
         single_hh = 100* (f1160310/ f0453010),
         blk_pc =  f0453810,
         his_pc =  f0454210,
         unemploy =   (f1451215/f1451015)* 100 ,
         unemploy_male  = ((f1451215 - f1451515)/(f1451015 - f1451315 )) * 100,
         work_in = f1460615, 
         extr =  f1462115,
         mnf =  f1458715, 
         median_income =f1434615) %>%
  dplyr:: select(state, cofips, coname, childpov10, rucc, single_hh, blk_pc, his_pc, unemploy, unemploy_male, work_in, extr, mnf, median_income)%>%
  filter(complete.cases(.))%>%
  filter(!state %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  as.data.frame()

pop_dat$metro_cut<-cut(as.numeric(pop_dat$rucc),
                    breaks =c(0,6,10),
                    include.lowest = T,
                    labels = c("Metropolitan", "Non-Metropolitan"))

head(pop_dat)
summary(pop_dat)

# import data from ACS 2015
library(tidycensus)
census_api_key('7db876e1da66a89dd43fa1040e25b844c058760e')
edu_table <- get_acs(geography = "county", year = 2015, geometry = T, output = "wide",
                        table = "B06009", cache_table = T)

edudat <- edu_table %>%
  mutate(less_high = ((B06009_002E + B06009_003E)/B06009_001E) * 100, 
         cofips = substr(GEOID, 1, 5))%>%
  dplyr:: select(less_high, cofips)

dat <- geo_join(edudat, pop_dat, by_sp = "cofips", by_df = "cofips", how = "inner")
summary(dat)

ggplot(dat) + 
  geom_sf(aes(fill = childpov10, color = childpov10)) +
  scale_fill_viridis_c()+
  scale_color_viridis_c() +
  ggtitle("Spatial Variation of County-Level Child Poverty Rate, 2019") +
  theme(plot.title = element_text(size = 10))

dat$pov_cut<-cut(as.numeric(dat$childpov10),
                       breaks =c(0,15,30,100),
                       include.lowest = T,
                       labels = c("Less than 15", "15 to 30", "Over 30"))

ggplot(dat) + 
  geom_sf(aes(fill = factor(pov_cut)), lwd = 0.1) +
  scale_fill_manual(values = c("white", "grey", "blue"), name= "% Child Poverty")+ 
  ggtitle("Spatial Variation of County-Level Child Poverty Rate, 2019") +
  theme(plot.title = element_text(size = 10), legend.text = element_text(size= 7), legend.position="bottom")

ggplot(dat) +
  geom_sf(aes(fill = factor(metro_cut)), lwd = 0.1) +
  scale_fill_manual(values = c("lightblue","white"), name= "County Type")

# linear regression model
fit <- lm(childpov10 ~ single_hh + blk_pc + his_pc + unemploy + unemploy_male + work_in + extr +  mnf + less_high + metro_cut+ log(median_income), data = dat)
summary(fit)

# Test of normality
ggplot(fit) + 
  stat_qq(aes(sample = .stdresid)) + 
  geom_abline() + 
  ggtitle(label = "Normal Q-Q plot for regression residuals") +  theme(plot.title = element_text(size = 10))

# the Shapiro-Wilk test
shapiro.test(fit$residuals)

# Test of spatial autocorrelation
# construct the k = 4 nearest neighbor list
knn <- knearneigh(coordinates(as_Spatial(dat)), k = 4)
knn <- knn2nb(knn)
plot(st_geometry(dat), main = "US Counties, K-Nearest Neighbors")
plot(knn, coordinates(as_Spatial(dat)), col="red", pch=".",lwd=0.5, add=T)

# spatial autocorrelation test
data.wts <- nb2listw(knn, style = "W")
lm.morantest(fit, listw = data.wts)

moran.test(x = dat$childpov10[!is.na(dat$childpov10)], listw = data.wts)
# local moran statistics
local <- localmoran(x = dat$childpov10[!is.na(dat$childpov10)], listw = data.wts, alternative = "two.sided")

data <- as_Spatial(dat)

data$local <- local[,1]
data$localp <- local[,5]
data$sinc <- scale(data$childpov10)
data$inc <- lag.listw(var = data$sinc, x = data.wts)
data_quad_sig <- NA
data_quad_sig[(data$sinc >= 0 & data$inc >= 0) & (data$localp <= 0.05)] <- 1
data_quad_sig[(data$sinc <= 0 & data$inc <= 0) & (data$localp <= 0.05)] <- 2
data_quad_sig[(data$sinc >= 0 & data$inc <= 0) & (data$localp <= 0.05)] <- 3
data_quad_sig[(data$sinc <= 0 & data$inc >= 0) & (data$localp <= 0.05)] <- 4
data_quad_sig[data$localp > 0.05] <- 5


bks <- seq(1, 5, 1)
labels <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not Clustered")
np <- findInterval(data_quad_sig, bks)


colors <- c("red", "blue", "lightpink", "lightblue", "white")

plot(data, col = colors[np], main = "Moran LISA Cluster Map - Child Poverty Rate")
legend("bottomleft", legend = labels,fill = colors, bty = "n")


# Spatial regression Model
library(spatialreg)
# spatial error model
fit.err <- errorsarlm(childpov10 ~ single_hh + blk_pc + his_pc + unemploy + unemploy_male + work_in + extr +  mnf + less_high + metro_cut+ log(median_income), 
                      data = dat, listw = data.wts)
summary(fit.err)

# Spatial lag model
fit.lag <- lagsarlm(childpov10 ~ single_hh + blk_pc + his_pc + unemploy + unemploy_male + work_in + extr +  mnf + less_high + metro_cut+ log(median_income), 
                    data = dat, listw = data.wts, type = "lag")
summary(fit.lag)


# Spatial Regime Model
efit<-errorsarlm(childpov10 ~ single_hh + blk_pc + his_pc + unemploy + unemploy_male + work_in + extr +  mnf + less_high + log(median_income), 
                 data = dat, listw = data.wts, method = "MC")

efit.0 <-errorsarlm(childpov10 ~ metro_cut/(single_hh + blk_pc + his_pc + unemploy + unemploy_male + work_in + extr +  mnf + less_high + log(median_income)), 
                 data = dat, listw = data.wts, method = "MC")
anova(efit, efit.0)

# neighbor list for metro/non-metro areas
#k = 4  nearest neighbors for the metropolitan counties

knn_m<-knearneigh(coordinates(as_Spatial(dat[dat$metro_cut =="Metropolitan",])), k=4)
knn_m<-knn2nb(knn_m, sym=T)
wts_m<-nb2listw(knn_m)

#k = 4 nearest neighbors for the nonmetropolitan counties
knn_n<-knearneigh(coordinates(as_Spatial(dat[dat$metro_cut=="Non-Metropolitan",])), k=4)
knn_n<-knn2nb(knn_n, sym=T)
wts_n<-nb2listw(knn_n)

efit.1<- errorsarlm(childpov10 ~ (single_hh + blk_pc + his_pc + unemploy + unemploy_male + work_in + extr +  mnf + less_high + log(median_income)), 
                       data = dat[dat$metro_cut == "Metropolitan",], 
                       listw = wts_m, method = "MC")
summary(efit.1, Nagelkerke=T)

efit.2 <- errorsarlm(childpov10 ~ (single_hh + blk_pc + his_pc + unemploy + unemploy_male + work_in + extr +  mnf + less_high + log(median_income)), 
                    data = dat[dat$metro_cut == "Non-Metropolitan",], 
                    listw = wts_n, method = "MC")
summary(efit.2, Nagelkerke=T)


# export result
library(stargazer)

# descriptive statistics
sum_table <- dat %>% 
  as.data.frame() %>%
  mutate(logincome = log(median_income), metro = ifelse(metro_cut == "Metropolitan", 1, 0)) %>%
  dplyr::select(childpov10, metro, single_hh, blk_pc, his_pc, unemploy, unemploy_male, work_in, extr, mnf, less_high, logincome)

stargazer(sum_table, omit.summary.stat = c("p25", "p75"))

sum_table1 <- dat %>%
  as.data.frame() %>%
  mutate(logincome = log(median_income)) %>%
  filter(metro_cut == "Metropolitan")

stargazer(fit, fit.err, fit.lag)
stargazer(efit.1, efit.2)
stargazer(efit.2)
