#setwd("~/custom_folder")

library(readxl)
require(tibble)
require(dplyr)
require(rstan)
require(rstanarm)
require(ggmcmc)
require(MASS)
require(ciTools) 
require(DescTools)
source('rank_by_distance2.R')   # Custom ranking function
source('functions.R')   # Custom function

set.seed(1098);
################################################################################
# This script is responsible for running Bayesian regression model for wastewater
# concentration as response variable and estimated prevalence as explaatory varialbe    
#
# Written By: Boseung Choi
# Last Edited: Nov 10, 2022
#
################################################################################

result.path = "results/"


site_names=c('MSD01', 'MSD02','MSD03','MSD04','MSD05'); zone = 'msd_agg' #'aggregate'
#site_names=c('MSD01'); zone = 'msd_01' #zone_1' 
#site_names=c('MSD02'); zone = 'msd_02' #zone_2'
#site_names=c('MSD03','MSD04','MSD05'); zone = 'msd_0345' #zone_345'


##### population data       
pop = as.data.frame(c(349850,295910,55928,32460,31269))
site = as.data.frame(c("MSD01","MSD02","MSD03","MSD04","MSD05"))
dat_pop =  cbind(pop,site); colnames(dat_pop) <-c("pop","site") 
dat_pop <- dat_pop %>%
  filter(site %in% site_names)
site_pop = sum(unique(dat_pop$pop))       


################ prevalence estimation result 

dat_pred.prev <- read.csv(paste0(result.path,'prevalence_estimate_',zone,'.csv'))%>%
  as_tibble() %>%
  transmute(date = as.Date(date),
            prev = as.numeric(med),
            prev.top = as.numeric(top),
            prev.bot = as.numeric(bottom),
            pop = site_pop,
            prev.cnt = round(as.numeric(med) * site_pop),
            prev.cnt.top = round(as.numeric(top) * site_pop),
            prev.cnt.bot = round(as.numeric(bottom) * site_pop)
            )%>%
  group_by(date) %>%
  summarize(prev = mean(prev), pop =pop, prev.cnt =round(mean(prev.cnt)),
            prev.top = mean(prev.top), prev.cnt.top =round(mean(prev.cnt.top)),
            prev.bot = mean(prev.bot), prev.cnt.bot =round(mean(prev.cnt.bot))
            ) %>%
  distinct() 

min_date = min(dat_pred.prev$date)
max_date = max(dat_pred.prev$date); 

time_point=as.numeric(dat_pred.prev$date-min(dat_pred.prev$date))
dat_pred.prev=add_column(dat_pred.prev,time_point)
week = round(dat_pred.prev$time_point/7)
dat_pred.prev = add_column(dat_pred.prev, week)



######################### Wastewater  
dat_water <-read.csv(file = "data/WW_data.csv") %>%
  transmute(date = as.Date(Date,format="%m/%d/%y"),
            l_avg_n1 = log(avg_n1), 
            sample_no = as.numeric(Sample_No), 
            site = as.character(Site),
            qdi = as.numeric(QDI_n1),
            avg_PMMoV = as.numeric(avg_PMMoV),
            qdi_PMMoV = as.numeric(QDI_PMMoV),
            avg_n1 = as.numeric(avg_n1),
            avg_n1_m = as.numeric(avg_n1)/as.numeric(avg_PMMoV),
            type="water")%>%
  filter(site %in% site_names) %>%
  filter(date<=max_date & date >=min_date) %>% 
  group_by(date) %>%
  summarize(avg_n1 = mean(avg_n1), l_avg_n1 = mean(l_avg_n1), avg_PMMoV = mean(avg_PMMoV),
            avg_n1_m = mean(avg_n1_m)) %>%
  distinct() 

dat_water<-dat_water[order(as.numeric(dat_water$date)),]

time_point=as.numeric(dat_water$date-min_date)
dat_water=add_column(dat_water,time_point)

week = round(dat_water$time_point/7)
dat_water = add_column(dat_water, week)

######################### hospitalization data
hospital <- read_excel("data/Hospitalization_data.xlsx")

dat_host <-hospital %>%
  mutate(date = as.Date(date)) %>%
  filter(date<=max_date & date >=min_date) %>% 
  distinct() 

time_point=as.numeric(dat_host$date-min_date)
dat_host=add_column(dat_host,time_point)
week = round(dat_host$time_point/7)
dat_host = add_column(dat_host, week)

#########################  sero prevalence testing data   
dat_obs <-read.csv(paste0(result.path,'data_',zone,'.csv'))%>%
  transmute(date = as.Date(date),
            num_antibody_positive = as.numeric(num_antibody_positive),
            num_antibody_total = as.numeric(num_antibody_total), 
            percent_positive = as.numeric(percent_positive)
  )%>%
  filter(date<=max_date & date >=min_date) 

time_point=as.numeric(dat_obs$date-min_date)
dat_obs=add_column(dat_obs,time_point)

week = round(dat_obs$time_point/7)
dat_obs = add_column(dat_obs, week)


####### observed incidence
dat_inc <- read.csv('data/Weekly_Inc.csv') %>%
  as_tibble() %>%
  transmute(date = as.Date(Date),
            inc = as.numeric(Weekly_Incidence),   # inc per 100 K 
            site  = as.character(Shed_Site), 
            pop = as.numeric(Pop) ) %>%
  filter(site %in% site_names) %>%
  filter(date<=max_date & date >=min_date)

dat_inc<-dat_inc[order(as.numeric(dat_inc$date)),]

####### Estimated incidence
dat_pred.inc <- read.csv(paste0(result.path,'incidence_estimate_',zone,'.csv'))%>%
  as_tibble() %>%
  transmute(date = as.Date(date),
            inc = as.numeric(med),
            inc.top = as.numeric(top),
            inc.bot = as.numeric(bottom),
            pop = site_pop ) %>%
  filter(date<=max_date & date >=min_date) %>%
  group_by(date) %>%
  summarize(inc = sum(inc), inc.top = sum(inc.top),inc.bot = sum(inc.bot),
            pop =pop) %>%
  distinct()  

time_point=as.numeric(dat_pred.inc$date-min(dat_pred.inc$date))
dat_pred.inc=add_column(dat_pred.inc,time_point)
week = round(dat_pred.inc$time_point/7)
dat_pred.inc = add_column(dat_pred.inc, week)



## week data generate 
dat_pred.prev.week <- dat_pred.prev %>% 
  as_tibble() %>%
  filter(date<=max_date & date >=min_date) %>% 
  group_by(week) %>%
  summarize(prev = mean(prev), pop =pop, prev.cnt =round(mean(prev.cnt)),
            prev.top = mean(prev.top), prev.cnt.top =round(mean(prev.cnt.top)),
            prev.bot = mean(prev.bot), prev.cnt.bot =round(mean(prev.cnt.bot)),
            date = max(date)) %>% 
  distinct() 

dat_water.week <- dat_water %>% 
  as_tibble() %>%
  group_by(week) %>%
  summarize(avg_n1 = mean(avg_n1), l_avg_n1 = mean(l_avg_n1), avg_PMMoV= mean(avg_PMMoV),
            avg_n1_m = mean(avg_n1_m)) %>% 
  distinct() 

dat_host.week <- dat_host %>% 
  as_tibble() %>%
  group_by(week) %>%
  summarize(admission_rate = mean(admission_rate)) %>% 
  distinct() 

dat_pred.inc.week <- dat_pred.inc %>% 
  as_tibble() %>%
  filter(date<=max_date & date >=min_date) %>% 
  group_by(week) %>%
  summarize(inc = mean(inc), pop =pop,
            inc.top = mean(inc.top), 
            inc.bot = mean(inc.bot), 
            date = max(date)
  ) %>% 
  distinct() 


###### Merge data
dat_reg <- merge( x= dat_pred.prev.week, y = dat_water.week, by = "week" ,all = TRUE)

dat_reg <-dat_reg %>% 
  filter(date <="2021-09-01") %>%
  filter(!is.na(avg_n1)) 

################################################################################
####### Fitting Bayesian linear regression 
lm.fit.Bayes <-stan_glm(avg_n1_m ~  prev,  data= dat_reg, family = gaussian(),
                        prior = cauchy(),
                        prior_intercept = cauchy(),
                        warmup =1000,
                        it = 3000,
                        chains = 2,
                        refresh = 1000,
                        seed=12345)
print(summary(lm.fit.Bayes, digits = 9, probs=c(0.025, 0.5, 0.975)))

################################################################################
## 
## Generating synthetic data accoding to vaccination status or delta variant status
## simul_competition() for generating prevalence and incidence for vaccination status
## simul_competition_delta() for generating prevalence and incidence for 
## without delta variant and vaccination status
##


####### synthetic prevalence without vaccination  
simul.para = read.csv(file = paste0(result.path,'parameters_',zone,'.csv'))
res = simul_competition(simul.al=0)
dat_pred.prev1 <- res$I1 %>%
  as_tibble() %>%
  transmute(date = as.Date(date),
            prev = as.numeric(med),
            prev.top = as.numeric(top),
            prev.bot = as.numeric(bottom),
            pop = site_pop,
            prev.cnt = round(as.numeric(med) * site_pop),
            prev.cnt.top = round(as.numeric(top) * site_pop),
            prev.cnt.bot = round(as.numeric(bottom) * site_pop)
  )%>%
  group_by(date) %>%
  summarize(prev = mean(prev), pop =pop, prev.cnt =round(mean(prev.cnt)),
            prev.top = mean(prev.top), prev.cnt.top =round(mean(prev.cnt.top)),
            prev.bot = mean(prev.bot), prev.cnt.bot =round(mean(prev.cnt.bot))
  ) %>%
  distinct() 

time_point=as.numeric(dat_pred.prev1$date-min(dat_pred.prev$date))
dat_pred.prev1=add_column(dat_pred.prev1,time_point)
week = round(dat_pred.prev1$time_point/7)
dat_pred.prev1 = add_column(dat_pred.prev1, week)

####### synthetic incidence without vaccination  

dat_pred.inc1 <- res$ds1  %>%
  as_tibble() %>%
  transmute(date = as.Date(date),
            inc = as.numeric(med),
            inc.top = as.numeric(top),
            inc.bot = as.numeric(bottom),
            pop = site_pop,
            inc.cnt = round(as.numeric(med) * site_pop),
            inc.cnt.top = round(as.numeric(top) * site_pop),
            inc.cnt.bot = round(as.numeric(bottom) * site_pop)
  )%>%
  group_by(date) %>%
  summarize(inc = sum(inc), pop =pop, inc.cnt =round(sum(inc.cnt)),
            inc.top = sum(inc.top), inc.cnt.top =round(sum(inc.cnt.top)),
            inc.bot = sum(inc.bot), inc.cnt.bot =round(sum(inc.cnt.bot))
  ) %>%
  distinct() 

time_point=as.numeric(dat_pred.inc1$date-min(dat_pred.inc1$date))
dat_pred.inc1=add_column(dat_pred.inc1,time_point)
week = round(dat_pred.inc1$time_point/7)
dat_pred.inc1 = add_column(dat_pred.inc1, week)

####### synthetic prevalence without delta variant   

simul.para = read.csv(file = paste0(result.path,'parameters_',zone,'.csv'))
data = simul_competition_delta()
dat_pred.prev1 <- data$I1 %>%
  #dat_pred.prev1 <- data.I1 %>%
  as_tibble() %>%
  transmute(date = as.Date(date),
            prev = as.numeric(med),
            prev.top = as.numeric(top),
            prev.bot = as.numeric(bottom),
            pop = site_pop,
            prev.cnt = round(as.numeric(med) * site_pop),
            prev.cnt.top = round(as.numeric(top) * site_pop),
            prev.cnt.bot = round(as.numeric(bottom) * site_pop)
  )%>%
  group_by(date) %>%
  summarize(prev = mean(prev), pop =pop, prev.cnt =round(mean(prev.cnt)),
            prev.top = mean(prev.top), prev.cnt.top =round(mean(prev.cnt.top)),
            prev.bot = mean(prev.bot), prev.cnt.bot =round(mean(prev.cnt.bot))
  ) %>%
  distinct() 

time_point=as.numeric(dat_pred.prev1$date-min(dat_pred.prev$date))
dat_pred.prev1=add_column(dat_pred.prev1,time_point)
week = round(dat_pred.prev1$time_point/7)
dat_pred.prev1 = add_column(dat_pred.prev1, week)


####### synthetic incidence without delta variant   

dat_pred.inc1 <- data$dS1%>%
  as_tibble() %>%
  transmute(date = as.Date(date),
            inc = as.numeric(med),
            inc.top = as.numeric(top),
            inc.bot = as.numeric(bottom),
            pop = site_pop,
            inc.cnt = round(as.numeric(med) * site_pop),
            inc.cnt.top = round(as.numeric(top) * site_pop),
            inc.cnt.bot = round(as.numeric(bottom) * site_pop)
  )%>%
  group_by(date) %>%
  summarize(inc = sum(inc), pop =pop, inc.cnt =round(sum(inc.cnt)),
            inc.top = sum(inc.top), inc.cnt.top =round(sum(inc.cnt.top)),
            inc.bot = sum(inc.bot), inc.cnt.bot =round(sum(inc.cnt.bot))
  ) %>%
  distinct() 

time_point=as.numeric(dat_pred.inc1$date-min(dat_pred.inc1$date))
dat_pred.inc1=add_column(dat_pred.inc1,time_point)
week = round(dat_pred.inc1$time_point/7)
dat_pred.inc1 = add_column(dat_pred.inc1, week)



