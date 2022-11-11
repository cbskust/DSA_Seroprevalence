#setwd("~/custom_folder")


################################################################################
# This code estimates all parameters of the SVIIRT model including specificity 
# and sensitivity using Hamiltonian Monte Carlo (HMC) method. 
# The input data is sero prevalence of COVID19 testing result. 
# This code loads Rstan code for HMC method. 
# The output results include posterior samples of parameters, estimated prevalence,
# estimated sero-prevalence, estimated incidence, summary of input data.  
# This code is for aggregated sewer zone.
#
# Written By: Boseung Choi (cbskust@gmail.com)
# Last Edited: November 10, 2022
#
################################################################################

### Needed libraries ###
source('rank_by_distance2.R')   # Custom ranking function
source('functions.R')   # Custom function

set.seed(12345); 

################################################################################

### Specify data paths ###

result.path = "results/"; 


### setting site_namces according to sewer zones; aggregated, msd1, msd2, msd3 to 5;
site_names=c('MSD01', 'MSD02','MSD03','MSD04','MSD05'); area = 'msd_agg' #'aggregate'

# Stan model path
stan.path <- paste0('stan_SVIIRT_',area,'.stan') ;  


# Start and end of relevant data
date.start <- 0
date.end <- 300
# Plot beyond last data pt
h.T <- 61
# Save paths
path.dS <- paste0(result.path,'incidence_estimate_',area,'.csv')
path.I <- paste0(result.path,'prevalence_estimate_',area,'.csv')
path.T <-paste0(result.path,'seroprevalence_estimate_',area,'.csv')
path.data <- paste0(result.path,'data_',area,'.csv')
path.parameters <- paste0(result.path,'parameters_',area,'.csv')

################################################################################

### Red in and format data ###

# Read in current data
data.current_zone.raw <- read_data(area)

# Subset data frame
data.current_zone.raw <- subset(data.current_zone.raw,
                                select=c('time_point',
                                         'n_result',
                                         'test_date'))

# Group by date, not observation
data.current_zone.formatted <- aggregate(n_result ~ time_point,
                                         data.current_zone.raw, sum)

# Rename columns
colnames(data.current_zone.formatted) <- c('time_point',
                                           'num_antibody_positive')

# Add column of test date
data.current_zone.formatted$date <- sort(unique(
  as.Date(data.current_zone.raw$test_date)))

origin.date <- data.current_zone.formatted$date[1]

# Add column of total tests by day
tests.by.day <- aggregate(n_result ~ time_point,
                          data.current_zone.raw, length)$n_result
data.current_zone.formatted$num_antibody_total <- tests.by.day


################################################################################

### Format data and helper variables for Stan ###

# Subset data for specified wave
data.current_zone.subset <- data.current_zone.formatted[
  (data.current_zone.formatted$time_point >= date.start) &
    (data.current_zone.formatted$time_point <= date.end),]

# Number of data points
k <- length(data.current_zone.subset$time_point)                

# Initial date
t.origin <- data.current_zone.subset$time_point[1]             

# Start from day one
data.current_zone.subset$time_point <- data.current_zone.subset$time_point - 
  t.origin + 1

print(summary(data.current_zone.subset$num_antibody_positive/data.current_zone.subset$num_antibody_total))

# Data as list for Stan
data.SVIRT <- list(k=k,
                  t0=0.5,
                  ti=data.current_zone.subset$time_point,
                  xi=data.current_zone.subset$num_antibody_positive,
                  ni=data.current_zone.subset$num_antibody_total
)

init_fun <-function(...) list(beta1 = runif(1, 0.1, 1), beta2 = runif(1,0.1,1))

################################################################################

### Run Stan code on data ###
tempdir()
fit <- stan(
  file=stan.path,           # Stan program
  data=data.SVIRT,           # named list of data
  chains=2,                 # number of Markov chains
  warmup=1000,              # number of warmup iterations per chain
  iter=2000,                # total number of iterations per chain
  cores=3,                  # number of cores
  refresh=1000,            # show progress every 'refresh' iterations
  control=list(adapt_delta=.9), 
  init=init_fun
)

print(summary(fit, digits = 9, probs=c(0.025, 0.5, 0.975)))
plot(fit, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit, plotfun = "trace", pars = "beta1", inc_warmup = F)
plot(fit, plotfun = "trace", pars = "rho", inc_warmup = F)

################################################################################

### Solve trajectories ###

# Percent positive column
data.current_zone.subset$percent_positive <- data.current_zone.subset$num_antibody_positive /
  data.current_zone.subset$num_antibody_total

# Rank parameters
parameters <- c('beta1', 'alpha','beta2',  'gamma', 'delta', 'rho', 'eps', 'psi','spec','sens')
parameters.frame <- as.data.frame((fit))[parameters]
ranked <- rank_by_distance2(parameters.frame, metric='median')
parameters.ranked <- ranked[[1]]
row <- ranked[[2]]

# Get 95% credible parameter draws
ci = 0.90
parameters.ranked.ci <- parameters.ranked[seq(1, as.integer(ci * nrow(parameters.ranked), 1)),]

# Sample trajectories
sample.size <- 500
parameter.draws <- parameters.ranked.ci[sample(nrow(parameters.ranked.ci), sample.size), ]

# Iteration to solve to 
time.point <- data.current_zone.subset$time_point[nrow(data.current_zone.subset)]+h.T
time.step <- 0.5

# Create matrices of trajectories
trajectory.matrix.dS <- matrix(nrow=sample.size,
                              ncol=time.point / time.step)
trajectory.matrix.I <- matrix(nrow=sample.size,
                             ncol=time.point / time.step)
trajectory.matrix.T <- matrix(nrow=sample.size,
                              ncol=time.point / time.step)

# Index through all parameter draws
for (i in 1:sample.size){
  # Get current row
  row <- parameter.draws[i, ]
  
  # get parameters
  beta1 <- row$beta1
  alpha <- row$alpha
  beta2 <- row$beta2
  gamma <- row$gamma
  delta <- row$delta
  rho <- row$rho
  eps <- row$eps
  psi <- row$psi
  
  # Solve DE
  sol <- euler(t=time.point,
               dt=time.step,
               fun=ode.sviirt,
               ic=c(1 - rho*1.1- eps - psi,0, rho, rho*0.01, eps, psi))
  
  # Save I and T only
  trajectory.matrix.dS[i, ] <- beta1*(sol[,1]*sol[,3]+ 1.20*sol[,1]*sol[,4]) + 
    beta2*(sol[,2]*sol[,3] + 1.20 *sol[,2] * sol[,4]);  # Incidence
  trajectory.matrix.I[i, ] <- sol[,3]+sol[,4]  # Prevalace 
  trajectory.matrix.T[i,] <- sol[, 6] # Seropositivity
}

bep=0.0
# Get confidence bounds
top.dS <- apply(trajectory.matrix.dS, 2, max)*(1+bep)
bottom.dS <- apply(trajectory.matrix.dS, 2, min)*(1-bep)
median.dS <-apply(trajectory.matrix.dS, 2, median)

top.I <- apply(trajectory.matrix.I, 2, max)*(1+bep)
bottom.I <- apply(trajectory.matrix.I, 2, min)*(1-bep)
median.I <-apply(trajectory.matrix.I, 2, median)

top.T <- apply(trajectory.matrix.T, 2, max)*(1+bep)
bottom.T <- apply(trajectory.matrix.T, 2, min)*(1-bep)
median.T <- apply(trajectory.matrix.T, 2, median)

# Format trajectory data
data.dS <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                     top=top.dS,
                     med=median.dS,
                     bottom=bottom.dS)
data.I <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                        top=top.I,
                        med=median.I,
                        bottom=bottom.I)
data.T <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                     top=top.T,
                     med=median.T,
                     bottom=bottom.T)

# Add psuedo median line
data.dS$med <- ((data.dS$top - data.dS$bottom) / 2) + data.dS$bottom
data.I$med <- ((data.I$top - data.I$bottom) / 2) + data.I$bottom
data.T$med <- ((data.T$top - data.T$bottom) / 2) + data.T$bottom


# Add date column
data.dS$date <- data.dS$time + as.Date(data.current_zone.subset$date[1])
data.I$date <- data.I$time + as.Date(data.current_zone.subset$date[1])
data.T$date <- data.T$time + as.Date(data.current_zone.subset$date[1])

# Save csvs
write.csv(data.dS, path.dS, row.names = F) 
write.csv(data.I, path.I, row.names = F)  
write.csv(data.T, path.T, row.names = F)
write.csv(data.current_zone.subset, path.data, row.names = F)
write.csv(as.data.frame(fit)[,-11], path.parameters, row.names = F)

################## End of code ###############################################