require(dplyr)
library(rstan)                           # Stan MCMC packaged
library(reshape)                         # Data modifications
require(tibble)
require(rstanarm)
require(ggmcmc)
require(MASS)
require(ciTools) 
require(DescTools)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')

###############################################################################


read_data <-function(area){
  data.path <- paste0('data/',area,'.csv')
  df_test = read.csv(data.path)
  df_test <-df_test %>%
    mutate(vc_date1 = as.Date(vc_date1),
           vc_date2 = as.Date(vc_date2),
           vc_date3 = as.Date(vc_date3),
           APDATETIME = as.POSIXct(APDATETIME,format="%m/%d/%y %H:%M"),
           test_date  = as.Date(test_date),
           vaccinated = vaccinated,
           s_result = s_result,
           n_result = n_result,
           site = site,
           time_point = time_point
    )
}

# General Euler ODE solvers
euler <- function (t, dt, fun, ic)
{ 
  p <- length(ic)
  n <- t/dt 
  xmat <- matrix(0,ncol=p,nrow=n)
  x <- ic 
  xmat[1,] <- x
  for (i in 2:n) { 
    x <- x + fun(x)*dt
    xmat[i,] <- x
  } 
  ts(xmat,start=0,deltat=dt) 
} 


#Stochastic DE 


ode.sviirt <- function(x, b1=beta1, a=alpha, b2 = beta2, g=gamma, d=delta)
{
  c(
    -b1 * x[1] * x[3] - b1* 1.20 * x[1] * x[4] - a * x[1] ,
    a * x[1] - b2 * x[2] * x[3] - b2*1.20 * x[2] * x[4],
    b1 * x[1] * x[3] +  b2 * x[2] * x[3] - g * x[3],
    b1 * 1.20 * x[1] * x[4] +  b2 * 1.20 *  x[2] * x[4] - g * x[4],
    g * (x[3] + x[4]) - d * x[5],
    d * x[5]
  )
}

ode.svirt <- function(x, b1=beta1, a=alpha, b2 = beta2, g=gamma, d=delta)
{
  c(
    -b1 * x[1] * x[3] - a * x[1],
    a * x[1] - b2 * x[2] * x[3],
    b1 * x[1] * x[3] + b2 * x[2] *x[3] - g * x[3],
    g * x[3] - d * x[4],
    d * x[4]
  )
}

simul_competition <- function(parameters.frame=simul.para, simul.al=0) {
  
  ode.sviirt <- function(x, b1=beta1, a=alpha, b2 = beta2, g=gamma, d=delta)
  {
    c(
      -b1 * x[1] * x[3] - b1* 1.20 * x[1] * x[4] - a * x[1] ,
      a * x[1] - b2 * x[2] * x[3] - b2*1.20 * x[2] * x[4],
      b1 * x[1] * x[3] +  b2 * x[2] * x[3] - g * x[3],
      b1 * 1.20 * x[1] * x[4] +  b2 * 1.20 *  x[2] * x[4] - g * x[4],
      g * (x[3] + x[4]) - d * x[5],
      d * x[5]
    )
  }
  
  h.T <- 30; 
  # Rank parameters
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
  time.point <- 131+h.T
  time.step <- 0.5
  
  # Create matrices of trajectories
  trajectory.matrix.dS1 <- matrix(nrow=sample.size,
                                  ncol=time.point / time.step)
  trajectory.matrix.I1 <- matrix(nrow=sample.size,
                                 ncol=time.point / time.step)
  trajectory.matrix.T1 <- matrix(nrow=sample.size,
                                 ncol=time.point / time.step)

  # Index through all parameter draws
  for (i in 1:sample.size){
    # Get current row
    row <- parameter.draws[i, ]
    # get parameters
    beta1 <- row$beta1
    alpha <- simul.al
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
                 ic=c(1 - rho - eps - psi,0, rho,rho*0.01, eps, psi))
    
    # Save I and T only
    trajectory.matrix.dS1[i, ] <- beta1*(sol[,1]*sol[,3]+ 1.20*sol[,1]*sol[,4]) + 
      beta2*(sol[,2]*sol[,3] + 1.20 *sol[,2] * sol[,4]);  # Incidence
    trajectory.matrix.I1[i, ] <- sol[,3]+sol[,4]  # Prevalence 
    trajectory.matrix.T1[i,] <- sol[, 6] # Seropositivity
  }
  
  bep=0.0
  # Get confidence bounds
  top.dS <- apply(trajectory.matrix.dS1, 2, max)*(1+bep)
  bottom.dS <- apply(trajectory.matrix.dS1, 2, min)*(1-bep)
  median.dS <-apply(trajectory.matrix.dS1, 2, median)
  
  top.I <- apply(trajectory.matrix.I1, 2, max)*(1+bep)
  bottom.I <- apply(trajectory.matrix.I1, 2, min)*(1-bep)
  median.I <-apply(trajectory.matrix.I1, 2, median)
  
  top.T <- apply(trajectory.matrix.T1, 2, max)*(1+bep)
  bottom.T <- apply(trajectory.matrix.T1, 2, min)*(1-bep)
  median.T <- apply(trajectory.matrix.T1, 2, median)
  
  # Format trajectory data
  data.dS1 <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                         top=top.dS,
                         med=median.dS,
                         bottom=bottom.dS)
  
  data.I1 <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                        top=top.I,
                        med=median.I,
                        bottom=bottom.I)
  
  data.T1 <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                        top=top.T,
                        med=median.T,
                        bottom=bottom.T)
  
  # Add psuedo median line
  data.dS1$med <- ((data.dS1$top - data.dS1$bottom) / 2) + data.dS1$bottom
  data.I1$med <- ((data.I1$top - data.I1$bottom) / 2) + data.I1$bottom
  data.T1$med <- ((data.T1$top - data.T1$bottom) / 2) + data.T1$bottom

  
  # Add date column
  data.dS1$date <- floor(data.dS1$time) + as.Date("2021-04-23")
  data.I1$date <- floor(data.I1$time) + as.Date("2021-04-23")
  data.T1$date <- floor(data.T1$time) + as.Date("2021-04-23")

  #res = list(ds1 = data.dS1, I1 = data.I1, T1=data.T1, V1=data.V1, tr.I1 = trajectory.matrix.I1, para = parameter.draws)
  res = list(ds1 = data.dS1, I1 = data.I1, T1=data.T1)
}

simul_competition_delta <- function(parameters.frame=simul.para, simul.al=NA) {
  
  ode.svirt <- function(x, b1=beta1, a=alpha, b2 = beta2, g=gamma, d=delta)
  {
    c(
      -b1 * x[1] * x[3] - a * x[1],
      a * x[1] - b2 * x[2] * x[3],
      b1 * x[1] * x[3] + b2 * x[2] *x[3] - g * x[3],
      g * x[3] - d * x[4],
      d * x[4]
    )
  }
  
  h.T <- 30; 
  # Rank parameters
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
  time.point <- 131+h.T
  time.step <- 0.5
  
  # Create matrices of trajectories
  trajectory.matrix.dS1 <- matrix(nrow=sample.size,
                                  ncol=time.point / time.step)
  trajectory.matrix.I1 <- matrix(nrow=sample.size,
                                 ncol=time.point / time.step)
  trajectory.matrix.T1 <- matrix(nrow=sample.size,
                                 ncol=time.point / time.step)

  # Index through all parameter draws
  for (i in 1:sample.size){
    # Get current row
    row <- parameter.draws[i, ]
    # get parameters
    beta1 <- row$beta1
    alpha <-ifelse(is.na(simul.al),row$alpha,simul.al)
    beta2 <- row$beta2
    gamma <- row$gamma
    delta <- row$delta
    rho <- row$rho
    eps <- row$eps
    psi <- row$psi
    
    # Solve DE
    sol <- euler(t=time.point,
                 dt=time.step,
                 fun=ode.svirt,
                 ic=c(1 - rho - eps - psi,0, rho, eps, psi))
    
    # Save I and T only
    trajectory.matrix.dS1[i, ] <- beta1*sol[,3]*sol[,1] + beta2*sol[,3]*sol[,2]  # Incidence
    trajectory.matrix.I1[i, ] <- sol[,3]  # Prevalence 
    trajectory.matrix.T1[i,] <- sol[, 5] # Seropositivity
  }
  
  bep=0.0
  # Get confidence bounds
  top.dS <- apply(trajectory.matrix.dS1, 2, max)*(1+bep)
  bottom.dS <- apply(trajectory.matrix.dS1, 2, min)*(1-bep)
  median.dS <-apply(trajectory.matrix.dS1, 2, median)
  
  top.I <- apply(trajectory.matrix.I1, 2, max)*(1+bep)
  bottom.I <- apply(trajectory.matrix.I1, 2, min)*(1-bep)
  median.I <-apply(trajectory.matrix.I1, 2, median)
  
  top.T <- apply(trajectory.matrix.T1, 2, max)*(1+bep)
  bottom.T <- apply(trajectory.matrix.T1, 2, min)*(1-bep)
  median.T <- apply(trajectory.matrix.T1, 2, median)
  
  # Format trajectory data
  data.dS1 <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                         top=top.dS,
                         med=median.dS,
                         bottom=bottom.dS)
  
  data.I1 <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                        top=top.I,
                        med=median.I,
                        bottom=bottom.I)
  
  data.T1 <- data.frame(time=seq(from=0, to=time.point-time.step, by=time.step),
                        top=top.T,
                        med=median.T,
                        bottom=bottom.T)
  
  # Add psuedo median line
  data.dS1$med <- ((data.dS1$top - data.dS1$bottom) / 2) + data.dS1$bottom
  data.I1$med <- ((data.I1$top - data.I1$bottom) / 2) + data.I1$bottom
  data.T1$med <- ((data.T1$top - data.T1$bottom) / 2) + data.T1$bottom

  # Add date column
  data.dS1$date <- floor(data.dS1$time) + as.Date("2021-04-23")
  data.I1$date <- floor(data.I1$time) + as.Date("2021-04-23")
  data.T1$date <- floor(data.T1$time) + as.Date("2021-04-23")

  res = list(dS1 = data.dS1, I1 = data.I1, T1=data.T1)
}

