
functions {
  
  real[] SVIRT(real t , real[] y , real[] parms , real[] rdata , int[] idata){
    
    real beta1 = parms[1];    // Infection rate from S to I 
    real alpha = parms[2];    // vaccination rate 
    real beta2 = parms[3];    // Infection rate from V to I 
    real gamma = parms[4];   // Recovery rate
    real delta = parms[5];   // Rate recovered produce antibodies
    
    real rho = parms[6];     // Initially infected
    real eps = parms[7];     // Initially recovered
    real psi = parms[8];     // Initial amount with antibodies

    real dydt[6];            // System of ODES
    
      dydt[1] = - beta1 *y[1] * y[3] - beta1*1.50*y[1]*y[4]- alpha * y[1] ;    // dS/dt
      dydt[2] = alpha * y[1] - beta2 *y[2] * y[3] - beta2* 1.50 *y[2]*y[4];      // dV/dt
      dydt[3] = beta1 *y[1] * y[3] +beta2 *y[2] * y[3] -gamma * y[3];  //dI1/dt
      dydt[4] = beta1 *1.50* y[1] * y[4] +beta2 *1.50*y[2] * y[4] -gamma * y[4];   // dI2/dt
      dydt[5] = gamma * (y[3]+y[4]) - delta*y[5];        // dR/dt
      dydt[6] = delta * y[5];                     // dT/dt
    
    return dydt;           
  }
  
  real PR(real pob, real spec, real sens ){
  
  real sen = sens;
  real spe =  spec; //.966
  real val;
  
  val = sen * pob + (1 - spe) * (1 - pob);
  return val; 
  }
  
}

data {
  
  int<lower=0> k;           //number of time points  
  real<lower=0> t0;         // initial time 
  real <lower=0> ti[k];     // test time (in days)
  real<lower=0> xi[k];      // positive tests  
  real<lower=0> ni[k];      // total tests 

}

transformed data{
  
  real x_r[0];
  int x_i[0];
  
}

parameters {
  
  
  real<lower=0.0> beta1;
  real<lower=0.0> alpha;
  real<lower=0.0> beta2; 
  real<lower=0.1> gamma;       
  real<lower=0.0> delta;    
  
  real<lower=0.0> rho;      
  real <lower=0.0> eps;
  real <lower=0.0> psi;     
  real <lower=0.0> spec;
  real <lower=0.0> sens;
  
}

transformed parameters{
  
}

model {
  
  real temp[k,6];     // Solution matrix
  real parms[8];      // Parameters vector
  real init[6];       // Initial condition vector

  parms[1] = beta1;
  parms[2] = alpha;
  parms[3] = beta2;
  parms[4] = gamma; 
  parms[5] = delta;

  parms[6] = rho;
  parms[7] = eps; 
  parms[8] = psi;

  init[1] = 1 - rho - rho*0.001 - eps - psi; // Initial value of susceptible
  init[2] = 0;                              // Initial value of vaccination
  init[3] = rho;                            // Initial value of infected. 
  init[4] = rho*0.001;                       // Initial value of delta variant infected
  init[5] = eps;                            // Initial value of recovered
  init[6] = psi;                            // Initial value of seroprevalent
  
  temp = integrate_ode_rk45(SVIRT,init,t0,ti,parms,x_r,x_i);  // Solution matrix

  for(i in 1:k){
    
    target += xi[i] * log(PR(temp[i,6],spec,sens));           // Likelihood of positives
    target += (ni[i]-xi[i]) * log(1-PR(temp[i,6],spec, sens)); // Likelihood of negatives
    
  }
  
  target += gamma_lpdf(beta1| 40.97*1, 92.32*1);  
  target += gamma_lpdf(alpha| 0.0100*10000000, 1000000);  // for msd01 
  target += gamma_lpdf(beta2| 40.97*1, 92.32*1);  
  target += gamma_lpdf(gamma| 21.80*1, 90.32*1);
  target += gamma_lpdf(delta| 24.29, 232.00);
  target += gamma_lpdf(rho| 0.0012*4648, 4648);
  target += gamma_lpdf(eps| 1.74, 1039.09);
  target += gamma_lpdf(psi| 112.50, 5035.15);
  target += beta_lpdf(spec| 21.7*1, 3.83*1); 
  target += beta_lpdf(sens| 7.1*10, 3.83*10); 
}
