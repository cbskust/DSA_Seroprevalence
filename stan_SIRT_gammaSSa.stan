
functions {
  
  real[] SIRT(real t , real[] y , real[] parms , real[] rdata , int[] idata){
    
    real beta = parms[1];    // Infection rate
    real gamma = parms[2];   // Recovery rate
    real delta = parms[3];   // Rate recovered produce antibodies
    
    real rho = parms[4];     // Initially infected
    real eps = parms[5];     // Initially recovered
    real psi = parms[6];     // Initial amount with antibodies

    real dydt[4];            // System of ODES
    
      dydt[1] = - beta *y[1] * y[2];                // dS/dt
      dydt[2] = beta*y[1] * y[2] - gamma*y[2];      // dI/dt
      dydt[3] = gamma*y[2] - delta*y[3];        // dR/dt
      dydt[4] = delta*y[3];                     // dT/dt
    
    return dydt;           
  }
  
  real PR(real pob , real spec){
  
  real sens = 1;
  real spe =  spec ; //.966;
  real val;
  
  val = sens * pob + (1 - spe) * (1 - pob);
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
  
  real<lower=0.0> beta;                    
  real<lower=0.1> gamma;       
  real <lower=0.0> delta;    
  
  real<lower=0.0> rho;      
  real <lower=0.0> eps;
  real <lower=0.0> psi;     
  real <lower=0.0> spec;  
  
}

transformed parameters{
  
}

model {
  
  real temp[k,4];     // Solution matrix
  real parms[6];      // Parameters vector
  real init[4];       // Initial condition vector

  parms[1] = beta;  
  parms[2] = gamma; 
  parms[3] = delta;

  parms[4] = rho;
  parms[5] = eps; 
  parms[6] = psi;

  init[1] = 1 - rho - eps - psi;      // Initial number susceptible
  init[2] = rho;                      // Initial number infected
  init[3] = eps;                      // Initial number recovered
  init[4] = psi;                      // Initial numbber seroprevalent
  
  
  
  temp = integrate_ode_rk45(SIRT,init,t0,ti,parms,x_r,x_i,  // Solution matrix
                            1.0E-6,1.0E-6,1.0E6); 
  
  for(i in 1:k){
    
    target += xi[i] * log(PR(temp[i,4],spec));           // Likelihood of positives
    target += (ni[i]-xi[i]) * log(1-PR(temp[i,4],spec)); // Likelihood of negatives
    
  }
  
 //target += gamma_lpdf(beta| 8.74, 30.88);
  target += gamma_lpdf(beta| 40, 28.88);  //14
  target += gamma_lpdf(gamma| 7.85, 2*30.70);
  target += gamma_lpdf(delta| 25.41, 239.25);
  target += gamma_lpdf(rho| 3.1, 11206.96);
  target += gamma_lpdf(eps| 1.79, 854.29);
  target += gamma_lpdf(psi| 112, 5112.49);
   target += beta_lpdf(spec| 109, 3.83);
}
