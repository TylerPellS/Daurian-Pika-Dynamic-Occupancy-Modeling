
functions {
  real do_missing(int prev, int nmiss, int end, vector phi, vector gam, int yi){
    real ret;
    real thisp = prev * 1.0; // make it numeric
    for(y in 1:nmiss){   // thisp is going above 1
      thisp = thisp * phi[yi+y-1] + (1-thisp) * gam[yi+y-1];  // prob of being occupied
    }
    if(end==1) ret = thisp*phi[yi+nmiss] + (1-thisp) * gam[yi+nmiss];   // prob of being a one at next observation
    if(end==0) ret = thisp*(1-phi[yi+nmiss]) + (1-thisp) * (1-gam[yi+nmiss]);
    return log(ret);  
  }
}

data {
   int<lower=0> nsite;
   int<lower=0> nsurv;  // number of surveyed years
   int<lower=0> nyear; // number of years between first and last survey
   int<lower=0> Kpsi, Kphi, Kgam;   // number of covariates
   
   matrix[nsite,Kpsi] Xpsi;   // model matrix
   array[nsite] matrix[nyear,Kphi] Xphi;   // model matrix
   array[nsite] matrix[nyear,Kgam] Xgam;   // model matrix
   
   array[nsurv-1] int<lower=0> interval ;  // interval between surveys
   array[nsite,nsurv]  int<lower=0,upper=1> Y;    // occupied or not
   
}

transformed data {
  
  array[nsurv] int year_ndx;
  year_ndx[1] = 1;
  for(s in 2:nsurv){
    year_ndx[s] = year_ndx[s-1] + interval[s-1];
  }

}

parameters {
   real<lower=0,upper=1> gamma_0;  // colonization probability
   real<lower=0,upper=1> phi_0;    // staying colonized
   real<lower=0,upper=1> psi1_0;  
   vector[Kpsi] beta_psi;  
   vector[Kphi] beta_phi;
   vector[Kgam] beta_gam;
   
}

transformed parameters {
  
    vector[nsite] psi1 = inv_logit(logit(psi1_0) + Xpsi * beta_psi );
  
    array[nsite] vector[nyear] phi;   // persistence for each site and each year
    array[nsite] vector[nyear] gamma;
    {
      vector[nyear] thismu;
      for(r in 1:nsite){
        thismu = logit(phi_0) + Xphi[r]*beta_phi;
        phi[r] = inv_logit(thismu);
        thismu = logit(gamma_0) + Xgam[r]*beta_gam;
        gamma[r] = inv_logit(thismu);
      }
    }
}

model {
  // priors
  beta_psi ~ normal(0,1);
  beta_phi ~ normal(0,1);
  beta_gam ~ normal(0,1);
  
   // likelihood
  for (r in 1:nsite){
    
       // likelihood of occupancy at first survey 
    target += Y[r,1]*log(psi1[r]) + (1-Y[r,1])*log((1-psi1[r]));
    
    for (t in 2:nsurv){
      
      if(interval[t-1]>1){
        target += do_missing(Y[r,t-1], (interval[t-1]-1), Y[r,t], phi[r], gamma[r], year_ndx[t-1]);  // interval - 1 is number of missing years
    
      } // end multi year interval
      
      if(interval[t-1]==1){
        if(Y[r,t-1]==1){
          target += Y[r,t] * log(phi[r][year_ndx[t]-1]) + (1-Y[r,t])*log((1-phi[r][year_ndx[t]-1]));
        }
      
        if(Y[r,t-1]==0){
          target += Y[r,t] * log(gamma[r][year_ndx[t]-1]) + (1-Y[r,t])*log((1-gamma[r][year_ndx[t]-1]));
        }
      }  // end single year interval
    }
  }
}



