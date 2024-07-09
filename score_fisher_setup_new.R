###multivariate amle for Weibull regression with both measurement error and misclassification
###(keep 4th order Taylor expansion)
###2024.4.6


###Score function of amle4 for AFT weibull regression with both measurement error and misclassification#################
###Parameter description
##p: number of variables measured with error, e.g., p=2
##betav: start points for parameter estimation, including shape parameters 
##m: conditional mean of x|w; n*p matrix
##tau2: conditional variance of x|w
##z: error free covariates matrix
##zo: observed binary misclassification variable
##c11, c12, c21, c22: estimation of \Delta^{-1}\Delta_e
##t: observed survival time
##delta: censored indicator
##lam: shape parameter in weibull distribution

score_amle4_mis = function(betav){
  b00 = betav[1]; b01 = betav[p+2] 
  beta0 = matrix(betav[2:(p+1)], ncol=1)
  beta1 = matrix(betav[(p+3):(2+2*p)], ncol=1)
  lam = betav[2*(1+p)+1] #shape parameter in Weibull distribution
  #####some notations
  zo0 = 1*(zo==0); zo1 = 1*(zo==1)
  Lam0 = Lam1 = t^lam  #similar to cumulative baseline hazard
  bh0 = bh1 = lam*(t^(lam-1)) #similar to baseline hazard
  eta0 = as.numeric(exp(b00+m%*%beta0)*Lam0)  #n*1 vector
  eta1 = as.numeric(exp(b01+m%*%beta1)*Lam1)  #n*1 vector
  A0 = (exp(b00+m%*%beta0)*bh0)^delta*exp(-eta0)
  A1 = (exp(b01+m%*%beta1)*bh1)^delta*exp(-eta1)
  #original g1(eta0),g1(eta1),g2(eta0),g2(eta1)
  g1_eta0 = delta-(2*delta+1)*eta0+eta0^2
  g1_eta1 = delta-(2*delta+1)*eta1+eta1^2
  g2_eta0 = delta-(14*delta+1)*eta0+(18*delta+7)*eta0^2-2*(2*delta+3)*eta0^3+eta0^4
  g2_eta1 = delta-(14*delta+1)*eta1+(18*delta+7)*eta1^2-2*(2*delta+3)*eta1^3+eta1^4
  
  gamma1_b0 = 0.5*as.numeric(t(beta0) %*% tau2 %*% beta0)
  gamma1_b1 = 0.5*as.numeric(t(beta1) %*% tau2 %*% beta1)
  ##the following terms are some constant terms about elements in tau2
  #for zo=0
  s1112 = 4*tau2[1,1]*tau2[1,2] 
  s2221 = 4*tau2[2,2]*tau2[2,1] 
  s1212 = 2*(tau2[1,1]*tau2[2,2] + 2*tau2[1,2]^2)
  
  b40t0 = tau2[1,1]^2*beta0[1]^4 + tau2[2,2]^2*beta0[2]^4
  b31t0 = s1112*beta0[2]*beta0[1]^3 + s2221*beta0[1]*beta0[2]^3
  b22t0 = s1212*beta0[1]^2*beta0[2]^2 
  gamma2_b0 = 1/8*(b40t0 + b31t0 + b22t0)
  #for zo=1
  b40t1 = tau2[1,1]^2*beta1[1]^4 + tau2[2,2]^2*beta1[2]^4
  b31t1 = (s1112*beta1[2])*beta1[1]^3 + (s2221*beta1[1])*beta1[2]^3 
  b22t1 = s1212*beta1[1]^2*beta1[2]^2
  gamma2_b1 = 1/8*(b40t1 + b31t1 + b22t1)
  
  B0 = 1 + g1_eta0*gamma1_b0 + g2_eta0*gamma2_b0
  B1 = 1 + g1_eta1*gamma1_b1 + g2_eta1*gamma2_b1
  A01 = c11*A0*B0; A02 = c12*A1*B1 
  A03 = c21*A0*B0; A04 = c22*A1*B1
  #first derivative of g1_eta0, g1_eta1, g2_eta0 and g2_eta1
  pg1_eta0 = -(2*delta+1)*eta0+2*eta0^2
  pg1_eta1 = -(2*delta+1)*eta1+2*eta1^2
  pg2_eta0 = -(14*delta+1)*eta0+2*(18*delta+7)*eta0^2-6*(2*delta+3)*eta0^3+4*eta0^4
  pg2_eta1 = -(14*delta+1)*eta1+2*(18*delta+7)*eta1^2-6*(2*delta+3)*eta1^3+4*eta1^4
  cc = diag(p); pgamma1_b0 = pgamma1_b1 = rep(0,p)
  for(i in 1:p){
    ci = matrix(cc[i,],nrow=1) #c1,c2,...,cp
    pgamma1_b0[i] = as.numeric(ci%*%tau2%*%beta0) #P\gamma_1(\beta0), 1*p vector
    pgamma1_b1[i] = as.numeric(ci%*%tau2%*%beta1) #p\gamma_1(\beta1), 1*p vector
  }
  
  #partial diervative of gamma2_eta0 and gamma2_eta1
  g31_0 = 4*tau2[1,1]^2*beta0[1]^3 + 3*(s1112*beta0[2])*beta0[1]^2 + 
          (s2221*beta0[2]^3) + 2*beta0[1]*(s1212*beta0[2]^2)
  g32_0 = 4*tau2[2,2]^2*beta0[2]^3 + 3*(s2221*beta0[1])*beta0[2]^2 + 
          (s1112*beta0[1]^3) + 2*beta0[2]*(s1212*beta0[1]^2)
  pgamma2_b0 = 1/8*c(g31_0,g32_0) #p \gamma_2(\beta0), 1*p vector
  #for zo=1
  g31_1 = 4*tau2[1,1]^2*beta1[1]^3 + 3*(s1112*beta1[2])*beta1[1]^2 + 
          (s2221*beta1[2]^3) + 2*beta1[1]*(s1212*beta1[2]^2)
  g32_1 = 4*tau2[2,2]^2*beta1[2]^3 + 3*(s2221*beta1[1])*beta1[2]^2 + 
          (s1112*beta1[1]^3) + 2*beta1[2]*(s1212*beta1[1]^2) 
  pgamma2_b1 = 1/8*c(g31_1,g32_1) #p \gamma_2(\beta_1), 1*p vector
  
  ###score function for b00 and b01
  #for zo=0
  PA01_0 = A01*(delta-eta0) + c11*A0*(pg1_eta0*gamma1_b0 + pg2_eta0*gamma2_b0)
  PA03_0 = A03*(delta-eta0) + c21*A0*(pg1_eta0*gamma1_b0 + pg2_eta0*gamma2_b0)
  #for zo=1
  PA02_0 = A02*(delta-eta1) + c12*A1*(pg1_eta1*gamma1_b1 + pg2_eta1*gamma2_b1)
  PA04_0 = A04*(delta-eta1) + c22*A1*(pg1_eta1*gamma1_b1 + pg2_eta1*gamma2_b1)
  #score for b00 and b01
  s00 = (PA01_0/(A01+A02))*zo0 + (PA03_0/(A03+A04))*zo1
  s01 = (PA02_0/(A01+A02))*zo0 + (PA04_0/(A03+A04))*zo1
  #score for each variable measured with error
  s00p = s01p = matrix(0,n,p)
  for(i in 1:p){
    PA01i = PA01_0*m[,i] + c11*A0*(pgamma1_b0[i]*g1_eta0 + pgamma2_b0[i]*g2_eta0)
    PA03i = PA03_0*m[,i] + c21*A0*(pgamma1_b0[i]*g1_eta0 + pgamma2_b0[i]*g2_eta0)
    PA02i = PA02_0*m[,i] + c12*A1*(pgamma1_b1[i]*g1_eta1 + pgamma2_b1[i]*g2_eta1)
    PA04i = PA04_0*m[,i] + c22*A1*(pgamma1_b1[i]*g1_eta1 + pgamma2_b1[i]*g2_eta1)
    s00p[,i] = (PA01i/(A01+A02))*zo0 + (PA03i/(A03+A04))*zo1
    s01p[,i] = (PA02i/(A01+A02))*zo0 + (PA04i/(A03+A04))*zo1
  }
  #score for a
  pg1_a_eta0 = log(t)*pg1_eta0;  pg1_a_eta1 = log(t)*pg1_eta1
  pg2_a_eta0 = log(t)*pg2_eta0;  pg2_a_eta1 = log(t)*pg2_eta1
  tempa_0 = delta/lam+(delta-eta0)*log(t); 
  tempa_1 = delta/lam+(delta-eta1)*log(t)
  tempaa0 = pg1_a_eta0*gamma1_b0+pg2_a_eta0*gamma2_b0
  tempaa1 = pg1_a_eta1*gamma1_b1+pg2_a_eta1*gamma2_b1
  PA01_a = A01*tempa_0+c11*A0*tempaa0;  PA03_a = A03*tempa_0+c21*A0*tempaa0;
  PA02_a = A02*tempa_1+c12*A1*tempaa1;  PA04_a = A04*tempa_1+c22*A1*tempaa1;
  s0_a = ((PA01_a+PA02_a)/(A01+A02))*zo0+((PA03_a+PA04_a)/(A03+A04))*zo1
  sfc = cbind(s00,s00p,s01,s01p,s0_a)
  sf = colMeans(sfc)
  return(sf)
}

###variance estimation, using sandwich variance
var_amle4_mis = function(betav, t, delta, m, zo, tau2, c11, c12, c21, c22, p, n){
  b00 = betav[1]; b01 = betav[p+2] 
  beta0 = matrix(betav[2:(p+1)], ncol=1)
  beta1 = matrix(betav[(p+3):(2+2*p)], ncol=1)
  lam = betav[2*(1+p)+1] #shape parameter in Weibull distribution
  #####some notations
  zo0 = 1*(zo==0); zo1 = 1*(zo==1)
  Lam0 = Lam1 = t^lam  #similar to cumulative baseline hazard
  bh0 = bh1 = lam*(t^(lam-1)) #similar to baseline hazard
  eta0 = as.numeric(exp(b00+m%*%beta0)*Lam0)  #n*1 vector
  eta1 = as.numeric(exp(b01+m%*%beta1)*Lam1)  #n*1 vector
  A0 = (exp(b00+m%*%beta0)*bh0)^delta*exp(-eta0)
  A1 = (exp(b01+m%*%beta1)*bh1)^delta*exp(-eta1)
  #original g1(eta0),g1(eta1),g2(eta0),g2(eta1)
  g1_eta0 = delta-(2*delta+1)*eta0+eta0^2
  g1_eta1 = delta-(2*delta+1)*eta1+eta1^2
  g2_eta0 = delta-(14*delta+1)*eta0+(18*delta+7)*eta0^2-2*(2*delta+3)*eta0^3+eta0^4
  g2_eta1 = delta-(14*delta+1)*eta1+(18*delta+7)*eta1^2-2*(2*delta+3)*eta1^3+eta1^4
  
  gamma1_b0 = 0.5*as.numeric(t(beta0) %*% tau2 %*% beta0)
  gamma1_b1 = 0.5*as.numeric(t(beta1) %*% tau2 %*% beta1)
  ##the following terms are some constant terms about elements in tau2
  #for zo=0
  s1112 = 4*tau2[1,1]*tau2[1,2] 
  s2221 = 4*tau2[2,2]*tau2[2,1]
  s1212 = 2*(tau2[1,1]*tau2[2,2] + 2*tau2[1,2]^2)
  
  b40t0 = tau2[1,1]^2*beta0[1]^4 + tau2[2,2]^2*beta0[2]^4
  b31t0 = s1112*beta0[2]*beta0[1]^3 + s2221*beta0[1]*beta0[2]^3
  b22t0 = s1212*beta0[1]^2*beta0[2]^2 
  gamma2_b0 = 1/8*(b40t0 + b31t0 + b22t0)
  #for zo=1
  b40t1 = tau2[1,1]^2*beta1[1]^4 + tau2[2,2]^2*beta1[2]^4
  b31t1 = (s1112*beta1[2])*beta1[1]^3 + (s2221*beta1[1])*beta1[2]^3 
  b22t1 = s1212*beta1[1]^2*beta1[2]^2
  gamma2_b1 = 1/8*(b40t1 + b31t1 + b22t1)
  
  B0 = 1 + g1_eta0*gamma1_b0 + g2_eta0*gamma2_b0
  B1 = 1 + g1_eta1*gamma1_b1 + g2_eta1*gamma2_b1
  A01 = c11*A0*B0; A02 = c12*A1*B1 
  A03 = c21*A0*B0; A04 = c22*A1*B1
  #first derivative of g1_eta0, g1_eta1, g2_eta0 and g2_eta1
  pg1_eta0 = -(2*delta+1)*eta0+2*eta0^2
  pg1_eta1 = -(2*delta+1)*eta1+2*eta1^2
  pg2_eta0 = -(14*delta+1)*eta0+2*(18*delta+7)*eta0^2-6*(2*delta+3)*eta0^3+4*eta0^4
  pg2_eta1 = -(14*delta+1)*eta1+2*(18*delta+7)*eta1^2-6*(2*delta+3)*eta1^3+4*eta1^4
  cc = diag(p); pgamma1_b0 = pgamma1_b1 = rep(0,p)
  for(i in 1:p){
    ci = matrix(cc[i,],nrow=1) #c1,c2,...,cp
    pgamma1_b0[i] = as.numeric(ci%*%tau2%*%beta0) #P\gamma_1(\beta0), 1*p vector
    pgamma1_b1[i] = as.numeric(ci%*%tau2%*%beta1) #p\gamma_1(\beta1), 1*p vector
  }
  
  #partial diervative of gamma2_eta0 and gamma2_eta1
  g31_0 = 4*tau2[1,1]^2*beta0[1]^3 + 3*(s1112*beta0[2])*beta0[1]^2 + 
          (s2221*beta0[2]^3) + 2*beta0[1]*(s1212*beta0[2]^2)
  g32_0 = 4*tau2[2,2]^2*beta0[2]^3 + 3*(s2221*beta0[1])*beta0[2]^2 + 
          (s1112*beta0[1]^3) + 2*beta0[2]*(s1212*beta0[1]^2)
  pgamma2_b0 = 1/8*c(g31_0,g32_0) #p \gamma_2(\beta0), 1*p vector
  #for zo=1
  g31_1 = 4*tau2[1,1]^2*beta1[1]^3 + 3*(s1112*beta1[2])*beta1[1]^2 + 
          (s2221*beta1[2]^3) + 2*beta1[1]*(s1212*beta1[2]^2)
  g32_1 = 4*tau2[2,2]^2*beta1[2]^3 + 3*(s2221*beta1[1])*beta1[2]^2 + 
          (s1112*beta1[1]^3) + 2*beta1[2]*(s1212*beta1[1]^2) 
  pgamma2_b1 = 1/8*c(g31_1,g32_1) #p \gamma_2(\beta_1), 1*p vector
  
  ###score function for b00 and b01
  #for zo=0
  PA01_0 = A01*(delta-eta0) + c11*A0*(pg1_eta0*gamma1_b0 + pg2_eta0*gamma2_b0)
  PA03_0 = A03*(delta-eta0) + c21*A0*(pg1_eta0*gamma1_b0 + pg2_eta0*gamma2_b0)
  #for zo=1
  PA02_0 = A02*(delta-eta1) + c12*A1*(pg1_eta1*gamma1_b1 + pg2_eta1*gamma2_b1)
  PA04_0 = A04*(delta-eta1) + c22*A1*(pg1_eta1*gamma1_b1 + pg2_eta1*gamma2_b1)
  #score for each variable measured with error
  PA01_p = PA03_p = matrix(0,n,p)
  PA02_p = PA04_p = matrix(0,n,p)
  for(i in 1:p){
    PA01_p[,i] = PA01_0*m[,i] + c11*A0*(pgamma1_b0[i]*g1_eta0 + pgamma2_b0[i]*g2_eta0)
    PA03_p[,i] = PA03_0*m[,i] + c21*A0*(pgamma1_b0[i]*g1_eta0 + pgamma2_b0[i]*g2_eta0)
    PA02_p[,i] = PA02_0*m[,i] + c12*A1*(pgamma1_b1[i]*g1_eta1 + pgamma2_b1[i]*g2_eta1)
    PA04_p[,i] = PA04_0*m[,i] + c22*A1*(pgamma1_b1[i]*g1_eta1 + pgamma2_b1[i]*g2_eta1)
  }
  #score for lam
  pg1_a_eta0 = log(t)*pg1_eta0;  pg1_a_eta1 = log(t)*pg1_eta1
  pg2_a_eta0 = log(t)*pg2_eta0;  pg2_a_eta1 = log(t)*pg2_eta1
  
  tempa_0 = delta/lam+(delta-eta0)*log(t)
  tempa_1 = delta/lam+(delta-eta1)*log(t)
  tempaa0 = pg1_a_eta0*gamma1_b0+pg2_a_eta0*gamma2_b0
  tempaa1 = pg1_a_eta1*gamma1_b1+pg2_a_eta1*gamma2_b1
  PA01_a = A01*tempa_0+c11*A0*tempaa0;  PA03_a = A03*tempa_0+c21*A0*tempaa0;
  PA02_a = A02*tempa_1+c12*A1*tempaa1;  PA04_a = A04*tempa_1+c22*A1*tempaa1;
 
  ###for the second derivative
  ppg1_eta0 = -(2*delta+1)*eta0+4*eta0^2
  ppg1_eta1 = -(2*delta+1)*eta1+4*eta1^2
  ppg2_eta0 = -(14*delta+1)*eta0+4*(18*delta+7)*eta0^2-18*(2*delta+3)*eta0^3+16*eta0^4
  ppg2_eta1 = -(14*delta+1)*eta1+4*(18*delta+7)*eta1^2-18*(2*delta+3)*eta1^3+16*eta1^4
  
  temp1_delta_pg1_eta0 = ((delta-eta0)*pg1_eta0+ppg1_eta0)*gamma1_b0
  temp2_delta_pg1_eta0 = ((delta-eta0)*pg2_eta0+ppg2_eta0)*gamma2_b0
  P0PA01_0 = PA01_0*(delta-eta0)-A01*eta0+c11*A0*(temp1_delta_pg1_eta0+temp2_delta_pg1_eta0)
  P0PA03_0 = PA03_0*(delta-eta0)-A03*eta0+c21*A0*(temp1_delta_pg1_eta0+temp2_delta_pg1_eta0)
  
  pt = 2*(1+p)+1
  fall = matrix(0,pt,pt)  #second derivative matrix
  #second derivative for the intercept term corresponding to zo=0
  Bijp1 = (P0PA01_0/(A01+A02) - (PA01_0/(A01+A02))^2)*zo0
  Bijp3 = (P0PA03_0/(A03+A04) - (PA03_0/(A03+A04))^2)*zo1 
  fall[1,1] = mean(Bijp1+Bijp3)
  #the first row for 1+p covariates measured with error
  temp1_0 = (delta-eta0)*g1_eta0 + pg1_eta0; 
  temp2_0 = (delta-eta0)*g2_eta0 + pg2_eta0;
  for(j in 2:(p+1)){
    tem0 = A0*(temp1_0*pgamma1_b0[j-1]+temp2_0*pgamma2_b0[j-1])
    PPPA01ij = P0PA01_0*m[,j-1] + c11*tem0
    PPPA03ij = P0PA03_0*m[,j-1] + c21*tem0 #A0*(temp1_0*pgamma1_b0[j-1] + temp2_0*pgamma2_b0[j-1])
    Bijp1 = (PPPA01ij/(A01+A02) - PA01_0*PA01_p[,j-1]/(A01+A02)^2)*zo0
    Bijp3 = (PPPA03ij/(A03+A04) - PA03_0*PA03_p[,j-1]/(A03+A04)^2)*zo1
    fall[1,j] = mean(Bijp1 + Bijp3)
  }
  fall[1,2+p] = mean(-(PA01_0*PA02_0/(A01+A02)^2)*zo0 - (PA03_0*PA04_0/(A03+A04)^2)*zo1)
  for(j in (3+p):(2+2*p)){
    Bijp1 = -(PA01_0*PA02_p[,j-2-p]/(A01+A02)^2)*zo0
    Bijp3 = -(PA03_0*PA04_p[,j-2-p]/(A03+A04)^2)*zo1
    fall[1,j] =  mean(Bijp1 + Bijp3)
  }
  
  ppg1_0a_eta0 = log(t)*ppg1_eta0
  ppg2_0a_eta0 = log(t)*ppg2_eta0
  tempa1_pg_eta0 = (tempa_0*pg1_eta0+ppg1_0a_eta0)*gamma1_b0
  tempa2_pg_eta0 = (tempa_0*pg2_eta0+ppg2_0a_eta0)*gamma2_b0
  PaPA01_0 = PA01_a*(delta-eta0)-A01*eta0*log(t)+c11*A0*(tempa1_pg_eta0+tempa2_pg_eta0)
  PaPA03_0 = PA03_a*(delta-eta0)-A03*eta0*log(t)+c21*A0*(tempa1_pg_eta0+tempa2_pg_eta0)

  Bijp1 = (PaPA01_0/(A01+A02)-PA01_0*(PA01_a+PA02_a)/(A01+A02)^2)*zo0
  Bijp3 = (PaPA03_0/(A03+A04)-PA03_0*(PA03_a+PA04_a)/(A03+A04)^2)*zo1
  fall[1,3+2*p] =  mean(Bijp1 + Bijp3)
  ###construct derivative of ppgamma1_b0, ppgamma2_b0, ppgamma1_b1 and ppgamma2_b1
  ppgamma1_b0 = ppgamma1_b1 = matrix(0,p,p)
  cc = diag(p)
  for(i in 1:p){
    ci = matrix(cc[i,],nrow=1)
    for(j in 1:p){
      cj = matrix(cc[j,],nrow=1)
      ppgamma1_b0[i,j] = as.numeric(ci%*%tau2%*%t(cj))
      ppgamma1_b1[i,j] = ppgamma1_b0[i,j]
    }
  }
  ####pp \gamma_2(\beta_0)
  g311_0 = 12*tau2[1,1]^2*beta0[1]^2+6*(s1112*beta0[2])*beta0[1]+2*(s1212*beta0[2]^2)
  g312_0 = 3*s1112*beta0[1]^2+3*s2221*beta0[2]^2+4*beta0[1]*s1212*beta0[2]
  g321_0 = 3*s2221*beta0[2]^2+3*s1112*beta0[1]^2+4*beta0[2]*s1212*beta0[1]
  g322_0 = 12*tau2[2,2]^2*beta0[2]^2+6*(s2221*beta0[1])*beta0[2]+2*(s1212*beta0[1]^2)
  ppgamma2_b0 = 1/8*matrix(c(g311_0,g312_0,g321_0,g322_0), byrow = T,p,p) 
  #for zo=1
  g311_1 = 12*tau2[1,1]^2*beta1[1]^2+6*(s1112*beta1[2])*beta1[1]+2*(s1212*beta1[2]^2)
  g312_1 = 3*s1112*beta1[1]^2+3*s2221*beta1[2]^2+4*beta1[1]*s1212*beta1[2]
  g321_1 = 3*s2221*beta1[2]^2+3*s1112*beta1[1]^2+4*beta1[2]*s1212*beta1[1]
  g322_1 = 12*tau2[2,2]^2*beta1[2]^2+6*(s2221*beta1[1])*beta1[2]+2*(s1212*beta1[1]^2)
  ppgamma2_b1 = 1/8*matrix(c(g311_1,g312_1,g321_1,g322_1), byrow = T,p,p) #pp \gamma_2(\beta_1)
  #i stands for row, j stands for column
  for(i in 2:(p+1)){
    tem1 = A0*(temp1_0*pgamma1_b0[i-1]+temp2_0*pgamma2_b0[i-1])
    for(j in 1:pt){
      if(j < i) fall[i,j] = fall[j,i]
      else if(j <= (p+1)){
        tem2 = A0*(temp1_0*pgamma1_b0[j-1]+temp2_0*pgamma2_b0[j-1])
        tem3 = A0*(ppgamma1_b0[i-1,j-1]*g1_eta0+ppgamma2_b0[i-1,j-1]*g2_eta0)
        PPPA01ij = P0PA01_0*m[,i-1]*m[,j-1]+c11*tem1*m[,j-1]+c11*tem2*m[,i-1]+c11*tem3
        PPPA03ij = P0PA03_0*m[,i-1]*m[,j-1]+c21*tem1*m[,j-1]+c21*tem2*m[,i-1]+c21*tem3
        Bijp1 = (PPPA01ij/(A01+A02) - PA01_p[,i-1]*PA01_p[,j-1]/(A01+A02)^2)*zo0 
        Bijp3 = (PPPA03ij/(A03+A04) - PA03_p[,i-1]*PA03_p[,j-1]/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp1+Bijp3)
      }
      else if(j == (2+p)){
        Bijp1 = -(PA01_p[,i-1]*PA02_0/(A01+A02)^2)*zo0 
        Bijp3 = -(PA03_p[,i-1]*PA04_0/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp1+Bijp3)
      }
      else if(j <= 2*(p+1)){
        Bijp1 = -(PA01_p[,i-1]*PA02_p[,j-2-p]/(A01+A02)^2)*zo0 
        Bijp3 = -(PA03_p[,i-1]*PA04_p[,j-2-p]/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp1+Bijp3)
      }
      else{
        PaPA01_i = PaPA01_0*m[,i-1]+c11*A0*((tempa_0*g1_eta0+pg1_a_eta0)*pgamma1_b0[i-1]+
                                            (tempa_0*g2_eta0+pg2_a_eta0)*pgamma2_b0[i-1])
        PaPA03_i = PaPA03_0*m[,i-1]+c21*A0*((tempa_0*g1_eta0+pg1_a_eta0)*pgamma1_b0[i-1]+
                                            (tempa_0*g2_eta0+pg2_a_eta0)*pgamma2_b0[i-1])
        Bijp1 = (PaPA01_i/(A01+A02) - PA01_p[,i-1]*(PA01_a+PA02_a)/(A01+A02)^2)*zo0 
        Bijp3 = (PaPA03_i/(A03+A04) - PA03_p[,i-1]*(PA03_a+PA04_a)/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp1+Bijp3)
      }
    }
  }
  #for row 2+p
  #temp1_0 = (delta-eta0)*g1_eta0 + pg1_eta0; temp2_0 = (delta-eta0)*g2_eta0 + pg2_eta0
  temp1_delta_pg1_eta1 = ((delta-eta1)*pg1_eta1+ppg1_eta1)*gamma1_b1
  temp2_delta_pg1_eta1 = ((delta-eta1)*pg2_eta1+ppg2_eta1)*gamma2_b1
  P0PA02_0 = PA02_0*(delta-eta1)-A02*eta1+c12*A1*(temp1_delta_pg1_eta1+temp2_delta_pg1_eta1)
  P0PA04_0 = PA04_0*(delta-eta1)-A04*eta1+c22*A1*(temp1_delta_pg1_eta1+temp2_delta_pg1_eta1)
  
  ppg1_0a_eta1 = log(t)*ppg1_eta1
  ppg2_0a_eta1 = log(t)*ppg2_eta1
  tempa1_pg_eta1 = (tempa_1*pg1_eta1+ppg1_0a_eta1)*gamma1_b1
  tempa2_pg_eta1 = (tempa_1*pg2_eta1+ppg2_0a_eta1)*gamma2_b1
  PaPA02_0 = PA02_a*(delta-eta1)-A02*eta1*log(t)+c12*A1*(tempa1_pg_eta1+tempa2_pg_eta1)
  PaPA04_0 = PA04_a*(delta-eta1)-A04*eta1*log(t)+c22*A1*(tempa1_pg_eta1+tempa2_pg_eta1)
  
  temp1_1 = (delta-eta1)*g1_eta1 + pg1_eta1
  temp2_1 = (delta-eta1)*g2_eta1 + pg2_eta1
  i = p+2
  for(j in 1:pt){
    if(j < i) fall[i,j] = fall[j,i]
    else if(j == (2+p)){
      Bijp2 = (P0PA02_0/(A01+A02) - (PA02_0/(A01+A02))^2)*zo0
      Bijp4 = (P0PA04_0/(A03+A04) - (PA04_0/(A03+A04))^2)*zo1 
      fall[i,j] = mean(Bijp2+Bijp4)
    }
    else if(j<=2+2*p){   
      tem1 = A1*(temp1_1*pgamma1_b1[j-2-p]+temp2_1*pgamma2_b1[j-2-p])
      PPPA02ij = P0PA02_0*m[,j-2-p]+c12*tem1
      PPPA04ij = P0PA04_0*m[,j-2-p]+c22*tem1
      Bijp2 = (PPPA02ij/(A01+A02) - PA02_0*PA02_p[,j-2-p]/(A01+A02)^2)*zo0 
      Bijp4 = (PPPA04ij/(A03+A04) - PA04_0*PA04_p[,j-2-p]/(A03+A04)^2)*zo1
      fall[i,j] = mean(Bijp2+Bijp4)
    }
    else{
      Bijp2 = (PaPA02_0/(A01+A02)-PA02_0*(PA01_a+PA02_a)/(A01+A02)^2)*zo0
      Bijp4 = (PaPA04_0/(A03+A04)-PA04_0*(PA03_a+PA04_a)/(A03+A04)^2)*zo1
      fall[i,pt] =  mean(Bijp2 + Bijp4) #pt=3+2*p
    }
  }
  #for row (3+p):(2+2*p)
  for(i in (3+p):(2+2*p)){
    tem1 = A1*(temp1_1*pgamma1_b1[i-2-p]+temp2_1*pgamma2_b1[i-2-p])
    for(j in 1:pt){
      if(j < i) fall[i,j] = fall[j,i]
      else if(j <= 2*(p+1)){
        tem2 = A1*(temp1_1*pgamma1_b1[j-2-p]+temp2_1*pgamma2_b1[j-2-p])
        tem3 = A1*(ppgamma1_b1[i-2-p,j-2-p]*g1_eta1+ppgamma2_b1[i-2-p,j-2-p]*g2_eta1)
        PPPA02ij = P0PA02_0*m[,i-2-p]*m[,j-2-p]+c12*tem1*m[,j-2-p]+c12*tem2*m[,i-2-p]+c12*tem3
        PPPA04ij = P0PA04_0*m[,i-2-p]*m[,j-2-p]+c22*tem1*m[,j-2-p]+c22*tem2*m[,i-2-p]+c22*tem3
        Bijp2 = (PPPA02ij/(A01+A02) - PA02_p[,i-2-p]*PA02_p[,j-2-p]/(A01+A02)^2)*zo0 
        Bijp4 = (PPPA04ij/(A03+A04) - PA04_p[,i-2-p]*PA04_p[,j-2-p]/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp2+Bijp4)
      }
      else{
        PaPA02_i = PaPA02_0*m[,i-2-p]+c12*A1*((tempa_1*g1_eta1+pg1_a_eta1)*pgamma1_b1[i-2-p]+
                                              (tempa_1*g2_eta1+pg2_a_eta1)*pgamma2_b1[i-2-p])
        PaPA04_i = PaPA04_0*m[,i-2-p]+c22*A1*((tempa_1*g1_eta1+pg1_a_eta1)*pgamma1_b1[i-2-p]+
                                              (tempa_1*g2_eta1+pg2_a_eta1)*pgamma2_b1[i-2-p])
        Bijp2 = (PaPA02_i/(A01+A02) - PA02_p[,i-2-p]*(PA01_a+PA02_a)/(A01+A02)^2)*zo0 
        Bijp4 = (PaPA04_i/(A03+A04) - PA04_p[,i-2-p]*(PA03_a+PA04_a)/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp2+Bijp4)
      }
    } 
  }
  ##for the 2*(p+1)+1 (or pt) 
  for(j in 1:(pt-1)){
    fall[pt,j] = fall[j,pt]
  }
  ppg1_a_eta0 = (log(t))^2*ppg1_eta0; ppg1_a_eta1 = (log(t))^2*ppg1_eta1
  ppg2_a_eta0 = (log(t))^2*ppg2_eta0; ppg2_a_eta1 = (log(t))^2*ppg2_eta1
  tempaa_0 = -delta/lam^2-eta0*(log(t))^2; tempaa_1 = -delta/lam^2-eta1*(log(t))^2

  tempaa0f = (tempa_0*pg1_a_eta0+ppg1_a_eta0)*gamma1_b0+(tempa_0*pg2_a_eta0+ppg2_a_eta0)*gamma2_b0;
  tempaa1f = (tempa_1*pg1_a_eta1+ppg1_a_eta1)*gamma1_b1+(tempa_1*pg2_a_eta1+ppg2_a_eta1)*gamma2_b1;
  PaPA01_a = PA01_a*tempa_0+A01*tempaa_0+c11*A0*tempaa0f
  PaPA03_a = PA03_a*tempa_0+A03*tempaa_0+c21*A0*tempaa0f
  PaPA02_a = PA02_a*tempa_1+A02*tempaa_1+c12*A1*tempaa1f
  PaPA04_a = PA04_a*tempa_1+A04*tempaa_1+c22*A1*tempaa1f
  Bijp12 = ((PaPA01_a+PaPA02_a)/(A01+A02) - (PA01_a+PA02_a)^2/(A01+A02)^2)*zo0 
  Bijp34 = ((PaPA03_a+PaPA04_a)/(A03+A04) - (PA03_a+PA04_a)^2/(A03+A04)^2)*zo1 
  fall[pt,pt] = mean(Bijp12+Bijp34)
  inv_fall = solve(fall)
  ###score function
  s00 = (PA01_0/(A01+A02))*zo0 + (PA03_0/(A03+A04))*zo1
  s01 = (PA02_0/(A01+A02))*zo0 + (PA04_0/(A03+A04))*zo1
  s00p = s01p = matrix(0,n,p)
  for(i in 1:p){
    s00p[,i] = (PA01_p[,i]/(A01+A02))*zo0 + (PA03_p[,i]/(A03+A04))*zo1
    s01p[,i] = (PA02_p[,i]/(A01+A02))*zo0 + (PA04_p[,i]/(A03+A04))*zo1
  }
  s0_a = ((PA01_a+PA02_a)/(A01+A02))*zo0+((PA03_a+PA04_a)/(A03+A04))*zo1
  ss = cbind(s00,s00p,s01,s01p,s0_a)
  jcob = t(ss)%*%ss/n  #pt*pt
  var_call = (1/n)*(inv_fall%*%jcob%*%inv_fall) ##sandwich variance
  return (var_call)
}