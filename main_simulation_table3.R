###This program can reproduce simulation results in Table 3 of main text
###author: Zhiqiang CAO
###date: 2024.4.9

rm(list = ls(all = TRUE))

library(nleqslv)
library(MASS)
library(mvtnorm)
library(survival)

#input program settings
setwd("C:/Users/82655/Dropbox/research/measurment_error/paper7_mem_mis_cox/R_code/github")
#input program settings
source("score_fisher_setup_new.R")

###generate simulated data
gend = function(n,c1,lam,trubet,covt,cove,pt,pi0,pi1){
  ##lam is scale parameter and a is shape parameter 
  dataxz = rmvnorm(n=n, mean=c(0,0), sigma=covt)
  x1 = dataxz[,1]; x2 = dataxz[,2] 
  datae = rmvnorm(n=n, mean=c(0,0), sigma=cove) 
  e1 = datae[,1]; e2 = datae[,2]
  #surrgorate variable for true variable x
  w1 = aa1 + bb1 * x1 + e1
  w2 = aa2 + bb2 * x2 + e2  
  z = rbinom(n, 1, prob = pt)
  b00 = trubet[1]; b10 = trubet[2]; b20 = trubet[3]; 
  b01 = trubet[4]; b11 = trubet[5]; b21 = trubet[6];
  
  x10 = x1[z==0]; x11 = x1[z==1] 
  x20 = x2[z==0]; x21 = x2[z==1]
  n0 = sum(z==0); n1 = sum(z==1)
  zo = numeric(n)
  zo[z==1] = rbinom(n1, 1, 1-pi1)
  zo[z==0] = rbinom(n0, 1, pi0)
  
  #generate survival time
  u = runif(n)
  t1 = rep(n)
  t1[z==0] = (-log(1-u[z==0])/(exp(b00+b10*x10+b20*x20)))^(1/lam)
  t1[z==1] = (-log(1-u[z==1])/(exp(b01+b11*x11+b21*x21)))^(1/lam) 
  cen1 = runif(n, min=0, max=c1)
  delta1 = 1*(t1<=cen1)    
  time1 = delta1*t1+(1-delta1)*cen1  #observed time
  
  datat = data.frame(
    t = time1,
    status = delta1,
    x1 = x1,
    x2 = x2,
    w1 = w1,
    w2 = w2,
    z = z,
    zo = zo
  )
  return(datat)
}


#sample size: 
n = 5000 
#simulation times
N = 1000 

#intercept terms in additive error model
aa1 = aa2 = 0  
#scale-bias terms in additive error model
bb1 = bb2 = 0.5  
#variacne-covariance of true variables
covt = matrix(c(1, 0, 0, 1), 2, 2) 
#correlation coefficient between error variables 
rho_uu = 0.75

qna = qnorm(0.975)
seed0 = seq(1, N) + 110007
varfac = 100
p = 2; 

cho_set = 1:12
scen_set = 1:6
final_res = NULL
for (choice in cho_set){ ###choice of beta
  if(choice == 1){
    trubet = c(0, 0.3933033, 0.3933033, 0, 0.2933033, 0.2933033) #0.1 
  }else if(choice == 2){
    trubet = c(0, 0.5562149, 0.5562149, 0, 0.4562149, 0.4562149) #0.2 
  }else if(choice == 3){
    trubet = c(0, 0.6812213, 0.6812213, 0, 0.5812213, 0.5812213) #0.3 
  }else if(choice == 4){
    trubet = c(0, 0.7866066, 0.7866066, 0, 0.6866066, 0.6866066) #0.4
  }else if(choice == 5){
    trubet = c(0, 0.879453, 0.879453, 0, 0.779453, 0.779453)     #0.5
  }else if(choice == 6){
    trubet = c(0, 0.9633924, 0.9633924, 0, 0.8633924, 0.8633924) #0.6
  }else if(choice == 7){
    trubet = c(0, 0.4890096, 0.4890096, 0, 0.3890096, 0.3890096)   #0.1
  }else if(choice == 8){
    trubet = c(0, 0.6915641, 0.6915641, 0, 0.5915641, 0.5915641)   #0.2
  }else if(choice == 9){
    trubet = c(0, 0.8469896, 0.8469896, 0, 0.7469896, 0.7469896)   #0.3
  }else if(choice == 10){
    trubet = c(0, 0.9780193, 0.9780193, 0, 0.8780193, 0.8780193)   #0.4
  }else if(choice == 11){
    trubet = c(0, 1.093459, 1.093459, 0, 0.993459, 0.993459)    #0.5
  }else if(choice == 12){
    trubet = c(0, 1.197824, 1.197824, 0, 1.097824, 1.097824)    #0.6
  }
  
  ##variacne-covariance of error variables
  if(choice<=6){  #for relative large \tau_m
    cove = matrix(c(0.8, rho_uu*sqrt(0.8*0.8), rho_uu*sqrt(0.8*0.8), 0.8), 2, 2)
  }else{          #for relative small \tau_m
    cove = matrix(c(0.25, rho_uu*sqrt(0.25*0.25), rho_uu*sqrt(0.25*0.25), 0.25), 2, 2)
  }
  
  choice_res = NULL
  for(case in scen_set){ #misclassification cases
    if(case == 1){
      po = 0.44
      pi0 = 0.15; pi1 = 0.15      #Pr(V=0|V_o=0)=0.889; Pr(V=1|V_o=1)=0.800
    }else if(case == 2){
      po = 0.44
      pi0 = 0.218; pi1 = 0.218    #Pr(V=0|V_o=0)=0.847; Pr(V=1|V_o=1)=0.700
    }else if(case == 3){
      po = 0.44
      pi0 = 0.2775; pi1 = 0.2775  #Pr(V=0|V_o=0)=0.819; Pr(V=1|V_o=1)=0.600
    }else if(case == 4){
      po = 0.44
      pi0 = 0.327; pi1 = 0.327    #Pr(V=0|V_o=0)=0.809; Pr(V=1|V_o=1)=0.500
    }else if(case == 5){
      po = 0.44
      pi0 = 0.3654; pi1 = 0.3654  #Pr(V=0|V_o=0)=0.819; Pr(V=1|V_o=1)=0.400
    }else if(case == 6){
      po = 0.44
      pi0 = 0.3937; pi1 = 0.3937  #Pr(V=0|V_o=0)=0.847; Pr(V=1|V_o=1)=0.300
    }
    
    pt = (po-pi0)/(1-pi0-pi1); 
    delta = matrix(c((1-pi0)*(1-pt)+pi1*pt,0,0,pi0*(1-pt)+(1-pi1)*pt),2,2)
    delta_e = matrix(c((1-pi0)*(1-pt),pi1*pt,pi0*(1-pt),(1-pi1)*pt),byrow = TRUE,2,2)
    ccm = solve(delta)%*%delta_e
    c11 = ccm[1,1]; c12 = ccm[1,2]
    c21 = ccm[2,1]; c22 = ccm[2,2]
    
    #measurement error
    sigwx = matrix(c(bb1 * covt[1, 1], bb1 * covt[1, 2], bb2 * covt[2, 1], bb2 *
                       covt[2, 2]), byrow = T, ncol = 2)
    sigxw = t(sigwx)
    msig_11 = bb1 ^ 2 * covt[1, 1] + cove[1, 1]
    msig_12 = bb1 * bb2 * covt[1, 2] + cove[1, 2]
    msig_22 = bb2 ^ 2 * covt[2, 2] + cove[2, 2]
    sigww = matrix(c(msig_11, msig_12, msig_12, msig_22), 2, 2)
    rho2 = solve(sigww)%*%sigwx
    tau2 = covt - sigxw %*% solve(sigww) %*% t(sigxw)
    rho2f = sigxw %*% solve(sigww)
    mt = c(0, 0)
    mmu = matrix(rep(mt, n), nrow = 2)
    mux = matrix(c(0,0),2,1);
    muw = matrix(c(aa1,aa2),2,1)
    mu1 = c(1,t(mux)-t(muw)%*%rho2)
    mu2 = cbind(matrix(rep(0,2),2,1),rho2)
    R0 = rbind(mu1,mu2)
    R = kronecker(R0,ccm)
    taub2 = max(abs(trubet))^2*max(tau2);
    
    extrfac = 10 * max(trubet)
    betae1 = betae2 = betae3 = matrix(0, N, 6)
    
    i = k = 1;
    while (i<=N) {
      cat("iter=",i,"\n")
      set.seed(seed0+k)
      dtaset = gend(n,c1=1.6,lam=1,trubet,covt,cove,pt,pi0,pi1) #consider lam=1 and c1=1.6, 50% censoring rate
      t = dtaset[,1]; delta = dtaset[,2]
      x1 = dtaset[,3]; x2 = dtaset[,4]
      w1 = dtaset[,5]; w2 = dtaset[,6]
      z = dtaset[,7]; zo = dtaset[,8]
      zo0 = 1*(zo==0); zo1 = 1*(zo==1)
      w10 = w1*zo0; w20 = w2*zo0
      w11 = w1*zo1; w21 = w2*zo1
      #naive estimation
      naive_weibull = survreg(Surv(t,delta)~zo0 + w10 + w20 + zo1 + w11 + w21 - 1,dist = "weibull")
      scale_naive = naive_weibull$scale
      est_naive = naive_weibull$coefficients
      est_naive = -est_naive/scale_naive  #transform AFT to weibull regression coefficient
      est_temp1 = c(est_naive,1/scale_naive)  #transform scale parameter to a
      cond1 = max(abs(est_temp1)) <= extrfac
      #RCE
      wm = cbind(zo0,zo1,w1*zo0,w1*zo1,w2*zo0,w2*zo1)
      mwm = wm%*%R
      zm0 = mwm[,1]; zm1 = mwm[,2]
      m10 = mwm[,3]; m20 = mwm[,5]; #when Zo = 0
      m11 = mwm[,4]; m21 = mwm[,6]; #when Zo = 1, length is n
      rc_weibull = survreg(Surv(t,delta)~zm0 + m10 + m20 + zm1 + m11 + m21 - 1,dist = "weibull")
      scale_rc = rc_weibull$scale
      est_rc = rc_weibull$coefficients
      est_rc = -est_rc/scale_rc  #transform from AFT to weibull regression coefficient
      est_temp2 = c(est_rc,1/scale_rc)
      cond2 = max(abs(est_temp2)) <= extrfac
      #AMLE with fourth Taylor expansion
      betav = est_temp1  #start points
      dataw = cbind(w1, w2)
      mw = apply(dataw, 2, mean)
      mmu1 = matrix(rep(mw, n), nrow = 2)
      m = rho2f %*% (t(dataw) - mmu1) + mmu
      m = t(m); m1 = m[, 1]; m2 = m[, 2]
      amle4 = nleqslv(betav, score_amle4_mis)
      flag = amle4$termcd == 1 #converged or not
      est_temp3 = amle4$x
      cond3 = flag & max(abs(est_temp3)) <= extrfac
      cc1 = cond1 & cond2 & cond3
      if(cc1==TRUE){
        #variance of naive 
        var_temp1 = diag(naive_weibull$var) #for gamma and log(scale)
        cond1f = max(var_temp1)>varfac
        #variance of RCE
        var_temp2 = diag(rc_weibull$var)     #for gamma and log(scale)
        cond2f = max(var_temp2)>varfac
        #variance of ALE
        var_temp3 = var_amle4_mis(est_temp3, t, delta, m, zo, tau2, c11, c12, c21, c22, p, n)
        var_temp3 = diag(var_temp3)
        cond3f = is.nan(var_temp3[1])|is.nan(var_temp3[2])|is.nan(var_temp3[3])|is.nan(var_temp3[4])|is.nan(var_temp3[5])|is.nan(var_temp3[6])|is.nan(var_temp3[7])|any(var_temp3<0)|max(var_temp3)>varfac
        cc2 = !cond1f & !cond2f & !cond3f 
        if(cc2 == TRUE){
          #naive estimation
          betae1[i, ] = est_temp1[1:6]   #consider \beta_{ij}, not \lambda
          #RCE
          betae2[i, ] = est_temp2[1:6]
          #AMLE
          betae3[i, ] = est_temp3[1:6]
          i = i+1; k = k+1
        }
        else {i = i; k = k+1}
      }
      else {i = i; k = k+1}
    }
    
    #naive estimate
    ebetae1 = apply(betae1, 2, mean)
    bias1 = ebetae1 - trubet
    rb1 = (bias1/trubet)*100  #relative bias (%)
    resu1 = rb1[c(2:3,5:6)] #consider \beta_{ij} (i = 1; 2, j = 0; 1)
    #RCE
    ebetae2 = apply(betae2, 2, mean)
    bias2 = ebetae2 - trubet
    rb2 = (bias2/trubet)*100  #relative bias (%)
    resu2 = rb2[c(2:3,5:6)] #consider \beta_{ij} (i = 1; 2, j = 0; 1)
    #AMLE
    ebetae3 = apply(betae3, 2, mean)
    bias3 = ebetae3 - trubet
    rb3 = (bias3/trubet)*100  #relative bias (%)
    resu3 = rb3[c(2:3,5:6)] #consider \beta_{ij} (i = 1; 2, j = 0; 1)
    resut = rbind(resu1, resu2, resu3)
    colnames(resut) = c("RB10", "RB20", "RB11", "RB21")
    rownames(resut) = c("naive","rce","amle")
    aarb_res = apply(abs(resut),1,mean)
    res = data.frame(resut,AARB=aarb_res)
    write.csv(round(res,3),paste0("cho_",choice,"case_",case,"n_",n,".csv"))
    cat("case=",case,"\n")
    choice_res = cbind(choice_res,aarb_res)
  }
  colnames(choice_res) = paste0("case",1:6)
  rownames(choice_res) = c("naive","rce","amle")
  final_res = rbind(final_res,choice_res)
  write.csv(round(choice_res,1),paste0("cho_",choice,"n_",n,".csv"))
}
write.csv(round(final_res,1),paste0("table3_res.csv"))


