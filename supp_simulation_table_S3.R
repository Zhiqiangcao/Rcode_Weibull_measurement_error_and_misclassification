###compare performances of true, naive, RCE and AMLE based on weibull regression with both measurement error and misclassification
###This program can reproduce simulation results in Table S3 in supplementary material
###author: Zhiqiang CAO
###R version: R version 4.0.5 (2021-03-31) -- "Shake and Throw"

rm(list = ls(all = TRUE))

library(nleqslv)
library(MASS)
library(mvtnorm)
library(survival)

#input source programs
#make sure the current working directory is the code folder, i.e.,  
#change directory where unzipped files are located
source("C:/Users/user/Dropbox/research/measurment_error/paper7_mem_mis_cox/R_code/github/score_fisher_setup_new.R")


###############data generation################################
gend = function(n,c1,lam,trubet,covt,cove,pt,pi0,pi1){
  ##lam is shape parameter 
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

#####################################main program#########################
#sample size: 5000
n = 5000
#simulation times
N = 1000 

#two choices of regression coefficients
trubet = c(0, 0.3, 0.3, 0, 0.2, 0.2) 
c_set = c(4.7,4.2,3.9);     #25% censored


#intercept terms in additive error model
aa1 = aa2 = 0  
#scale-bias terms in additive error model
bb1 = 0.5; bb2 = 0.5  
#variacne-covariance of true variables
covt = matrix(c(1.0, 0, 0, 1.0), 2, 2)  
rho_uu = 0.75
cove = matrix(c(0.8, rho_uu*sqrt(0.8*0.8), rho_uu*sqrt(0.8*0.8), 0.8), 2, 2)
qna = qnorm(0.975)
#gamma choice
lambda_set = c(5/7,1,1.4);    #three choices for shape parameter

final_res = NULL
for(j in 1:3){
  lam = lambda_set[j]
  c1 = c_set[j]
  for(scenario in 1:3){ #three cases of misclassification rates 
    if(scenario == 1){
      po = 0.22
      pi0 = 0.1; pi1 = 0.1
    }else if(scenario == 2){
      po = 0.22
      pi0 = 0.14; pi1 = 0.23
    }else if(scenario == 3){
      po = 0.44
      pi0 = 0.2; pi1 = 0.2
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
    p = 2; taub2 = max(trubet)^2*max(tau2)
    
    #keep all estimation results 
    betae1 = betae2  = betae3 = betae4 = matrix(0, N, 7)
    var1 = var2 = var3 = var3 = var4 = matrix(0, N, 7)
    cova1 = cova2 = cova3 = cova4 = matrix(0, N, 7)
    
    trutheta = c(trubet,lam)  #true values of parameter
    extrfac = 10 * max(trutheta)
    seed0 = seq(1, N) + 123456
    varfac = 100
    i = k = 1;
    
    while (i<=N) {
      cat("gamma=",lam, "case=",scenario,"iter=",i,"\n")
      set.seed(seed0+k)
      dtaset = gend(n,c1,lam,trubet,covt,cove,pt,pi0,pi1)
      t = dtaset[,1]; delta = dtaset[,2]
      x1 = dtaset[,3]; x2 = dtaset[,4]
      w1 = dtaset[,5]; w2 = dtaset[,6]
      z = dtaset[,7]; zo = dtaset[,8]
      #naive estimation
      zo0 = 1*(zo==0); zo1 = 1*(zo==1)
      w10 = w1*zo0; w20 = w2*zo0
      w11 = w1*zo1; w21 = w2*zo1
      naive_weibull = survreg(Surv(t,delta)~zo0 + w10 + w20 + zo1 + w11 + w21 - 1,dist = "weibull")
      scale_naive = naive_weibull$scale
      est_naive = naive_weibull$coefficients
      est_naive = -est_naive/scale_naive  #transform AFT to weibull regression coefficient
      est_temp1 = c(est_naive,1/scale_naive)  #transform scale parameter to a
      cond1 = max(abs(est_temp1)) <= extrfac
      
      #true estimation
      z0 = 1*(z==0); z1 = 1*(z==1)
      x10 = x1*z0; x20 = x2*z0
      x11 = x1*z1; x21 = x2*z1
      true_weibull = survreg(Surv(t,delta)~z0 + x10 + x20 + z1 + x11 + x21 - 1,dist = "weibull")
      scale_true = true_weibull$scale
      est_true = true_weibull$coefficients
      est_true = -est_true/scale_true  #transform from AFT to weibull regression coefficient
      est_temp2 = c(est_true,1/scale_true)
      cond2 = max(abs(est_temp2)) <= extrfac
      
      #RCE estimation
      wm = cbind(zo0,zo1,w1*zo0,w1*zo1,w2*zo0,w2*zo1)
      mwm = wm%*%R
      zm0 = mwm[,1]; zm1 = mwm[,2]
      m10 = mwm[,3]; m20 = mwm[,5]; 
      m11 = mwm[,4]; m21 = mwm[,6]; 
      rc_weibull = survreg(Surv(t,delta)~zm0 + m10 + m20 + zm1 + m11 + m21 - 1,dist = "weibull")
      scale_rc = rc_weibull$scale
      est_rc = rc_weibull$coefficients
      est_rc = -est_rc/scale_rc  #transform from AFT to weibull regression coefficient
      est_temp3 = c(est_rc,1/scale_rc)
      cond3 = max(abs(est_temp3)) <= extrfac
      
      #AMLE with fourth order Taylor expansion
      betav = est_temp1  #start points
      dataw = cbind(w1, w2)
      mw = apply(dataw, 2, mean)
      mmu1 = matrix(rep(mw, n), nrow = 2)
      m = rho2f %*% (t(dataw) - mmu1) + mmu
      m = t(m); m1 = m[, 1]; m2 = m[, 2]
      amle4 = nleqslv(betav, score_amle4_mis)
      flag = amle4$termcd == 1 
      est_temp4 = amle4$x
      cond4 = flag & max(abs(est_temp4)) <= extrfac
      cc1 = cond1 & cond2 & cond3 & cond4
      if(cc1==TRUE){
        #variance of naive 
        var_temp1 = diag(naive_weibull$var) #for gamma and log(scale)
        cond1 = max(var_temp1)>varfac
        #variance of true
        var_temp2 = diag(true_weibull$var)  #for gamma and log(scale)
        cond2 = max(var_temp2)>varfac
        #variance of rc
        var_temp3 = diag(rc_weibull$var)     #for gamma and log(scale)
        cond3 = max(var_temp3)>varfac
        #variance of amle4
        var_temp4 = var_amle4_mis(est_temp4, t, delta, m, zo, tau2, c11, c12, c21, c22, p, n)
        var_temp4 = diag(var_temp4)
        cond4 = is.nan(var_temp4[1])|is.nan(var_temp4[2])|is.nan(var_temp4[3])|is.nan(var_temp4[4])|is.nan(var_temp4[5])|is.nan(var_temp4[6])|is.nan(var_temp4[7])|any(var_temp4<0)|max(var_temp4)>varfac
        cc2 = !cond1 & !cond2 & !cond3 & !cond4
        if(cc2 == TRUE){
          #naive 
          betae1[i, ] = est_temp1
          var_est1 = var_temp1[1:6]/scale_naive^2
          var_scale1 = scale_naive^2*var_temp1[7]  #delta mathod
          var_a1 = (1/scale_naive^4)*var_scale1
          var1[i,] = c(var_est1,var_a1)
          low1 = betae1[i,]-qna*sqrt(var1[i,])
          high1 = betae1[i,]+qna*sqrt(var1[i,])
          cova1[i,] = as.numeric(low1 <= trutheta & trutheta <= high1)
          #true
          betae2[i, ] = est_temp2
          var_est2 = var_temp2[1:6]/scale_true^2
          var_scale2 = scale_true^2*var_temp2[7]  #delta mathod
          var_a2 = (1/scale_true^4)*var_scale2
          var2[i,] = c(var_est2,var_a2)
          low2 = betae2[i, ] - qna*sqrt(var2[i, ])
          high2 = betae2[i, ] + qna*sqrt(var2[i, ])
          cova2[i,] = as.numeric(low2 <= trutheta & trutheta <= high2)
          #rc
          betae3[i, ] = est_temp3
          var_est3 = var_temp3[1:6]/scale_rc^2
          var_scale3 = scale_rc^2*var_temp3[7]  #delta mathod
          var_a3 = (1/scale_rc^4)*var_scale3
          var3[i,] = c(var_est3,var_a3)
          low3 = betae3[i, ] - qna*sqrt(var3[i, ])
          high3 = betae3[i, ] + qna*sqrt(var3[i, ])
          cova3[i,] = as.numeric(low3 <= trutheta & trutheta <= high3)
          #amle
          betae4[i, ] = est_temp4
          var4[i, ] = var_temp4
          low4 = betae4[i, ] - qna*sqrt(var4[i, ])
          high4 = betae4[i, ] + qna*sqrt(var4[i, ])
          cova4[i,] = as.numeric(low4 <= trutheta & trutheta <= high4)
          i = i+1; k = k+1
        }
        else {i = i; k = k+1}
      }
      else {i = i; k = k+1}
    }
    
    #summary results
    #naive estimate
    ebetae1 = apply(betae1, 2, mean)
    bias1 = ebetae1 - trutheta
    se1 = apply(sqrt(var1), 2, mean)
    rb1 = (bias1/trutheta)*100  #relative bias (%)
    ecova1 = apply(cova1, 2, mean)
    #Keep beta_{10}, beta_{20}, beta_{11} and beta_{21} and lambda
    resu1 = c(rb1[2], se1[2], ecova1[2], 
              rb1[3], se1[3], ecova1[3],
              rb1[5], se1[5], ecova1[5], 
              rb1[6], se1[6], ecova1[6],
              rb1[7], se1[7], ecova1[7])
    #true
    ebetae2 = apply(betae2, 2, mean)
    bias2 = ebetae2 - trutheta
    se2 = apply(sqrt(var2), 2, mean)
    rb2 = (bias2/trutheta)*100  #relative bias (%)
    ecova2 = apply(cova2, 2, mean)
    #Keep beta_{10}, beta_{20}, beta_{11} and beta_{21} and lambda
    resu2 = c(rb2[2], se2[2], ecova2[2],
              rb2[3], se2[3], ecova2[3],
              rb2[5], se2[5], ecova2[5],
              rb2[6], se2[6], ecova2[6],
              rb2[7], se2[7], ecova2[7])
    #RCE
    ebetae3 = apply(betae3, 2, mean)
    bias3 = ebetae3 - trutheta
    se3 = apply(sqrt(var3), 2, mean)
    rb3 = (bias3/trutheta)*100  #relative bias (%)
    ecova3 = apply(cova3, 2, mean)
    #Keep beta_{10}, beta_{20}, beta_{11} and beta_{21} and lambda
    resu3 = c(rb3[2], se3[2], ecova3[2],
              rb3[3], se3[3], ecova3[3],
              rb3[5], se3[5], ecova3[5],
              rb3[6], se3[6], ecova3[6],
              rb3[7], se3[7], ecova3[7])
    #amle
    ebetae4 = apply(betae4, 2, mean)
    bias4 = ebetae4 - trutheta
    se4 = apply(sqrt(var4), 2, mean)
    rb4 = (bias4/trutheta)*100  #relative bias (%)
    ecova4 = apply(cova4, 2, mean)
    #Keep beta_{10}, beta_{20}, beta_{11} and beta_{21} and lambda
    resu4 = c(rb4[2], se4[2], ecova4[2],
              rb4[3], se4[3], ecova4[3],
              rb4[5], se4[5], ecova4[5],
              rb4[6], se4[6], ecova4[6],
              rb4[7], se4[7], ecova4[7])
    
    resut = rbind(resu2, resu1, resu3, resu4)
    colnames(resut) = c("RB10", "se10","CP10",
                        "RB20", "se20","CP20",
                        "RB11", "se11","CP11",
                        "RB21", "se21","CP21",
                        "RB_lam", "se_lam", "CP_lam")
    rownames(resut) = c("true","naive","rc","amle")
    final_res = rbind(final_res,resut)
  }
}
#final results in Table S3 
print(round(final_res,3))

