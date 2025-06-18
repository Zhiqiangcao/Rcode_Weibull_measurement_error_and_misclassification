###naive, RCE, AMLE4 of Weibull regression for male's data from EPIC-InterAct Study
###This program can produce similar results in the right panel of Table 6 of main text
###note:  real data has been masked (they are different to real data used in paper) since data is confidential
###R version 4.0.5 (2021-03-31) -- "Shake and Throw"

rm(list=ls()) #clean environment 
library(MASS)
library(rootSolve) 
library(nleqslv)
library(survival)


#input source programs
#make sure the current working directory is the code folder, i.e.,  
#change directory where unzipped files are located
source("C:/Users/82655/Dropbox/research/measurment_error/paper7_mem_mis_cox/R_code/github/score_fisher_setup_real_data_ana.R")
load("C:/Users/82655/Dropbox/research/measurment_error/paper7_mem_mis_cox/R_code/github/epic_mask.RData")

epic = epic_mask
sex = 2 #for female
if(sex == 1) epic = epic[epic$sex == 1,]  #male
if(sex == 2) epic = epic[epic$sex == 2,]  #female
n = dim(epic)[1]

##three error-prone covariates 
fv = scale(sqrt(epic$qge02))+scale(sqrt(epic$qge0401)) #Fv (fruit+vegetable)
prot = sqrt(epic$qe_prot)  #Protein
cho = sqrt(epic$qe_cho)    #Carbohydrate

#misclassification category variable
PA = 1*(epic$pa_index>=4) #physical activity

#standardization 
fv = scale(fv)
prot = scale(prot)
cho = scale(cho)

#err-free covariates
bmi = as.numeric(scale(epic$bmi_adj))  
edu = 1*(epic$l_school>=3)
smoke = 1*(epic$smoke_stat==3)
#center 
is11 = 1*epic$centre==11; is12 = 1*epic$centre==12; is13 = 1*epic$centre==13
is14 = 1*epic$centre==14; is15 = 1*epic$centre==15; is16 = 1*epic$centre==16
is21 = 1*epic$centre==21; is22 = 1*epic$centre==22; is23 = 1*epic$centre==23
is24 = 1*epic$centre==24; is25 = 1*epic$centre==25;
is31 = 1*epic$centre==31; is32 = 1*epic$centre==32; is33 = 1*epic$centre==33
is34 = 1*epic$centre==34; is35 = 1*epic$centre==35;
is41 = 1*epic$centre==41; is42 = 1*epic$centre==42; is51 = 1*epic$centre==51
is52 = 1*epic$centre==52; is71 = 1*epic$centre==71; is72 = 1*epic$centre==72


###############################################Naive############################
PA0 = 1*(PA==0); PA1= 1*(PA==1)
fv0 = fv*PA0; fv1 = fv*PA1; 
prot0 = prot*PA0; prot1 = prot*PA1; 
cho0 = cho*PA0; cho1 = cho*PA1; 

###naive
naive = survreg(Surv(st,event)~PA0+fv0+prot0+cho0+PA1+fv1+prot1+cho1+bmi+edu+smoke+factor(centre)-1,data = epic,dist = "weibull")
naive_s = summary(naive)
res_naive = naive_s$table 
np = dim(res_naive)[1]
res_naive_temp = round(res_naive[1:11,],4)
scale_naive = naive$scale
naive_est = -res_naive_temp[,1]/scale_naive  #transform from aft to weibull regression
shape_est = 1/scale_naive
naive_se = res_naive_temp[,2]/scale_naive
shape_se = sqrt(diag(naive$var)[np])/scale_naive
naive_estval = c(naive_est,shape_est)
naive_std = c(naive_se,shape_se)
zvalue = naive_estval/naive_std  #z test, norm
np = length(naive_estval)
pv = rep(0,np)
for(i in 1:np){
  t_temp = zvalue[i]
  if(t_temp>0) pv[i] = (1-pnorm(t_temp))*2
  else pv[i] = pnorm(t_temp)*2
}
res_naive1 = data.frame(Estimate = naive_estval,
                        Std = naive_std,
                        z_value = zvalue,
                        p_value = pv)
#naive estimation results
naive_res = round(res_naive1,3)
rownames(naive_res) = c("Intercept0","FV0","prot0","cho0","Intercept1","FV1","prot1","cho1","bmi","edu","smoke","gamma")
naive_res

############################################RCE####################################
#infromation of misclassication
po = 0.18        
pi0 = 0.14; pi1 = 0.48  
pt = (po-pi0)/(1-pi0-pi1)
delta = matrix(c((1-pi0)*(1-pt)+pi1*pt,0,0,pi0*(1-pt)+(1-pi1)*pt),2,2)
delta_e = matrix(c((1-pi0)*(1-pt),pi1*pt,pi0*(1-pt),(1-pi1)*pt),byrow = TRUE,2,2)
ccm = solve(delta)%*%delta_e
c11 = ccm[1,1]; c12 = ccm[1,2]
c21 = ccm[2,1]; c22 = ccm[2,2]

#infromation of measurement error
sigxx = matrix(c(0.4199,0.0595,-0.0770,
                 0.0595,1.0106,0.0348,
                 -0.0770,0.0348,0.0576),byrow=T,3,3)
siguu = matrix(c(0.7620,0.3147,0.3005,
                 0.3147,0.8988,0.4970,
                 0.3005,0.4970,0.5228),byrow=T,3,3)
sigxw = matrix(c(0.3161,0.0188,-0.1278,
                 0.0448,0.3198, 0.1002,
                 -0.0580,0.0104, 0.1658),byrow=T,3,3)

#compute E(X|W)
w = cbind(fv,prot,cho)
p = dim(w)[2]
sigww = cov(w)  #using sigww from validation study
umw = apply(w,2,mean) #mu_w
mux = umw  
exw = matrix(rep(mux,n),byrow=TRUE,ncol=p)+t(sigxw%*%solve(sigww)%*%t(w-umw)) 


sigwx = t(sigxw)
rho2 = solve(sigww)%*%sigwx
mt = numeric(p)
mmu = matrix(rep(mt, n), nrow = p)
mux = matrix(numeric(p),p,1);
muw = matrix(numeric(p),p,1);
mu1 = c(1,t(mux)-t(muw)%*%rho2)
mu2 = cbind(matrix(numeric(p),p,1),rho2)
R0 = rbind(mu1,mu2)
R = kronecker(R0,ccm)

wm = cbind(PA0,PA1,fv0,fv1,prot0,prot1,cho0,cho1)
mwm = wm%*%R
PAm0 = mwm[,1]; PAm1 = mwm[,2] 
fvm0 = mwm[,3]; fvm1 = mwm[,4] 
protm0 = mwm[,5]; protm1 = mwm[,6]
chom0 = mwm[,7]; chom1 = mwm[,8]

rc = survreg(Surv(st,event)~PAm0+fvm0+protm0+chom0+PAm1+fvm1+protm1+chom1+bmi+edu+smoke+
               factor(centre)-1,data = epic,dist = "weibull")
rc_s = summary(rc)
res_rc = rc_s$table
np = dim(res_rc)[1]
res_rc_temp = round(res_rc[1:11,],4)
scale_rc = rc$scale
rc_est = -res_rc_temp[,1]/scale_rc  #transform from aft to weibull regression
shape_est = 1/scale_rc
rc_se = res_rc_temp[,2]/scale_rc
shape_se = sqrt(diag(rc$var)[np])/scale_rc
rc_estval = c(rc_est,shape_est)
rc_std = c(rc_se,shape_se)
zvalue = rc_estval/rc_std  #z test, norm
np = length(rc_estval)
pv = rep(0,np)
for(i in 1:np){
  t_temp = zvalue[i]
  if(t_temp>0) pv[i] = (1-pnorm(t_temp))*2
  else pv[i] = pnorm(t_temp)*2
}
res_rce = data.frame(Estimate = rc_estval,
                     Std = rc_std,
                     z_value = zvalue,
                     p_value = pv)
rc_res = round(res_rce,3)
rownames(rc_res) = c("Intercept0","FV0","prot0","cho0","Intercept1","FV1","prot1","cho1","bmi","edu","smoke","gamma")
rc_res


#final results summary
inter0_n = naive_res[1,c(1,2,4)]
fv0_n = naive_res[2,c(1,2,4)]
prot0_n = naive_res[3,c(1,2,4)]
cho0_n = naive_res[4,c(1,2,4)]
inter1_n = naive_res[5,c(1,2,4)]
fv1_n = naive_res[6,c(1,2,4)]
prot1_n = naive_res[7,c(1,2,4)]
cho1_n = naive_res[8,c(1,2,4)]
bmi_n = naive_res[9,c(1,2,4)]
edu_n = naive_res[10,c(1,2,4)]
smoke_n = naive_res[11,c(1,2,4)]
Naive = as.numeric(c(inter0_n,fv0_n,prot0_n,cho0_n,inter1_n,fv1_n,prot1_n,cho1_n,
                     bmi_n,edu_n,smoke_n))

inter0_r = rc_res[1,c(1,2,4)]
fv0_r = rc_res[2,c(1,2,4)]
prot0_r = rc_res[3,c(1,2,4)]
cho0_r = rc_res[4,c(1,2,4)]
inter1_r = rc_res[5,c(1,2,4)]
fv1_r = rc_res[6,c(1,2,4)]
prot1_r = rc_res[7,c(1,2,4)]
cho1_r = rc_res[8,c(1,2,4)]
bmi_r = rc_res[9,c(1,2,4)]
edu_r = rc_res[10,c(1,2,4)]
smoke_r = rc_res[11,c(1,2,4)]
RCE = as.numeric(c(inter0_r,fv0_r,prot0_r,cho0_r,inter1_r,fv1_r,prot1_r,cho1_r,
                   bmi_r,edu_r,smoke_r))


res = cbind(Naive,RCE)
second_row = rep(c("Estimate","SE","P-value"),11)
first_row = c(rep(rep(c("Intercept","FV","Protein","Carbohydrate"),each=3),2),rep(c("BMI","Edu","Smoke"),each=3))
final_res = data.frame(first_row,second_row,res)
final_res
