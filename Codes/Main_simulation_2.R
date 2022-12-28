#############In this file: how to run the simulations######################
############and apply the simulation file. Several dynamic effects can be included. 
########### Now: re-do the simulation function together with the leads and lags. 
## check results against time and group specific coefficients###
## it is really important ot put lots of weight on the beta_alpha, otherwise hte ##
library("plm")
library("dplyr")
library("ggplot2")
library(data.table)
library(here)
library(tictoc)
library(tidyr)
library(Rfast)
library(pracma)
library(readr)
source(here("Codes","functions.R"))
setwd(here("Codes","Data_Simulation"))#####specifies the working directory where the simulation runs are saved: here in the data simulation folder
##########################################################################################################################
######------------------------------Generate the data--------------#######################################################
######------------------Specify global parameters of the simulations-------###############################################
##########################################################################################################################
t=10###number of time periods#############################################################################################
c<-c(-0.1,0.1)###Coefficients for additional time-constant regressors###################################################
beta1<-sample(c,t,replace=T)
beta=0#########################The coefficient for abortion/parameter of interest --------------##########################
lead=5##Number of leads, maximum is 5############################################
lag=5###Number of lags, maximum is 5
beta_alpha<-5#####Turns time-varying unobserved heterogeneity on/off. If beta_alpha=1: time-varying unobserved heterogeneity is present
n=100#############Number of observations##########
lpm=0##If lpm=1: linear probability model, else: numeric
####the betas for the leads and lags###Extension for dynamic effects##
ll<-0######Indicator for lags: if ll=1: leads/lags present############
ifelse(ll==1,beta_lead<-seq(0,0,le=lead),beta_lead<-rep(0,lead))##Coefficient for leads
ifelse(ll==1,beta_lag<-seq(0,0,le=lag),beta_lag<-rep(0,lag))##Coefficients for lags
beta_age<-1# determines the way that alpha is influenced by the age of the abortion
###Make time dummies#########################################################
time=rep(1:t,n)
ai<-factor(time)
dummies_t = model.matrix(~ai-1)
beta_time=rep(0,t)#optional time-trend coefficients
###for each individual draw an initial age and then age them forward#########
########optional age trends, not currently in the simulations in the paper##
range=(t-1)
age_init=20
age_end=29
age=round(runif(n,20,29),0)
age=rep(age,each=t)
age=age+rep(seq(0,(t-1),le=t),n)
age.ai<-factor(age)
dummies.age = model.matrix(~age.ai-1)
###The coefficient for the age trend##########################################
beta_age_2=(seq(0,0,le=dim(dummies.age)[2]))^2
##############################################################################
########Finish global parameters##############################################
##----------------------------Run the simulations: sim times------------------
##############################################################################
sim=20
x<-seq(1:10000)####Vector of random seeds: replace with custom seeds if desired####
seed.vec=sample(x,sim)###############
sim.data.fun<-function(sim,n){
  for(i in 1:sim)
  {
    seed=seed.vec[i]
    result<-group.fun(beta_age,lead,lag,ll)
    GM=result$GM
    alpha_gi=result$alpha_gi
    dummies=result$dummies
    alpha=result$alpha
   # age_new=result$age_new
    df<-result$df
    tic("sleeping")
    print("falling asleep...")
    #simu_fun(n,seed,beta,beta1,beta_alpha,t,G,GM,alpha_gi,alpha,age_new,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
    simu_fun_2(n,seed,beta,beta1,beta_alpha,t,G,GM,dummies, alpha_gi,alpha,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
    print("...waking up")
    toc()
  }
}
######Run the simlations, data files will be saved in the data_simulations folder###############
n=1000
sim.data.fun(sim,n)
####always have to re-run the age and time dummies, otherwise dimension problems; make a gunction 
n=1500
sim.data.fun(sim,n)
n=2000
sim.data.fun(sim,n)
n=2500
sim.data.fun(sim,n)
n=3000
sim.data.fun(sim,n)

#########----------The results from the calculations for the parameter of interest, other covariates, standard errors, unobserved 
#########----------heterogeneity profiles, both BIC's are saved on disk----------#########
main.sim<-function(n,b,ll)
{
  
  filenames <-list.files(pattern=paste0(n,"_","data_"))#,"_","data_",beta==0,"_",ll,"_",lpm,sep="")) 
  filenames<-sample(filenames,length(filenames))
  ######Read in the G_max, needs to be calculated beforehand
  sigma2<- readRDS(paste0(n,"sigma2"))
  run=length(filenames)###Number of calculations for the BIC estimator, maximum is the number of simulated files#####
  for(j in 1:run)
  {
    DAtA<-read.csv(filenames[j])
    delta_mat<-matrix(NA,t,1)
    theta_1<-c()
    #theta<-matrix(NA,11,1)
    se<-c()
    G_star<-1:7
    BIC<-c()
    BIC1<-c()
    obj<-c()
    for(i in 1:length(G_star))
    {
      X0 = cbind(DAtA$V4 ,DAtA[,9:18])
      #X0 = cbind(DAtA$value,DAtA[,5])
      # X0<-DAtA$value
      Y0 = DAtA$Y
      XX = X0#; % Useful if X0 contains additional variables. Change to number of variables used in regression
      YY = Y0
      #X = as.matrix(cbind(rep(1,t*n),XX));
      X = as.matrix(XX);
      Y = as.numeric(YY);
      K = ifelse(dim(X)[2]>1,dim(X)[2],1)#
      n = length(Y)/t##number of observations
      results_tr<-BM_fun(G_star[i],sim=1,t,X,Y)
      gi_class=results_tr$gi_class
      theta_1[i]<-results_tr$beta[1]
      theta_all<-results_tr$beta
      delta<-results_tr$delta
      theta<-cbind(theta_1,theta_all)
      delta_mat<-cbind(delta_mat,delta)
      se[i]<-results_tr$se[1]
      ob_group<-opt_group(Y,X,G_star[i],theta_all,sigma2,n,t,K,delta,gi_class)
      BIC[i]=ob_group$BIC
      BIC1[i]=ob_group$BIC1
      obj[i]=ob_group$obj
      print(BIC[i])
    }
    #########----------The results from the calculations for the parameter of interest, other covariates, standard errors, unobserved 
    #########----------heterogeneity profiles, both BIC's are saved on disk----------#########
    write.table(theta_1,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_theta_1_rep",sep=""),col.names = F,append=TRUE)
    write.table(theta,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_theta_rep",sep=""),col.names = F,append=TRUE)
    write.table(se,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_se_rep",sep=""),col.names = F,append=TRUE)
    write.table(delta_mat,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_delta_mat_rep",sep=""),col.names = F,append=TRUE)
    write.table(BIC,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_BIC",sep=""),col.names = F,append=TRUE)
    write.table(BIC1,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_BIC1",sep=""),col.names = F,append=TRUE)
    print(j)
  }
}
b=beta######is parameter of interest equal to zero
#ll=0###lags yes/no
main.sim(100,b=beta,ll)














library("plm")
library("dplyr")
library("ggplot2")
library(data.table)
library(here)
library(tictoc)
library(tidyr)
library(Rfast)
library(pracma)
library(readr)
#library(fixest)
source(here("Codes","functions_2.R"))
setwd(here("Codes","Data_Simulation"))
library("plm")

#########################################################################################################
#sim=10
#x<-seq(1:10000)####Vector of random seeds: replace with custom seeds if desired####
#seed.vec=sample(x,sim)###################################################################################
#seed=seed.vec[1]
#result<-group.fun(beta_age,lead,lag,ll)
#GM=result$GM
#alpha_gi=result$alpha_gi
#dummies=result$dummies
#alpha=result$alpha
#age_new=age.ai#result$age_new
#df<-result$df
#df<-as.matrix(data.frame(apply(df, 2, function(x) as.numeric(as.character(x)))))###transform to numeric
################# ################################
#####Heterogeneous treatment effects #############
##### What kind of treatment effect heterogeneity do we want? 
##### This version has treatment effect heterogeneity by group#######
#delta<-c(0.1,0.2,0.3)#previous version  beta=0.5###The coefficient for abortion/parameter of interest --------------##########################
#################################
### If we wanted to have treatment effect heterogeneity by timing, we need a specific 
### beta associated with the treatment time/age, not with the groups#### 
### One other thing to look at: if the coefficients are group specific, would reweighting help?## 
### Get the estimated group sizes and assignments and check for that ###


####first: look at the distribution of treatments, i.e. how many treated in each period: actually
#### looks pretty good (about 3% each year and evenly distributed)#### 

#### 


#################################
#####df: this is the data frame with the index and the treatment period dummy and the leads and lags dummies##
#result_df<-simu_fun(n,seed,beta=delta,beta1,beta_alpha,t,G,GM, alpha_gi,alpha,age_new,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
#y<-result_df$Y
#y<-as.vector(y)#as.matrix(data.frame(apply(y, 2, function(x) as.numeric(as.character(x)))))
######in X included should be all the lag indicators, time dummies and age dummies. For now: only value and leads and lags
#X<-cbind(result_df[,4:(4+lead+lag)],dummies_t)###full matrix
#X<-X[,-(2+(lead-2))]# X$mid_o_lead_1 ####drop first lead before treatment, only necessary when all time periods are included as leads or lags leads and lags are included#X[,-(2+lead)]####drop first lead before treatment
#X<-as.matrix(data.frame(apply(X, 2, function(x) as.numeric(as.character(x)))))
#delta_hat<-solve(t(X)%*%X)%*%t(X)%*%y### in this way, the estimator estimates the unweighted average. 
#delta_hat
#delta_w1<-(table(GM)[1]*delta[1]+table(GM)[2]*delta[2]+table(GM)[3]*delta[3])/(n)
#delta_w2<-mean(delta)

###below: I should also include the time fixed effects at some point#####
####Implement heterogeneous treatment effects by timing, i.e. by time period## 

######----------------------------Run the simulations: sim times------------------#######
#########################################################################################
sim=10
x<-seq(1:10000)####Vector of random seeds: replace with custom seeds if desired####
seed.vec=sample(x,sim)###############
sim.data.fun<-function(sim,n,beta,het){
  theta_all<-matrix(NA,sim,(length(beta1)+1))##this should have the dimension of the covariate vector
  delta_mat<-matrix(NA,t,1)
  gi_class<-matrix(NA,n,sim)
  ab_prob<-matrix(NA,t,sim)
  fe<-c()
  coeff_ols<-c()
  for(i in 1:sim)
  {
    #####Old simulation set up:
    seed=seed.vec[i]
    result<-group.fun(beta_age,lead,lag,ll)
    GM=result$GM
    alpha_gi=result$alpha_gi
    dummies=result$dummies
    alpha=result$alpha
    age_new=result$age_new
    df<-result$df
    result_df_old<-simu_fun_old(n,seed,beta,beta1,beta_alpha,t,G,GM,alpha_gi,alpha,age_new,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
    Y<-as.numeric(result_df_old$Y)
    #X<-cbind(result_df[,4:(4+lead+lag)])###full matrix
    X<-cbind(result_df_old[,4:4:dim(result_df_old)[2]])###full matrix
    X_test<-as_tibble(cbind(X$value,rep(1:t,n)), name_repair=T)
    X<-X[,-2]
    X<-X[,-2]
    #simu_fun_old(n,seed,beta,beta1,beta_alpha,t,G,GM,alpha_gi,alpha,age_new,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
    X<-as.matrix(data.frame(apply(X, 2, function(x) as.numeric(as.character(x)))))
    K = ifelse(dim(X)[2]>1,dim(X)[2],1)#
    n = length(Y)/t##number of observations
    results_tr<-BM_fun(G=3,sim=20,t,X,Y,K)
    gi_class[,i]=results_tr$gi_class
    #theta_1[i]<-results_tr$beta[1]
    theta_all[i,]<-results_tr$beta
    df<-cbind()
    fe_estim <- plm(Y~value,data=result_df_old,index = c("ID","year"), model = "within",effect="twoways")
    fe[i]<-fe_estim$coeff# std0.10814 
    result_0<-lm(Y~X)
    coeff_ols[i]<-result_0$coeff[2]
    print(i)
  }  
    
    
    se<-c()
    seed=seed.vec[i]
    result<-group.fun(beta_age,lead,lag,ll)
    GM=result$GM
    alpha_gi=result$alpha_gi
    dummies=result$dummies
    alpha=result$alpha
    age_new=result$age_new
    df<-result$df
    ######################################################
    tic("sleeping")
    print("falling asleep...")
    #beta=delta
    result_df<-simu_fun(n,seed,beta=beta,beta1,beta_alpha,t,G,GM, alpha_gi,alpha,age_new,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm=0,het)
    Y<-as.numeric(result_df$Y)
    #X<-cbind(result_df[,4:(4+lead+lag)])###full matrix
    X<-cbind(result_df[,4:4:dim(result_df)[2]])###full matrix
    X_test<-as_tibble(cbind(X$value,rep(1:t,n)), name_repair=T)
    ab_prob[,i]<- as.matrix(X_test%>%group_by(V2)%>%summarise(mean(V1)))[,2]
    X<-X[,-2]
    X<-X[,-2]
    X<-as.matrix(data.frame(apply(X, 2, function(x) as.numeric(as.character(x)))))
    
    K = ifelse(dim(X)[2]>1,dim(X)[2],1)#
    n = length(Y)/t##number of observations
    results_tr<-BM_fun(G=3,sim=1,t,X,Y,K)
    gi_class[,i]=results_tr$gi_class
    #theta_1[i]<-results_tr$beta[1]
    theta_all[i,]<-results_tr$beta
    delta<-results_tr$delta
  #theta<-cbind(theta,theta_all)
    delta_mat<-cbind(delta_mat,delta)
    #delta_mat1[,i]<-results_tr$delta[,1]
    #delta_mat2[,i]<-results_tr$delta[,2]
    #delta_mat3[,i]<-results_tr$delta[,3]
    
    se[i]<-results_tr$se[1]
    print(i)
     print("...waking up")
    toc()
  }
  #########----------The results from the calculations for the parameter of interest, other covariates, standard errors, unobserved 
  #########----------heterogeneity profiles, both BIC's are saved on disk----------#########
 #beta100<-beta*100
  write.table(theta_all,file=paste0(n,"_",var(beta)==0,"_",beta[1]<beta[2],"_",ll,"_","_theta_all_rep_",het,sep=""),col.names = F,append=TRUE)
  write.table(gi_class,file=paste0(n,"_",var(beta)==0,"_",beta[1]<beta[2],"_",ll,"_","_gi_class_",het,sep=""),col.names = F,append=TRUE)
  write.table(se,file=paste0(n,"_",var(beta)==0,"_",beta[1]<beta[2],"_",ll,"_","_se_rep_",het,sep=""),col.names = F,append=TRUE)
  write.table(delta_mat,file=paste0(n,"_",var(beta)==0,"_",beta[1]<beta[2],"_",ll,"_",lpm,"_delta_mat_rep_",het,sep=""),col.names = F,append=TRUE)
  write.table(ab_prob,file=paste0(n,"_",var(beta)==0,"_",beta[1]<beta[2],"_",ll,"_",lpm,"_ab_prob_",het,sep=""),col.names = F,append=TRUE)
  write.table(rowMeans(ab_prob),file=paste0(n,"_",var(beta)==0,"_",beta[1]<beta[2],"_",ll,"_",lpm,"_ab_prob_mean_",het,sep=""),col.names = F,append=TRUE)
  write.table( colMeans(theta_all)[1]-mean(beta),file=paste0(n,"_",var(beta)==0,"_",beta[1]<beta[2],"_",ll,"_",lpm,"_mean_dev_",het,sep=""),col.names = F,append=TRUE)
  print(colMeans(theta_all)[1]-mean(beta))
  return(list(theta_all=theta_all,delta_mat=delta_mat,gi_class=gi_class,thet_mean=colMeans(theta_all)[1]-mean(beta)))
}


####Generating n and t specific inputs #######
#####evaluate the group specific uh ##########
t=10###number of time periods, right now the unobserved time profiles are only valid for ten time periods. #############################################################################################
n=2000
c<-c(-0.2,0.2)###Coefficients for additional time-constant regressors###################################################
beta1<-sample(c,t,replace=T)
beta=0.5
lead=1##Number of leads, maximum is 5############################################
lag=1##Number of lags, maximum is 5
beta_alpha<-1#####Turns time-varying unobserved heterogeneity on/off. If beta_alpha=1: time-varying unobserved heterogeneity is present
#n=80#############Number of observations##########
lpm=0##If lpm=1: linear probability model, else: numeric
####the betas for the leads and lags###Extension for dynamic effects##
ll<-0######Indicator for lags: if ll=1: leads/lags present############
ifelse(ll==1,beta_lead<-seq(0.0,0.0,le=lead),beta_lead<-rep(0,lead))##Coefficient for leads
ifelse(ll==1,beta_lag<-seq(0.15,0.1,le=lag),beta_lag<-rep(0,lag))##Coefficients for lags
beta_age<-1# determines the way that alpha is influenced by the age of the abortion
###Make time dummies#####################
time=rep(1:t,n)
ai<-factor(time)
dummies_t = model.matrix(~ai-1)
beta_time=rep(0,t)#optional time-trend coefficients
###for each individual draw an initial age and then age them forward#########
########optional age trends, not currently in the simulations in the paper##
range=(t-1)
age_init=20
age_end=29
age=round(runif(n,20,29),0)
age=rep(age,each=t)
age=age+rep(seq(0,(t-1),le=t),n)
age.ai<-factor(age)
dummies.age = model.matrix(~age.ai-1)
###The coefficient for the age trend#####################################################################
beta_age_2=(seq(0-2,-0.2,le=dim(dummies.age)[2]))^2



beta=c(0.2,0.2,0.2)
theta_all1<-sim.data.fun(sim,n,beta,het="gs")
beta=c(0.1,0.2,0.3)
theta_all2<-sim.data.fun(sim,n,beta,het="gs")
beta=c(0.3,0.2,0.1)###look again; what is group high risk####
theta_all3<-sim.data.fun(sim,n,beta,het="gs")
#####Now evaluate the timing/period specific heterogeneity#######
#####homogeneous effect for reference ###########################
####change the naming of the table files ########################
beta=rep(0.2,le=t)
theta_all4<-sim.data.fun(sim,n,beta,het="ts")
z<-seq(0.1,0.3,le=t)
beta<-sample(z,t)
#ATE=mean(delta)###unweighted ATE
theta_all5<-sim.data.fun(sim,n,beta,het="ts")





#######When comparing the unobserved heterogeneity trends, these should be ordered by mean value##
table(theta_all$gi_class)






