###--------------In this file: how to run the simulations--------##################
###--------------The simulation set- up for time specific effects,#################
###--------------time heterogeneity and the group heterogeneity------##############
###--------------Below is also the set-up for the group specific coefficients######
###--------------at the bottom: the baseline for homogenous coefficients ############
library(furrr)
library(purrr)
library(tidyverse)
library(stringr)
library(progressr)
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
n=1000 #############Number of observations##########
t=10  ###number of time periods#############################################################################################
c<-c(-0.05,0.05) ###Coefficients for additional time-constant, additional regressors###################################################
beta1<-sample(c,t,replace=T)
beta=seq(-0.2,0.2,le=t)##The coefficient for abortion/parameter of interest --------------##########################
lead=5##Number of leads, maximum is 5############################################
lag=5###Number of lags, maximum is 5
beta_alpha<-5#####Turns time-varying unobserved heterogeneity on/off. If beta_alpha=1: time-varying unobserved heterogeneity is present
alpha_i<-rep(0,n)
lpm=0##If lpm=1: linear probability model, else: numeric
####the betas for the leads and lags###Extension for dynamic effects##
ll<-0######Indicator for lags: if ll=1: leads/lags present############
ifelse(ll==1,beta_lead<-seq(0,0,le=lead),beta_lead<-rep(0,lead))##Coefficient for leads
ifelse(ll==1,beta_lag<-seq(0.1,0.1,le=lag),beta_lag<-rep(0,lag))##Coefficients for lags
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
sim=100
x<-seq(1:10000)####Vector of random seeds: replace with custom seeds if desired####
seed.vec=sample(x,sim)###
######Run the simulations, data files will be saved in the data_simulations folder###############
sim.data.fun<-function(seed.vec){
  seed=seed.vec
  result<-group.fun(beta_age,lead,lag,ll)
  GM=result$GM
  alpha_gi=result$alpha_gi
  dummies=result$dummies
  alpha=result$alpha
  df<-result$df
  tic("sleeping")
  print("falling asleep...")
  #simu_fun(n,seed,beta,beta1,beta_alpha,t,G,GM,alpha_gi,alpha,age_new,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
  simu_fun_hs(n,seed,beta,beta1,beta_alpha,t,G,GM,dummies, alpha_gi,alpha_i,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm,hs=1)
  print("...waking up")
  toc()
}
#########generate the samples#########map(seed.vec,sim.data.fun)###
plan(multisession, workers = 4)#adapts the code in the map function to be distributed 
future_map(seed.vec,sim.data.fun,.progress = TRUE)###run the simulation function 
#####################
#####################read in the data files "hs" for heterogenous effects ####
filenames <-list.files(pattern=paste0(n,"_data_","hs_"))
filenames<-sample(filenames,length(filenames))#filenames[1:5]
######Read in the G_max, needs to be calculated beforehand
sigma2<- readRDS(paste0(n,"sigma2"))
##### do the main simulations 
main_sim<-function(filenames)
{
  b=beta
  DAtA<-as_tibble(read.csv(filenames))
  G_star<-1:7
  g_func<-function(G_star,DAtA){
    X0 = DAtA%>%select(Tr,starts_with("X")) 
    Y0 = DAtA$Y
    XX = X0#; % Useful if X0 contains additional variables. Change to number of variables used in regression
    YY = Y0
    X = as.matrix(XX);
    Y = as.numeric(YY);
    K = ifelse(dim(X)[2]>1,dim(X)[2],1)#
    n = length(Y)/t##number of observations
    results_tr<-BM_fun(G_star,sim=1,t,X,Y)
    gi_class=results_tr$gi_class
    theta_1<-results_tr$beta[1]
    theta_all<-results_tr$beta
    delta<-results_tr$delta
    theta<-cbind(theta_1,theta_all)
    delta_mat<-cbind(delta)
    se<-results_tr$se[1]
    ob_group<-opt_group(Y,X,G_star,theta_all,sigma2,n,t,K,delta,gi_class)
    BIC=ob_group$BIC
    BIC1=ob_group$BIC1
    write.table(cbind(G_star,theta_1),file=paste0(n,"_","hs","_",ll,"_",lpm,"_theta_1_rep",sep=""),col.names = F,row.names=F,append=TRUE)
    write.table(theta,file=paste0(n,"_","hs","_",ll,"_",lpm,"_theta_rep",sep=""),col.names = F,append=TRUE)
    write.table(se,file=paste0(n,"_","hs","_",ll,"_",lpm,"_se_rep",sep=""),col.names = F,append=TRUE)
    write.table(delta_mat,file=paste0(n,"_","hs","_",ll,"_",lpm,"_delta_mat_rep",sep=""),col.names = F,append=TRUE)
    write.table(BIC,file=paste0(n,"_","hs","_",ll,"_",lpm,"_BIC",sep=""),col.names = F,append=TRUE)
    write.table(BIC1,file=paste0(n,"_","hs","_",ll,"_",lpm,"_BIC1",sep=""),col.names = F,append=TRUE)
     }
  map(G_star, g_func,DAtA)
  ab_prob<-DAtA%>%group_by(t)%>%summarise(sum(Tr))
  write.table(cbind(1:t,ab_prob),file=paste0(n,"_","hs","_",ll,"_",lpm,"_abprob",sep=""),col.names = F,append=TRUE)
}

#####Distribute the functions to workers and run the estimation for each simulation##########
##### sample. The results are stored in the Data_simulation folder as 
plan(multisession, workers = 4)
future_map(filenames,main_sim,.progress = TRUE)
####calculate the OLS and the TWFE estimates ###
filenames <-list.files(pattern=paste0(n,"_data_","hs_"))
for(j in 1:length(filenames))
{
  DAtA<-read.csv(filenames[j])
  X0 = DAtA%>%select(Tr,starts_with("X")) 
  Y0 = DAtA$Y
  XX = X0#; % Useful if X0 contains additional variables. Change to number of variables used in regression
  YY = Y0
  X = as.matrix(XX);
  #X1<-as.matrix(X[,1])
  Y = YY;
  K = ifelse(dim(X)[2]>1,dim(X)[2],1)# 1 in our case, because only one covariate size, otherwise dim(X,2);##
  N = length(Y)/t##number of observations
  df<-data.frame(Y,X)
  result_0<-lm(Y~X,data=df)
  coeff_0<-result_0$coeff[2]
  print( coeff_0)
  se.ols<-coef(summary(result_0))[2, "Std. Error"]
  ols<-t(c(coeff_0,se.ols))
  df<-data.frame(df,DAtA$i,DAtA$t)
  # ####Estimate the fixed effect model###########################################################
  fe_estim <- plm(Y~Tr,data=df,index = c("DAtA.i","DAtA.t"), model = "within",effect="twoways")###fixed effects 
  print(fe_estim)
  fixef(fe_estim)
  fe_coef<-fe_estim$coeff
  wi_summary=summary(fe_estim)
  pval <- wi_summary[["coefficients"]][ , "Pr(>|t|)"]
  se.fe<-coef(summary(fe_estim))[, "Std. Error"]
  fe<-t(c(fe_coef,se.fe,pval))
  write.table(fe,file=paste0(n,"_","hs","_",ll,"_",lpm,"_fe",sep=""),col.names = F,append=TRUE)
  write.table(ols,file=paste0(n,"_","hs","_",ll,"_",lpm,"_ols",sep=""),col.names = F,append=TRUE)
  print(j)
  ab_prob<-DAtA%>%group_by(t)%>%summarise(sum(Tr))
  write.table(cbind(1:t,ab_prob),file=paste0(n,"_","hs","_",ll,"_",lpm,"_abprob",sep=""),col.names = F,append=TRUE)
}


#####Analyze mean results from simulations #############

theta_1<-as.data.frame(read.table((here("Codes","Data_Simulation",(paste0(n,"_","hs","_",ll,"_",lpm,"_theta_1_rep"))))))
OLS<-colMeans(read.table(here("Codes","Data_Simulation",(paste0(n,"_","hs","_",ll,"_",lpm,"_ols",sep=""))))) 
FE<-colMeans(read.table((here("Codes","Data_Simulation",(paste0(n,"_","hs","_",ll,"_",lpm,"_fe",sep=""))))))
ab_prob<-as.data.frame(read.table((here("Codes","Data_Simulation",(paste0(n,"_","hs","_",ll,"_",lpm,"_abprob"))))))

#####calculate the abortion probabilities over all simulation run in each period### 


help_ab<-ab_prob%>%group_by(V3)%>%summarize(mean_size = mean(V4, na.rm = TRUE))
res_ab<-round(help_ab[,2],0)/n

write.table(res_ab,file=paste0(n,"_","hs","_",ll,"_","ab_prob",sep=""),col.names = F,append=TRUE)

help<-theta_1%>%
  group_by(V1) %>%
  summarize(mean_size = mean(V2, na.rm = TRUE))
mn = c(OLS[2],FE[2],help$mean_size[1:7])
mn
write.table(res_coef,file=paste0(n,"_","hs","_",ll,"_","coef",sep=""),col.names = F,append=TRUE)


############Group Specific coefficients #################
library(furrr)
library(purrr)
library(tidyverse)
library(stringr)
#############In this file: how to run the simulations######################
############and apply the simulation file. Several dynamic effects can be included. 
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
n=1000#############Number of observations##########
t=10###number of time periods#############################################################################################
G=3 ###number of groups ###
c<-c(-0.05,0.05)###Coefficients for additional time-constant regressors###################################################
beta1<-sample(c,t,replace=T)
beta=c(0,-0.2,0.2)#########################The coefficient for abortion/parameter of interest, by group.--------------##########################
lead=5##Number of leads, maximum is 5############################################
lag=5###Number of lags, maximum is 5
beta_alpha<-5#####Turns time-varying unobserved heterogeneity on/off. If beta_alpha=1: time-varying unobserved heterogeneity is present
alpha_i<-rep(0,n)
lpm=0##If lpm=1: linear probability model, else: numeric
####the betas for the leads and lags###Extension for dynamic effects##
ll<-0######Indicator for lags: if ll=1: leads/lags present############
ifelse(ll==1,beta_lead<-seq(0,0,le=lead),beta_lead<-rep(0,lead))##Coefficient for leads
ifelse(ll==1,beta_lag<-seq(0.1,0.1,le=lag),beta_lag<-rep(0,lag))##Coefficients for lags
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
sim=100
x<-seq(1:10000)####Vector of random seeds: replace with custom seeds if desired####
seed.vec=sample(x,sim)###############
####need to specifiy n beforehand, I think this is not true, can pass other arguments afterwards #########
sim.data.fun<-function(seed.vec){
   seed=seed.vec
    result<-group.fun(beta_age,lead,lag,ll)
    GM=result$GM
    alpha_gi=result$alpha_gi
    dummies=result$dummies
    alpha=result$alpha
    df<-result$df
    tic("sleeping")
    print("falling asleep...")
    #simu_fun(n,seed,beta,beta1,beta_alpha,t,G,GM,alpha_gi,alpha,age_new,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
  simu_fun_new(n,seed,beta,beta1,beta_alpha,t,G,GM,dummies, alpha_gi,alpha_i,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm,gs=1)
     print("...waking up")
    toc()
  }
plan(multisession, workers = 4)
future_map(seed.vec,sim.data.fun)
#####################
#####################gs indicates group-specific##########
filenames <-list.files(pattern=paste0(n,"_data_","gs_"))
filenames<-sample(filenames,length(filenames))
######Read in the G_max, needs to be calculated beforehand
sigma2<- readRDS(paste0(n,"sigma2"))
###global parameters 
main_sim<-function(filenames)
{
  b=beta
  DAtA<-as_tibble(read.csv(filenames))
  G_star<-1:7
  g_func<-function(G_star,DAtA){
    X0 = DAtA%>%select(Tr,starts_with("X")) #
    Y0 = DAtA$Y
    XX = X0#; % Useful if X0 contains additional variables. Change to number of variables used in regression
    YY = Y0
    X = as.matrix(XX);
    Y = as.numeric(YY);
    K = ifelse(dim(X)[2]>1,dim(X)[2],1)#
    n = length(Y)/t##number of observations
    results_tr<-BM_fun(G_star,sim=1,t,X,Y)
    gi_class=results_tr$gi_class
    theta_1<-results_tr$beta[1]
    theta_all<-results_tr$beta
    delta<-results_tr$delta
    theta<-cbind(theta_1,theta_all)
    delta_mat<-cbind(delta)
    se<-results_tr$se[1]
    ob_group<-opt_group(Y,X,G_star,theta_all,sigma2,n,t,K,delta,gi_class)
    BIC=ob_group$BIC
    BIC1=ob_group$BIC1
    write.table(cbind(G_star,theta_1),file=paste0(n,"_","gs","_",ll,"_",lpm,"_theta_1_rep",sep=""),col.names = F,row.names=F,append=TRUE)
    write.table(theta,file=paste0(n,"_","gs","_",ll,"_",lpm,"_theta_rep",sep=""),col.names = F,append=TRUE)
    write.table(se,file=paste0(n,"_","gs","_",ll,"_",lpm,"_se_rep",sep=""),col.names = F,append=TRUE)
    write.table(delta_mat,file=paste0(n,"_","gs","_",ll,"_",lpm,"_delta_mat_rep",sep=""),col.names = F,append=TRUE)
    write.table(BIC,file=paste0(n,"_","gs","_",ll,"_",lpm,"_BIC",sep=""),col.names = F,append=TRUE)
    write.table(BIC1,file=paste0(n,"_","gs","_",ll,"_",lpm,"_BIC1",sep=""),col.names = F,append=TRUE)
  }
  map(G_star, g_func,DAtA)
  ab_prob<-DAtA%>%group_by(gi)%>%summarise(sum(Tr))
  write.table(cbind(1:3,ab_prob),file=paste0(n,"_","gs","_",ll,"_",lpm,"_abprob",sep=""),col.names = F,row.names=F,append=TRUE)
}
######run the simulations 
plan(multisession, workers = 4)
future_map(filenames,main_sim)
#######Calculate OLS and TWFE estimates #######
filenames <-list.files(pattern=paste0(n,"_data_","gs_"))
for(j in 1:length(filenames))
{
  DAtA<-read.csv(filenames[j])
  X0 = DAtA%>%select(Tr,starts_with("X")) 
  Y0 = DAtA$Y
  XX = X0#; % Useful if X0 contains additional variables. Change to number of variables used in regression
  YY = Y0
  X = as.matrix(XX);
  #X1<-as.matrix(X[,1])
  Y = YY;
  K = ifelse(dim(X)[2]>1,dim(X)[2],1)# 1 in our case, because only one covariate size, otherwise dim(X,2);##
  N = length(Y)/t##number of observations
  df<-data.frame(Y,X)
  result_0<-lm(Y~X,data=df)
  coeff_0<-result_0$coeff[2]
  print( coeff_0)
  # p_val<-summary(result_0)$coefficients[2,4] 
  se.ols<-coef(summary(result_0))[2, "Std. Error"]
  ols<-t(c(coeff_0,se.ols))
  df<-data.frame(df,DAtA$i,DAtA$t)
  # ####Estimate the fixed effect model###########################################################
  fe_estim <- plm(Y~Tr,data=df,index = c("DAtA.i","DAtA.t"), model = "within",effect="twoways")###fixed effects 
  print(fe_estim)
  # fe_estim <- plm(Y~X,data=df,index = c("ID","year"), model = "within",effect="twoways")###fixed effects 
  fixef(fe_estim)
  fe_coef<-fe_estim$coeff
  wi_summary=summary(fe_estim)
  pval <- wi_summary[["coefficients"]][ , "Pr(>|t|)"]
  se.fe<-coef(summary(fe_estim))[, "Std. Error"]
  fe<-t(c(fe_coef,se.fe,pval))
  write.table(fe,file=paste0(n,"_","gs","_",ll,"_",lpm,"_fe",sep=""),col.names = F,append=TRUE)
  write.table(ols,file=paste0(n,"_","gs","_",ll,"_",lpm,"_ols",sep=""),col.names = F,append=TRUE)
  print(j)
  # colSums(summary(fixef(fe_estim))[ , c("Estimate", "Pr(>|t|)")]) # only estimates and p-values
}

#######################Read in simulation results #######
theta_1<-as.data.frame(read.table((here("Codes","Data_Simulation",(paste0(n,"_","gs","_",ll,"_",lpm,"_theta_1_rep"))))))
OLS<-colMeans(read.table(here("Codes","Data_Simulation",(paste0(n,"_","gs","_",ll,"_",lpm,"_ols",sep=""))))) 
FE<-colMeans(read.table((here("Codes","Data_Simulation",(paste0(n,"_","gs","_",ll,"_",lpm,"_fe",sep=""))))))


help_ab<-ab_prob%>%group_by(V3)%>%summarize(mean_size = mean(V4, na.rm = TRUE))
res_ab<-round(help_ab[,2],0)/n

write.table(res_ab,file=paste0(n,"_","gs","_",ll,"_","ab_prob",sep=""),col.names = F,append=TRUE)

help<-theta_1%>%
  group_by(V1) %>%
  summarize(mean_size = mean(V2, na.rm = TRUE))
mn = c(OLS[2],FE[2],help$mean_size[1:7])
mn
write.table(res_coef,file=paste0(n,"_","gs","_",ll,"_","coef",sep=""),col.names = F,append=TRUE)


######Baseline: Homogenous treatment effects #############
library(furrr)
library(purrr)
library(tidyverse)
library(stringr)

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
n=1000#############Number of observations##########
t=10###number of time periods#############################################################################################
c<-c(-0.05,0.05)###Coefficients for additional time-constant regressors###################################################
beta1<-sample(c,t,replace=T)
beta=0#########################The coefficient for abortion/parameter of interest --------------##########################
lead=5##Number of leads, maximum is 5############################################
lag=5###Number of lags, maximum is 5
beta_alpha<-5#####Turns time-varying unobserved heterogeneity on/off. If beta_alpha=1: time-varying unobserved heterogeneity is present
alpha_i<-rep(0,n)
lpm=0##If lpm=1: linear probability model, else: numeric
####the betas for the leads and lags###Extension for dynamic effects##
ll<-0######Indicator for lags: if ll=1: leads/lags present############
ifelse(ll==1,beta_lead<-seq(0,0,le=lead),beta_lead<-rep(0,lead))##Coefficient for leads
ifelse(ll==1,beta_lag<-seq(0.1,0.1,le=lag),beta_lag<-rep(0,lag))##Coefficients for lags
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
sim=100
x<-seq(1:10000)####Vector of random seeds: replace with custom seeds if desired####
seed.vec=sample(x,sim)###############
####need to specifiy n beforehand, I think this is not true, can pass other arguments afterwards #########
sim.data.fun<-function(seed.vec){
  seed=seed.vec
  result<-group.fun(beta_age,lead,lag,ll)
  GM=result$GM
  alpha_gi=result$alpha_gi
  dummies=result$dummies
  alpha=result$alpha
  df<-result$df
  tic("sleeping")
  print("falling asleep...")
  simu_fun_2(n,seed,beta,beta1,beta_alpha,t,G,GM,dummies, alpha_gi,alpha_i,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
  print("...waking up")
  toc()
}
######Run the simlations, data files will be saved in the data_simulations folder###############
plan(multisession, workers = 4)
future_map(seed.vec,sim.data.fun)
#####################
####################Read in simulated data ###########
filenames <-list.files(pattern=paste0(n,"_data_"))
filenames<-sample(filenames,length(filenames))#
######Read in the G_max, needs to be calculated beforehand
sigma2<- readRDS(paste0(n,"sigma2"))
###global parameters 
main_sim<-function(filenames)
{
  b=beta
  DAtA<-as_tibble(read.csv(filenames[1]))
  G_star<-1:7
  g_func<-function(G_star,DAtA){
    X0 = DAtA%>%select(Tr,starts_with("X")) #for now include the X_p's and the treatment, not the leads and lags
    Y0 = DAtA$Y
    XX = X0#; % Useful if X0 contains additional variables. Change to number of variables used in regression
    YY = Y0
    X = as.matrix(XX);
    Y = as.numeric(YY);
    K = ifelse(dim(X)[2]>1,dim(X)[2],1)#
    n = length(Y)/t##number of observations
    results_tr<-BM_fun(G_star,sim=1,t,X,Y)
    gi_class=results_tr$gi_class
    theta_1<-results_tr$beta[1]
    theta_all<-results_tr$beta
    delta<-results_tr$delta
    theta<-cbind(theta_1,theta_all)
    delta_mat<-cbind(delta)
    se<-results_tr$se[1]
    ob_group<-opt_group(Y,X,G_star,theta_all,sigma2,n,t,K,delta,gi_class)
    BIC=ob_group$BIC
    BIC1=ob_group$BIC1
    write.table(cbind(G_star,theta_1),file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_theta_1_rep",sep=""),col.names = F,row.names=F,append=TRUE)
    write.table(theta,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_theta_rep",sep=""),col.names = F,append=TRUE)
    write.table(se,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_se_rep",sep=""),col.names = F,append=TRUE)
    write.table(delta_mat,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_delta_mat_rep",sep=""),col.names = F,append=TRUE)
    write.table(BIC,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_BIC",sep=""),col.names = F,append=TRUE)
    write.table(BIC1,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_BIC1",sep=""),col.names = F,append=TRUE)
  }
  map(G_star, g_func,DAtA)
}

######Now implement this with furrr ###
plan(multisession, workers = 4)
future_map(filenames,main_sim)
b=beta
filenames <-list.files(pattern=paste0(n,"_data_",beta==b))
for(j in 1:length(filenames))
{
  DAtA<-read.csv(filenames[j])
  X0 = DAtA%>%select(Tr,starts_with("X")) 
  Y0 = DAtA$Y
  XX = X0#; % Useful if X0 contains additional variables. Change to number of variables used in regression
  YY = Y0
  X = as.matrix(XX);
  #X1<-as.matrix(X[,1])
  Y = YY;
  K = ifelse(dim(X)[2]>1,dim(X)[2],1)# 1 in our case, because only one covariate size, otherwise dim(X,2);##
  N = length(Y)/t##number of observations
  df<-data.frame(Y,X)
  result_0<-lm(Y~X,data=df)
  coeff_0<-result_0$coeff[2]
  print( coeff_0)
  # p_val<-summary(result_0)$coefficients[2,4] 
  se.ols<-coef(summary(result_0))[2, "Std. Error"]
  ols<-t(c(coeff_0,se.ols))
  df<-data.frame(df,DAtA$i,DAtA$t)
  # ####Estimate the fixed effect model###########################################################
  fe_estim <- plm(Y~Tr,data=df,index = c("DAtA.i","DAtA.t"), model = "within",effect="twoways")###fixed effects 
  print(fe_estim)
  # fe_estim <- plm(Y~X,data=df,index = c("ID","year"), model = "within",effect="twoways")###fixed effects 
  fixef(fe_estim)
  fe_coef<-fe_estim$coeff
  wi_summary=summary(fe_estim)
  pval <- wi_summary[["coefficients"]][ , "Pr(>|t|)"]
  se.fe<-coef(summary(fe_estim))[, "Std. Error"]
  fe<-t(c(fe_coef,se.fe,pval))
  write.table(fe,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_fe",sep=""),col.names = F,append=TRUE)
  write.table(ols,file=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_ols",sep=""),col.names = F,append=TRUE)
  print(j)
  # colSums(summary(fixef(fe_estim))[ , c("Estimate", "Pr(>|t|)")]) # only estimates and p-values
}

theta_1<-as.data.frame(read.table((here("Codes","Data_Simulation",(paste0(n,"_",beta==b,"_",ll,"_",lpm,"_theta_1_rep"))))))
OLS<-colMeans(read.table(here("Codes","Data_Simulation",(paste0(n,"_",beta==b,"_",ll,"_",lpm,"_ols",sep=""))))) 
FE<-colMeans(read.table((here("Codes","Data_Simulation",(paste0(n,"_",beta==b,"_",ll,"_",lpm,"_fe",sep=""))))))


help_ab<-ab_prob%>%group_by(V3)%>%summarize(mean_size = mean(V4, na.rm = TRUE))
res_ab<-round(help_ab[,2],0)/n

write.table(res_ab,file=paste0(n,"_","hm","_",ll,"_","ab_prob",sep=""),col.names = F,append=TRUE)

help<-theta_1%>%
  group_by(V1) %>%
  summarize(mean_size = mean(V2, na.rm = TRUE))
mn = c(OLS[2],FE[2],help$mean_size[1:7])
mn
write.table(res_coef,file=paste0(n,"_","hm","_",ll,"_","coef",sep=""),col.names = F,append=TRUE)

