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
t=10###number of time periods#############################################################################################
c<-c(-0.05,0.05)###Coefficients for additional time-constant regressors###################################################
beta1<-sample(c,t,replace=T)
beta=0#########################The coefficient for abortion/parameter of interest --------------##########################
lead=5##Number of leads, maximum is 5############################################
lag=5###Number of lags, maximum is 5
beta_alpha<-0#####Turns time-varying unobserved heterogeneity on/off. If beta_alpha=1: time-varying unobserved heterogeneity is present
n=1000#############Number of observations##########
lpm=0##If lpm=1: linear probability model, else: numeric
####the betas for the leads and lags###Extension for dynamic effects##
ll<-1######Indicator for lags: if ll=1: leads/lags present############
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
##############################################################################
########Finish global parameters##############################################
##----------------------------Run the simulations: sim times------------------
##############################################################################
sim=10
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
    age_new=result$age_new
    df<-result$df
    tic("sleeping")
    print("falling asleep...")
    simu_fun(n,seed,beta,beta1,beta_alpha,t,G,GM,alpha_gi,alpha,age_new,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
    print("...waking up")
    toc()
  }
}
######Run the simlations, data files will be saved in the data_simulations folder###############
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
  filenames <- list.files(pattern=paste0(n,"_",beta==b,"_",ll,"_",lpm,"_", "_data_lag.csv"))
  filenames<-sample(filenames,length(filenames))
  ######Read in the G_max, needs to be calculated beforehand
  sigma2<- readRDS(paste0(n,"sigma2"))
  run=length(filenames)###Number of calculations for the BIC estimator, maximum is the number of simulated files#####
  for(j in 1:run)
  {
    DAtA<-read.csv(filenames[j])
    delta_mat<-matrix(NA,t,1)
    theta_1<-c()
    theta<-matrix(NA,11,1)
    se<-c()
    G_star<-1:5
    BIC<-c()
    BIC1<-c()
    obj<-c()
    for(i in 1:length(G_star))
    {
      X0 = cbind(DAtA$V4 ,DAtA[,9:18])
      Y0 = DAtA$Y
      XX = X0#; % Useful if X0 contains additional variables. Change to number of variables used in regression
      YY = Y0
      X = as.matrix(XX);
      Y = YY;
      K = ifelse(dim(X)[2]>1,dim(X)[2],1)#
      n = length(Y)/t##number of observations
      results_tr<-BM_fun(G_star[i],sim=1,10,X,Y)
      gi_class=results_tr$gi_class
      theta_1[i]<-results_tr$beta[1]
      theta_all<-results_tr$beta
      delta<-results_tr$delta
      theta<-cbind(theta,theta_all)
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
b=0######is parameter of interest equal to zero
ll=0###lags yes/no

main.sim(1000,b,ll)
main.sim(1500,b,ll)
main.sim(2000,b,ll)
main.sim(2500,b,ll)
main.sim(3000,b,ll)


