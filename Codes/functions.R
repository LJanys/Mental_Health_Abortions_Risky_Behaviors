
BM_fun<-function(G,sim=3,t,X,Y)#The main estimator, including standard errors. 
{
  #######This function code was adapted to R from the Matlab code by Bonhomme and Manresa (2015) grouped fixed-effects estimator###
  K = ifelse(dim(X)[2]>1,dim(X)[2],1)# Number of covariates
  N = length(Y)/t##number of observations
  Resbeta=matrix(0,sim,K)#;empty matrices to store results of the iterative estimator
  Resdelta=matrix(0,sim,G*t)#;
  ResQ=matrix(0,sim,1)#;
  ####Replace missing information with zeros in both the response and the independent variable(s). 
  y.miss<-matrix(Y,nrow=N,ncol=t,byrow = T)
  U_miss = matrix(0,N,t)
  for (i in 1:N)
  {
    for (j in 1:t)
    {
      ifelse(is.na(y.miss[i,j])==1, (U_miss[i,j] = 1), (U_miss[i,j] = 0))
    }
  } 
  Y_old<-Y
  X_old<-X
  Y[is.na(Y)==1]=0
  X[is.na(X)==1]=0
  DAtA_nonmiss=DAtA
  DAtA_nonmiss[is.na(DAtA_nonmiss)==1]=0
  gi_auxaux=matrix(0,N,G)
  gi=matrix(0,N,G)
  gitot=matrix(0,N*t,G*t)
  N_nonmiss=N###von mir
  for (jsim in 1:sim)
  {
    ####Random starting values are drawn from a std. normal distribution
    beta_init=rnorm(K)
    ###Calculate residuals for either multiple regressors or one regressor####
    ifelse((dim(X)[2]>1),  (W=Y-X%*%beta_init),  (W=Y-X*beta_init))
    # select randomly G "centers" among the individuals 
    V0=sample(N,N_nonmiss) 
    V=V0[1:G] 
    delta_init=matrix(0,t,G)
    for (g in 1:G)
    {
      a=(V[g]-1)*t+1
      b= V[g]*t
      delta_init[,g]=W[a:b]  ###draws the residuals of the randomly chosen centers
    }
    par_init= c(as.vector(delta_init),beta_init)
    # Iterations
    deltapar=1
    while(deltapar>0)
    {
      # Step 1: assignment (classification into groups)
      ###also needs adjustment when X has more than one element##
      for (g in 1:G)
      {
        ifelse(dim(X)[2]>1,(U=(as.vector(Y-X%*%beta_init)-kronecker(rep(1,N),delta_init[,g]))^2),(U = (as.vector(Y-X*beta_init)-kronecker(rep(1,N),delta_init[,g]))^2))#;###;## )
        RU = t(matrix(U,t,N))
        gi_auxaux[,g] = t(rowSums((RU)))###sum residuals for each individual over all time periods
      }
      gi_class = colMins(t(gi_auxaux))## Integer, column indicator for the class with the smallest residual. 
      
      for (g in 1:G)
      {
        gi[,g] = gi_class==g
      }
      
      # Step 2: update parameters (ols)
      for (i in 1:N)
      {
        giaux=t(kronecker(diag(t),gi[i,]))
        gitot[((i-1)*t+1):(t*i),] = giaux
      }
      
      Xtot<-cbind(gitot, X)
      par_new=pinv(Xtot)%*%Y
      A=par_new-par_init
      deltapar=norm(A,type="2")#one value###
      par_init=par_new
      delta_init_vect=par_init[1:(G*t)]
      delta_init=t(matrix(delta_init_vect,G,t))
      ############################################################################################################        
      beta_init=par_init[(G*t+1):(G*t+K)]
      print(beta_init)
    }
    Resbeta[jsim,]=beta_init
    Resdelta[jsim,]=delta_init_vect
    ResQ[jsim]=t(Y-Xtot%*%par_init)%*%(Y-Xtot%*%par_init)
  }
  #Pick the parameters that minimize the sum of squared residuals.
  Resbeta_final = Resbeta####Matrix containing all 
  Resbetadelta_final = Resdelta#
  ResQfinal = ResQ#
  
  temp<-cbind(ResQfinal,rep(1:sim))
  temp<-as.data.frame(temp)
  foo=temp
  foo.sorted=foo[order(foo[,1]),]
  x1<-foo.sorted[,1]
  x2<-foo.sorted[,2]
  Resbeta[x2,]
  beta_final=Resbeta_final[x2[1],]
  delta_final=Resbetadelta_final[x2[1],]
  delta_final_r=t(matrix(delta_final,G,t))
  par_final=c(delta_final,beta_final)
  ###assignment to groups######################################                                                      
  for(g in 1:G)
  {
    ifelse(dim(X)[2]>1,(U=(as.vector(Y-X%*%beta_final)-kronecker(rep(1,N),delta_final_r[,g]))^2),(U = (as.vector(Y-X*beta_final)-kronecker(rep(1,N),delta_final_r[,g]))^2))
    RU = t(matrix(U,t,N))
    gi_auxaux[,g] = (rowSums(RU))###sum residuals for each individual over all time periods
  }
  gi_class = colMins(t(gi_auxaux))
  for (g in 1:G)
  {
    gi[,g] = gi_class==g
  }
  # update
  for (i in 1:N)
  {
    giaux=t(kronecker(diag(t),gi[i,]))
    gitot[((i-1)*t+1):(t*i),] = giaux
  }
  Xtot_final<-cbind(gitot, X)
  par_final=pinv(Xtot_final)%*%Y
  delta_final_vect=par_final[1:(G*t)]
  delta_final=t(matrix(delta_final_vect,G,t))
  ####################################     
  beta_final=par_final[(G*t+1):(G*t+K)]
  DAtA_new = DAtA_nonmiss
  ginumb<-gi_class
  ## Standard error  estimation %%
  Ybar_gt<-matrix(0,N*t,1)
  #Ybar_gt=zeros(N*T,1);
  Xbar_gt<-matrix(0,N*t,K)
  #Xbar_gt=zeros(N*T,K);
  MYbar_gt=matrix(0,G*t,1)
  #MYbar_gt=zeros(G*T,1);
  MXbar_gt=matrix(0,G*t,K)
  #=zeros(G*T,K);
  #gi=zeros(N,G);
  gi<-matrix(0,N,G)
  for(g in 1:G)
  {
    gi[,g]=ginumb==g
  }
  gisum=colSums(gi)
  for(i in 1:N)
  {
    
    for (j in 1:t){ Yt=Y[seq(j,N*t,by=t)]
    Ybar_gt[(i-1)*t+j]=mean(Yt[ginumb==ginumb[i]])
    Xt=as.matrix(X[seq(j,N*t,by=t),],N*t,K)
    Xbar_gt[(i-1)*t+j,]=mean(Xt[ginumb==ginumb[i],])
    }
  }
  # Attention: This is the large - T variance: We have validated coverage rates for our example in our simulation set-up. 
  ###delta_final###
  theta_par<-beta_final
  ei=Y-Ybar_gt-(X-Xbar_gt)%*%theta_par
  Rei = t(matrix(ei,t,N)) 
  Omega =matrix(0,N*t,N*t)
  Mi=matrix(0,t,t)
  for(i in 1:N)
  {
    Mi =(Rei[i,])%*%t(Rei[i,])
    Omega[((i-1)*t+1):(i*t),((i-1)*t+1):(i*t)] = Mi
  }
  gitot=kronecker(gi,diag(t))
  Xtot<-cbind(gitot, X)
  V= solve(t(Xtot)%*%Xtot)%*%t(Xtot)%*%Omega%*%Xtot%*%solve(t(Xtot)%*%Xtot)
  V=V*N*t/(N*t-G*t-K)
  std_cluster = sqrt(diag(V))
  # standard errors
  std_cluster2=std_cluster[(G*t+1):(G*t+K)]
  std_cluster[(G*t+1):(G*t+K)]
  return(list(beta=beta_final,delta=delta_final, se=std_cluster2,Ybar_gt=Ybar_gt,Xbar_gt=Xbar_gt,Y=Y,X=X,gi_class=gi_class))
}

opt_group_G_max<-function(Y,X,G,theta_par,N,t,K,delta)
{
  obj1<-matrix(NA,N,t)
  for(i in 1:N)
  {
    for(j in 1:t)
    {
      obj1[i,j]<-(Y[(i-1)*t+j]-X[(i-1)*t+j,]%*%theta_par-delta[j,gi_class[i]])^2 
    }
  }
  obj<-sum(obj1)
  sigma2<-obj/(N*t-G*t-N-K)
  ###
  sigma22 = mean(sigma2)
  BIC =obj/(N*t)+sigma2*log(N*t)*G(t+N-G)/(N*t)
  BIC = obj/(N*t)+sigma22*log(N*t)*G(t+N-G)/(N*t)
  return(list(BIC=BIC,sigma2=sigma2,obj=obj))
}

#####################Function where the the groups and the demeaned X and Y are calculated using the training deltas and betas. 
gi_fun<-function(G, theta, delta, X,Y)
{
  K = ifelse(dim(X)[2]>1,dim(X)[2],1)# 1 in our case, because only one covariate size, otherwise dim(X,2);##
  N = length(Y)/t##number of observations
  gi_auxaux=matrix(0,N,G)
  gi=matrix(0,N,G)
  gitot=matrix(0,N*t,G*t)
  beta_final=theta#c(theta_1[2,2],rep(0.5,le=10))
  delta_final<-as.matrix(delta)#cbind(delta_rep$V3,delta_rep$V4)
  par_final=c(delta_final,beta_final)
  gi_auxaux=matrix(0,N,G)
  gi=matrix(0,N,G)
  gitot=matrix(0,N*t,G*t)
  ###assignment to groups######################################                                                      
  for(g in 1:G)
  {
    ifelse(dim(X)[2]>1,(U=(as.vector(Y-X%*%beta_final)-kronecker(rep(1,N),delta_final[,g]))^2),(U = (as.vector(Y-X*beta_final)-kronecker(rep(1,N),delta_final[,g]))^2))
    # U = (Y-X%*%beta_final-kronecker(rep(1,N),delta_final_r[,g]))^2#;##
    RU = t(matrix(U,t,N))####dimensions checken###
    gi_auxaux[,g] = (rowSums(RU))###sum residuals for each individual over all time periods, dim(N*g)
  }
  
  gi_class = colMins(t(gi_auxaux))
  for (g in 1:G)
  {
    gi[,g] = gi_class==g
  }
  for (i in 1:N)
  {
    giaux=t(kronecker(diag(t),gi[i,]))
    gitot[((i-1)*t+1):(t*i),] = giaux
  }
  Xtot_final<-cbind(gitot, X)
  par_final=pinv(Xtot_final)%*%Y
  delta_final_vect=par_final[1:(G*t)]
  delta_final=t(matrix(delta_final_vect,G,t))
  ####################################     
  beta_final=par_final[(G*t+1):(G*t+K)]
  ##Only for example purposes##
  ginumb<-gi_class
  ## Standard error  estimation %%
  Ybar_gt<-matrix(0,N*t,1)
  #Ybar_gt=zeros(N*T,1);
  Xbar_gt<-matrix(0,N*t,K)
  #Xbar_gt=zeros(N*T,K);
  MYbar_gt=matrix(0,G*t,1)
  #MYbar_gt=zeros(G*T,1);
  MXbar_gt=matrix(0,G*t,K)
  #=zeros(G*T,K);
  #gi=zeros(N,G);
  gi<-matrix(0,N,G)
  for(g in 1:G)
  {
    gi[,g]=ginumb==g
  }
  
  gisum=colSums(gi)
  for(i in 1:N)
  {
    #if(gisum[ginumb[i]]>1)
    #{
    for (j in 1:t){ Yt=Y[seq(j,N*t,by=t)]
    Ybar_gt[(i-1)*t+j]=mean(Yt[ginumb==ginumb[i]])
    Xt=as.matrix(X[seq(j,N*t,by=t),],N*t,K)
    Xbar_gt[(i-1)*t+j,]=mean(Xt[ginumb==ginumb[i],])
    }
  }
  return(list( Ybar_gt=Ybar_gt,Xbar_gt=Xbar_gt,gi_class=gi_class))
}

################Function calculating both BIC's using sigma2_max from the G_max as input###########
opt_group<-function(Y,X,G,theta_par,sigma2,N,t,K,delta,gi_class)####
{
  obj1<-matrix(NA,N,t)
  delta=as.matrix(delta)
  for(i in 1:N)
  {
    for(j in 1:t)
    {
      obj1[i,j]<-(Y[(i-1)*t+j]-X[(i-1)*t+j,]%*%theta_par-delta[j,gi_class[i]])^2 
    }
  }
  obj<-sum(obj1)
  BIC1 = obj/(N*t)+sigma2*log(N*t)*(G*t+N+K)/(N*t)##less steep penalty
  BIC = obj/(N*t)+sigma2*log(N*t)*G*(t+N-G+K)/(N*t)#steeper penalty
  return(list(BIC=BIC,BIC1=BIC1,obj=obj))
}

##############################################################################################
####Function: Calculates the group function, that is an input 
group.fun<-function(beta_age,lead,lag,ll){
  G=3
  #########group membership equation############
  age<-round(rnorm(n,28,10),0)
  age[age<16]=16
  age[age>40]=40
  sum(age<16)
  sum(age>40)
  tt<-c()
  age_new<-rep(age,each=t)
  alpha<-25+age*beta_age+rnorm(n,0,5)
  alpha.grid<-quantile(alpha, c(0,.5,0.9,1))###determines the share of the groups. 
  GM=c()
  ####Assigns group status based on alpha##########
  for(i in 1:n)
  {
    if(alpha[i]<alpha.grid[2])
    {
      GM[i]<-1
    }
    if(alpha[i]>alpha.grid[2]&alpha[i]<alpha.grid[3])
    {
      GM[i]<-2
    }
    if(alpha[i]>alpha.grid[3])
    {
      GM[i]<-3
    }
  }
  #####Assigns the time of abortion based on group membership##########
  #####Could be a more flexible specification: make the age and the abortion probability be determined proportionally to the unobserved 
  #####mental health curves##############
  for(i in 1:n)
  {
    if(alpha[i]<alpha.grid[2])
    {
      x<-sort(c(1:(5),1:(5)))
      tt[i]<-sample(x,1,replace=T)
    }
    if(alpha[i]>alpha.grid[2]&alpha[i]<alpha.grid[3])
    {
      x<-sort(c(5:(t),1:t))
      tt[i]<-sample(x,1,replace=T)
    }
    if(alpha[i]>alpha.grid[3])
    {
      x<-sort(c(5:(t),3:(t)))
      tt[i]<-sample(x,1,replace=T)
    }
  }
  X<-factor(tt)
  dummies = model.matrix(~X-1)
  ####Now we can add another routine that determines whether or not someone actually did have an abortion
  ####based on their group membership#####################
  A=c()
  for(i in 1:n)
  {
    if (GM[i]==1)
    {
      pi=0.2
      A[i]<-rbinom(1,1,pi)
    }
    if (GM[i]==2)
    {
      pi=0.2
      A[i]<-rbinom(1,1,pi)
    }
    if (GM[i]==3)
    {
      pi=0.6
      A[i]<-rbinom(1,1,pi)
    }
  }
  dummies = model.matrix(~X-1)
  dummies<-dummies*A
  ####assign individual id's###########################################
  dummies1<-data.frame(cbind(1:n,dummies))
  names(dummies1)<- c("ID",sprintf("%s",seq(1:t)))
  ######re-shape into long format to get the lags and leads##########
  long<- dummies1 %>%
    pivot_longer(!ID, names_to = "year")
  df<-data.table(long)
  df[, sprintf("mid_o_lead_%0d", 1:lead) := shift(value, c(1:lead), type = 'lead'),by=ID]
  df[, sprintf("mid_o_lag_%0d", 1:lag) := shift(value, c(1:lag), type = 'lag'),by=ID]
  ####I replace the NA's with zeros.#########################################
  df[is.na(df)==T]=0
  ###########################################
  alpha_gi=matrix(0,t,G)
  a=0
  alpha_gi[,1]=0.01
  alpha_gi[,2]=-1+exp((1:t)/10)^{1/2}
  alpha_gi[,3]=-1+exp((1:t)/10)^{1.2}
  #####  #####  #####  #####  #####  #####  #####  #####  #####  #####  #####  #####  #####
  ###Plot the curves of the unobserved heterogeneity###
  plot(alpha_gi[,1],type="l",lty=1,lwd=2,ylim=c(0,1),ylab="P(MH)",xlab="Time")
  lines(alpha_gi[,2],lty=2,lwd=2,col="red")
  lines(alpha_gi[,3],lty=3,lwd=2,col="blue")
  ##This is the group function ending: return: dummies, group assignment. 
  return(list(GM=GM,dummies=dummies, alpha_gi=alpha_gi,df=df))
}


simu_fun<-function(n,seed,beta,beta1,beta_alpha,t,G,GM,alpha_gi,alpha,age_new,df,beta_lead,beta_lag,lead,lag,dummies_t,lpm)
{
  X_t<-matrix(rep(rbinom(n*t,1,0.5),t),n,t)
  X_p<-matrix(NA,nrow=n*t,ncol=t)##This generates some exogenous, time-constant regressors.
  for(i in 1:t)
  {
    X_p[,i]<-rep(X_t[,i],each=t)
  }
  ######make long versions of the individual fixed effects and the group####  
  alpha=matrix(rep(t(alpha),n),ncol=ncol(alpha),byrow=TRUE)
  GM=rep(GM,each=t)
  Y<-c()
  k<-n*t
  pi_x=c()####empty vector generate the predicted probabilities that are bounded by 0 and 1######################
  Y_star<-c()##Potentially latent  model/predicted values, otherwise actual dependent var value####  
  ###transform the 
  for(i in 1:k)
  {
    Y_star[i]=dummies.age[i,]%*%beta_age_2+dummies_t[i,]%*%beta_time+beta*df[i,3]+beta_alpha*alpha[i,GM[i]]+t(X_p[i,])%*%beta1+t(df[i,4:(4+lead-1)])%*%beta_lead+t(df[i,(4+lead):(4+lead+lag-1)])%*%beta_lag+lpm*runif(1,-0.5,0.5)+(1-lpm)*rnorm(1,0,1)
    #######age_dummies+ time dummies+main effect+
    pi_x[i]=Y_star[i]
    if(Y_star[i]<0)
    {
      pi_x[i]=0
    }
    if(Y_star[i]>1)
    {
      pi_x[i]=1
    }
    Y[i]=rbinom(1, size=1, prob=pi_x[i])
  }
  print(sum(Y==1))
  
  #print(cor(X1[,3],X1[,6]))####this gives the correlation between the event (value) and the unobserved 
  # plot(density(Y_star))
  ##############without lpm, Y is the actual values, with lpm=1: Y is a binary value, generated from the latent model##################
  ifelse(lpm==1,Y<-Y, Y<-Y_star)
  df_result<-as.data.frame(cbind(Y,df))
  ######save the data csv forma###################################################################################################################################
  #write.csv(df,file=here("Data_Simulation",paste0(seed,n,"_",beta==0,"_",ll,"_",lpm,"_", "_data_lag.csv",sep="")))###Comment out if not want to save the data set
  return(df_result=df_result)
}
