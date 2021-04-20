############---------------------------Main analysis file------------------------------------############################
###########---------------------------for the simulations--------------------------------------------------------########
##########This file generates all plots from the simulations-----------------------------------------------------########
library("plm")
library("dplyr")
library("ggplot2")
library(data.table)
library(tictoc)
library(tidyr)
library(Rfast)
library(pracma)
library(patchwork)
library(here)
library(readr)
library(dotwhisker)
library(broom)
setwd(here("Codes","Data_Simulation"))
####Summarizing the changes in the coefficients####
###First for N=1000###############################
main.ana<-function(n,b,ll)
{
  theta_1<-as.data.frame(read.table((here("Codes","Data_Simulation",(paste0(n,"_",beta==b,"_",ll,"_",lpm,"_theta_1_rep"))))))
  OLS<-colMeans(read.table(here("Codes","Data_Simulation",(paste0(n,"_",beta==b,"_",ll,"_",lpm,"_ols",sep=""))))) 
  FE<-colMeans(read.table((here("Codes","Data_Simulation",(paste0(n,"_",beta==b,"_",ll,"_",lpm,"_fe",sep=""))))))
  help<-theta_1%>%
    group_by(V1) %>%
    summarize(mean_size = mean(V2, na.rm = TRUE))
  mn = c(OLS[2],FE[2],help$mean_size[2:5])
  #######We can make the analogous plot with Bettinas code#####################
  se<-as.data.frame(read.table(here("Codes","Data_Simulation",((paste0(n,"_",beta==b,"_",ll,"_",lpm,"_se_rep"))))))
  help<-se%>%
    group_by(V1) %>%
    summarize(mean_size = mean(V2, na.rm = TRUE))
  se = c(OLS[3],FE[3],help$mean_size[2:5])
  name = c("OLS", "OLS FE", "GFE G2", "GFE G3", "GFE G4", "GFE G5")
  ci95_low = mn - 1.96*se
  ci95_up = mn + 1.96*se
  dat2 <-data.frame(name, mn,se,ci95_low, ci95_up)
  dat <- data.frame(name,mn,se)
  str(dat)
  pdf(file=here("Abortion_Draft",paste0(n,"_lag","_coeff_plot_sim.pdf")))
  print(ggplot(dat, aes(name, mn)) + 
          theme_bw()+
          geom_hline(yintercept=0, lty=2, lwd=1, colour="black")+
          geom_errorbar(aes(ymin=mn - 1.96*se, ymax=mn + 1.96*se), lwd=1, colour="red", width=0)+
          geom_point(size=4, pch=21, fill="black") +
          labs(x="model coefficients", y="magnitude") +
          scale_x_discrete(limits = c("OLS", "OLS FE", "GFE G2", "GFE G3", "GFE G4","GFE G5")))
  dev.off()
  
  ########BIC####################################################################
  BIC<-as.data.frame(read.table((here("Codes","Data_Simulation",(paste0(n,"BIC"))))))
  help<-BIC%>%
    group_by(V1) %>%
    summarize(mean_size = mean(V2, na.rm = TRUE))
  highlight.gene <- 3
  help$highlight <- c("normal","normal","highlight","normal","normal")#ifelse(a$GeneName == highlight.gene, "highlight", "normal")
  mycolours <- c("highlight" = "red", "normal" = "grey50")
  textdf <- help[help$V1==3,]
  
  pdf(file=here("Abortion_Draft",paste0(n,"_lag","sim_BIC.pdf")))
  print(ggplot(data = help, aes(x = V1, y = mean_size)) +
          geom_point(size = 4, aes(colour = highlight)) +
          scale_color_manual("Status", values = mycolours) +
          geom_text(data = textdf, aes(x = V1 * 1.2, y = mean_size, label = "True G")) +
          theme(legend.position = "none") +
          theme()+
          labs(x="Number of Groups", y="BIC")  +
          ylim(2,14))
  # ylim((min(help$mean_size)-1.5),(max(help$mean_size)+1.5)))
  dev.off()
  
  ########BIC 1####################################################################
  BIC1<-as.data.frame(read.table((here("Codes","Data_Simulation",(paste0(n,"BIC1"))))))
  help<-BIC1%>%
    group_by(V1) %>%
    summarize(mean_size = mean(V2, na.rm = TRUE))
  highlight.gene <- 3
  help$highlight <- c("normal","normal","highlight","normal","normal")#ifelse(a$GeneName == highlight.gene, "highlight", "normal")
  mycolours <- c("highlight" = "red", "normal" = "black")
  textdf <- help[help$V1==3,]
  
  pdf(file=here("Abortion_Draft",paste0(n,"_lag","sim_BIC1.pdf")))
  print(ggplot(data = help, aes(x = V1, y = mean_size)) +
          geom_point(size = 4, aes(colour = highlight)) +
          scale_color_manual("Status", values = mycolours) +
          geom_text(data = textdf, aes(x = V1 * 1.1, y = mean_size, label = "True G")) +
          theme(legend.position = "none") +
          theme()+
          ylim(2,14)+
          labs(x="Number of Groups", y="BIC") )
  dev.off()
  ############################################################################
  # #2000_FALSE_1_delta_mat_rep
  # delta_mat1<-as.data.frame(read.table((here("Codes","Data_Simulation",(paste0(n,"_",beta==b,"_",ll,"_","delta_mat_rep"))))))
  # 
  # delta_mat2<-delta_mat1[,3:dim(delta_mat1)[2]]/10
  # delta_mat<-cbind(delta_mat1[,1],rep(0,dim(delta_mat2)[1]),delta_mat2)
  # names(delta_mat[,1])<-c("t")
  # ##Means of each column by V2
  # 
  # ####If we wanted to analyze the averages of the different groups/how this tracks the real curves and what happens 
  # ####if we introduce other curves, we can first sort them by "magnitude", i.e. some kind of norm.
  # 
  # #for one group######
  # g1<-delta_mat[1:10,3]
  # pdf(file=here("Abortion_Draft",paste0(n,"g1.pdf")))
  # plot(g1,ylim=c(min(delta_mat2),max(delta_mat2)),type="l",ylab=expression(alpha[g]),xlab=expression(t),cex.lab=1.1,lwd=2)
  # dev.off()
  # ###for two groups
  # g2<-as.data.frame(delta_mat[1:10,4:5])
  # g2<-g2[, order(colSums(g2))]
  # pdf(file=here("Abortion_Draft",paste0(n,"g2.pdf")))
  # plot(g2[1:10,1],ylim=c(min(delta_mat2),max(delta_mat2)),type="l",ylab=expression(alpha[g]),xlab=expression(t),cex.lab=1.1,lwd=2)
  # lines(g2[1:10,2],lty=2,col="blue",lwd=2)
  # dev.off()
  # ###for three groups
  # g3<-as.data.frame(delta_mat[1:10,6:8])
  # g3<-g3[, order(colSums(g3))]####sorts these into the right risk order groups
  # pdf(file=here("Abortion_Draft",paste0(n,"g3.pdf")))
  # plot(g3[1:10,1],ylim=c(min(delta_mat2),max(delta_mat2)),type="l",ylab=expression(alpha[g]),xlab=expression(t),cex.lab=1.1,lwd=2)
  # lines(g3[1:10,2],lty=2,col="red",lwd=2)
  # lines(g3[1:10,3],lty=2,col="blue",lwd=2)
  # dev.off()
  # ###for four groups
  # 
  # g4<-as.data.frame(delta_mat[1:10,9:12])
  # g4<-g4[, order(colSums(g4))]
  # pdf(file=here("Abortion_Draft",paste0(n,"g4.pdf")))
  # plot(g4[1:10,1],ylim=c(min(delta_mat2),max(delta_mat2)),type="l",ylab=expression(alpha[g]),xlab=expression(t),cex.lab=1.1,lwd=2)
  # lines(g4[1:10,2],lty=2,col="green",lwd=2)
  # lines(g4[1:10,3],lty=2,col="red",lwd=2)
  # lines(g4[1:10,4],lty=2,col="blue",lwd=2)
  # dev.off()
  # ###for five groups
  # 
  # g5<-as.data.frame(delta_mat[1:10,13:17])
  # g5<-g5[, order(colSums(g5))]
  # pdf(file=here("Abortion_Draft",paste0(n,"g5.pdf")))
  # plot(g5[1:10,1],ylim=c(min(delta_mat2),max(delta_mat2)),type="l",ylab=expression(alpha[g]),xlab=expression(t),cex.lab=1.1,lwd=2)
  # lines(g5[1:10,2],lty=2,col="green",lwd=2)
  # lines(g5[1:10,3],lty=2,col="pink",lwd=2)
  # lines(g5[1:10,4],lty=2,col="red",lwd=2)
  # lines(g5[1:10,5],lty=2,col="blue",lwd=2)
  # dev.off()
  
}

main.ana(1000,0,ll=1)
main.ana(1500,0,ll=1)
main.ana(2000,0,ll=1)
main.ana(2500,0,ll=1)
main.ana(3000,0,ll=1)


