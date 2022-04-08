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
 
}

main.ana(1000,0,ll=1)
main.ana(1500,0,ll=1)
main.ana(2000,0,ll=1)
main.ana(2500,0,ll=1)
main.ana(3000,0,ll=1)


