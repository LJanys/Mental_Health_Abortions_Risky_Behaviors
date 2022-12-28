# Install HonestDiD package from https://github.com/asheshrambachan/HonestDiD, if not installed 
 # install.packages('remotes') 
 # Sys.setenv('R_REMOTES_NO_ERRORS_FROM_WARNINGS' = 'true') 
 # remotes::install_github('asheshrambachan/HonestDiD') 
library(furrr)
library(purrr)
library(tidyverse)
library(stringr)
library(progressr)
library(mvtnorm)
library(Rglpk)
#############In this file: how to run the simulations######################
############and apply the simulation file. Several dynamic effects can be included. 
## 
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
library(HonestDiD) 
library(readxl) 

setwd(here("Results_Event_study"))

options(warn = -1) 

# Define a function to read in multiple sheets of an excel file simultaneously 
multiplesheets <- function(fname) { 
  sheets <- readxl::excel_sheets(fname) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
  names(data_frame) <- sheets 
  return(data_frame) 
} 

# Read in Excel Data 
# Complete Your Data Path if 'SenAnlysDataFormat (2) copy.xlsx' is not in Your Working Directory 
data <- multiplesheets('SenAnlysDataFormat (1).xlsx') 

# Reformat Data 
data <- c(list(sigma=as.matrix(data$sigma)), 
          lapply(data[c('beta', 'tVec', 'referencePeriod', 'prePeriodIndices', 'postPeriodIndices')], unlist)) 

beta_hat_full<-read.csv("coef_se_staggered_cs.csv")
beta_hat<-beta_hat_full$estimate
sigma_hat=read.csv("vcv_staggered_cs.csv")[,-1]

PrePeriods=6
PostPeriods=7

t_vec=seq(-6,6,le=13)#data$tVec
#length(data$postPeriodIndices),
#length(data$prePeriodIndices), 
# Create Event Study Plot 
event_study_plot <- createEventStudyPlot( 
  betahat = beta_hat, 
  sigma = as.matrix(sigma_hat), 
  numPrePeriods = PrePeriods ,
  numPostPeriods = PostPeriods,
  alpha =  0.05 , 
  timeVec = t_vec, 
  referencePeriod = -1,#data$referencePeriod, 
  useRelativeEventTime = TRUE 
) 
pdf(file="event_study_roth.pdf")
event_study_plot 
dev.off()


BC_l_vec = c(1,1,1,1,1,1,1)#BC_l_vec,
BC_DeltaSDRM_RobustResults = createSensitivityResults_relativeMagnitudes(betahat = beta_hat, 
                                                                        sigma = as.matrix(sigma_hat), 
                                                                         bound = "deviation from linear trend",
                                                                          l_vec =BC_l_vec,
                                                                         numPrePeriods = PrePeriods, 
                                                                         numPostPeriods = PostPeriods,
                                                                         Mbarvec = seq(from = 0.1, to = 0.5, by = 0.02),
                                                                         gridPoints = 100, grid.lb = -1, grid.ub = 1)



head(BC_DeltaSDRM_RobustResults)

BC_OriginalResults = constructOriginalCS(betahat = beta_hat, 
                                         sigma = as.matrix(sigma_hat), 
                                         l_vec =BC_l_vec,
                                         numPrePeriods = PrePeriods, 
                                         numPostPeriods = PostPeriods
                                         )

# Construct sensitivity plot.
robustResults<-BC_DeltaSDRM_RobustResults

originalResults<-BC_OriginalResults



createSensitivityPlot_relativeMagnitudes <- function(robustResults, originalResults, rescaleFactor = 1, maxMbar = Inf, add_xAxis = TRUE) {
  # Set Mbar for OLS to be the min Mbar in robust results minus the gap between Mbars in robust
  Mbargap <- min( diff( sort( robustResults$Mbar) ) )
  Mbarmin <- min( robustResults$Mbar)
  
  originalResults$Mbar <- Mbarmin - Mbargap
  df <- bind_rows(originalResults, robustResults)
  
  # Rescale all the units by rescaleFactor
  df <- df %>% mutate_at( c("Mbar", "ub", "lb"), ~ .x * rescaleFactor)
  
  # Filter out observations above maxM (after rescaling)
  df <- df %>% filter(Mbar <= maxMbar)
  
  p <- ggplot(data = df, aes(x=Mbar)) +
    geom_errorbar(aes(ymin = lb, ymax = ub, color = factor(method)),
                  width = Mbargap * rescaleFactor / 2) +
    scale_color_manual(values = c("red", '#01a2d9')) +
    theme(legend.title=element_blank(), legend.position="bottom") +
    labs(x = latex2exp::TeX("$\\bar{M}$"), y = "")
  
  if (add_xAxis) {
    p <- p + geom_hline(yintercept = 0)
  }
  return(p)
}

BC_DeltaSDRM_SensitivityPlot<-createSensitivityPlot_relativeMagnitudes(robustResults, originalResults, rescaleFactor = 1, maxMbar = Inf, add_xAxis = TRUE) 



pdf(file="sens_plot_mbar.pdf")                                                                     
BC_DeltaSDRM_SensitivityPlot
dev.off()


