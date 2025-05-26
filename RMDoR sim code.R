###################
####           ####
#### Project 3 ####
####           ####
###################
# Calculation of the Average duration of response

######################
# Simulations part 1 #
######################

library(survival)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Small check 
shapeD = 3; scaleD = 6; shapeP = 0.715037; scaleP = 12; shapeR = 1.7; scaleR = 13
par = c(shapeD, scaleD, shapeP, scaleP, shapeR, scaleR, cL, cR)
Prog  <- rweibull(n,shape=par[3],scale=par[4])
Death <-  Prog + rweibull(n,shape=par[1],scale=par[2]) #rep(10000,n)
Resp  <- rweibull(n,shape=par[5],scale=par[6])
sum(( Resp < Prog) & ( Prog <= Resp + 1 ))
sum(( Resp < Prog))
sum(Resp - Prog)







# parameters for the Weibull distributions for time to progression, 
# response, time between progression and death, time to censoring
# Computes the DoR in case the progression and response times are set equal
# to the right and to the mid end-points of the intervals.

# time to progression
shapeP = 2
scaleP = 25
# time to response
shapeR = 2
scaleR = 25  
# time between progression and death
shapeD = 3       
scaleD = 6
# time to censoring
cL = 5                 # censoring times are simulated from unif[cL,CR]
cR = 70

par = c(shapeD, scaleD, shapeP, scaleP, shapeR, scaleR, cL, cR)
# note: if you change one of the parameters, you need to update   par

# Scenario's, note: I included 0
S0 = seq(0,100,by=0.01)     # as a check (DoR should be equal to DoRInt below)
S1 = c(0, 2, 4,  6,  8, 10, 12, 18, 24, 30, 36, 52, 72, 1000)
S2 = c(0, 2, 4,  8, 10, 12, 18, 24, 30, 36, 52, 72, 1000)
S3 = c(0, 2, 4,  8, 12, 18, 24, 30, 36, 52, 72, 1000)
S4 = c(0, 2, 4,  8, 12, 18, 30, 36, 52, 72, 1000)
S5 = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72)
S6 = c(0, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 21, 27, 33, 39, 45, 51, 57, 63, 69)


# remaining assumptions
n = 500          # sample size
M = 1000 #1000      # number of iterations for simulation studies
Scenario <- S1   # choice of the scenario
tau = 40        # area between curves is computed at interval [0,tau]

# This function simulates observations for n patients.
SimulateData <- function(n,Scenario,par)
{ #First data for progression, death and response are simulated
  Prog  <- rweibull(n,shape=par[3],scale=par[4])
  Death <-  Prog + rweibull(n,shape=par[1],scale=par[2]) #rep(10000,n)
  Resp  <- rweibull(n,shape=par[5],scale=par[6])
  #Cens  <- rep(10000,n) 
  Cens  <- runif(n,par[7],par[8]) 
  
  # No interval censoring
  # Observations are calculated from the simulated data
  # PFS
  FUpfs <- apply(cbind(Prog,Death),FUN=min,MARGIN=1)
  statuspfs <- as.numeric(FUpfs < Cens)
  Tpfs <- statuspfs*FUpfs + (1-statuspfs)*Cens 
  
  # RPFS
  FUrpfs <- apply(cbind(Prog,Resp,Death),FUN=min,MARGIN=1)
  statusrpfs <- as.numeric(FUrpfs < Cens) 
  Trpfs <- statusrpfs*FUrpfs + (1-statusrpfs)*Cens 
  
  
  # interval censoring, upper  
  # progression and response dates are calculated based on chosen scenario
  ProgInt <- RespInt <- ProgInt_r_check <- ProgInt_m_check <- ProgInt_ORR <- rep(0,n) 
  for (i in 1:n){
    ProgInt[i] <- min(Scenario[Prog[i] <= Scenario])  
    RespInt[i] <- min(Scenario[Resp[i] <= Scenario])
    
    ProgInt_r_check[i] <- ProgInt[i]
    ProgInt_m_check[i] <- ProgInt[i]
    ProgInt_ORR[i] <- ProgInt[i]
    
    if((RespInt[i] < ProgInt[i]) & ((RespInt[i] + 1) > Prog[i]) == T)
    {
      # Right censored check RR
      ProgInt_r_check[i] <- RespInt[i] + 1
      
      # Mid point censored check RR
      ProgInt_m_check[i] <- RespInt[i] + 0.5
      
      # ORR progression check
      ProgInt_ORR[i] <- RespInt[i]
    }
  }
  
  # Observations are calculated
  FUpfsInt <- apply(cbind(ProgInt,Death),FUN=min,MARGIN=1)
  statuspfsInt <- as.numeric(FUpfsInt < Cens)
  TpfsInt <- statuspfsInt*FUpfsInt + (1-statuspfsInt)*Cens 
  
  FUrpfsInt <- apply(cbind(ProgInt,RespInt,Death),FUN=min,MARGIN=1)
  statusrpfsInt <- as.numeric(FUrpfsInt < Cens) 
  TrpfsInt <- statusrpfsInt*FUrpfsInt + (1-statusrpfsInt)*Cens 
  
  
  # interval censoring, middle point of interval  
  # progression and response dates are calculated based on chosen scenario
  TpfsIntM <- TpfsInt 
  TrpfsIntM <- TrpfsInt
  s1 <- s2 <- rep(0,n) 
  for (i in 1:n){
    s1[i] <- max(Scenario[Prog[i] >= Scenario])    # compute left end point
    s2[i] <- max(Scenario[Resp[i] >= Scenario])
  }
  
  # Observations are calculated
  for (i in 1:n) 
  { if (statuspfsInt[i] == 1)
    if (ProgInt[i] < Death[i])
      TpfsIntM[i] <- (s1[i]+TpfsInt[i])/2
  }  
  
  for (i in 1:n) 
  { if (statusrpfsInt[i] == 1)
    if (min(ProgInt[i],RespInt[i]) < Death[i])
      TrpfsIntM[i] <- (min(s1[i],s2[i])+TrpfsInt[i])/2
  }  
  
  
  # With 4 weeks check
  ###################################################
  # interval censoring (right) Prog in the right point
  FUpfsInt_r_conf <- apply(cbind(ProgInt_r_check,Death),FUN=min,MARGIN=1)
  statuspfsInt_r_conf <- as.numeric(FUpfsInt_r_conf < Cens)
  TpfsInt_r_conf <- statuspfsInt_r_conf*FUpfsInt_r_conf + (1-statuspfsInt_r_conf)*Cens 
  
  FUrpfsInt_r_conf <- apply(cbind(ProgInt_r_check,RespInt,Death),FUN=min,MARGIN=1)
  statusrpfsInt_r_conf <- as.numeric(FUrpfsInt_r_conf < Cens) 
  TrpfsInt_r_conf <- statusrpfsInt_r_conf*FUrpfsInt_r_conf + (1-statusrpfsInt_r_conf)*Cens 
  
  
  # With 4 weeks check
  ###################################################
  # interval censoring (mid) Prog in the mid point
  TpfsIntM_conf <- TpfsInt_r_conf 
  TrpfsIntM_conf <- TrpfsInt_r_conf
  s1 <- s2 <- rep(0,n) 
  for (i in 1:n){
    s1[i] <- max(Scenario[Prog[i] >= Scenario])    # compute left end point
    s2[i] <- max(Scenario[Resp[i] >= Scenario])
  }
  
  # Observations are calculated
  for (i in 1:n) 
  { if (statuspfsInt_r_conf[i] == 1)
    if (ProgInt_r_check[i] < Death[i])
      TpfsIntM_conf[i] <- (s1[i]+TpfsInt_r_conf[i])/2
  }  
  
  for (i in 1:n) 
  { if (statusrpfsInt_r_conf[i] == 1)
    if (min(ProgInt_r_check[i],RespInt[i]) < Death[i])
      TrpfsIntM_conf[i] <- (min(s1[i],s2[i])+TrpfsInt_r_conf[i])/2
  }  
  
  
    # With 4 weeks check ORR 
  ###################################################
  # ORR with right censoring
  FUpfsInt_ORR <- apply(cbind(ProgInt_ORR,Death),FUN=min,MARGIN=1)
  statuspfsInt_ORR <- as.numeric(FUpfsInt_ORR < Cens)
  TpfsInt_ORR <- statuspfsInt_ORR*FUpfsInt_ORR + (1-statuspfsInt_ORR)*Cens 
  
  FUrpfsInt_ORR <- apply(cbind(ProgInt_ORR,RespInt,Death),FUN=min,MARGIN=1)
  statusrpfsInt_ORR <- as.numeric(FUrpfsInt_ORR < Cens) 
  TrpfsInt_ORR <- statusrpfsInt_ORR*FUrpfsInt_ORR + (1-statusrpfsInt_ORR)*Cens 
  
  # With 4 weeks check ORR 
  ###################################################
  # ORR with mid point
  TpfsIntM_ORR <- TpfsInt_ORR 
  TrpfsIntM_ORR <- TrpfsInt_ORR
  s1 <- s2 <- rep(0,n) 
  for (i in 1:n)
  {
    s1[i] <- max(Scenario[Prog[i] >= Scenario])    # compute left end point
    s2[i] <- max(Scenario[Resp[i] >= Scenario])
  }
  
  # Observations are calculated
  for (i in 1:n) 
  {
    if (statuspfsInt_ORR[i] == 1)
    if (ProgInt_ORR[i] < Death[i])
      TpfsIntM_ORR[i] <- (s1[i]+TpfsInt_ORR[i])/2
  }  
  
  for (i in 1:n) 
  { 
    if (statusrpfsInt_ORR[i] == 1)
    if (min(ProgInt_ORR[i],RespInt[i]) < Death[i])
      TrpfsIntM_ORR[i] <- (min(s1[i],s2[i])+TrpfsInt_ORR[i])/2
  }  
   

  
  list(timePFS=Tpfs,StatusPFS=statuspfs,timeRPFS=Trpfs,StatusRPFS=statusrpfs,
       timePFSInt_no_conf=TpfsInt,StatusPFSInt_no_conf=statuspfsInt,timeRPFSInt_no_conf=TrpfsInt,StatusRPFSInt_no_conf=statusrpfsInt,
       Prog=Prog,Resp=Resp,
       timePFSIntM_no_conf=TpfsIntM,StatusPFSIntM_no_conf=statuspfsInt,timeRPFSIntM_no_conf=TrpfsIntM,StatusRPFSIntM_no_conf=statusrpfsInt,
       timePFSInt_conf=TpfsInt_r_conf,StatusPFSInt_conf=statuspfsInt_r_conf,timeRPFSInt_conf=TrpfsInt_r_conf,StatusRPFSInt_conf=statusrpfsInt_r_conf,
       timePFSIntM_conf=TpfsIntM_conf,StatusPFSIntM_conf=statuspfsInt_r_conf,timeRPFSIntM_conf=TrpfsIntM_conf,StatusRPFSIntM_conf=statusrpfsInt_r_conf,
       timePFSInt_ORR=TpfsInt_ORR,StatusPFSInt_ORR=statuspfsInt_ORR,timeRPFSInt_ORR=TrpfsInt_ORR,StatusRPFSInt_ORR=statusrpfsInt_ORR,
       timePFSIntM_ORR=TpfsIntM_ORR,StatusPFSIntM_ORR=statuspfsInt_ORR,timeRPFSIntM_ORR=TrpfsIntM_ORR,StatusRPFSIntM_ORR=statusrpfsInt_ORR)
}

SimulateData2 <- function(n,Scenario,par)
{ #First data for progression, death and response are simulated
  Prog  <- rweibull(n,shape=par[3],scale=par[4])
  Death <-  Prog + rweibull(n,shape=par[1],scale=par[2]) #rep(10000,n)
  Resp  <- rweibull(n,shape=par[5],scale=par[6])
  Cens  <- rep(10000,n) 
  #Cens  <- runif(n,par[7],par[8]) 
  
  # No interval censoring
  # Observations are calculated from the simulated data
  # PFS
  FUpfs <- apply(cbind(Prog,Death),FUN=min,MARGIN=1)
  statuspfs <- as.numeric(FUpfs < Cens)
  Tpfs <- statuspfs*FUpfs + (1-statuspfs)*Cens 
  
  # RPFS
  FUrpfs <- apply(cbind(Prog,Resp,Death),FUN=min,MARGIN=1)
  statusrpfs <- as.numeric(FUrpfs < Cens) 
  Trpfs <- statusrpfs*FUrpfs + (1-statusrpfs)*Cens 
  
  
  # interval censoring, upper  
  # progression and response dates are calculated based on chosen scenario
  ProgInt <- RespInt <- ProgInt_r_check <- ProgInt_m_check <- ProgInt_ORR <- rep(0,n) 
  for (i in 1:n){
    ProgInt[i] <- min(Scenario[Prog[i] <= Scenario])  
    RespInt[i] <- min(Scenario[Resp[i] <= Scenario])
    
    ProgInt_r_check[i] <- ProgInt[i]
    ProgInt_m_check[i] <- ProgInt[i]
    ProgInt_ORR[i] <- ProgInt[i]
    
    if((RespInt[i] < ProgInt[i]) & ((RespInt[i] + 1) > Prog[i]) == T)
    {
      # Right censored check RR
      ProgInt_r_check[i] <- RespInt[i] + 1
      
      # Mid point censored check RR
      ProgInt_m_check[i] <- RespInt[i] + 0.5
      
      # ORR progression check
      ProgInt_ORR[i] <- RespInt[i]
    }
  }
  
  # Observations are calculated
  FUpfsInt <- apply(cbind(ProgInt,Death),FUN=min,MARGIN=1)
  statuspfsInt <- as.numeric(FUpfsInt < Cens)
  TpfsInt <- statuspfsInt*FUpfsInt + (1-statuspfsInt)*Cens 
  
  FUrpfsInt <- apply(cbind(ProgInt,RespInt,Death),FUN=min,MARGIN=1)
  statusrpfsInt <- as.numeric(FUrpfsInt < Cens) 
  TrpfsInt <- statusrpfsInt*FUrpfsInt + (1-statusrpfsInt)*Cens 
  
  
  # interval censoring, middle point of interval  
  # progression and response dates are calculated based on chosen scenario
  TpfsIntM <- TpfsInt 
  TrpfsIntM <- TrpfsInt
  s1 <- s2 <- rep(0,n) 
  for (i in 1:n){
    s1[i] <- max(Scenario[Prog[i] >= Scenario])    # compute left end point
    s2[i] <- max(Scenario[Resp[i] >= Scenario])
  }
  
  # Observations are calculated
  for (i in 1:n) 
  { if (statuspfsInt[i] == 1)
    if (ProgInt[i] < Death[i])
      TpfsIntM[i] <- (s1[i]+TpfsInt[i])/2
  }  
  
  for (i in 1:n) 
  { if (statusrpfsInt[i] == 1)
    if (min(ProgInt[i],RespInt[i]) < Death[i])
      TrpfsIntM[i] <- (min(s1[i],s2[i])+TrpfsInt[i])/2
  }  
  
  
  # With 4 weeks check
  ###################################################
  # interval censoring (right) Prog in the right point
  FUpfsInt_r_conf <- apply(cbind(ProgInt_r_check,Death),FUN=min,MARGIN=1)
  statuspfsInt_r_conf <- as.numeric(FUpfsInt_r_conf < Cens)
  TpfsInt_r_conf <- statuspfsInt_r_conf*FUpfsInt_r_conf + (1-statuspfsInt_r_conf)*Cens 
  
  FUrpfsInt_r_conf <- apply(cbind(ProgInt_r_check,RespInt,Death),FUN=min,MARGIN=1)
  statusrpfsInt_r_conf <- as.numeric(FUrpfsInt_r_conf < Cens) 
  TrpfsInt_r_conf <- statusrpfsInt_r_conf*FUrpfsInt_r_conf + (1-statusrpfsInt_r_conf)*Cens 
  
  
  # With 4 weeks check
  ###################################################
  # interval censoring (mid) Prog in the mid point
  TpfsIntM_conf <- TpfsInt_r_conf 
  TrpfsIntM_conf <- TrpfsInt_r_conf
  s1 <- s2 <- rep(0,n) 
  for (i in 1:n){
    s1[i] <- max(Scenario[Prog[i] >= Scenario])    # compute left end point
    s2[i] <- max(Scenario[Resp[i] >= Scenario])
  }
  
  # Observations are calculated
  for (i in 1:n) 
  { if (statuspfsInt_r_conf[i] == 1)
    if (ProgInt_r_check[i] < Death[i])
      TpfsIntM_conf[i] <- (s1[i]+TpfsInt_r_conf[i])/2
  }  
  
  for (i in 1:n) 
  { if (statusrpfsInt_r_conf[i] == 1)
    if (min(ProgInt_r_check[i],RespInt[i]) < Death[i])
      TrpfsIntM_conf[i] <- (min(s1[i],s2[i])+TrpfsInt_r_conf[i])/2
  }  
  
  
  # With 4 weeks check ORR 
  ###################################################
  # ORR with right censoring
  FUpfsInt_ORR <- apply(cbind(ProgInt_ORR,Death),FUN=min,MARGIN=1)
  statuspfsInt_ORR <- as.numeric(FUpfsInt_ORR < Cens)
  TpfsInt_ORR <- statuspfsInt_ORR*FUpfsInt_ORR + (1-statuspfsInt_ORR)*Cens 
  
  FUrpfsInt_ORR <- apply(cbind(ProgInt_ORR,RespInt,Death),FUN=min,MARGIN=1)
  statusrpfsInt_ORR <- as.numeric(FUrpfsInt_ORR < Cens) 
  TrpfsInt_ORR <- statusrpfsInt_ORR*FUrpfsInt_ORR + (1-statusrpfsInt_ORR)*Cens 
  
  # With 4 weeks check ORR 
  ###################################################
  # ORR with mid point
  TpfsIntM_ORR <- TpfsInt_ORR 
  TrpfsIntM_ORR <- TrpfsInt_ORR
  s1 <- s2 <- rep(0,n) 
  for (i in 1:n)
  {
    s1[i] <- max(Scenario[Prog[i] >= Scenario])    # compute left end point
    s2[i] <- max(Scenario[Resp[i] >= Scenario])
  }
  
  # Observations are calculated
  for (i in 1:n) 
  {
    if (statuspfsInt_ORR[i] == 1)
      if (ProgInt_ORR[i] < Death[i])
        TpfsIntM_ORR[i] <- (s1[i]+TpfsInt_ORR[i])/2
  }  
  
  for (i in 1:n) 
  { 
    if (statusrpfsInt_ORR[i] == 1)
      if (min(ProgInt_ORR[i],RespInt[i]) < Death[i])
        TrpfsIntM_ORR[i] <- (min(s1[i],s2[i])+TrpfsInt_ORR[i])/2
  }  
  
  
  
  list(timePFS=Tpfs,StatusPFS=statuspfs,timeRPFS=Trpfs,StatusRPFS=statusrpfs,
       timePFSInt_no_conf=TpfsInt,StatusPFSInt_no_conf=statuspfsInt,timeRPFSInt_no_conf=TrpfsInt,StatusRPFSInt_no_conf=statusrpfsInt,
       Prog=Prog,Resp=Resp,
       timePFSIntM_no_conf=TpfsIntM,StatusPFSIntM_no_conf=statuspfsInt,timeRPFSIntM_no_conf=TrpfsIntM,StatusRPFSIntM_no_conf=statusrpfsInt,
       timePFSInt_conf=TpfsInt_r_conf,StatusPFSInt_conf=statuspfsInt_r_conf,timeRPFSInt_conf=TrpfsInt_r_conf,StatusRPFSInt_conf=statusrpfsInt_r_conf,
       timePFSIntM_conf=TpfsIntM_conf,StatusPFSIntM_conf=statuspfsInt_r_conf,timeRPFSIntM_conf=TrpfsIntM_conf,StatusRPFSIntM_conf=statusrpfsInt_r_conf,
       timePFSInt_ORR=TpfsInt_ORR,StatusPFSInt_ORR=statuspfsInt_ORR,timeRPFSInt_ORR=TrpfsInt_ORR,StatusRPFSInt_ORR=statusrpfsInt_ORR,
       timePFSIntM_ORR=TpfsIntM_ORR,StatusPFSIntM_ORR=statuspfsInt_ORR,timeRPFSIntM_ORR=TrpfsIntM_ORR,StatusRPFSIntM_ORR=statusrpfsInt_ORR)
}


# this function computes the area under the KM curve at the interval [0,tau]
ComputeAreaUnderCurve <- function(time,status,tau)
{ KMout <- summary(survfit(Surv(time, status)~1))
tjumps <- c(0,KMout$time,tau)
tjumps[tjumps>tau] <- tau
Svalue <- c(1,KMout$surv)
sum((tjumps[2:length(tjumps)]-tjumps[1:(length(tjumps)-1)])*Svalue)
}


# Main program
sim_func_mDoR_calc <- function(shapeD, scaleD, shapeP, scaleP, shapeR, scaleR, cL, cR, tau, n, iter = 1000, scan_scenario, CI_low = 0.05, CI_upper = 0.95)
{
  par = c(shapeD, scaleD, shapeP, scaleP, shapeR, scaleR, cL, cR)
  M <- iter
  
  DoRInt_no_conf <- DoRIntM_no_conf <- DoRInt_conf <- DoRIntM_conf <- DoRInt_ORR <- DoRIntM_ORR <- DoR <- rep(0,M)
  AUCPFSIntM_no_conf <- AUCRPFSIntM_no_conf <- AUCPFSInt_no_conf <- AUCRPFSInt_no_conf <- AUCPFS <- AUCRPFS <- rep(0,M)  
  AUCPFSInt_conf <- AUCRPFSInt_conf <- AUCPFSIntM_conf <- AUCRPFSIntM_conf <- AUCPFSInt_ORR <- AUCRPFSInt_ORR <- AUCPFSIntM_ORR <- AUCRPFSIntM_ORR <- rep(0,M)  
  atrisk_daily          <- NULL
  atrisk_interval       <- NULL
  atrisk_interval_check <- NULL
  atrisk_interval_ORR   <- NULL
  atrisk_midpoint       <- NULL
  atrisk_midpoint_check <- NULL
  atrisk_midpoint_ORR   <- NULL
  
  # AUC stands for area under curve
  
  for (i in 1:M){
    data <- SimulateData(n, scan_scenario, par)
    
    # Without 4 weeks check
      # interval censoring (right)
      AUCPFSInt_no_conf[i]  <- ComputeAreaUnderCurve(data$timePFSInt_no_conf,data$StatusPFSInt_no_conf,tau)
      AUCRPFSInt_no_conf[i] <- ComputeAreaUnderCurve(data$timeRPFSInt_no_conf,data$StatusRPFSInt_no_conf,tau) 
      DoRInt_no_conf[i]     <- AUCPFSInt_no_conf[i] - AUCRPFSInt_no_conf[i]     # area between curves
      
      # interval censoring Midpoint
      AUCPFSIntM_no_conf[i]  <- ComputeAreaUnderCurve(data$timePFSIntM_no_conf,data$StatusPFSIntM_no_conf,tau)
      AUCRPFSIntM_no_conf[i] <- ComputeAreaUnderCurve(data$timeRPFSIntM_no_conf,data$StatusRPFSIntM_no_conf,tau) 
      DoRIntM_no_conf[i]     <- AUCPFSIntM_no_conf[i] - AUCRPFSIntM_no_conf[i]     # area between curves
    
    # With 4 weeks check
    ###################################################
      # interval censoring (right) Prog in the right point
      AUCPFSInt_conf[i]  <- ComputeAreaUnderCurve(data$timePFSInt_conf,data$StatusPFSInt_conf,tau)
      AUCRPFSInt_conf[i] <- ComputeAreaUnderCurve(data$timeRPFSInt_conf,data$StatusRPFSInt_conf,tau) 
      DoRInt_conf[i]     <- AUCPFSInt_conf[i] - AUCRPFSInt_conf[i]     # area between curves
      
      # interval censoring Midpoint Prog in the mid point
      AUCPFSIntM_conf[i]  <- ComputeAreaUnderCurve(data$timePFSIntM_conf,data$StatusPFSIntM_conf,tau)
      AUCRPFSIntM_conf[i] <- ComputeAreaUnderCurve(data$timeRPFSIntM_conf,data$StatusRPFSIntM_conf,tau) 
      DoRIntM_conf[i]     <- AUCPFSIntM_conf[i] - AUCRPFSIntM_conf[i]     # area between curves
    ###################################################
      # With 4 weeks check ORR
      ###################################################
      # interval censoring (right) Prog ORR
      AUCPFSInt_ORR[i]  <- ComputeAreaUnderCurve(data$timePFSInt_ORR,data$StatusPFSInt_ORR,tau)
      AUCRPFSInt_ORR[i] <- ComputeAreaUnderCurve(data$timeRPFSInt_ORR,data$StatusRPFSInt_ORR,tau) 
      DoRInt_ORR[i]     <- AUCPFSInt_ORR[i] - AUCRPFSInt_ORR[i]     # area between curves
      
      # interval censoring Midpoint Prog ORR
      AUCPFSIntM_ORR[i]  <- ComputeAreaUnderCurve(data$timePFSIntM_ORR,data$StatusPFSIntM_ORR,tau)
      AUCRPFSIntM_ORR[i] <- ComputeAreaUnderCurve(data$timeRPFSIntM_ORR,data$StatusRPFSIntM_ORR,tau) 
      DoRIntM_ORR[i]     <- AUCPFSIntM_ORR[i] - AUCRPFSIntM_ORR[i]     # area between curves
    ####################################################
      
    # no interval-censoring, events exactly observed
    AUCPFS[i]  <- ComputeAreaUnderCurve(data$timePFS,data$StatusPFS,tau)
    AUCRPFS[i] <- ComputeAreaUnderCurve(data$timeRPFS,data$StatusRPFS,tau) 
    DoR[i]     <- AUCPFS[i]-AUCRPFS[i]                # area between curves
    
    # Count the number of pts at risk in each time tau
    atrisk_daily[i]          <- sum(data$timePFS >= tau)
    atrisk_interval[i]       <- sum(data$timePFSInt_no_conf >= tau)
    atrisk_interval_check[i] <- sum(data$timePFSInt_conf >= tau)
    atrisk_interval_ORR[i]   <- sum(data$timePFSInt_ORR >= tau)
    atrisk_midpoint[i]       <- sum(data$timePFSIntM_no_conf >= tau)
    atrisk_midpoint_check[i] <- sum(data$timePFSIntM_conf >= tau)
    atrisk_midpoint_ORR[i]   <- sum(data$timePFSIntM_ORR >= tau)
    
  }
  
  return(list(Average_DoR = mean(DoR),
              DoR_CI_low  = quantile(DoR, probs = c(CI_low, CI_upper))[1],
              DoR_CI_upper= quantile(DoR, probs = c(CI_low, CI_upper))[2],
              Average_DoR_Int_no_conf = mean(DoRInt_no_conf), 
              DoR_Int_no_conf_CI_low  = quantile(DoRInt_no_conf, probs = c(CI_low, CI_upper))[1],
              DoR_Int_no_conf_CI_upper= quantile(DoRInt_no_conf, probs = c(CI_low, CI_upper))[2],
              Average_DoR_IntMid_no_conf = mean(DoRIntM_no_conf),
              DoR_IntMid_no_conf_CI_low  = quantile(DoRIntM_no_conf, probs = c(CI_low, CI_upper))[1],
              DoR_IntMid_no_conf_CI_upper= quantile(DoRIntM_no_conf, probs = c(CI_low, CI_upper))[2],
              Average_DoR_Int_conf_right = mean(DoRInt_conf),
              DoR_Int_conf_right_CI_low  = quantile(DoRInt_conf, probs = c(CI_low, CI_upper))[1],
              DoR_Int_conf_right_CI_upper= quantile(DoRInt_conf, probs = c(CI_low, CI_upper))[2],
              Average_DoR_IntMid_conf = mean(DoRIntM_conf),
              DoR_IntMid_conf_CI_low  = quantile(DoRIntM_conf, probs = c(CI_low, CI_upper))[1],
              DoR_IntMid_conf_CI_upper= quantile(DoRIntM_conf, probs = c(CI_low, CI_upper))[2],
              Average_DoR_Int_ORR = mean(DoRInt_ORR), 
              DoR_Int_ORR_CI_low  = quantile(DoRInt_ORR, probs = c(CI_low, CI_upper))[1],
              DoR_Int_ORR_CI_upper= quantile(DoRInt_ORR, probs = c(CI_low, CI_upper))[2],
              Average_DoR_IntMid_ORR = mean(DoRIntM_ORR),
              DoR_IntMid_ORR_CI_low  = quantile(DoRIntM_ORR, probs = c(CI_low, CI_upper))[1],
              DoR_IntMid_ORR_CI_upper= quantile(DoRIntM_ORR, probs = c(CI_low, CI_upper))[2],
              atrisk_daily = mean(atrisk_daily),
              atrisk_interval = mean(atrisk_interval),
              atrisk_interval_check = mean(atrisk_interval_check),
              atrisk_interval_ORR = mean(atrisk_interval_ORR),
              atrisk_midpoint = mean(atrisk_midpoint),
              atrisk_midpoint_check = mean(atrisk_midpoint_check),
              atrisk_midpoint_ORR = mean(atrisk_midpoint_ORR),
              DoR = DoR,
              DoRInt_no_conf = DoRInt_no_conf,
              DoRIntM_no_conf = DoRIntM_no_conf,
              DoRInt_conf = DoRInt_conf,
              DoRIntM_conf = DoRIntM_conf,
              DoRInt_ORR = DoRInt_ORR,
              DoRIntM_ORR = DoRIntM_ORR))
}

# Main program
sim_func_mDoR_calc2 <- function(shapeD, scaleD, shapeP, scaleP, shapeR, scaleR, cL, cR, tau, n, iter = 1000, scan_scenario, CI_low = 0.05, CI_upper = 0.95)
{
  par = c(shapeD, scaleD, shapeP, scaleP, shapeR, scaleR, cL, cR)
  M <- iter
  
  DoRInt_no_conf <- DoRIntM_no_conf <- DoRInt_conf <- DoRIntM_conf <- DoRInt_ORR <- DoRIntM_ORR <- DoR <- rep(0,M)
  AUCPFSIntM_no_conf <- AUCRPFSIntM_no_conf <- AUCPFSInt_no_conf <- AUCRPFSInt_no_conf <- AUCPFS <- AUCRPFS <- rep(0,M)  
  AUCPFSInt_conf <- AUCRPFSInt_conf <- AUCPFSIntM_conf <- AUCRPFSIntM_conf <- AUCPFSInt_ORR <- AUCRPFSInt_ORR <- AUCPFSIntM_ORR <- AUCRPFSIntM_ORR <- rep(0,M)  
  atrisk_daily          <- NULL
  atrisk_interval       <- NULL
  atrisk_interval_check <- NULL
  atrisk_interval_ORR   <- NULL
  atrisk_midpoint       <- NULL
  atrisk_midpoint_check <- NULL
  atrisk_midpoint_ORR   <- NULL
  
  # AUC stands for area under curve
  
  for (i in 1:M){
    data <- SimulateData2(n, scan_scenario, par)
    
    # Without 4 weeks check
    # interval censoring (right)
    AUCPFSInt_no_conf[i]  <- ComputeAreaUnderCurve(data$timePFSInt_no_conf,data$StatusPFSInt_no_conf,tau)
    AUCRPFSInt_no_conf[i] <- ComputeAreaUnderCurve(data$timeRPFSInt_no_conf,data$StatusRPFSInt_no_conf,tau) 
    DoRInt_no_conf[i]     <- AUCPFSInt_no_conf[i] - AUCRPFSInt_no_conf[i]     # area between curves
    
    # interval censoring Midpoint
    AUCPFSIntM_no_conf[i]  <- ComputeAreaUnderCurve(data$timePFSIntM_no_conf,data$StatusPFSIntM_no_conf,tau)
    AUCRPFSIntM_no_conf[i] <- ComputeAreaUnderCurve(data$timeRPFSIntM_no_conf,data$StatusRPFSIntM_no_conf,tau) 
    DoRIntM_no_conf[i]     <- AUCPFSIntM_no_conf[i] - AUCRPFSIntM_no_conf[i]     # area between curves
    
    # With 4 weeks check
    ###################################################
    # interval censoring (right) Prog in the right point
    AUCPFSInt_conf[i]  <- ComputeAreaUnderCurve(data$timePFSInt_conf,data$StatusPFSInt_conf,tau)
    AUCRPFSInt_conf[i] <- ComputeAreaUnderCurve(data$timeRPFSInt_conf,data$StatusRPFSInt_conf,tau) 
    DoRInt_conf[i]     <- AUCPFSInt_conf[i] - AUCRPFSInt_conf[i]     # area between curves
    
    # interval censoring Midpoint Prog in the mid point
    AUCPFSIntM_conf[i]  <- ComputeAreaUnderCurve(data$timePFSIntM_conf,data$StatusPFSIntM_conf,tau)
    AUCRPFSIntM_conf[i] <- ComputeAreaUnderCurve(data$timeRPFSIntM_conf,data$StatusRPFSIntM_conf,tau) 
    DoRIntM_conf[i]     <- AUCPFSIntM_conf[i] - AUCRPFSIntM_conf[i]     # area between curves
    ###################################################
    # With 4 weeks check ORR
    ###################################################
    # interval censoring (right) Prog ORR
    AUCPFSInt_ORR[i]  <- ComputeAreaUnderCurve(data$timePFSInt_ORR,data$StatusPFSInt_ORR,tau)
    AUCRPFSInt_ORR[i] <- ComputeAreaUnderCurve(data$timeRPFSInt_ORR,data$StatusRPFSInt_ORR,tau) 
    DoRInt_ORR[i]     <- AUCPFSInt_ORR[i] - AUCRPFSInt_ORR[i]     # area between curves
    
    # interval censoring Midpoint Prog ORR
    AUCPFSIntM_ORR[i]  <- ComputeAreaUnderCurve(data$timePFSIntM_ORR,data$StatusPFSIntM_ORR,tau)
    AUCRPFSIntM_ORR[i] <- ComputeAreaUnderCurve(data$timeRPFSIntM_ORR,data$StatusRPFSIntM_ORR,tau) 
    DoRIntM_ORR[i]     <- AUCPFSIntM_ORR[i] - AUCRPFSIntM_ORR[i]     # area between curves
    ####################################################
    
    # no interval-censoring, events exactly observed
    AUCPFS[i]  <- ComputeAreaUnderCurve(data$timePFS,data$StatusPFS,tau)
    AUCRPFS[i] <- ComputeAreaUnderCurve(data$timeRPFS,data$StatusRPFS,tau) 
    DoR[i]     <- AUCPFS[i]-AUCRPFS[i]                # area between curves
    
    # Count the number of pts at risk in each time tau
    atrisk_daily[i]          <- sum(data$timePFS >= tau)
    atrisk_interval[i]       <- sum(data$timePFSInt_no_conf >= tau)
    atrisk_interval_check[i] <- sum(data$timePFSInt_conf >= tau)
    atrisk_interval_ORR[i]   <- sum(data$timePFSInt_ORR >= tau)
    atrisk_midpoint[i]       <- sum(data$timePFSIntM_no_conf >= tau)
    atrisk_midpoint_check[i] <- sum(data$timePFSIntM_conf >= tau)
    atrisk_midpoint_ORR[i]   <- sum(data$timePFSIntM_ORR >= tau)
    
  }
  
  return(list(Average_DoR = mean(DoR),
              DoR_CI_low  = quantile(DoR, probs = c(CI_low, CI_upper))[1],
              DoR_CI_upper= quantile(DoR, probs = c(CI_low, CI_upper))[2],
              Average_DoR_Int_no_conf = mean(DoRInt_no_conf), 
              DoR_Int_no_conf_CI_low  = quantile(DoRInt_no_conf, probs = c(CI_low, CI_upper))[1],
              DoR_Int_no_conf_CI_upper= quantile(DoRInt_no_conf, probs = c(CI_low, CI_upper))[2],
              Average_DoR_IntMid_no_conf = mean(DoRIntM_no_conf),
              DoR_IntMid_no_conf_CI_low  = quantile(DoRIntM_no_conf, probs = c(CI_low, CI_upper))[1],
              DoR_IntMid_no_conf_CI_upper= quantile(DoRIntM_no_conf, probs = c(CI_low, CI_upper))[2],
              Average_DoR_Int_conf_right = mean(DoRInt_conf),
              DoR_Int_conf_right_CI_low  = quantile(DoRInt_conf, probs = c(CI_low, CI_upper))[1],
              DoR_Int_conf_right_CI_upper= quantile(DoRInt_conf, probs = c(CI_low, CI_upper))[2],
              Average_DoR_IntMid_conf = mean(DoRIntM_conf),
              DoR_IntMid_conf_CI_low  = quantile(DoRIntM_conf, probs = c(CI_low, CI_upper))[1],
              DoR_IntMid_conf_CI_upper= quantile(DoRIntM_conf, probs = c(CI_low, CI_upper))[2],
              Average_DoR_Int_ORR = mean(DoRInt_ORR), 
              DoR_Int_ORR_CI_low  = quantile(DoRInt_ORR, probs = c(CI_low, CI_upper))[1],
              DoR_Int_ORR_CI_upper= quantile(DoRInt_ORR, probs = c(CI_low, CI_upper))[2],
              Average_DoR_IntMid_ORR = mean(DoRIntM_ORR),
              DoR_IntMid_ORR_CI_low  = quantile(DoRIntM_ORR, probs = c(CI_low, CI_upper))[1],
              DoR_IntMid_ORR_CI_upper= quantile(DoRIntM_ORR, probs = c(CI_low, CI_upper))[2],
              atrisk_daily = mean(atrisk_daily),
              atrisk_interval = mean(atrisk_interval),
              atrisk_interval_check = mean(atrisk_interval_check),
              atrisk_interval_ORR = mean(atrisk_interval_ORR),
              atrisk_midpoint = mean(atrisk_midpoint),
              atrisk_midpoint_check = mean(atrisk_midpoint_check),
              atrisk_midpoint_ORR = mean(atrisk_midpoint_ORR),
              DoR = DoR,
              DoRInt_no_conf = DoRInt_no_conf,
              DoRIntM_no_conf = DoRIntM_no_conf,
              DoRInt_conf = DoRInt_conf,
              DoRIntM_conf = DoRIntM_conf,
              DoRInt_ORR = DoRInt_ORR,
              DoRIntM_ORR = DoRIntM_ORR))
}

# Run simulation setting to replicate the experimental arm in the paper.
# We want to compare and see differences in the average duration of response 
# in three levels. Censoring (right, midpoint, no-censoring), schedule (1 month, 2 months, 3 months) 
# for different tau values. 
#########################################################################################################
# Experimental arm
###################
set.seed(2405)
S1 = 0:100
S2 = seq(0,100,2)
S3 = seq(0,100,3)
S4 = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72)
S5 = c(0, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 21, 27, 33, 39, 45, 51, 57, 63, 69)

tau.val <- 50
iter <- 1000

# Strategy S1
#############

database_exp_400_scheme_S1 <- as.data.frame(matrix(NA, tau.val, 4))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S1)
  database_exp_400_scheme_S1[z,] <- c(k$Average_DoR, k$Average_DoR_Int, k$Average_DoR_IntMid, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S1) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "tau")

##############

# Strategy S2
#############
set.seed(2405)

database_exp_400_scheme_S2 <- as.data.frame(matrix(NA, tau.val, 4))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S2)
  database_exp_400_scheme_S2[z,] <- c(k$Average_DoR, k$Average_DoR_Int, k$Average_DoR_IntMid, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S2) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "tau")

##############

# Strategy S3
#############
set.seed(2405)

database_exp_400_scheme_S3 <- as.data.frame(matrix(NA, tau.val, 4))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S3)
  database_exp_400_scheme_S3[z,] <- c(k$Average_DoR, k$Average_DoR_Int, k$Average_DoR_IntMid, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S3) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "tau")

##############

# Strategy S4
#############
set.seed(2405)

database_exp_400_scheme_S4 <- as.data.frame(matrix(NA, tau.val, 4))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S4)
  database_exp_400_scheme_S4[z,] <- c(k$Average_DoR, k$Average_DoR_Int, k$Average_DoR_IntMid, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S4) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "tau")

##############
# Strategy S5
#############
set.seed(2405)

database_exp_400_scheme_S5 <- as.data.frame(matrix(NA, tau.val, 4))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S5)
  database_exp_400_scheme_S5[z,] <- c(k$Average_DoR, k$Average_DoR_Int, k$Average_DoR_IntMid, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S5) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "tau")

##############
plot(database_exp_400_scheme_S1$tau, database_exp_400_scheme_S1$Average_DoR, type = 'l', ylim = c(0,17), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, database_exp_400_scheme_S1$Average_DoR_Int,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))
lines(database_exp_400_scheme_S1$tau, database_exp_400_scheme_S1$Average_DoR_Int_Mid,  lty=3  , col = 'darkgreen', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))

lines(database_exp_400_scheme_S2$tau, database_exp_400_scheme_S2$Average_DoR, pch="*",type = 'l', ylim = c(0,17), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S2$tau, database_exp_400_scheme_S2$Average_DoR_Int, pch="*", lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))
lines(database_exp_400_scheme_S2$tau, database_exp_400_scheme_S2$Average_DoR_Int_Mid, pch="*", lty=3  , col = 'darkgreen', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))

lines(database_exp_400_scheme_S3$tau, database_exp_400_scheme_S3$Average_DoR,pch="+", type = 'l', ylim = c(0,17), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S3$tau, database_exp_400_scheme_S3$Average_DoR_Int, pch="+", lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))
lines(database_exp_400_scheme_S3$tau, database_exp_400_scheme_S3$Average_DoR_Int_Mid, pch="+", lty=3  , col = 'darkgreen', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))

lines(database_exp_400_scheme_S4$tau, database_exp_400_scheme_S4$Average_DoR,pch="+", type = 'l', ylim = c(0,17), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S4$tau, database_exp_400_scheme_S4$Average_DoR_Int, pch="+", lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))
lines(database_exp_400_scheme_S4$tau, database_exp_400_scheme_S4$Average_DoR_Int_Mid, pch="+", lty=3  , col = 'darkgreen', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))


par(new=T)

# Control arm
#############
# Strategy S1
#############
set.seed(2405)

database_ctl_400_scheme_S1 <- as.data.frame(matrix(NA, tau.val, 8))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S1)
  database_ctl_400_scheme_S1[z,] <- c(k$Average_DoR, k$Average_DoR_Int_no_conf, k$Average_DoR_IntMid_no_conf, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S1) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############

# Strategy S2
#############
set.seed(2405)

database_ctl_400_scheme_S2 <- as.data.frame(matrix(NA, tau.val, 8))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S2)
  database_ctl_400_scheme_S2[z,] <- c(k$Average_DoR, k$Average_DoR_Int_no_conf, k$Average_DoR_IntMid_no_conf, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S2) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############

# Strategy S3
#############
set.seed(2405)

database_ctl_400_scheme_S3 <- as.data.frame(matrix(NA, tau.val, 4))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S3)
  database_ctl_400_scheme_S3[z,] <- c(k$Average_DoR, k$Average_DoR_Int, k$Average_DoR_IntMid, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S3) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "tau")

##############

# Strategy S4
#############
set.seed(2405)

database_ctl_400_scheme_S4 <- as.data.frame(matrix(NA, tau.val, 4))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S4)
  database_ctl_400_scheme_S4[z,] <- c(k$Average_DoR, k$Average_DoR_Int, k$Average_DoR_IntMid, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S4) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "tau")

##############
# Strategy S5
#############
set.seed(2405)

database_ctl_400_scheme_S5 <- as.data.frame(matrix(NA, tau.val, 4))
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S5)
  database_ctl_400_scheme_S5[z,] <- c(k$Average_DoR, k$Average_DoR_Int, k$Average_DoR_IntMid, tau)
  z <- z + 1 
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S5) <- c("Average_DoR", "Average_DoR_Int", "Average_DoR_Int_Mid", "tau")

##############
plot(database_ctl_400_scheme_S1$tau, database_ctl_400_scheme_S1$Average_DoR, type = 'l', ylim = c(0,17), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, database_ctl_400_scheme_S1$Average_DoR_Int,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))
lines(database_ctl_400_scheme_S1$tau, database_ctl_400_scheme_S1$Average_DoR_Int_Mid,  lty=3  , col = 'darkgreen', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))

lines(database_ctl_400_scheme_S2$tau, database_ctl_400_scheme_S2$Average_DoR, pch="*",type = 'l', ylim = c(0,17), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S2$tau, database_ctl_400_scheme_S2$Average_DoR_Int, pch="*", lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))
lines(database_ctl_400_scheme_S2$tau, database_ctl_400_scheme_S2$Average_DoR_Int_Mid, pch="*", lty=3  , col = 'darkgreen', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))

lines(database_ctl_400_scheme_S3$tau, database_ctl_400_scheme_S3$Average_DoR,pch="+", type = 'l', ylim = c(0,17), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S3$tau, database_ctl_400_scheme_S3$Average_DoR_Int, pch="+", lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))
lines(database_ctl_400_scheme_S3$tau, database_ctl_400_scheme_S3$Average_DoR_Int_Mid, pch="+", lty=3  , col = 'darkgreen', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))

lines(database_ctl_400_scheme_S4$tau, database_ctl_400_scheme_S4$Average_DoR,pch="+", type = 'l', ylim = c(0,17), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S4$tau, database_ctl_400_scheme_S4$Average_DoR_Int, pch="+", lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))
lines(database_ctl_400_scheme_S4$tau, database_ctl_400_scheme_S4$Average_DoR_Int_Mid, pch="+", lty=3  , col = 'darkgreen', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,17))

# Experimental square difference
par(mfrow = c(2,2))
plot(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.1), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.1))

plot(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.3), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.3))

plot(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.7), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.7))

plot(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int)^2, type = 'l', ylim = c(0,2.2), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,2.2))

plot(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int)^2, type = 'l', ylim = c(0,2.2), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,2.2))

# Control square difference
par(mfrow = c(2,2))
plot(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.05), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.05))

plot(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.1), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.3))

plot(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.7), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.7))

plot(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int)^2, type = 'l', ylim = c(0,2.2), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,2.2))

plot(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int)^2, type = 'l', ylim = c(0,2.2), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,2.2))

# Control and experimental square difference
par(mfrow = c(3,2))
plot(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.1), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.1))
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.1), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.1))
abline(v = S1,  cex=0.01, lwd=0.1, lty="dotted")

plot(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.3), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.3))
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.3), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.3))
abline(v = S2,  cex=0.01, lwd=0.1, lty="dotted")

plot(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.7), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.7))
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int)^2, type = 'l', ylim = c(0,0.7), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,0.7))
abline(v = S3,  cex=0.01, lwd=0.1, lty="dotted")

plot(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int)^2, type = 'l', ylim = c(0,2.2), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,2.2))
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int)^2, type = 'l', ylim = c(0,2.2), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid)^2,  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,2.2))
abline(v = S4,  cex=0.01, lwd=0.1, lty="dotted")

plot(database_exp_400_scheme_S1$tau, ((database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int)^2), type = 'l', ylim = c(0,1), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_exp_400_scheme_S1$tau, ((database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid)^2),  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,1))
lines(database_ctl_400_scheme_S1$tau, ((database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int)^2), type = 'l', ylim = c(0,1), xlab = "Tau", ylab = "Average DoR daily vs interval")
lines(database_ctl_400_scheme_S1$tau, ((database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid)^2),  lty=2  , col = 'red', xlab = "Tau", ylab = "Average DoR daily vs interval", ylim = c(0,1))
abline(v = S5,  cex=0.01, lwd=0.1, lty="dotted")

dataset <- as.data.frame(cbind(Tau = rep(database_exp_400_scheme_S1$tau, 20),
                               Trial = rep(rep(c("Experimental","Control"), each = 100), 5),
                               Censoring = rep(rep(c("Interval","Midpoint"), each = 50), 10),
                               Scenario = rep(c("Scenario1","Scenario2","Scenario3","Scenario4","Scenario5"), each = 200),
                               DoR_diff_square = round(c((database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid)^2,
                                            (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid)^2,
                                            (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid)^2,
                                            (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid)^2,
                                            (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid)^2,
                                            (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid)^2,
                                            (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid)^2,
                                            (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid)^2,
                                            (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid)^2,
                                            (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid)^2),4),
                           DoR_diff_abs = round(c(abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int), abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid),
                                            abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int), abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid),
                                            abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int), abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid),
                                            abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int), abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid),
                                            abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int), abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid),
                                            abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int), abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid),
                                            abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int), abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid),
                                            abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int), abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid),
                                            abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int), abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid),
                                            abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int), abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid)),4),
                           DoR_diff = round(c((database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid),
                                                  (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid),
                                                  (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid),
                                                  (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid),
                                                  (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid),
                                                  (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid),
                                                  (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid),
                                                  (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid),
                                                  (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid),
                                                  (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid)),4)))


g1 <- ggplot(dataset, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff_square), color = Censoring)) + geom_line(aes(linetype = Trial)) + facet_grid(~Scenario) + theme_light()
g2 <- ggplot(dataset, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff_abs), color = Censoring)) + geom_line(aes(linetype = Trial)) + facet_grid(~Scenario) + theme_light()
g3 <- ggplot(dataset, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff), color = Censoring)) + geom_line(aes(linetype = Trial)) + facet_grid(~Scenario) + theme_light()








#plot(database_ctl_400_scheme_S2$tau, database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int, type = 'l', ylim = c(-0.2,0.4))
#lines(database_ctl_400_scheme_S2$tau, database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid, type = 'o', col = 1, ylim = c(-0.2,0.4))
#lines(database_ctl_400_scheme_S2$tau, database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_conf, type = 'l', col = 2, ylim = c(-0.2,0.4))
#lines(database_ctl_400_scheme_S2$tau, database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf, type = 'o', col = 2, ylim = c(-0.2,0.4))
#lines(database_ctl_400_scheme_S2$tau, database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR, type = 'l', col = 3, ylim = c(-0.2,0.4))
#lines(database_ctl_400_scheme_S2$tau, database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR_Mid, type = 'o', col = 3, ylim = c(-0.2,0.4))

# Scenarios including all estimated scenarios 



















# Run simulation setting to replicate the experimental arm in the paper.
# We want to compare and see differences in the average duration of response 
# in three levels. Censoring (right, midpoint, no-censoring), schedule (1 month, 2 months, 3 months) 
# for different tau values. 
#########################################################################################################
# Experimental arm
###################
set.seed(2405)
S1 = 0:100
S2 = seq(0,100,2)
S3 = seq(0,100,3)
S4 = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72)
S5 = c(0, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 21, 27, 33, 39, 45, 51, 57, 63, 69)

tau.val <- 50
iter <- 1000

# Strategy S1
#############

database_exp_400_scheme_S1 <- as.data.frame(matrix(NA, tau.val, 14))
dor_exp_1 <- NULL
dor_exp_1_no_conf <- NULL
dor_exp_1_mid_no_conf <- NULL
dor_exp_1_conf <- NULL
dor_exp_1_mid_conf <- NULL
dor_exp_1_ORR <- NULL
dor_exp_1_mid_ORR <- NULL
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S1)
  database_exp_400_scheme_S1[z,] <- c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_exp_1[[tau]] <- k$DoR
  dor_exp_1_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_exp_1_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_exp_1_conf[[tau]]        <- k$DoRInt_conf
  dor_exp_1_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_exp_1_ORR[[tau]]         <- k$DoRInt_ORR
  dor_exp_1_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S1) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############

# Strategy S2
#############
set.seed(2405)

database_exp_400_scheme_S2 <- as.data.frame(matrix(NA, tau.val, 14))
z <- 1
dor_exp_2 <- NULL
dor_exp_2_no_conf <- NULL
dor_exp_2_mid_no_conf <- NULL
dor_exp_2_conf <- NULL
dor_exp_2_mid_conf <- NULL
dor_exp_2_ORR <- NULL
dor_exp_2_mid_ORR <- NULL
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S2)
  database_exp_400_scheme_S2[z,] <-  c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_exp_2[[tau]] <- k$DoR
  dor_exp_2_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_exp_2_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_exp_2_conf[[tau]]        <- k$DoRInt_conf
  dor_exp_2_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_exp_2_ORR[[tau]]         <- k$DoRInt_ORR
  dor_exp_2_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S2) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############

# Strategy S3
#############
set.seed(2405)

database_exp_400_scheme_S3 <- as.data.frame(matrix(NA, tau.val, 14))
z <- 1
dor_exp_3 <- NULL
dor_exp_3_no_conf <- NULL
dor_exp_3_mid_no_conf <- NULL
dor_exp_3_conf <- NULL
dor_exp_3_mid_conf <- NULL
dor_exp_3_ORR <- NULL
dor_exp_3_mid_ORR <- NULL
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S3)
  database_exp_400_scheme_S3[z,] <-  c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_exp_3[[tau]] <- k$DoR
  dor_exp_3_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_exp_3_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_exp_3_conf[[tau]]        <- k$DoRInt_conf
  dor_exp_3_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_exp_3_ORR[[tau]]         <- k$DoRInt_ORR
  dor_exp_3_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S3) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############

# Strategy S4
#############
set.seed(2405)

database_exp_400_scheme_S4 <- as.data.frame(matrix(NA, tau.val, 14))
z <- 1
dor_exp_4 <- NULL
dor_exp_4_no_conf <- NULL
dor_exp_4_mid_no_conf <- NULL
dor_exp_4_conf <- NULL
dor_exp_4_mid_conf <- NULL
dor_exp_4_ORR <- NULL
dor_exp_4_mid_ORR <- NULL
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S4)
  database_exp_400_scheme_S4[z,] <-  c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_exp_4[[tau]] <- k$DoR
  dor_exp_4_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_exp_4_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_exp_4_conf[[tau]]        <- k$DoRInt_conf
  dor_exp_4_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_exp_4_ORR[[tau]]         <- k$DoRInt_ORR
  dor_exp_4_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S4) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############
# Strategy S5
#############
set.seed(2405)

database_exp_400_scheme_S5 <- as.data.frame(matrix(NA, tau.val, 14))
z <- 1
dor_exp_5 <- NULL
dor_exp_5_no_conf <- NULL
dor_exp_5_mid_no_conf <- NULL
dor_exp_5_conf <- NULL
dor_exp_5_mid_conf <- NULL
dor_exp_5_ORR <- NULL
dor_exp_5_mid_ORR <- NULL
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S5)
  database_exp_400_scheme_S5[z,] <-  c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_exp_5[[tau]] <- k$DoR
  dor_exp_5_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_exp_5_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_exp_5_conf[[tau]]        <- k$DoRInt_conf
  dor_exp_5_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_exp_5_ORR[[tau]]         <- k$DoRInt_ORR
  dor_exp_5_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_exp_400_scheme_S5) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

#############
# Control arm
#############
# Strategy S1
#############
set.seed(2405)

database_ctl_400_scheme_S1 <- as.data.frame(matrix(NA, tau.val, 14))
z <- 1
dor_ctl_1 <- NULL
dor_ctl_1_no_conf <- NULL
dor_ctl_1_mid_no_conf <- NULL
dor_ctl_1_conf <- NULL
dor_ctl_1_mid_conf <- NULL
dor_ctl_1_ORR <- NULL
dor_ctl_1_mid_ORR <- NULL
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S1)
  database_ctl_400_scheme_S1[z,] <-  c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_ctl_1[[tau]] <- k$DoR
  dor_ctl_1_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_ctl_1_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_ctl_1_conf[[tau]]        <- k$DoRInt_conf
  dor_ctl_1_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_ctl_1_ORR[[tau]]         <- k$DoRInt_ORR
  dor_ctl_1_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S1) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############

# Strategy S2
#############
set.seed(2405)

database_ctl_400_scheme_S2 <- as.data.frame(matrix(NA, tau.val, 14))
z <- 1
dor_ctl_2 <- NULL
dor_ctl_2_no_conf <- NULL
dor_ctl_2_mid_no_conf <- NULL
dor_ctl_2_conf <- NULL
dor_ctl_2_mid_conf <- NULL
dor_ctl_2_ORR <- NULL
dor_ctl_2_mid_ORR <- NULL
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S2)
  database_ctl_400_scheme_S2[z,] <-  c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_ctl_2[[tau]] <- k$DoR
  dor_ctl_2_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_ctl_2_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_ctl_2_conf[[tau]]        <- k$DoRInt_conf
  dor_ctl_2_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_ctl_2_ORR[[tau]]         <- k$DoRInt_ORR
  dor_ctl_2_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S2) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############

# Strategy S3
#############
set.seed(2405)

database_ctl_400_scheme_S3 <- as.data.frame(matrix(NA, tau.val, 14))
z <- 1
dor_ctl_3 <- NULL
dor_ctl_3_no_conf <- NULL
dor_ctl_3_mid_no_conf <- NULL
dor_ctl_3_conf <- NULL
dor_ctl_3_mid_conf <- NULL
dor_ctl_3_ORR <- NULL
dor_ctl_3_mid_ORR <- NULL
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S3)
  database_ctl_400_scheme_S3[z,] <-  c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_ctl_3[[tau]] <- k$DoR
  dor_ctl_3_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_ctl_3_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_ctl_3_conf[[tau]]        <- k$DoRInt_conf
  dor_ctl_3_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_ctl_3_ORR[[tau]]         <- k$DoRInt_ORR
  dor_ctl_3_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S3) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############

# Strategy S4
#############
set.seed(2405)

database_ctl_400_scheme_S4 <- as.data.frame(matrix(NA, tau.val, 14))
z <- 1
dor_ctl_4 <- NULL
dor_ctl_4_no_conf <- NULL
dor_ctl_4_mid_no_conf <- NULL
dor_ctl_4_conf <- NULL
dor_ctl_4_mid_conf <- NULL
dor_ctl_4_ORR <- NULL
dor_ctl_4_mid_ORR <- NULL
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S4)
  database_ctl_400_scheme_S4[z,] <-  c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_ctl_4[[tau]] <- k$DoR
  dor_ctl_4_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_ctl_4_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_ctl_4_conf[[tau]]        <- k$DoRInt_conf
  dor_ctl_4_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_ctl_4_ORR[[tau]]         <- k$DoRInt_ORR
  dor_ctl_4_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S4) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")

##############
# Strategy S5
#############
set.seed(2405)

database_ctl_400_scheme_S5 <- as.data.frame(matrix(NA, tau.val, 14))
z <- 1
dor_ctl_5 <- NULL
dor_ctl_5_no_conf <- NULL
dor_ctl_5_mid_no_conf <- NULL
dor_ctl_5_conf <- NULL
dor_ctl_5_mid_conf <- NULL
dor_ctl_5_ORR <- NULL
dor_ctl_5_mid_ORR <- NULL
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.715037, scaleP = 13.35673, shapeR = 1.7, scaleR = 10, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S5)
  database_ctl_400_scheme_S5[z,] <-  c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_ctl_5[[tau]] <- k$DoR
  dor_ctl_5_no_conf[[tau]]     <- k$DoRInt_no_conf
  dor_ctl_5_mid_no_conf[[tau]] <- k$DoRIntM_no_conf
  dor_ctl_5_conf[[tau]]        <- k$DoRInt_conf
  dor_ctl_5_mid_conf[[tau]]    <- k$DoRIntM_conf
  dor_ctl_5_ORR[[tau]]         <- k$DoRInt_ORR
  dor_ctl_5_mid_ORR[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(database_ctl_400_scheme_S5) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")
##############

##############################################
# Calculate the CI for the difference of DoR #
##############################################
dor_1_diff <- mapply('-', dor_exp_1, dor_ctl_1, SIMPLIFY = FALSE)
dor_2_diff <- mapply('-', dor_exp_2, dor_ctl_2, SIMPLIFY = FALSE)
dor_3_diff <- mapply('-', dor_exp_3, dor_ctl_3, SIMPLIFY = FALSE)
dor_4_diff <- mapply('-', dor_exp_4, dor_ctl_4, SIMPLIFY = FALSE)
dor_5_diff <- mapply('-', dor_exp_5, dor_ctl_5, SIMPLIFY = FALSE)

dor_1_diff_low <- as.numeric(unlist(lapply(dor_1_diff, quantile, probs = 0.05)))
dor_1_diff_hi  <- as.numeric(unlist(lapply(dor_1_diff, quantile, probs = 0.95)))

dor_2_diff_low <- as.numeric(unlist(lapply(dor_2_diff, quantile, probs = 0.05)))
dor_2_diff_hi  <- as.numeric(unlist(lapply(dor_2_diff, quantile, probs = 0.95)))

dor_3_diff_low <- as.numeric(unlist(lapply(dor_3_diff, quantile, probs = 0.05)))
dor_3_diff_hi  <- as.numeric(unlist(lapply(dor_3_diff, quantile, probs = 0.95)))

dor_4_diff_low <- as.numeric(unlist(lapply(dor_4_diff, quantile, probs = 0.05)))
dor_4_diff_hi  <- as.numeric(unlist(lapply(dor_4_diff, quantile, probs = 0.95)))

dor_5_diff_low <- as.numeric(unlist(lapply(dor_5_diff, quantile, probs = 0.05)))
dor_5_diff_hi  <- as.numeric(unlist(lapply(dor_5_diff, quantile, probs = 0.95)))
################################
dor_1_ratio <- mapply('/', dor_ctl_1, dor_exp_1, SIMPLIFY = FALSE)
dor_2_ratio <- mapply('/', dor_ctl_2, dor_exp_2, SIMPLIFY = FALSE)
dor_3_ratio <- mapply('/', dor_ctl_3, dor_exp_3, SIMPLIFY = FALSE)
dor_4_ratio <- mapply('/', dor_ctl_4, dor_exp_4, SIMPLIFY = FALSE)
dor_5_ratio <- mapply('/', dor_ctl_5, dor_exp_5, SIMPLIFY = FALSE)

dor_1_ratio_low <- as.numeric(unlist(lapply(dor_1_ratio, quantile, probs = 0.05)))
dor_1_ratio_hi  <- as.numeric(unlist(lapply(dor_1_ratio, quantile, probs = 0.95)))

dor_2_ratio_low <- as.numeric(unlist(lapply(dor_2_ratio, quantile, probs = 0.05)))
dor_2_ratio_hi  <- as.numeric(unlist(lapply(dor_2_ratio, quantile, probs = 0.95)))

dor_3_ratio_low <- as.numeric(unlist(lapply(dor_3_ratio, quantile, probs = 0.05)))
dor_3_ratio_hi  <- as.numeric(unlist(lapply(dor_3_ratio, quantile, probs = 0.95)))

dor_4_ratio_low <- as.numeric(unlist(lapply(dor_4_ratio, quantile, probs = 0.05)))
dor_4_ratio_hi  <- as.numeric(unlist(lapply(dor_4_ratio, quantile, probs = 0.95)))

dor_5_ratio_low <- as.numeric(unlist(lapply(dor_5_ratio, quantile, probs = 0.05)))
dor_5_ratio_hi  <- as.numeric(unlist(lapply(dor_5_ratio, quantile, probs = 0.95)))




dor_1_right_noconf_diff <- mapply('-', dor_exp_1_no_conf, dor_ctl_1_no_conf, SIMPLIFY = FALSE)
dor_2_right_noconf_diff <- mapply('-', dor_exp_2_no_conf, dor_ctl_2_no_conf, SIMPLIFY = FALSE)
dor_3_right_noconf_diff <- mapply('-', dor_exp_3_no_conf, dor_ctl_3_no_conf, SIMPLIFY = FALSE)
dor_4_right_noconf_diff <- mapply('-', dor_exp_4_no_conf, dor_ctl_4_no_conf, SIMPLIFY = FALSE)
dor_5_right_noconf_diff <- mapply('-', dor_exp_5_no_conf, dor_ctl_5_no_conf, SIMPLIFY = FALSE)

dor_1_right_noconf_diff_low <- as.numeric(unlist(lapply(dor_1_right_noconf_diff, quantile, probs = 0.05)))
dor_1_right_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_1_right_noconf_diff, quantile, probs = 0.95)))

dor_2_right_noconf_diff_low <- as.numeric(unlist(lapply(dor_2_right_noconf_diff, quantile, probs = 0.05)))
dor_2_right_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_2_right_noconf_diff, quantile, probs = 0.95)))

dor_3_right_noconf_diff_low <- as.numeric(unlist(lapply(dor_3_right_noconf_diff, quantile, probs = 0.05)))
dor_3_right_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_3_right_noconf_diff, quantile, probs = 0.95)))

dor_4_right_noconf_diff_low <- as.numeric(unlist(lapply(dor_4_right_noconf_diff, quantile, probs = 0.05)))
dor_4_right_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_4_right_noconf_diff, quantile, probs = 0.95)))

dor_5_right_noconf_diff_low <- as.numeric(unlist(lapply(dor_5_right_noconf_diff, quantile, probs = 0.05)))
dor_5_right_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_5_right_noconf_diff, quantile, probs = 0.95)))
#######################################################

dor_1_right_noconf_ratio <- mapply('/', dor_ctl_1_no_conf, dor_exp_1_no_conf, SIMPLIFY = FALSE)
dor_2_right_noconf_ratio <- mapply('/', dor_ctl_2_no_conf, dor_exp_2_no_conf, SIMPLIFY = FALSE)
dor_3_right_noconf_ratio <- mapply('/', dor_ctl_3_no_conf, dor_exp_3_no_conf, SIMPLIFY = FALSE)
dor_4_right_noconf_ratio <- mapply('/', dor_ctl_4_no_conf, dor_exp_4_no_conf, SIMPLIFY = FALSE)
dor_5_right_noconf_ratio <- mapply('/', dor_ctl_5_no_conf, dor_exp_5_no_conf, SIMPLIFY = FALSE)

dor_1_right_noconf_ratio[[1]] <- rep(0, 1000)
dor_2_right_noconf_ratio[[1]] <- rep(0, 1000)
dor_2_right_noconf_ratio[[2]] <- rep(0, 1000)
dor_3_right_noconf_ratio[[1]] <- rep(0, 1000)
dor_3_right_noconf_ratio[[2]] <- rep(0, 1000)
dor_3_right_noconf_ratio[[3]] <- rep(0, 1000)
dor_4_right_noconf_ratio[[1]] <- rep(0, 1000)
dor_5_right_noconf_ratio[[1]] <- rep(0, 1000)

dor_1_right_noconf_ratio_low <- as.numeric(unlist(lapply(dor_1_right_noconf_ratio, quantile, probs = 0.05)))
dor_1_right_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_1_right_noconf_ratio, quantile, probs = 0.95)))

dor_2_right_noconf_ratio_low <- as.numeric(unlist(lapply(dor_2_right_noconf_ratio, quantile, probs = 0.05)))
dor_2_right_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_2_right_noconf_ratio, quantile, probs = 0.95)))

dor_3_right_noconf_ratio_low <- as.numeric(unlist(lapply(dor_3_right_noconf_ratio, quantile, probs = 0.05)))
dor_3_right_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_3_right_noconf_ratio, quantile, probs = 0.95)))

dor_4_right_noconf_ratio_low <- as.numeric(unlist(lapply(dor_4_right_noconf_ratio, quantile, probs = 0.05)))
dor_4_right_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_4_right_noconf_ratio, quantile, probs = 0.95)))

dor_5_right_noconf_ratio_low <- as.numeric(unlist(lapply(dor_5_right_noconf_ratio, quantile, probs = 0.05)))
dor_5_right_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_5_right_noconf_ratio, quantile, probs = 0.95)))








dor_1_mid_noconf_diff <- mapply('-', dor_exp_1_mid_no_conf, dor_ctl_1_mid_no_conf, SIMPLIFY = FALSE)
dor_2_mid_noconf_diff <- mapply('-', dor_exp_2_mid_no_conf, dor_ctl_2_mid_no_conf, SIMPLIFY = FALSE)
dor_3_mid_noconf_diff <- mapply('-', dor_exp_3_mid_no_conf, dor_ctl_3_mid_no_conf, SIMPLIFY = FALSE)
dor_4_mid_noconf_diff <- mapply('-', dor_exp_4_mid_no_conf, dor_ctl_4_mid_no_conf, SIMPLIFY = FALSE)
dor_5_mid_noconf_diff <- mapply('-', dor_exp_5_mid_no_conf, dor_ctl_5_mid_no_conf, SIMPLIFY = FALSE)


dor_1_mid_noconf_diff_low <- as.numeric(unlist(lapply(dor_1_mid_noconf_diff, quantile, probs = 0.05)))
dor_1_mid_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_1_mid_noconf_diff, quantile, probs = 0.95)))

dor_2_mid_noconf_diff_low <- as.numeric(unlist(lapply(dor_2_mid_noconf_diff, quantile, probs = 0.05)))
dor_2_mid_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_2_mid_noconf_diff, quantile, probs = 0.95)))

dor_3_mid_noconf_diff_low <- as.numeric(unlist(lapply(dor_3_mid_noconf_diff, quantile, probs = 0.05)))
dor_3_mid_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_3_mid_noconf_diff, quantile, probs = 0.95)))

dor_4_mid_noconf_diff_low <- as.numeric(unlist(lapply(dor_4_mid_noconf_diff, quantile, probs = 0.05)))
dor_4_mid_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_4_mid_noconf_diff, quantile, probs = 0.95)))

dor_5_mid_noconf_diff_low <- as.numeric(unlist(lapply(dor_5_mid_noconf_diff, quantile, probs = 0.05)))
dor_5_mid_noconf_diff_hi  <- as.numeric(unlist(lapply(dor_5_mid_noconf_diff, quantile, probs = 0.95)))

######################################################

dor_1_mid_noconf_ratio <- mapply('/', dor_ctl_1_mid_no_conf, dor_exp_1_mid_no_conf, SIMPLIFY = FALSE)
dor_2_mid_noconf_ratio <- mapply('/', dor_ctl_2_mid_no_conf, dor_exp_2_mid_no_conf, SIMPLIFY = FALSE)
dor_3_mid_noconf_ratio <- mapply('/', dor_ctl_3_mid_no_conf, dor_exp_3_mid_no_conf, SIMPLIFY = FALSE)
dor_4_mid_noconf_ratio <- mapply('/', dor_ctl_4_mid_no_conf, dor_exp_4_mid_no_conf, SIMPLIFY = FALSE)
dor_5_mid_noconf_ratio <- mapply('/', dor_ctl_5_mid_no_conf, dor_exp_5_mid_no_conf, SIMPLIFY = FALSE)

dor_1_mid_noconf_ratio[[1]] <- rep(0, 1000)
dor_2_mid_noconf_ratio[[1]] <- rep(0, 1000)
dor_2_mid_noconf_ratio[[2]] <- rep(0, 1000)
dor_3_mid_noconf_ratio[[1]] <- rep(0, 1000)
dor_3_mid_noconf_ratio[[2]] <- rep(0, 1000)
dor_3_mid_noconf_ratio[[3]] <- rep(0, 1000)
dor_4_mid_noconf_ratio[[1]] <- rep(0, 1000)
dor_5_mid_noconf_ratio[[1]] <- rep(0, 1000)

dor_1_mid_noconf_ratio_low <- as.numeric(unlist(lapply(dor_1_mid_noconf_ratio, quantile, probs = 0.05)))
dor_1_mid_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_1_mid_noconf_ratio, quantile, probs = 0.95)))

dor_2_mid_noconf_ratio_low <- as.numeric(unlist(lapply(dor_2_mid_noconf_ratio, quantile, probs = 0.05)))
dor_2_mid_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_2_mid_noconf_ratio, quantile, probs = 0.95)))

dor_3_mid_noconf_ratio_low <- as.numeric(unlist(lapply(dor_3_mid_noconf_ratio, quantile, probs = 0.05)))
dor_3_mid_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_3_mid_noconf_ratio, quantile, probs = 0.95)))

dor_4_mid_noconf_ratio_low <- as.numeric(unlist(lapply(dor_4_mid_noconf_ratio, quantile, probs = 0.05)))
dor_4_mid_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_4_mid_noconf_ratio, quantile, probs = 0.95)))

dor_5_mid_noconf_ratio_low <- as.numeric(unlist(lapply(dor_5_mid_noconf_ratio, quantile, probs = 0.05)))
dor_5_mid_noconf_ratio_hi  <- as.numeric(unlist(lapply(dor_5_mid_noconf_ratio, quantile, probs = 0.95)))


#dataset <- as.data.frame(cbind(Tau = rep(database_exp_400_scheme_S1$tau, 60),
#                               Trial = rep(rep(c("Experimental","Control"), each = 300), 5),
#                               Censoring = rep(rep(c("Interval","Midpoint","Interval_check","Midpoint_ckeck","Interval_ORR","Midpoint_ORR"), each = 50), 10),
#                               Scenario = rep(c("Scenario1","Scenario2","Scenario3","Scenario4","Scenario5"), each = 600),
#                               DoR_diff_square = round(c((database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR_Mid)^2,
#                                                         (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR_Mid)^2,
#                                                         (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR_Mid)^2,
#                                                         (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR_Mid)^2,
#                                                         (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR_Mid)^2,
#                                                         (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR_Mid)^2,
#                                                         (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR_Mid)^2,
#                                                         (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR_Mid)^2,
#                                                         (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR_Mid)^2,
#                                                         (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR_Mid)^2),4),
#                               DoR_diff_abs = round(c(abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int), abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid), abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_conf), abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf), abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR), abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR_Mid),
#                                                      abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int), abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid), abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_conf), abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf), abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR), abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR_Mid),
#                                                      abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int), abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid), abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_conf), abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid_conf), abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR), abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR_Mid),
#                                                      abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int), abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid), abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_conf), abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf), abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR), abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR_Mid),
#                                                      abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int), abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid), abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_conf), abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid_conf), abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR), abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR_Mid),
#                                                      abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int), abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid), abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_conf), abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid_conf), abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR), abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR_Mid),
#                                                      abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int), abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid), abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_conf), abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid_conf), abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR), abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR_Mid),
#                                                      abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int), abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid), abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_conf), abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid_conf), abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR), abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR_Mid),
#                                                      abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int), abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid), abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_conf), abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid_conf), abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR), abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR_Mid), 
#                                                      abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int), abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid), abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_conf), abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid_conf), abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR), abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR_Mid)),4),
#                               DoR_diff = (-1)*round(c((database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_conf), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR_Mid),
#                                                  (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_conf), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR_Mid),
#                                                  (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_conf), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR_Mid),
#                                                  (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_conf), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR_Mid),
#                                                  (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_conf), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR_Mid),
#                                                  (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_conf), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR_Mid),
#                                                  (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_conf), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR_Mid),
#                                                  (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_conf), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR_Mid),
#                                                  (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_conf), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR_Mid), 
#                                                  (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_conf), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR_Mid)),4)))

dataset <- as.data.frame(cbind(Tau = rep(database_exp_400_scheme_S1$tau, 60),
                               Trial = rep(rep(c("Setting A","Setting B"), each = 300), 5),
                               Censoring = rep(rep(c("E1","E2","D1","D2","Interval_ORR","Midpoint_ORR"), each = 50), 10),
                               Scenario = rep(c("Scenario1","Scenario2","Scenario3","Scenario4","Scenario5"), each = 600),
                               DoR_diff_square = round(c((database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR)^2, (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR_Mid)^2,
                                                         (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR)^2, (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR_Mid)^2,
                                                         (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR)^2, (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR_Mid)^2,
                                                         (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR)^2, (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR_Mid)^2,
                                                         (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR)^2, (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR_Mid)^2,
                                                         (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR)^2, (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR_Mid)^2,
                                                         (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR)^2, (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR_Mid)^2,
                                                         (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR)^2, (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR_Mid)^2,
                                                         (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_conf)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid_conf)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR)^2, (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR_Mid)^2,
                                                         (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_conf)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid_conf)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR)^2, (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR_Mid)^2),4),
                               DoR_diff_relative = round(c((database_exp_400_scheme_S1$Average_DoR_Int - database_exp_400_scheme_S1$Average_DoR)/database_exp_400_scheme_S1$Average_DoR_Int, abs(database_exp_400_scheme_S1$Average_DoR_Int_Mid - database_exp_400_scheme_S1$Average_DoR)/database_exp_400_scheme_S1$Average_DoR_Int_Mid, (database_exp_400_scheme_S1$Average_DoR_Int_conf - database_exp_400_scheme_S1$Average_DoR)/database_exp_400_scheme_S1$Average_DoR_Int_conf, (database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf - database_exp_400_scheme_S1$Average_DoR)/database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf, abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR), abs(database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR_Mid),
                                                           (database_ctl_400_scheme_S1$Average_DoR_Int - database_ctl_400_scheme_S1$Average_DoR)/database_exp_400_scheme_S1$Average_DoR_Int, abs(database_ctl_400_scheme_S1$Average_DoR_Int_Mid - database_ctl_400_scheme_S1$Average_DoR)/database_ctl_400_scheme_S1$Average_DoR_Int_Mid, (database_ctl_400_scheme_S1$Average_DoR_Int_conf - database_ctl_400_scheme_S1$Average_DoR)/database_ctl_400_scheme_S1$Average_DoR_Int_conf, (database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S1$Average_DoR)/database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf, abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR), abs(database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR_Mid),
                                                           (database_exp_400_scheme_S2$Average_DoR_Int - database_exp_400_scheme_S2$Average_DoR)/database_exp_400_scheme_S2$Average_DoR_Int, abs(database_exp_400_scheme_S2$Average_DoR_Int_Mid - database_exp_400_scheme_S2$Average_DoR)/database_exp_400_scheme_S2$Average_DoR_Int_Mid, (database_exp_400_scheme_S2$Average_DoR_Int_conf - database_exp_400_scheme_S2$Average_DoR)/database_exp_400_scheme_S2$Average_DoR_Int_conf, (database_exp_400_scheme_S2$Average_DoR_Int_Mid_conf - database_exp_400_scheme_S2$Average_DoR)/database_exp_400_scheme_S2$Average_DoR_Int_Mid_conf, abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR), abs(database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR_Mid),
                                                           (database_ctl_400_scheme_S2$Average_DoR_Int - database_ctl_400_scheme_S2$Average_DoR)/database_ctl_400_scheme_S2$Average_DoR_Int, abs(database_ctl_400_scheme_S2$Average_DoR_Int_Mid - database_ctl_400_scheme_S2$Average_DoR)/database_ctl_400_scheme_S2$Average_DoR_Int_Mid, (database_ctl_400_scheme_S2$Average_DoR_Int_conf - database_ctl_400_scheme_S2$Average_DoR)/database_ctl_400_scheme_S2$Average_DoR_Int_conf, (database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S2$Average_DoR)/database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf, abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR), abs(database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR_Mid),
                                                           (database_exp_400_scheme_S3$Average_DoR_Int - database_exp_400_scheme_S3$Average_DoR)/database_exp_400_scheme_S3$Average_DoR_Int, abs(database_exp_400_scheme_S3$Average_DoR_Int_Mid - database_exp_400_scheme_S3$Average_DoR)/database_exp_400_scheme_S3$Average_DoR_Int_Mid, (database_exp_400_scheme_S3$Average_DoR_Int_conf - database_exp_400_scheme_S3$Average_DoR)/database_exp_400_scheme_S3$Average_DoR_Int_conf, (database_exp_400_scheme_S3$Average_DoR_Int_Mid_conf - database_exp_400_scheme_S3$Average_DoR)/database_exp_400_scheme_S3$Average_DoR_Int_Mid_conf, abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR), abs(database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR_Mid),
                                                           (database_ctl_400_scheme_S3$Average_DoR_Int - database_ctl_400_scheme_S3$Average_DoR)/database_ctl_400_scheme_S3$Average_DoR_Int, abs(database_ctl_400_scheme_S3$Average_DoR_Int_Mid - database_ctl_400_scheme_S3$Average_DoR)/database_ctl_400_scheme_S3$Average_DoR_Int_Mid, (database_ctl_400_scheme_S3$Average_DoR_Int_conf - database_ctl_400_scheme_S3$Average_DoR)/database_ctl_400_scheme_S3$Average_DoR_Int_conf, (database_ctl_400_scheme_S3$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S3$Average_DoR)/database_ctl_400_scheme_S3$Average_DoR_Int_Mid_conf, abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR), abs(database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR_Mid),
                                                           (database_exp_400_scheme_S4$Average_DoR_Int - database_exp_400_scheme_S4$Average_DoR)/database_exp_400_scheme_S4$Average_DoR_Int, abs(database_exp_400_scheme_S4$Average_DoR_Int_Mid - database_exp_400_scheme_S4$Average_DoR)/database_exp_400_scheme_S4$Average_DoR_Int_Mid, (database_exp_400_scheme_S4$Average_DoR_Int_conf - database_exp_400_scheme_S4$Average_DoR)/database_exp_400_scheme_S4$Average_DoR_Int_conf, (database_exp_400_scheme_S4$Average_DoR_Int_Mid_conf - database_exp_400_scheme_S4$Average_DoR)/database_exp_400_scheme_S4$Average_DoR_Int_Mid_conf, abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR), abs(database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR_Mid),
                                                           (database_ctl_400_scheme_S4$Average_DoR_Int - database_ctl_400_scheme_S4$Average_DoR)/database_ctl_400_scheme_S4$Average_DoR_Int, abs(database_ctl_400_scheme_S4$Average_DoR_Int_Mid - database_ctl_400_scheme_S4$Average_DoR)/database_ctl_400_scheme_S4$Average_DoR_Int_Mid, (database_ctl_400_scheme_S4$Average_DoR_Int_conf - database_ctl_400_scheme_S4$Average_DoR)/database_ctl_400_scheme_S4$Average_DoR_Int_conf, (database_ctl_400_scheme_S4$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S4$Average_DoR)/database_ctl_400_scheme_S4$Average_DoR_Int_Mid_conf, abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR), abs(database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR_Mid),
                                                           (database_exp_400_scheme_S5$Average_DoR_Int - database_exp_400_scheme_S5$Average_DoR)/database_exp_400_scheme_S5$Average_DoR_Int, abs(database_exp_400_scheme_S5$Average_DoR_Int_Mid - database_exp_400_scheme_S5$Average_DoR)/database_exp_400_scheme_S5$Average_DoR_Int_Mid, (database_exp_400_scheme_S5$Average_DoR_Int_conf - database_exp_400_scheme_S5$Average_DoR)/database_exp_400_scheme_S5$Average_DoR_Int_conf, (database_exp_400_scheme_S5$Average_DoR_Int_Mid_conf - database_exp_400_scheme_S5$Average_DoR)/database_exp_400_scheme_S5$Average_DoR_Int_Mid_conf, abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR), abs(database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR_Mid), 
                                                           (database_ctl_400_scheme_S5$Average_DoR_Int - database_ctl_400_scheme_S5$Average_DoR)/database_ctl_400_scheme_S5$Average_DoR_Int, abs(database_ctl_400_scheme_S5$Average_DoR_Int_Mid - database_ctl_400_scheme_S5$Average_DoR)/database_ctl_400_scheme_S5$Average_DoR_Int_Mid, (database_ctl_400_scheme_S5$Average_DoR_Int_conf - database_ctl_400_scheme_S5$Average_DoR)/database_ctl_400_scheme_S5$Average_DoR_Int_conf, (database_ctl_400_scheme_S5$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S5$Average_DoR)/database_ctl_400_scheme_S5$Average_DoR_Int_Mid_conf, abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR), abs(database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR_Mid)),4),
                               DoR_diff = (-1)*round(c((database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_conf), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR), (database_exp_400_scheme_S1$Average_DoR - database_exp_400_scheme_S1$Average_DoR_ORR_Mid),
                                                       (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_conf), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR), (database_ctl_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR_ORR_Mid),
                                                       (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_conf), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR), (database_exp_400_scheme_S2$Average_DoR - database_exp_400_scheme_S2$Average_DoR_ORR_Mid),
                                                       (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_conf), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR), (database_ctl_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR_ORR_Mid),
                                                       (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_conf), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR), (database_exp_400_scheme_S3$Average_DoR - database_exp_400_scheme_S3$Average_DoR_ORR_Mid),
                                                       (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_conf), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR), (database_ctl_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR_ORR_Mid),
                                                       (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_conf), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR), (database_exp_400_scheme_S4$Average_DoR - database_exp_400_scheme_S4$Average_DoR_ORR_Mid),
                                                       (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_conf), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR), (database_ctl_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR_ORR_Mid),
                                                       (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_conf), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_Int_Mid_conf), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR), (database_exp_400_scheme_S5$Average_DoR - database_exp_400_scheme_S5$Average_DoR_ORR_Mid), 
                                                       (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_conf), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_Int_Mid_conf), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR), (database_ctl_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR_ORR_Mid)),4)))

dataset %>%
  filter(Scenario != "Scenario5") -> dataset

dataset %>%
  filter(Censoring != "Interval_ORR") -> dataset
dataset %>%
  filter(Censoring != "Midpoint_ORR") -> dataset


g1 <- ggplot(dataset, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff_square), color = Censoring)) + geom_line(aes(linetype = Trial)) + facet_wrap(~ Scenario,  nrow = 2) + theme_light() + xlim(0,40) + xlab("Tau") + ylab("Relative difference(exact - estimated)")
g2 <- ggplot(dataset, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff_relative), color = Censoring)) + geom_line(aes(linetype = Trial)) + facet_wrap(~ Scenario,  nrow = 2)  + theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0,40) + xlab("Tau") + ylab("Relative bias(estimated - exact)")
g3 <- ggplot(dataset, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff), color = Censoring)) + geom_line(aes(linetype = Trial)) + facet_wrap(~ Scenario,  nrow = 2) + theme_light() + xlim(0,40) + xlab("Tau") + ylab("RMDoR difference(exact - estimated)")

line <- 0.1
g_test <- ggplot(dataset, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff), color = Censoring)) + geom_line(aes(linetype = Trial)) + facet_wrap(~ Scenario,  nrow = 2) + theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0,40) + xlab("Tau") + ylab("RMDoR difference(estimated - exact)")

g_test <- g_test +  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 1), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 2), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 3),colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 4), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 6), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 7), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 8), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 9), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 10), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 11), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 12),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 13), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 14), colour="grey", linewidth = line, linetype = "dashed")  + geom_hline(yintercept = 0)
g_test <- g_test + geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 15), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 16), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 17),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 18),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 19),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 20),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 21),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 22),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 23), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 24),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 25),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 26),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 27),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 28), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 29),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 30),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 31), colour="grey", linewidth = line, linetype = "dashed")+geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 32), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 33), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 34), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 35), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 36), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 37), colour="grey", linewidth = line, linetype = "dashed") + geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 38), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 39), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario1"), aes(xintercept= 40), colour="grey", linewidth = line, linetype = "dashed")


g_test <- g_test +  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 2), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 4),colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 6), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 8), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 10), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 12), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 14), colour="grey", linewidth = line, linetype = "dashed")
g_test <- g_test + geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 16), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 18),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 20),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 22),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 24), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 26),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 28),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 30),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 32), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 34), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 36), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 38), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario2"), aes(xintercept= 40), colour="grey", linewidth = line, linetype = "dashed")

g_test <- g_test +  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 3), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 6), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 9), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 12), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 15), colour="grey", linewidth = line, linetype = "dashed")+ geom_hline(yintercept = 0)
g_test <- g_test +  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 18), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 21),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 24),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 27),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 30), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 33), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 36), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 39), colour="grey", linewidth = line, linetype = "dashed") +  geom_vline(data=filter(dataset, Scenario == "Scenario3"), aes(xintercept= 42), colour="grey", linewidth = line, linetype = "dashed")


g_test <- g_test +  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 1.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 3),colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 4.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 6), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 7.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 9), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 10.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 12),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 13.5), colour="grey", linewidth = line, linetype = "dashed")
g_test <- g_test + geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 15), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 16.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 18),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 21),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 24),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 27),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 30),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 33), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 36), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset, Scenario == "Scenario4"), aes(xintercept= 39), colour="grey", linewidth = line, linetype = "dashed")

g_test <- g_test  + guides(linetype=guide_legend(""))

plot((database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_Int - database_ctl_400_scheme_S1$Average_DoR_Int), ylim = c(-0.05,0.4), xlim = c(0,40), ylab = "Two arm differnce", xlab = "Tau") 
lines((database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_Int_Mid - database_ctl_400_scheme_S1$Average_DoR_Int_Mid), col = 2)
lines((database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_Int_conf - database_ctl_400_scheme_S1$Average_DoR_Int_conf), col = 3)
lines((database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf), col = 4)
lines((database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_ORR - database_ctl_400_scheme_S1$Average_DoR_ORR), col = 5)
lines((database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_ORR_Mid - database_ctl_400_scheme_S1$Average_DoR_ORR_Mid), col = 6)


dataset_CI <- as.data.frame(cbind(Tau = rep(database_exp_400_scheme_S1$tau, 72),
                                  Trial = rep(rep(c("Setting B","Setting A"), each = 450), 4),
                                  Type  = rep(rep(c("Est","CI_low","CI_hi"), each = 50), 24),
                                  Censoring = rep(rep(c("Daily","Daily_CI_low","Daily_CI_hi","E1","E1_CI_low","E1_CI_hi","E2","E2_CI_low","E2_CI_hi"), each = 50), 8),
                                  Scenario = rep(c("Scenario1","Scenario2","Scenario3","Scenario4"), each = 900),
                                  DoR = c(database_ctl_400_scheme_S1$Average_DoR, database_ctl_400_scheme_S1$DoR_lower_daily , database_ctl_400_scheme_S1$DoR_upper_daily,
                                          database_ctl_400_scheme_S1$Average_DoR_Int, database_ctl_400_scheme_S1$DoR_lower_int , database_ctl_400_scheme_S1$DoR_upper_int,
                                          database_ctl_400_scheme_S1$Average_DoR_Int_Mid, database_ctl_400_scheme_S1$DoR_lower_mid , database_ctl_400_scheme_S1$DoR_upper_daily,
                                          
                                          database_exp_400_scheme_S1$Average_DoR, database_exp_400_scheme_S1$DoR_lower_daily , database_exp_400_scheme_S1$DoR_upper_daily,
                                          database_exp_400_scheme_S1$Average_DoR_Int, database_exp_400_scheme_S1$DoR_lower_int , database_exp_400_scheme_S1$DoR_upper_int,
                                          database_exp_400_scheme_S1$Average_DoR_Int_Mid, database_exp_400_scheme_S1$DoR_lower_mid , database_exp_400_scheme_S1$DoR_upper_daily,
                                          
                                          database_ctl_400_scheme_S2$Average_DoR, database_ctl_400_scheme_S2$DoR_lower_daily , database_ctl_400_scheme_S2$DoR_upper_daily,
                                          database_ctl_400_scheme_S2$Average_DoR_Int, database_ctl_400_scheme_S2$DoR_lower_int , database_ctl_400_scheme_S2$DoR_upper_int,
                                          database_ctl_400_scheme_S2$Average_DoR_Int_Mid, database_ctl_400_scheme_S2$DoR_lower_mid , database_ctl_400_scheme_S2$DoR_upper_daily,
                                          
                                          database_exp_400_scheme_S2$Average_DoR, database_exp_400_scheme_S2$DoR_lower_daily , database_exp_400_scheme_S2$DoR_upper_daily,
                                          database_exp_400_scheme_S2$Average_DoR_Int, database_exp_400_scheme_S2$DoR_lower_int , database_exp_400_scheme_S2$DoR_upper_int,
                                          database_exp_400_scheme_S2$Average_DoR_Int_Mid, database_exp_400_scheme_S2$DoR_lower_mid , database_exp_400_scheme_S2$DoR_upper_daily,
                                          
                                          database_ctl_400_scheme_S3$Average_DoR, database_ctl_400_scheme_S3$DoR_lower_daily , database_ctl_400_scheme_S3$DoR_upper_daily,
                                          database_ctl_400_scheme_S3$Average_DoR_Int, database_ctl_400_scheme_S3$DoR_lower_int , database_ctl_400_scheme_S3$DoR_upper_int,
                                          database_ctl_400_scheme_S3$Average_DoR_Int_Mid, database_ctl_400_scheme_S3$DoR_lower_mid , database_ctl_400_scheme_S3$DoR_upper_daily,
                                          
                                          database_exp_400_scheme_S3$Average_DoR, database_exp_400_scheme_S3$DoR_lower_daily , database_exp_400_scheme_S3$DoR_upper_daily,
                                          database_exp_400_scheme_S3$Average_DoR_Int, database_exp_400_scheme_S3$DoR_lower_int , database_exp_400_scheme_S3$DoR_upper_int,
                                          database_exp_400_scheme_S3$Average_DoR_Int_Mid, database_exp_400_scheme_S3$DoR_lower_mid , database_exp_400_scheme_S3$DoR_upper_daily,
                                          
                                          database_ctl_400_scheme_S4$Average_DoR, database_ctl_400_scheme_S4$DoR_lower_daily , database_ctl_400_scheme_S4$DoR_upper_daily,
                                          database_ctl_400_scheme_S4$Average_DoR_Int, database_ctl_400_scheme_S4$DoR_lower_int , database_ctl_400_scheme_S4$DoR_upper_int,
                                          database_ctl_400_scheme_S4$Average_DoR_Int_Mid, database_ctl_400_scheme_S4$DoR_lower_mid , database_ctl_400_scheme_S4$DoR_upper_daily,
                                          
                                          database_exp_400_scheme_S4$Average_DoR, database_exp_400_scheme_S4$DoR_lower_daily , database_exp_400_scheme_S4$DoR_upper_daily,
                                          database_exp_400_scheme_S4$Average_DoR_Int, database_exp_400_scheme_S4$DoR_lower_int , database_exp_400_scheme_S4$DoR_upper_int,
                                          database_exp_400_scheme_S4$Average_DoR_Int_Mid, database_exp_400_scheme_S4$DoR_lower_mid , database_exp_400_scheme_S4$DoR_upper_daily)))

my_palette <- c(brewer.pal(4, "Set1")[rep(1:4, each = 3)], brewer.pal(3, "Set2"))

g_CI <- ggplot(dataset_CI, aes(x = as.numeric(Tau), y = as.numeric(DoR), color = Censoring)) + geom_line(aes(linetype = Type)) + facet_grid(Trial~Scenario) + theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0,40) + xlab("Tau") + ylab("RMDoR") + scale_linetype_manual(values=c(1,1,2)) + scale_size_manual(values =c(4,4,1)) +  scale_color_manual(values=my_palette)

# Two-arm analysis

diff_1_int          <- (database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_Int - database_ctl_400_scheme_S1$Average_DoR_Int)
diff_1_int_mid      <- (database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_Int_Mid - database_ctl_400_scheme_S1$Average_DoR_Int_Mid)
diff_1_int_conf     <- (database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_Int_conf - database_ctl_400_scheme_S1$Average_DoR_Int_conf)
diff_1_int_conf_mid <- (database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf)
diff_1_int_ORR      <- (database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_ORR - database_ctl_400_scheme_S1$Average_DoR_ORR)
diff_1_int_ORR_mid  <- (database_exp_400_scheme_S1$Average_DoR - database_ctl_400_scheme_S1$Average_DoR) - (database_exp_400_scheme_S1$Average_DoR_ORR_Mid - database_ctl_400_scheme_S1$Average_DoR_ORR_Mid)

diff_1_int_rela          <- ((-1)*diff_1_int)/(database_exp_400_scheme_S1$Average_DoR_Int - database_ctl_400_scheme_S1$Average_DoR_Int)
diff_1_int_mid_rela      <- ((-1)*diff_1_int_mid)/(database_exp_400_scheme_S1$Average_DoR_Int_Mid - database_ctl_400_scheme_S1$Average_DoR_Int_Mid)
diff_1_int_conf_rela     <- ((-1)*diff_1_int_conf)/(database_exp_400_scheme_S1$Average_DoR_Int_conf - database_ctl_400_scheme_S1$Average_DoR_Int_conf)
diff_1_int_conf_mid_rela <- ((-1)*diff_1_int_conf_mid)/(database_exp_400_scheme_S1$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S1$Average_DoR_Int_Mid_conf)
diff_1_int_ORR_rela      <- ((-1)*diff_1_int_ORR)/(database_exp_400_scheme_S1$Average_DoR_ORR - database_ctl_400_scheme_S1$Average_DoR_ORR)
diff_1_int_ORR_mid_rela  <- ((-1)*diff_1_int_ORR_mid)/(database_exp_400_scheme_S1$Average_DoR_ORR_Mid - database_ctl_400_scheme_S1$Average_DoR_ORR_Mid)

diff_2_int          <- (database_exp_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR) - (database_exp_400_scheme_S2$Average_DoR_Int - database_ctl_400_scheme_S2$Average_DoR_Int)
diff_2_int_mid      <- (database_exp_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR) - (database_exp_400_scheme_S2$Average_DoR_Int_Mid - database_ctl_400_scheme_S2$Average_DoR_Int_Mid)
diff_2_int_conf     <- (database_exp_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR) - (database_exp_400_scheme_S2$Average_DoR_Int_conf - database_ctl_400_scheme_S2$Average_DoR_Int_conf)
diff_2_int_conf_mid <- (database_exp_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR) - (database_exp_400_scheme_S2$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf)
diff_2_int_ORR      <- (database_exp_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR) - (database_exp_400_scheme_S2$Average_DoR_ORR - database_ctl_400_scheme_S2$Average_DoR_ORR)
diff_2_int_ORR_mid  <- (database_exp_400_scheme_S2$Average_DoR - database_ctl_400_scheme_S2$Average_DoR) - (database_exp_400_scheme_S2$Average_DoR_ORR_Mid - database_ctl_400_scheme_S2$Average_DoR_ORR_Mid)

diff_2_int_rela          <- ((-1)*diff_2_int)/(database_exp_400_scheme_S2$Average_DoR_Int - database_ctl_400_scheme_S2$Average_DoR_Int)
diff_2_int_mid_rela      <- ((-1)*diff_2_int_mid)/(database_exp_400_scheme_S2$Average_DoR_Int_Mid - database_ctl_400_scheme_S2$Average_DoR_Int_Mid)
diff_2_int_conf_rela     <- ((-1)*diff_2_int_conf)/(database_exp_400_scheme_S2$Average_DoR_Int_conf - database_ctl_400_scheme_S2$Average_DoR_Int_conf)
diff_2_int_conf_mid_rela <- ((-1)*diff_2_int_conf_mid)/(database_exp_400_scheme_S2$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S2$Average_DoR_Int_Mid_conf)
diff_2_int_ORR_rela      <- ((-1)*diff_2_int_ORR)/(database_exp_400_scheme_S2$Average_DoR_ORR - database_ctl_400_scheme_S2$Average_DoR_ORR)
diff_2_int_ORR_mid_rela  <- ((-1)*diff_2_int_ORR_mid)/(database_exp_400_scheme_S2$Average_DoR_ORR_Mid - database_ctl_400_scheme_S2$Average_DoR_ORR_Mid)

diff_3_int          <- (database_exp_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR) - (database_exp_400_scheme_S3$Average_DoR_Int - database_ctl_400_scheme_S3$Average_DoR_Int)
diff_3_int_mid      <- (database_exp_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR) - (database_exp_400_scheme_S3$Average_DoR_Int_Mid - database_ctl_400_scheme_S3$Average_DoR_Int_Mid)
diff_3_int_conf     <- (database_exp_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR) - (database_exp_400_scheme_S3$Average_DoR_Int_conf - database_ctl_400_scheme_S3$Average_DoR_Int_conf)
diff_3_int_conf_mid <- (database_exp_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR) - (database_exp_400_scheme_S3$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S3$Average_DoR_Int_Mid_conf)
diff_3_int_ORR      <- (database_exp_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR) - (database_exp_400_scheme_S3$Average_DoR_ORR - database_ctl_400_scheme_S3$Average_DoR_ORR)
diff_3_int_ORR_mid  <- (database_exp_400_scheme_S3$Average_DoR - database_ctl_400_scheme_S3$Average_DoR) - (database_exp_400_scheme_S3$Average_DoR_ORR_Mid - database_ctl_400_scheme_S3$Average_DoR_ORR_Mid)

diff_3_int_rela          <- ((-1)*diff_3_int)/(database_exp_400_scheme_S3$Average_DoR_Int - database_ctl_400_scheme_S3$Average_DoR_Int)
diff_3_int_mid_rela      <- ((-1)*diff_3_int_mid)/(database_exp_400_scheme_S3$Average_DoR_Int_Mid - database_ctl_400_scheme_S3$Average_DoR_Int_Mid)
diff_3_int_conf_rela     <- ((-1)*diff_3_int_conf)/(database_exp_400_scheme_S3$Average_DoR_Int_conf - database_ctl_400_scheme_S3$Average_DoR_Int_conf)
diff_3_int_conf_mid_rela <- ((-1)*diff_3_int_conf_mid)/(database_exp_400_scheme_S3$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S3$Average_DoR_Int_Mid_conf)
diff_3_int_ORR_rela      <- ((-1)*diff_3_int_ORR)/(database_exp_400_scheme_S3$Average_DoR_ORR - database_ctl_400_scheme_S3$Average_DoR_ORR)
diff_3_int_ORR_mid_rela  <- ((-1)*diff_3_int_ORR_mid)/(database_exp_400_scheme_S3$Average_DoR_ORR_Mid - database_ctl_400_scheme_S3$Average_DoR_ORR_Mid)

diff_4_int          <- (database_exp_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR) - (database_exp_400_scheme_S4$Average_DoR_Int - database_ctl_400_scheme_S4$Average_DoR_Int)
diff_4_int_mid      <- (database_exp_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR) - (database_exp_400_scheme_S4$Average_DoR_Int_Mid - database_ctl_400_scheme_S4$Average_DoR_Int_Mid)
diff_4_int_conf     <- (database_exp_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR) - (database_exp_400_scheme_S4$Average_DoR_Int_conf - database_ctl_400_scheme_S4$Average_DoR_Int_conf)
diff_4_int_conf_mid <- (database_exp_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR) - (database_exp_400_scheme_S4$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S4$Average_DoR_Int_Mid_conf)
diff_4_int_ORR      <- (database_exp_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR) - (database_exp_400_scheme_S4$Average_DoR_ORR - database_ctl_400_scheme_S4$Average_DoR_ORR)
diff_4_int_ORR_mid  <- (database_exp_400_scheme_S4$Average_DoR - database_ctl_400_scheme_S4$Average_DoR) - (database_exp_400_scheme_S4$Average_DoR_ORR_Mid - database_ctl_400_scheme_S4$Average_DoR_ORR_Mid)

diff_4_int_rela          <- ((-1)*diff_4_int)/(database_exp_400_scheme_S4$Average_DoR_Int - database_ctl_400_scheme_S4$Average_DoR_Int)
diff_4_int_mid_rela      <- ((-1)*diff_4_int_mid)/(database_exp_400_scheme_S4$Average_DoR_Int_Mid - database_ctl_400_scheme_S4$Average_DoR_Int_Mid)
diff_4_int_conf_rela     <- ((-1)*diff_4_int_conf)/(database_exp_400_scheme_S4$Average_DoR_Int_conf - database_ctl_400_scheme_S4$Average_DoR_Int_conf)
diff_4_int_conf_mid_rela <- ((-1)*diff_4_int_conf_mid)/(database_exp_400_scheme_S4$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S4$Average_DoR_Int_Mid_conf)
diff_4_int_ORR_rela      <- ((-1)*diff_4_int_ORR)/(database_exp_400_scheme_S4$Average_DoR_ORR - database_ctl_400_scheme_S4$Average_DoR_ORR)
diff_4_int_ORR_mid_rela  <- ((-1)*diff_4_int_ORR_mid)/(database_exp_400_scheme_S4$Average_DoR_ORR_Mid - database_ctl_400_scheme_S4$Average_DoR_ORR_Mid)

diff_5_int          <- (database_exp_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR) - (database_exp_400_scheme_S5$Average_DoR_Int - database_ctl_400_scheme_S5$Average_DoR_Int)
diff_5_int_mid      <- (database_exp_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR) - (database_exp_400_scheme_S5$Average_DoR_Int_Mid - database_ctl_400_scheme_S5$Average_DoR_Int_Mid)
diff_5_int_conf     <- (database_exp_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR) - (database_exp_400_scheme_S5$Average_DoR_Int_conf - database_ctl_400_scheme_S5$Average_DoR_Int_conf)
diff_5_int_conf_mid <- (database_exp_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR) - (database_exp_400_scheme_S5$Average_DoR_Int_Mid_conf - database_ctl_400_scheme_S5$Average_DoR_Int_Mid_conf)
diff_5_int_ORR      <- (database_exp_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR) - (database_exp_400_scheme_S5$Average_DoR_ORR - database_ctl_400_scheme_S5$Average_DoR_ORR)
diff_5_int_ORR_mid  <- (database_exp_400_scheme_S5$Average_DoR - database_ctl_400_scheme_S5$Average_DoR) - (database_exp_400_scheme_S5$Average_DoR_ORR_Mid - database_ctl_400_scheme_S5$Average_DoR_ORR_Mid)










#dataset2 <- as.data.frame(cbind(Tau = rep(database_exp_400_scheme_S1$tau, 30),
#                                Censoring = rep(rep(c("Interval","Midpoint","Interval_check","Midpoint_ckeck","Interval_ORR","Midpoint_ORR"), each = 50), 5),
#                                Scenario = rep(c("Scenario1","Scenario2","Scenario3","Scenario4","Scenario5"), each = 300),
#                                DoR_diff_square = round(c(diff_1_int^2, diff_1_int_mid^2, diff_1_int_conf^2, diff_1_int_conf_mid^2, diff_1_int_ORR^2, diff_1_int_ORR_mid^2,
#                                                          diff_2_int^2, diff_2_int_mid^2, diff_2_int_conf^2, diff_2_int_conf_mid^2, diff_2_int_ORR^2, diff_2_int_ORR_mid^2,
#                                                          diff_3_int^2, diff_3_int_mid^2, diff_3_int_conf^2, diff_3_int_conf_mid^2, diff_3_int_ORR^2, diff_3_int_ORR_mid^2,
#                                                          diff_4_int^2, diff_4_int_mid^2, diff_4_int_conf^2, diff_4_int_conf_mid^2, diff_4_int_ORR^2, diff_4_int_ORR_mid^2,
#                                                          diff_5_int^2, diff_5_int_mid^2, diff_5_int_conf^2, diff_5_int_conf_mid^2, diff_5_int_ORR^2, diff_5_int_ORR_mid^2), 4),
#                                DoR_diff_abs = round(c(abs(diff_1_int), abs(diff_1_int_mid), abs(diff_1_int_conf), abs(diff_1_int_conf_mid), abs(diff_1_int_ORR), abs(diff_1_int_ORR_mid),
#                                                       abs(diff_2_int), abs(diff_2_int_mid), abs(diff_2_int_conf), abs(diff_2_int_conf_mid), abs(diff_2_int_ORR), abs(diff_2_int_ORR_mid),
#                                                       abs(diff_3_int), abs(diff_3_int_mid), abs(diff_3_int_conf), abs(diff_3_int_conf_mid), abs(diff_3_int_ORR), abs(diff_3_int_ORR_mid),
#                                                       abs(diff_4_int), abs(diff_4_int_mid), abs(diff_4_int_conf), abs(diff_4_int_conf_mid), abs(diff_4_int_ORR), abs(diff_4_int_ORR_mid),
#                                                       abs(diff_5_int), abs(diff_5_int_mid), abs(diff_5_int_conf), abs(diff_5_int_conf_mid), abs(diff_5_int_ORR), abs(diff_5_int_ORR_mid)), 4),
#                                
#                                DoR_diff = (-1)*round(c((diff_1_int), (diff_1_int_mid), (diff_1_int_conf), (diff_1_int_conf_mid), (diff_1_int_ORR), (diff_1_int_ORR_mid),
#                                                   (diff_2_int), (diff_2_int_mid), (diff_2_int_conf), (diff_2_int_conf_mid), (diff_2_int_ORR), (diff_2_int_ORR_mid),
#                                                   (diff_3_int), (diff_3_int_mid), (diff_3_int_conf), (diff_3_int_conf_mid), (diff_3_int_ORR), (diff_3_int_ORR_mid),
#                                                   (diff_4_int), (diff_4_int_mid), (diff_4_int_conf), (diff_4_int_conf_mid), (diff_4_int_ORR), (diff_4_int_ORR_mid),
#                                                   (diff_5_int), (diff_5_int_mid), (diff_5_int_conf), (diff_5_int_conf_mid), (diff_5_int_ORR), (diff_5_int_ORR_mid)), 4)))

dataset2 <- as.data.frame(cbind(Tau = rep(database_exp_400_scheme_S1$tau, 30),
                                Censoring = rep(rep(c("E1","E2","D1","D2","Interval_ORR","Midpoint_ORR"), each = 50), 5),
                                Scenario = rep(c("Scenario1","Scenario2","Scenario3","Scenario4","Scenario5"), each = 300),
                                DoR_diff_square = round(c(diff_1_int^2, diff_1_int_mid^2, diff_1_int_conf^2, diff_1_int_conf_mid^2, diff_1_int_ORR^2, diff_1_int_ORR_mid^2,
                                                          diff_2_int^2, diff_2_int_mid^2, diff_2_int_conf^2, diff_2_int_conf_mid^2, diff_2_int_ORR^2, diff_2_int_ORR_mid^2,
                                                          diff_3_int^2, diff_3_int_mid^2, diff_3_int_conf^2, diff_3_int_conf_mid^2, diff_3_int_ORR^2, diff_3_int_ORR_mid^2,
                                                          diff_4_int^2, diff_4_int_mid^2, diff_4_int_conf^2, diff_4_int_conf_mid^2, diff_4_int_ORR^2, diff_4_int_ORR_mid^2,
                                                          diff_5_int^2, diff_5_int_mid^2, diff_5_int_conf^2, diff_5_int_conf_mid^2, diff_5_int_ORR^2, diff_5_int_ORR_mid^2), 4),
                                DoR_diff_relative = round(c((diff_1_int_rela), (diff_1_int_mid_rela), (diff_1_int_conf_rela), (diff_1_int_conf_mid_rela), (diff_1_int_ORR_rela), (diff_1_int_ORR_mid_rela),
                                                            (diff_2_int_rela), (diff_2_int_mid_rela), (diff_2_int_conf_rela), (diff_2_int_conf_mid_rela), (diff_2_int_ORR_rela), (diff_2_int_ORR_mid_rela),
                                                            (diff_3_int_rela), (diff_3_int_mid_rela), (diff_3_int_conf_rela), (diff_3_int_conf_mid_rela), (diff_3_int_ORR_rela), (diff_3_int_ORR_mid_rela),
                                                            (diff_4_int_rela), (diff_4_int_mid_rela), (diff_4_int_conf_rela), (diff_4_int_conf_mid_rela), (diff_4_int_ORR_rela), (diff_4_int_ORR_mid_rela),
                                                            (diff_5_int), (diff_5_int_mid), (diff_5_int_conf), (diff_5_int_conf_mid), (diff_5_int_ORR), (diff_5_int_ORR_mid)), 4),
                                
                                DoR_diff = (-1)*round(c((diff_1_int), (diff_1_int_mid), (diff_1_int_conf), (diff_1_int_conf_mid), (diff_1_int_ORR), (diff_1_int_ORR_mid),
                                                        (diff_2_int), (diff_2_int_mid), (diff_2_int_conf), (diff_2_int_conf_mid), (diff_2_int_ORR), (diff_2_int_ORR_mid),
                                                        (diff_3_int), (diff_3_int_mid), (diff_3_int_conf), (diff_3_int_conf_mid), (diff_3_int_ORR), (diff_3_int_ORR_mid),
                                                        (diff_4_int), (diff_4_int_mid), (diff_4_int_conf), (diff_4_int_conf_mid), (diff_4_int_ORR), (diff_4_int_ORR_mid),
                                                        (diff_5_int), (diff_5_int_mid), (diff_5_int_conf), (diff_5_int_conf_mid), (diff_5_int_ORR), (diff_5_int_ORR_mid)), 4)))

dataset2 %>%
  filter(Scenario != "Scenario5") -> dataset2

dataset2 %>%
  filter(Censoring != "Interval_ORR") -> dataset2
dataset2 %>%
  filter(Censoring != "Midpoint_ORR") -> dataset2


g1_two_arm <- ggplot(dataset2, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff_square), color = Censoring)) + geom_line()  + facet_wrap(~ Scenario,  nrow = 2) + theme_classic() + xlim(0,40)  + xlab("Tau") + ylab("Relative difference two arms(exact - estimated)")
g2_two_arm <- ggplot(dataset2, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff_relative), color = Censoring)) + geom_line()  + facet_wrap(~ Scenario,  nrow = 2)  + theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ xlim(0,40)+ xlab("Tau") + ylab("Relative difference two arms(estimated - exact)")
g3_two_arm <- ggplot(dataset2, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff), color = Censoring)) + geom_line()  + facet_wrap(~ Scenario,  nrow = 2) + theme_light() + xlim(0,40) + xlab("Tau") + ylab("RMDoR difference(exact - estimated)")

g_test2 <- ggplot(dataset2, aes(x = as.numeric(Tau), y = as.numeric(DoR_diff), color = Censoring)) + geom_line() + facet_wrap(~ Scenario,  nrow = 2) + theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ xlim(0,40) + xlab("Tau") + ylab("RMDoR difference two-arms(estimated - exact)")

g_test2 <- g_test2 +  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 1), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 2), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 3),colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 4), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 6), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 7), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 8), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 9), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 10), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 11), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 12),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 13), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 14), colour="grey", linewidth = line, linetype = "dashed")  + geom_hline(yintercept = 0)
g_test2 <- g_test2 + geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 15), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 16), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 17),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 18),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 19),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 20),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 21),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 22),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 23), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 24),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 25),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 26),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 27),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 28), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 29),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 30),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 31), colour="grey", linewidth = line, linetype = "dashed")+geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 32), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 33), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 34), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 35), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 36), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 37), colour="grey", linewidth = line, linetype = "dashed") + geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 38), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 39), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario1"), aes(xintercept= 40), colour="grey", linewidth = line, linetype = "dashed")


g_test2 <- g_test2 +  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 1), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 3),colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 7), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 9), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 11), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 13), colour="grey", linewidth = line, linetype = "dashed")
g_test2 <- g_test2 + geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 15), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 17),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 19),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 21),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 23), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 25),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 27),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 29),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 31), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 33), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 35), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 37), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario2"), aes(xintercept= 39), colour="grey", linewidth = line, linetype = "dashed")

g_test2 <- g_test2 +  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 1), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 4), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 7), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 10), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 13), colour="grey", linewidth = line, linetype = "dashed")+ geom_hline(yintercept = 0)
g_test2 <- g_test2 +  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 16), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 19),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 22),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 25),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 28), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 31), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 34), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 37), colour="grey", linewidth = line, linetype = "dashed") +  geom_vline(data=filter(dataset2, Scenario == "Scenario3"), aes(xintercept= 40), colour="grey", linewidth = line, linetype = "dashed")


g_test2 <- g_test2 +  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 1.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 3),colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 4.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 6), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 7.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 9), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 10.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 12),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 13.5), colour="grey", linewidth = line, linetype = "dashed")
g_test2 <- g_test2 + geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 15), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 16.5), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 18),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 21),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 24),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 27),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 30),  colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 33), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 36), colour="grey", linewidth = line, linetype = "dashed")+  geom_vline(data=filter(dataset2, Scenario == "Scenario4"), aes(xintercept= 39), colour="grey", linewidth = line, linetype = "dashed")



dataset_CI_two_arms <- as.data.frame(cbind(Tau = rep(database_exp_400_scheme_S1$tau, 36),
                                  Type  = rep(rep(c("Est","CI_low","CI_hi"), each = 50), 12),
                                  Censoring = c(rep(rep(c("Daily","Daily_CI_low","Daily_CI_hi"), each = 50),4),rep(rep(c("E1","E1_CI_low","E1_CI_hi"), each = 50),4),rep(rep(c("E2","E2_CI_low","E2_CI_hi"), each = 50),4)),
                                  Scenario = rep(rep(c("Scenario1","Scenario2","Scenario3","Scenario4"), each = 150),3),
                                  DoR = c(as.numeric(unlist(lapply(dor_1_diff, mean))), dor_1_diff_low, dor_1_diff_hi,
                                          as.numeric(unlist(lapply(dor_2_diff, mean))), dor_2_diff_low, dor_2_diff_hi,
                                          as.numeric(unlist(lapply(dor_3_diff, mean))), dor_3_diff_low, dor_3_diff_hi,
                                          as.numeric(unlist(lapply(dor_4_diff, mean))), dor_4_diff_low, dor_4_diff_hi,
                                          
                                          as.numeric(unlist(lapply(dor_1_right_noconf_diff, mean))), dor_1_right_noconf_diff_low, dor_1_right_noconf_diff_hi,
                                          as.numeric(unlist(lapply(dor_2_right_noconf_diff, mean))), dor_2_right_noconf_diff_low, dor_2_right_noconf_diff_hi,
                                          as.numeric(unlist(lapply(dor_3_right_noconf_diff, mean))), dor_3_right_noconf_diff_low, dor_3_right_noconf_diff_hi,
                                          as.numeric(unlist(lapply(dor_4_right_noconf_diff, mean))), dor_4_right_noconf_diff_low, dor_4_right_noconf_diff_hi,
                                          
                                          as.numeric(unlist(lapply(dor_1_mid_noconf_diff, mean))), dor_1_mid_noconf_diff_low, dor_1_mid_noconf_diff_hi,
                                          as.numeric(unlist(lapply(dor_2_mid_noconf_diff, mean))), dor_2_mid_noconf_diff_low, dor_2_mid_noconf_diff_hi,
                                          as.numeric(unlist(lapply(dor_3_mid_noconf_diff, mean))), dor_3_mid_noconf_diff_low, dor_3_mid_noconf_diff_hi,
                                          as.numeric(unlist(lapply(dor_4_mid_noconf_diff, mean))), dor_4_mid_noconf_diff_low, dor_4_mid_noconf_diff_hi)))

my_palette2 <- c(brewer.pal(4, "Set1")[rep(1:4, each = 3)], brewer.pal(3, "Set2"))

g_CI_two_arms <- ggplot(dataset_CI_two_arms, aes(x = as.numeric(Tau), y = as.numeric(DoR), color = Censoring)) + geom_line(aes(linetype = Type)) + facet_wrap(~ Scenario,  nrow = 2) + theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlim(0,40) + xlab("Tau") + ylab("RMDoR difference")+ scale_linetype_manual(values=c(1,1,2)) + scale_size_manual(values =c(4,4,1)) +   scale_color_manual(values=my_palette2)

scale_color_manual(values=c("red","red","red","blue","blue","blue","yellow","yellow","yellow"))





dataset_CI_two_arms_raio <- as.data.frame(cbind(Tau = rep(database_exp_400_scheme_S1$tau, 36),
                                           Type  = rep(rep(c("Est","CI_low","CI_hi"), each = 50), 12),
                                           Censoring = c(rep(rep(c("Daily","Daily_CI_low","Daily_CI_hi"), each = 50),4),rep(rep(c("E1","E1_CI_low","E1_CI_hi"), each = 50),4),rep(rep(c("E2","E2_CI_low","E2_CI_hi"), each = 50),4)),
                                           Scenario = rep(rep(c("Scenario1","Scenario2","Scenario3","Scenario4"), each = 150),3),
                                           DoR = c(as.numeric(unlist(lapply(dor_1_ratio, mean))), dor_1_ratio_low, dor_1_ratio_hi,
                                                   as.numeric(unlist(lapply(dor_2_ratio, mean))), dor_2_ratio_low, dor_2_ratio_hi,
                                                   as.numeric(unlist(lapply(dor_3_ratio, mean))), dor_3_ratio_low, dor_3_ratio_hi,
                                                   as.numeric(unlist(lapply(dor_4_ratio, mean))), dor_4_ratio_low, dor_4_ratio_hi,
                                                   
                                                   as.numeric(unlist(lapply(dor_1_right_noconf_ratio, mean))), dor_1_right_noconf_ratio_low, dor_1_right_noconf_ratio_hi,
                                                   as.numeric(unlist(lapply(dor_2_right_noconf_ratio, mean))), dor_2_right_noconf_ratio_low, dor_2_right_noconf_ratio_hi,
                                                   as.numeric(unlist(lapply(dor_3_right_noconf_ratio, mean))), dor_3_right_noconf_ratio_low, dor_3_right_noconf_ratio_hi,
                                                   as.numeric(unlist(lapply(dor_4_right_noconf_ratio, mean))), dor_4_right_noconf_ratio_low, dor_4_right_noconf_ratio_hi,
                                                   
                                                   as.numeric(unlist(lapply(dor_1_mid_noconf_ratio, mean))), dor_1_mid_noconf_ratio_low, dor_1_mid_noconf_ratio_hi,
                                                   as.numeric(unlist(lapply(dor_2_mid_noconf_ratio, mean))), dor_2_mid_noconf_ratio_low, dor_2_mid_noconf_ratio_hi,
                                                   as.numeric(unlist(lapply(dor_3_mid_noconf_ratio, mean))), dor_3_mid_noconf_ratio_low, dor_3_mid_noconf_ratio_hi,
                                                   as.numeric(unlist(lapply(dor_4_mid_noconf_ratio, mean))), dor_4_mid_noconf_ratio_low, dor_4_mid_noconf_ratio_hi)))

my_palette2 <- c(brewer.pal(4, "Set1")[rep(1:4, each = 3)], brewer.pal(3, "Set2"))

g_CI_two_arms_ratio <- ggplot(dataset_CI_two_arms_raio, aes(x = as.numeric(Tau), y = as.numeric(DoR), color = Censoring)) + geom_line(aes(linetype = Type)) + facet_wrap(~ Scenario,  nrow = 2) + theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Tau") + ylab("RMDoR ratio")+ scale_linetype_manual(values=c(1,1,2)) + scale_size_manual(values =c(4,4,1)) +   scale_color_manual(values=my_palette2) +   scale_x_continuous(breaks = c(5,15,25,35,45), limits = c(5,40))

scale_color_manual(values=c("red","red","red","blue","blue","blue","yellow","yellow","yellow"))
























# Check if the  Censoring or not affects the results 

set.seed(2405)
S1 = 0:100
S2 = seq(0,100,2)
S3 = seq(0,100,3)
S4 = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72)
S5 = c(0, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 21, 27, 33, 39, 45, 51, 57, 63, 69)

tau.val <- 40
iter <- 1000

# Strategy S1
#############

test <- as.data.frame(matrix(NA, tau.val, 14))
dor_exp_1test <- NULL
dor_exp_1_no_conftest <- NULL
dor_exp_1_mid_no_conftest <- NULL
dor_exp_1_conftest <- NULL
dor_exp_1_mid_conftest <- NULL
dor_exp_1_ORRtest <- NULL
dor_exp_1_mid_ORRtest <- NULL
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S3)
  test[z,] <- c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_exp_1test[[tau]] <- k$DoR
  dor_exp_1_no_conftest[[tau]]     <- k$DoRInt_no_conf
  dor_exp_1_mid_no_conftest[[tau]] <- k$DoRIntM_no_conf
  dor_exp_1_conftest[[tau]]        <- k$DoRInt_conf
  dor_exp_1_mid_conftest[[tau]]    <- k$DoRIntM_conf
  dor_exp_1_ORRtest[[tau]]         <- k$DoRInt_ORR
  dor_exp_1_mid_ORRtest[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(test) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")


set.seed(2405)
S1 = 0:100
S2 = seq(0,100,2)
S3 = seq(0,100,3)
S4 = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72)
S5 = c(0, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 21, 27, 33, 39, 45, 51, 57, 63, 69)

tau.val <- 40
iter <- 1000

# Strategy S1
#############

test2 <- as.data.frame(matrix(NA, tau.val, 14))
dor_exp_1test2 <- NULL
dor_exp_1_no_conftest2 <- NULL
dor_exp_1_mid_no_conftest2 <- NULL
dor_exp_1_conftest2 <- NULL
dor_exp_1_mid_conftest2 <- NULL
dor_exp_1_ORRtest2 <- NULL
dor_exp_1_mid_ORRtest2 <- NULL
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc2(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S3)
  test2[z,] <- c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_exp_1test2[[tau]] <- k$DoR
  dor_exp_1_no_conftest2[[tau]]     <- k$DoRInt_no_conf
  dor_exp_1_mid_no_conftest2[[tau]] <- k$DoRIntM_no_conf
  dor_exp_1_conftest2[[tau]]        <- k$DoRInt_conf
  dor_exp_1_mid_conftest2[[tau]]    <- k$DoRIntM_conf
  dor_exp_1_ORRtest2[[tau]]         <- k$DoRInt_ORR
  dor_exp_1_mid_ORRtest2[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(test2) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")





set.seed(2405)
S1 = 0:100
S2 = seq(0,100,2)
S3 = seq(0,100,3)
S4 = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72)
S5 = c(0, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 21, 27, 33, 39, 45, 51, 57, 63, 69)

tau.val <- 40
iter <- 1000

# Strategy S1
#############

test_s2 <- as.data.frame(matrix(NA, tau.val, 14))
dor_exp_1test_s2 <- NULL
dor_exp_1_no_conftest_s2 <- NULL
dor_exp_1_mid_no_conftest_s2 <- NULL
dor_exp_1_conftest_s2 <- NULL
dor_exp_1_mid_conftest_s2 <- NULL
dor_exp_1_ORRtest_s2 <- NULL
dor_exp_1_mid_ORRtest_s2 <- NULL
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S2)
  test_s2[z,] <- c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_exp_1test_s2[[tau]] <- k$DoR
  dor_exp_1_no_conftest_s2[[tau]]     <- k$DoRInt_no_conf
  dor_exp_1_mid_no_conftest_s2[[tau]] <- k$DoRIntM_no_conf
  dor_exp_1_conftest_s2[[tau]]        <- k$DoRInt_conf
  dor_exp_1_mid_conftest_s2[[tau]]    <- k$DoRIntM_conf
  dor_exp_1_ORRtest_s2[[tau]]         <- k$DoRInt_ORR
  dor_exp_1_mid_ORRtest_s2[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(test_s2) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")


set.seed(2405)
S1 = 0:100
S2 = seq(0,100,2)
S3 = seq(0,100,3)
S4 = c(0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72)
S5 = c(0, 1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 21, 27, 33, 39, 45, 51, 57, 63, 69)

tau.val <- 40
iter <- 1000

# Strategy S1
#############

test2_s2 <- as.data.frame(matrix(NA, tau.val, 14))
dor_exp_1test2_s2 <- NULL
dor_exp_1_no_conftest2_s2 <- NULL
dor_exp_1_mid_no_conftest2_s2 <- NULL
dor_exp_1_conftest2_s2 <- NULL
dor_exp_1_mid_conftest2_s2 <- NULL
dor_exp_1_ORRtest2_s2 <- NULL
dor_exp_1_mid_ORRtest2_s2 <- NULL
z <- 1
for (tau in 1:tau.val) 
{
  k <- sim_func_mDoR_calc2(shapeD = 3, scaleD = 6, shapeP = 0.760043, scaleP = 21.5417, shapeR = 1.5, scaleR = 6, cL = 5 , cR = 70, tau = tau, n = 400, iter = iter, scan_scenario = S2)
  test2_s2[z,] <- c(k$Average_DoR, k$DoR_CI_low, k$DoR_CI_upper, k$Average_DoR_Int_no_conf, k$DoR_Int_no_conf_CI_low, k$DoR_Int_no_conf_CI_upper, k$Average_DoR_IntMid_no_conf, k$DoR_IntMid_no_conf_CI_low, k$DoR_IntMid_no_conf_CI_upper, k$Average_DoR_Int_conf_right, k$Average_DoR_IntMid_conf, k$Average_DoR_Int_ORR, k$Average_DoR_IntMid_ORR, tau)
  z <- z + 1 
  dor_exp_1test2_s2[[tau]] <- k$DoR
  dor_exp_1_no_conftest2_s2[[tau]]     <- k$DoRInt_no_conf
  dor_exp_1_mid_no_conftest2_s2[[tau]] <- k$DoRIntM_no_conf
  dor_exp_1_conftest2_s2[[tau]]        <- k$DoRInt_conf
  dor_exp_1_mid_conftest2_s2[[tau]]    <- k$DoRIntM_conf
  dor_exp_1_ORRtest2_s2[[tau]]         <- k$DoRInt_ORR
  dor_exp_1_mid_ORRtest2_s2[[tau]]     <- k$DoRIntM_ORR
  print(round(z/tau.val,2))
}

colnames(test2_s2) <- c("Average_DoR", "DoR_lower_daily", "DoR_upper_daily","Average_DoR_Int", "DoR_lower_int", "DoR_upper_int", "Average_DoR_Int_Mid", "DoR_lower_mid", "DoR_upper_mid", "Average_DoR_Int_conf", "Average_DoR_Int_Mid_conf", "Average_DoR_ORR", "Average_DoR_ORR_Mid", "tau")


