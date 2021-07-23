###############
# COVID MODEL #
###############

#packages requred
library(readr)
library(ggplot2)
library(nimble)
library(tidyverse)
library(reshape2)
library(mgcv)
library(lubridate)
library(coda)
library(doParallel)
# Set working directory
setwd("~/Dropbox/Spatial modelling/Covid/Covid_rep")

# READ IN DATA: all data downloaded on 2021-07-22
# cases and deaths
England_CasesDeaths <- read_csv("England_CasesDeaths.csv")
# tests
England_Tests <- read_csv("England_newTest_pillars.csv")
# UK test pillar capacity 
UK_plannedCapacity <- read_csv("UK_plannedCapacity.csv")
# Vaccines Doses
England_vaccines <- read_csv("England_cumVac_2021-07-22.csv")
England_vaccines$date<-as.Date(England_vaccines$date)

# plot data 
colour_cases<-c("Cases (Specimen Date)"= "dark green", "Cases (Publish Date)"="green")
ggplot()+
  geom_line(data=England_CasesDeaths, aes(x=date, y=newCasesByPublishDate, colour="Cases (Publish Date)"))+
  geom_line(data=England_CasesDeaths, aes(x=date, y=newCasesBySpecimenDate, colour="Cases (Specimen Date)"))+
  scale_color_manual(name='Daily Numbers', values=colour_cases)
ggplot()+
  geom_line(data=England_CasesDeaths, aes(x=date, y=log(newCasesByPublishDate), colour="Cases (Publish Date)"))+
  geom_line(data=England_CasesDeaths, aes(x=date, y=log(newCasesBySpecimenDate), colour="Cases (Specimen Date)"))+
  scale_color_manual(name='Daily Numbers', values=colour_cases)

colour_deaths<-c("Deaths (Death Date-14)"= "dark blue", "Deaths (Publish Date)"="blue")
ggplot()+
  geom_line(data=England_CasesDeaths, aes(x=date, y=newDeaths28DaysByPublishDate, colour="Deaths (Publish Date)"))+
  geom_line(data=England_CasesDeaths, aes(x=as.Date(date-14), y=newDeaths28DaysByDeathDate, colour="Deaths (Death Date-14)"))+
  scale_color_manual(name='Daily Numbers', values=colour_deaths)

ggplot()+
  geom_line(data=England_CasesDeaths, aes(x=date, y=log(newDeaths28DaysByPublishDate), colour="Deaths (Publish Date)"))+
  geom_line(data=England_CasesDeaths, aes(x=as.Date(date-14), y=log(newDeaths28DaysByDeathDate), colour="Deaths (Death Date-14)"))+
  scale_color_manual(name='Daily Numbers', values=colour_deaths)

colour_tests<-c("Cases (Specimen Date)"= "dark green", "Cases (Publish Date)"="green", "Virus Test"="red", "UK Lab Test Capacity"="purple")
ggplot()+
  geom_line(data=England_CasesDeaths, aes(x=date, y=newCasesByPublishDate, colour="Cases (Publish Date)"))+
  geom_line(data=England_CasesDeaths, aes(x=date, y=newCasesBySpecimenDate, colour="Cases (Specimen Date)"))+
  geom_line(data=England_Tests, aes(x=date, y=newPillarOneTwoTestsByPublishDate, colour="Virus Test"))+
  geom_line(data=England_Tests, aes(x=date, y=newPillarOneTwoTestsByPublishDate, colour="Virus Test"))+
  geom_line(data=UK_plannedCapacity, aes(x=date, y=c(capacityPillarOne+capacityPillarTwo), colour="UK Lab Test Capacity"))+
  scale_color_manual(name='Daily Numbers', values=colour_tests)

colour_deathcases<-c("Cases (Specimen Date)"= "dark green", "Deaths (Death Date-14)"= "dark blue")
ggplot()+
  geom_line(data=England_CasesDeaths, aes(x=date, y=newCasesBySpecimenDate, colour="Cases (Specimen Date)"))+
  geom_line(data=England_CasesDeaths, aes(x=as.Date(date-14), y=newDeaths28DaysByDeathDate, colour="Deaths (Death Date-14)"))+
  scale_color_manual(name='Daily Numbers', values=colour_deathcases)
# Cases by publish date have spikes due to reporting technical difficulties 
# therefore it makes sense to use cases by specimen date
# and therefore Deaths by Death Date - adjusted by 14 days
# Tests are only available by publish Dates but most lab results are available
# within 72 hours and lateral flow tests are available within the hour
# Lab testing capacity in the UK is less than the tests in england as it doesn't
# take into account lateral flow capacity 

#  adjusted deaths
England_DeathDateADJUSTED<-England_CasesDeaths[,c("date","newDeaths28DaysByDeathDate")]
England_DeathDateADJUSTED$date<-as.Date(England_CasesDeaths$date-14)
colnames(England_DeathDateADJUSTED)[ncol(England_DeathDateADJUSTED)]<-c("newDeathsAdjusted")

# cases by specimen date
England_cases<-England_CasesDeaths[,c("date","newCasesBySpecimenDate")]

# Relevant data for dates where everything is observed 
England_CasesDeaths2<-full_join( England_cases,England_DeathDateADJUSTED, by=c("date"), all.x=TRUE)
All_data<-full_join(England_Tests[,c("date","newPillarOneTwoTestsByPublishDate")], England_CasesDeaths2, by=c("date"), all.x=TRUE)
All_data<-All_data[order(All_data$date),]
All_data_obs<-na.omit(All_data)

# add vaccination data
All_data_obs<-left_join(All_data_obs, England_vaccines[,4:6], by=c("date"))
All_data_obs$t<-as.numeric(as.factor(All_data_obs$date))
colnames(All_data_obs)<-c("date","TestsPublish","CasesSpecimen","DeathsAdj","cumFirstDose","cumSecondDose","t")
All_data_obs$cumFirstDose[is.na(All_data_obs$cumFirstDose)]<-0
All_data_obs$cumSecondDose[is.na(All_data_obs$cumSecondDose)]<-0

All_data_obs2<-All_data_obs[247:nrow(All_data_obs),]
#### MODEL ####

# NIMBLE model code. 
covid_code <- nimbleCode({
    for(t in 1:N){
      # Model total cases (latent).
     # C[t]~ dpois(P*IR[t])
      # Model for Confirmed Cases (observed).
      CC[t]~dbinom(size=Te[t],prob=(1-(1-IR[t])^phi))   
      # Model for total deaths (observed).
      D[t] ~ dpois(lambda=((P)*IFR[t]*IR[t])) #
      logit(IFR[t])<- alpha0 + alpha1*dose1[t] + alpha2*dose2[t]  #+time trend for improved treatment
      
    }
  # Priors:
  # RW1 for infection rate IR.
  IR[1] ~ T(dnorm(0.036,sd=sigma),0,1) 
  for(t in 1:c(N-1)){
    IR[t+1] ~ T(dnorm(IR[t],sd=sigma),0,1)
  }
 # for(t in 1:N){
    phi~dunif(1,1+gamma)
    gamma~dexp(0.1)
 # }
  sigma ~ T(dnorm(0,1),0,)
  # IFR paramters.
  alpha0 ~ dnorm(0, sd=1)
  alpha1 ~ T(dnorm(-1, sd=1),,0)
  alpha2 ~ T(dnorm(-1, sd=1),,0)
 
})

day_diff<-as.numeric(max(England_cases$date)-min(England_cases$date))
England_Population_2019<-56286961
data<-All_data_obs
# Constants
covid_constants <- list(N=nrow(data), P=England_Population_2019)
# Data 
covid_data <- list(D=data$DeathsAdj, CC=data$CasesSpecimen, 
                   Te=data$TestsPublish, 
                  dose1=(data$cumFirstDose-mean(data$cumFirstDose))/sd(data$cumFirstDose),
                   dose2=(data$cumSecondDose-mean(data$cumSecondDose))/sd(data$cumSecondDose))

# Set up the NIMBLE model for each chain.
n_chains<-2
covid_inits <- covid_model <- covid_compiled_model <- covid_mcmc_config <- covid_mcmc <- covid_compiled_mcmc <- list()
for(i in 1:n_chains){
  # Initial values
  covid_inits[[i]] <- list(alpha0=rnorm(1,0,sd=1),
                           alpha1=-rexp(1,0.01),
                           alpha2=-rexp(1,0.01),
                           sigma=rexp(1,0.01),
                           phi=c(1+rexp(1,0.5)),   #nrow(All_data_obs)
                           gamma=rexp(1,0.05),
                           IR=runif(nrow(data),0,1),
                           IFR=runif(nrow(data),0,1))
  # Build the model.
  covid_model[[i]] <- nimbleModel(covid_code,covid_constants,covid_data,covid_inits[[i]])
  # Compile the model.
  covid_compiled_model[[i]] <- compileNimble(covid_model[[i]])
  # Set up the MCMC.
  covid_mcmc_config[[i]] <- configureMCMC(covid_model[[i]],monitors=c('alpha0','IFR','IR','phi', 'sigma','phi','alpha1', 'alpha2', 'gamma')) #
  covid_mcmc[[i]] <- buildMCMC(covid_mcmc_config[[i]])
  # Compile the MCMC.
  covid_compiled_mcmc[[i]] <- compileNimble(covid_mcmc[[i]],project=covid_model[[i]])
}

#help(modelInitialization)
covid_compiled_model[[1]]$initializeInfo()
covid_compiled_model[[1]]$getLogProb("phi")
covid_compiled_model[[1]]$getLogProb("IR")
covid_compiled_model[[1]]$getLogProb("IFR")
covid_compiled_model[[1]]$getLogProb("alpha0")

#check intialisation

# Run the model in parallel. Change 'dopar' to 'do' to run the model on only one core.
registerDoParallel(cores=n_chains)
system.time({
  covid_samples <- as.mcmc.list(foreach(i=1:n_chains)%dopar%{
    runMCMC(covid_compiled_mcmc[[i]],niter=40000,nburnin=20000,inits=covid_inits[[i]],nchains=1,samplesAsCodaMCMC = TRUE,thin=10)
  })
})

# check trace plots for convergence

plot(covid_samples[,"IFR[1]"], main=c("IFR[1]"))
plot(covid_samples[,"alpha0"], main=c("alpha0"))
plot(covid_samples[,"alpha1"], main=c("alpha1"))
plot(covid_samples[,"alpha2"], main=c("alpha2"))
plot(covid_samples[,"IR[1]"], main=c("IR[1]"))

# Combine all MCMC chains.
covid_combined_samples <- as_tibble(do.call('rbind',covid_samples))

# Overall temporal effect on covid fatality.
covid_alpha0 <- select(covid_combined_samples,contains('alpha0'))%>%as.matrix()
median(covid_alpha0)
covid_alpha1 <- select(covid_combined_samples,contains('alpha1'))%>%as.matrix()
median(covid_alpha1)
covid_alpha2 <- select(covid_combined_samples,contains('alpha2'))%>%as.matrix()
median(covid_alpha2)

#preferential testing parameter
covid_phi <- select(covid_combined_samples,contains('phi'))%>%as.matrix()
median(covid_phi)
# max preferential testing parameter -1
covid_gamma <- select(covid_combined_samples,contains('gamma'))%>%as.matrix()
median(covid_gamma)

# SD of IR
covid_sigma <- select(covid_combined_samples,contains('sigma'))%>%as.matrix()
median(covid_sigma)

# IR
covid_IR<- select(covid_combined_samples,contains('IR'))%>%as.matrix()
covid_IR_median<-apply(covid_IR,2,median)
plot(1:nrow(data),covid_IR_median)

# IFR
covid_IFR<- select(covid_combined_samples,contains('IFR'))%>%as.matrix()
covid_IFR_median<-apply(covid_IFR,2,median)
plot(1:nrow(data),covid_IFR_median)

# Simulate cases
sim_Cases=t(apply(England_Population_2019*covid_IR*(1-covid_IFR),1,function(x)rpois(nrow(data),x)+data$DeathsAdj))

data$simCases0.5<-apply(sim_Cases,2,quantile,0.5)
data$simCases0.95<-apply(sim_Cases,2,quantile,0.95)
data$simCases0.05<-apply(sim_Cases,2,quantile,0.05)

ggplot(data)+
  geom_line(aes(x=date, y=CasesSpecimen, colour="Confirmed Cases"))+
  geom_line(aes(x=date, y=simCases0.5, colour="Simulated Cases"))+
  geom_hline(yintercept = England_Population_2019)



