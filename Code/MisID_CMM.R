#MisID models for Common Minke Whales
#Separated out from others to improve convergence
#WAIC and diagnostics were completed as described in Code/MisIDmodelsStan.R
#Zoe Rand
#Last updated 7/5/24

library(tidyverse)
library(cmdstanr)
library(bayesplot) #for diagnostics 
library(posterior) #for manipulating stan objects
#library(parallel) #run models in parallel--part of baseR
library(loo) #WAIC 

# Reading in Data ---------------------------------------------------------

catchdf<-read_csv("Data/MotherLengthFetalsex.csv")
cmminke_bin<- catchdf %>%
  filter(SpName == "Common minke") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1))

# Stan Model --------------------------------------------------------------
#for parallelization
options(mc.cores = parallel::detectCores())

file<-"Code/Stan Code/MisID_CMM.stan" #limited L50 and delta to improve convergence 
mod<-cmdstan_model(file) 




# Run Model ---------------------------------------------------------------

runStanMisID<-function(df){
  #no correction
  the_data1<-list(N = nrow(df), DatBin = df$Male,
                  fetalLengths = df$F_dec_ft, 
                  correction = 0, correctionMale = 0, correctionFemale = 0, 
                  maxL = 50)
  inits<-function(){
    list(a10 = rnorm(1), 
         b = rnorm(1)
    )
  }
  
  fit1<-mod$sample(data = the_data1, seed = 123, refresh = 200,
                   iter_warmup = 500, iter_sampling = 500, chains = 4, init=inits, adapt_delta = 0.98, max_treedepth = 15)
  parsofint<-c("a10", "b")
  #print(mcmc_trace(fit1$draws(parsofint)))
  fit1$save_object(file = paste0("Results/StanMCMC/", df$SpName[1], "_mod1.RDS"))
  
  #with correction Female
  the_data2<-list(N = nrow(df), DatBin = df$Male,
                  fetalLengths = df$F_dec_ft, 
                  correction = 1, correctionMale = 0, correctionFemale = 1, #females are mistaken for males
                  maxL = 50)
  inits2<-function(){
    list(a10 = rnorm(1), 
         b = rnorm(1), 
         L50 = as.array(rnorm(1, 0, 1)),
         delta = as.array(runif(1, 0, 2))
    )
  }
  
  fit2<-mod$sample(data = the_data2, seed = 123, refresh = 200,
                   iter_warmup = 1000, iter_sampling = 1000, chains = 4, init=inits2, adapt_delta = 0.98, max_treedepth = 15)
  
  parsofint2<-c("a10", "b", "L50", "delta")
  #print(mcmc_trace(fit2$draws(parsofint2), facet_args = list(ncol = 2)))
  fit2$save_object(file = paste0("Results/StanMCMC/", df$SpName[1], "_mod2.RDS"))
  
  #with correction Male
  the_data3<-list(N = nrow(df), DatBin = df$Male,
                  fetalLengths = df$F_dec_ft, 
                  correction = 1, correctionMale = 1, correctionFemale = 0, #males are mistaken for females
                  maxL = 50)
  inits2<-function(){
    list(a10 = rnorm(1), 
         b = rnorm(1), 
         L50 = as.array(rnorm(1, 0, 1)),
         delta = as.array(runif(1, 0, 2))
    )
  }
  
  fit3<-mod$sample(data = the_data3, seed = 123, refresh = 200,
                   iter_warmup = 1000, iter_sampling = 1000, chains = 4, init=inits2, adapt_delta = 0.98, max_treedepth = 15)
  
  parsofint3<-c("a10", "b", "L50", "delta")
  #print(mcmc_trace(fit3$draws(parsofint3), facet_args = list(ncol = 2)))
  
  fit3$save_object(file = paste0("Results/StanMCMC/", df$SpName[1], "_mod3.RDS"))
}

runStanMisID(cmminke_bin)

#checking convergence
Cm_mod1<-readRDS("Results/StanMCMC/Common minke_mod1.RDS")
Cm_mod2<-readRDS("Results/StanMCMC/Common minke_mod2.RDS")
Cm_mod3<-readRDS("Results/StanMCMC/Common minke_mod3.RDS")

parsofint<-c("a10", "b")
print(mcmc_trace(Cm_mod1$draws(parsofint), facet_args = list(ncol = 2)))

parsofint2<-c("a10", "b", "L50", "delta")
print(mcmc_trace(Cm_mod2$draws(parsofint2), facet_args = list(ncol = 2)))
print(mcmc_trace(Cm_mod3$draws(parsofint2), facet_args = list(ncol = 2)))
