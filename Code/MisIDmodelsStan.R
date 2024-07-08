#Fetal Mis-ID models using Stan
#This code is for all species but Common minke and Bryde's were completed
#separately because a few constraints were placed on those paramters to 
#improve convergence. See additional code (MisID_Brydes.R, MisID_CMM.R) and paper
#for more details.
#last updated 7/5/24


# Load Libraries ----------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(bayesplot) #for diagnostics 
library(posterior) #for manipulating stan objects
library(loo) #WAIC 
library(bayestestR) #for HDI
library(patchwork)#for plotting

# Reading in Data ---------------------------------------------------------

catchdf<-read_csv("Data/MotherLengthFetalsex.csv")


# Data Wrangling ----------------------------------------------------------
antbw_bin<- catchdf %>%
  filter(SpName == "Antarctic blue") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1))

fin_bin<- catchdf %>%
  filter(SpName == "Fin") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1))

sperm_bin<- catchdf %>%
  filter(SpName == "Sperm") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1))

hump_bin<- catchdf %>%
  filter(SpName == "Humpback") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1))

sei_bin<- catchdf %>%
  filter(SpName == "Sei") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1))

sei_bin<- sei_bin %>% filter(F_dec_ft <= 16.5) #remove fetuses that are too large
summary(sei_bin)

#see CMM specific code for model fit
cmminke_bin<- catchdf %>%
  filter(SpName == "Common minke") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1))


antminke_bin<- catchdf %>%
  filter(SpName == "Antarctic minke") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1))

brydes_bin<- catchdf %>%
  filter(SpName == "Brydes") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1))

brydes_bin<- brydes_bin %>% filter(F_dec_ft <= 15) #remove fetuses that are too large
summary(brydes_bin)

#list of dataframes
allDF<-list(antbw_bin, fin_bin, sperm_bin, hump_bin, sei_bin, cmminke_bin, antminke_bin, brydes_bin)

# Stan Model --------------------------------------------------------------
#for parallelization
options(mc.cores = parallel::detectCores())

file<-"Code/Stan Code/MisID.stan"
mod<-cmdstan_model(file) 

# Functions to run models -------------------------------------------------
#model 1 no misID
#model 2 female correction
#model 3 male correction
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
  fit1$save_object(file = paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod1.RDS"))
  
  #with correction Female
  the_data2<-list(N = nrow(df), DatBin = df$Male,
                  fetalLengths = df$F_dec_ft, 
                  correction = 1, correctionMale = 0, correctionFemale = 1, #females are mistaken for males
                  maxL = 50)
  inits2<-function(){
    list(a10 = rnorm(1), 
         b = rnorm(1), 
         L50 = as.array(rnorm(1)),
         delta = as.array(runif(1, 0, 10))
    )
  }
  
  fit2<-mod$sample(data = the_data2, seed = 222, refresh = 200,
                   iter_warmup = 1000, iter_sampling = 1000, chains = 4, init=inits2, adapt_delta = 0.98, max_treedepth = 15)
  
  parsofint2<-c("a10", "b", "L50", "delta")
  #print(mcmc_trace(fit2$draws(parsofint2), facet_args = list(ncol = 2)))
  fit2$save_object(file = paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod2.RDS"))
  
  #with correction Male
  the_data3<-list(N = nrow(df), DatBin = df$Male,
                  fetalLengths = df$F_dec_ft, 
                  correction = 1, correctionMale = 1, correctionFemale = 0, #males are mistaken for females
                  maxL = 50)
  inits2<-function(){
    list(a10 = rnorm(1), 
         b = rnorm(1), 
         L50 = as.array(rnorm(1)),
         delta = as.array(runif(1, 0, 10))
    )
  }
  
  fit3<-mod$sample(data = the_data3, seed = 222, refresh = 200,
                   iter_warmup = 1000, iter_sampling = 1000, chains = 4, init=inits2, adapt_delta = 0.98, max_treedepth = 15)
  
  parsofint3<-c("a10", "b", "L50", "delta")
  #print(mcmc_trace(fit3$draws(parsofint3), facet_args = list(ncol = 2)))
  
  fit3$save_object(file = paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod3.RDS"))
}


#run in parallel

mclapply(allDF, runStanMisID ) #all draws are saved in files, no output in environment from this




# WAIC comparison ---------------------------------------------------------

TraceandWAIC<-function(df){
  mod1<-readRDS(paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod1.RDS"))
  mod2<-readRDS(paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod2.RDS"))
  mod3<-readRDS(paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod3.RDS"))
  
  parsofint<-c("a10", "b")
  print(mcmc_trace(mod1$draws(parsofint), facet_args = list(ncol = 2)))
  
  parsofint2<-c("a10", "b", "L50", "delta")
  print(mcmc_trace(mod2$draws(parsofint2), facet_args = list(ncol = 2)))
  print(mcmc_trace(mod3$draws(parsofint2), facet_args = list(ncol = 2)))
  
  #loo1<-mod1$loo(variables = "loglike")
  #loo2<-mod2$loo(variables = "loglike")
  #loo3<-mod3$loo(variables = "loglike")
  WAIC1<-waic(mod1$draws("loglike"))
  WAIC2<-waic(mod2$draws("loglike"))
  WAIC3<-waic(mod3$draws("loglike"))
  
  return(loo_compare(WAIC1, WAIC2, WAIC3))
}

#test<-TraceandWAIC(allDF[[2]])
WAIC_res<-lapply(allDF, TraceandWAIC)



# WAIC results ------------------------------------------------------------


names(WAIC_res)<-c("Antarctic blue", "Fin", "Sperm", "Humpback", "Sei", "Common minke", "Antarctic minke", "Brydes")

WAIC_res
WAIC_res_tab<-lapply(WAIC_res, as_tibble, rownames = "model")

WAIC_tab<-bind_rows(WAIC_res_tab, .id = "species")

WAIC_tab_for_paper<-WAIC_tab %>% dplyr::select(species, model, waic, elpd_diff, se_diff)

#write_csv(WAIC_tab_for_paper, file = "Code/Results/WAIC_table_for_paper.csv")


# Diagnostics -------------------------------------------------------------
getRhatESS<-function(df){
  mod1<-readRDS(paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod1.RDS"))
  mod2<-readRDS(paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod2.RDS"))
  mod3<-readRDS(paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod3.RDS"))
  
  rhats1<-summary(bayesplot::rhat(mod1))[c(1,6)]
  rhats2<-summary(bayesplot::rhat(mod2))[c(1,6)]
  rhats3<-summary(bayesplot::rhat(mod3))[c(1,6)]
  
  neff1<-summary(bayesplot::neff_ratio(mod1))[c(1,6)]
  neff2<-summary(bayesplot::neff_ratio(mod2))[c(1,6)]
  neff3<-summary(bayesplot::neff_ratio(mod3))[c(1,6)]
  
  return_tab<-tibble("model" = c(1,2,3), rhat_min = c(rhats1[1], rhats2[1], rhats3[1]), 
                     rhat_max = c(rhats1[2], rhats2[2], rhats3[2]), 
                     neff_min = c(neff1[1], neff2[1], neff3[1]), 
                     neff_max = c(neff1[2], neff2[2], neff3[2]))
  
  return(return_tab)
}


diag_res<-lapply(allDF, getRhatESS)
names(diag_res)<-c("Antarctic blue", "Fin", "Sperm", "Humpback", "Sei", "Common minke", "Antarctic minke", "Brydes")
diag_tab<-bind_rows(diag_res, .id = "Species")
#write_csv(diag_tab, "Code/Results/misID_bayesiandiag.csv")



# Parameter table ---------------------------------------------------------
SummTable<-function(df){
  mod1<-readRDS(paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod1.RDS"))
  mod2<-readRDS(paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod2.RDS"))
  mod3<-readRDS(paste0("Code/Results/StanMCMC/", df$SpName[1], "_mod3.RDS"))
  
  parsofint<-c("a10", "b")
  sumstats1<-mod1$summary(parsofint) %>% 
    left_join(mod1$summary(parsofint, ~quantile(.x, probs = c(0.025, 0.975)))) %>%
    dplyr::select(variable, median, `2.5%`, `97.5%`) %>% add_column("model" = "model1")
  #sumstats1
  
  parsofint2<-c("a10", "b", "L50", "delta")
  sumstats2<-mod2$summary(parsofint2) %>% 
    left_join(mod2$summary(parsofint2, ~quantile(.x, probs = c(0.025, 0.975)))) %>%
    dplyr::select(variable, median, `2.5%`, `97.5%`) %>% add_column("model" = "model2")
  sumstats3<-mod3$summary(parsofint2) %>% 
    left_join(mod3$summary(parsofint2, ~quantile(.x, probs = c(0.025, 0.975)))) %>%
    dplyr::select(variable, median, `2.5%`, `97.5%`) %>% add_column("model" = "model3")
  
  sum_tab<-bind_rows(sumstats1, sumstats2, sumstats3)
  return(sum_tab)
}
sum_tabs<-lapply(allDF, SummTable)
names(sum_tabs)<-c("Antarctic blue", "Fin", "Sperm", "Humpback", "Sei", "Common minke", "Antarctic minke", "Brydes")
summary_table<-bind_rows(sum_tabs, .id = "Species") %>% pivot_wider(names_from = variable, values_from = c(median, `2.5%`, `97.5%`)) %>%
  select(Species, model, median_a10, `2.5%_a10`, `97.5%_a10`, median_b, `2.5%_b`, `97.5%_b`, 
         `median_L50[1]`, `2.5%_L50[1]`, `97.5%_L50[1]`, 
         `median_delta[1]`, `2.5%_delta[1]`, `97.5%_delta[1]`)
#fixing the order so it matches the WAIC table
summary_table_ordered<-WAIC_tab_for_paper %>% dplyr::select(species, model) %>% left_join(summary_table, by = c("species" = "Species", "model"))
write_csv(summary_table_ordered, "Code/Results/FetalMisID_summaryresults.csv")


# Results -----------------------------------------------------------------

bestmods<-c("2", "2", "2", "2", "2", "1", "3", "1")
readBestMods<-function(df, bestmod){
  bestMod<-readRDS(paste0("Results/StanMCMC/", df$SpName[1], "_mod", bestmod,".RDS"))
  return(bestMod)
}
modlist<-mapply(readBestMods, allDF, bestmods)
names(modlist)<-c("Antarctic blue", "Fin", "Sperm", "Humpback", "Sei", "Common minke", "Antarctic minke", "Brydes")

#get predictions
predm_list<-lapply(modlist, function(x){subset_draws(x$draws(), variable = "predm")})
sum_predm<-lapply(predm_list, function(x){summarise_draws(x, "median", ~quantile(.x, probs = c(0.025, 0.975)))[c("median", "2.5%", "97.5%")]})

sum_pred_tib<-mapply(function(x,y){a<-x %>% as_tibble %>% add_column(F_dec_ft = y$F_dec_ft); print(y$SpName[1]); return(a)},
                     sum_predm, allDF, SIMPLIFY = FALSE)
names(sum_pred_tib)<-c("Antarctic blue", "Fin", "Sperm", "Humpback", "Sei", "Common minke", "Antarctic minke", "Brydes")
pred_tab<-bind_rows(sum_pred_tib, .id = "Species")

#modes and highest density intervals

sum_predm2<-lapply(predm_list, function(x){apply(x, 3, ggdist::mode_hdci) %>% bind_rows(.id = "variable") %>% dplyr::select("y", "ymin", "ymax") %>%
    rename("mode" = "y", "2.5%" = "ymin", "97.5%" = "ymax")})

sum_pred_tib2<-mapply(function(x,y){a<-x %>% as_tibble %>% add_column(F_dec_ft = y$F_dec_ft); print(y$SpName[1]); return(a)},
                     sum_predm2, allDF, SIMPLIFY = FALSE)
names(sum_pred_tib2)<-c("Antarctic blue", "Fin", "Sperm", "Humpback", "Sei", "Common minke", "Antarctic minke", "Brydes")
pred_tab2<-bind_rows(sum_pred_tib2, .id = "Species")


# 99% correct -------------------------------------------------------------

#get L99 
#remove common minke whales and Bryde's whales because no correction
modlist1<-within(modlist, rm("Common minke", "Brydes"))

L99_list<-lapply(modlist1, function(x){subset_draws(x$draws(), variable = "L99")})

#check posteriors
for(i in 1:6){
  print(mcmc_dens(L99_list[[i]]))
}

#summarise posterior of L99
mode_hdci_L99<-lapply(L99_list,function(x){x %>% tidy_draws() %>% mode_hdci(`L99[1]`)})


L99_tab<-bind_rows(mode_hdci_L99, .id = "Species") 
L99_tab

#write_csv(L99_tab, "Code/Results/L99_table.csv")

# Probability of a negative slope -----------------------------------------
slope_list<-lapply(modlist, function(x){subset_draws(x$draws(), variable = "b")})


prob_neg<-function(x){
  a<-ifelse(x<0, 1, 0) #is slope negative
  b<-sum(a)
  c<-b/length(x) #divide by total iterations
}

prob_slope<-lapply(slope_list, prob_neg)

names(prob_slope)

prob_tab<-bind_rows(prob_slope, .id = "Species") %>% pivot_longer(everything(), names_to = "Species", values_to = "Prob")
prob_tab


#write_csv(prob_tab, "Code/Results/Prob_neg_fetal_length.csv")



# Number of data points ---------------------------------------------------

Num_fetus<-bind_rows(allDF) %>% group_by(SpName) %>% summarise(n= n())

sum(Num_fetus$n)

# Figures -----------------------------------------------------------------
#colors
pal<-c("#3696ad", "#72DC9A","#5e8ed5","#9482d6", "#c670c3", "#eb5e9f", "#fc576f", "#f56638", "#ffbd59")#add extra color for sperm whales
SpList<-c("Antarctic blue", "Sperm", "Pygmy blue","Fin","Humpback", "Sei",  "Brydes", "Antarctic minke", "Common minke")
col <-setNames(pal, SpList)


#datasets for plotting
#Antarctic blue
antbw_bin_forplot<- catchdf %>%
  filter(SpName == "Antarctic blue") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1)) %>%
  mutate(bins = cut(F_dec_ft, breaks = seq(0, max(F_dec_ft), by = 0.0833), ordered_result = TRUE)) %>%
  mutate(x_tmp = str_sub(bins, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  #group_by(bins) %>%
  group_by(bins, min, max) %>%
  summarise(n = n(), nmale = sum(Male), prop_male = nmale/n)


#Fin
fin_bin_forplot<- catchdf %>%
  filter(SpName == "Fin") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1)) %>%
  mutate(bins = cut(F_dec_ft, breaks = seq(0, max(F_dec_ft), by = 0.0833), ordered_result = TRUE)) %>%
  mutate(x_tmp = str_sub(bins, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  #group_by(bins) %>%
  group_by(bins, min, max) %>%
  summarise(n = n(), nmale = sum(Male), prop_male = nmale/n)


#Sperm
sperm_bin_forplot<- catchdf %>%
  filter(SpName == "Sperm") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1)) %>%
  mutate(bins = cut(F_dec_ft, breaks = seq(0, max(F_dec_ft), by = 0.0833), ordered_result = TRUE)) %>%
  mutate(x_tmp = str_sub(bins, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  #group_by(bins) %>%
  group_by(bins, min, max) %>%
  summarise(n = n(), nmale = sum(Male), prop_male = nmale/n)


#Humpback
hump_bin_forplot<- catchdf %>%
  filter(SpName == "Humpback") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1)) %>%
  mutate(bins = cut(F_dec_ft, breaks = seq(0, max(F_dec_ft), by = 0.0833), ordered_result = TRUE)) %>%
  mutate(x_tmp = str_sub(bins, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  #group_by(bins) %>%
  group_by(bins, min, max) %>%
  summarise(n = n(), nmale = sum(Male), prop_male = nmale/n)

#Sei

sei_bin_forplot<- catchdf %>%
  filter(SpName == "Sei") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1)) %>%
  filter(F_dec_ft <= 16.5) %>%
  mutate(bins = cut(F_dec_ft, breaks = seq(0, max(F_dec_ft), by = 0.0833), ordered_result = TRUE)) %>%
  mutate(x_tmp = str_sub(bins, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  #group_by(bins) %>%
  group_by(bins, min, max) %>%
  summarise(n = n(), nmale = sum(Male), prop_male = nmale/n)


#Common minke
cmminke_bin_forplot<- catchdf %>%
  filter(SpName == "Common minke") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1)) %>%
  mutate(bins = cut(F_dec_ft, breaks = seq(0, max(F_dec_ft), by = 0.0833), ordered_result = TRUE)) %>%
  mutate(x_tmp = str_sub(bins, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  #group_by(bins) %>%
  group_by(bins, min, max) %>%
  summarise(n = n(), nmale = sum(Male), prop_male = nmale/n)


#Antarctic minke
antminke_bin_forplot<- catchdf %>%
  filter(SpName == "Antarctic minke") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1)) %>%
  mutate(bins = cut(F_dec_ft, breaks = seq(0, max(F_dec_ft), by = 0.0833), ordered_result = TRUE)) %>%
  mutate(x_tmp = str_sub(bins, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  #group_by(bins) %>%
  group_by(bins, min, max) %>%
  summarise(n = n(), nmale = sum(Male), prop_male = nmale/n)


#Brydes
brydes_bin_forplot<- catchdf %>%
  filter(SpName == "Brydes") %>%
  mutate(Male = ifelse(F_S == 2, 0, 1)) %>%
  filter(F_dec_ft <= 15) %>%
  mutate(bins = cut(F_dec_ft, breaks = seq(0, max(F_dec_ft), by = 0.0833), ordered_result = TRUE)) %>%
  mutate(x_tmp = str_sub(bins, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  #group_by(bins) %>%
  group_by(bins, min, max) %>%
  summarise(n = n(), nmale = sum(Male), prop_male = nmale/n)

#combining dataframes together
plot_dfs<-list(antbw_bin_forplot, fin_bin_forplot, sperm_bin_forplot, hump_bin_forplot, sei_bin_forplot, cmminke_bin_forplot, antminke_bin_forplot, brydes_bin_forplot)
names(plot_dfs)<-c("Antarctic blue", "Fin", "Sperm", "Humpback", "Sei", "Common minke", "Antarctic minke", "Brydes")

plot_df_tab<-lapply(plot_dfs, function(x) { x$bins<-as.character(x$bins); return(x)}) %>% bind_rows(.id = "Species")


#plot species labels including sample size and probability of negative slope
Sp_labeller<-c(`Brydes` = paste0("Bryde's \n n = ", formatC(Num_fetus$n[Num_fetus$SpName == "Brydes"], big.mark = ","),", P(neg.) = ", sprintf("%.3f", prob_tab$Prob[prob_tab$Species == "Brydes"])), 
               `Antarctic minke` = paste0("Antarctic minke \n n = ", formatC(Num_fetus$n[Num_fetus$SpName == "Antarctic minke"], big.mark = ","), ", P(neg.) = ", sprintf("%.3f", prob_tab$Prob[prob_tab$Species == "Antarctic minke"])),
              `Common minke` = paste0("Common minke \n n = ", formatC(Num_fetus$n[Num_fetus$SpName == "Common minke"], big.mark = ","), ", P(neg.) = ", sprintf("%.3f", prob_tab$Prob[prob_tab$Species == "Common minke"])),
              `Sperm` = paste0("Sperm \n n = ", formatC(Num_fetus$n[Num_fetus$SpName == "Sperm"], big.mark = ","), ", P(neg.) = ", sprintf("%.3f", prob_tab$Prob[prob_tab$Species == "Sperm"])),
              `Antarctic blue` = paste0("Antarctic blue \n n = ", formatC(Num_fetus$n[Num_fetus$SpName == "Antarctic blue"], big.mark = ","), ", P(neg.) = ", sprintf("%.3f", prob_tab$Prob[prob_tab$Species == "Antarctic blue"])),
              `Fin` = paste0("Fin \n n = ", formatC(Num_fetus$n[Num_fetus$SpName == "Fin"], big.mark = ","), ", P(neg.) = ", sprintf("%.3f", prob_tab$Prob[prob_tab$Species == "Fin"])),
              `Humpback` = paste0("Humpback \n n = ", formatC(Num_fetus$n[Num_fetus$SpName == "Humpback"], big.mark = ","), ", P(neg.) = ", sprintf("%.3f", prob_tab$Prob[prob_tab$Species == "Humpback"])), 
              `Sei` = paste0("Sei \n n = ", formatC(Num_fetus$n[Num_fetus$SpName == "Sei"], big.mark = ","), ", P(neg.) = ", sprintf("%.3f", prob_tab$Prob[prob_tab$Species == "Sei"])))


#plot
FetalMisID<-ggplot() + geom_point(data = plot_df_tab, aes(x = min, y = prop_male, size = n, color = Species), alpha = 0.8) + 
  geom_ribbon(data = pred_tab2, aes(x = F_dec_ft, ymin = `2.5%`, ymax = `97.5%`, group = Species), alpha = 0, linetype = "dashed", linewidth = 0.4, color = "gray30") + 
  geom_line(data = pred_tab2, aes(x = F_dec_ft, y = mode), linewidth = 0.7, color = "gray30") + 
  theme_minimal(base_size = 12, base_family = "Arial") + 
  scale_color_manual(values = col, guide = "none") +
  scale_size_area(breaks = c(100, 500, 1000, 3000, 6000), max_size = 8) + 
  labs(x = "Fetal length (ft)", y = expression(paste("Sex ratio (", p[male], ")")), size = "Number of\n fetuses") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0.1)) + 
  facet_wrap(~factor(Species, levels = c("Common minke", "Antarctic minke", "Brydes", 
                                         "Sei", "Humpback", "Fin", "Sperm", "Antarctic blue")), ncol = 2, scales = "free_x", 
             labeller = as_labeller(Sp_labeller)) + 
  theme(axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        panel.spacing = unit(0.8, "lines"))

FetalMisID


#ggsave("Figures/Fig2FetalMisID.png", FetalMisID, height = 8, width = 7, units = "in", dpi = 600)



