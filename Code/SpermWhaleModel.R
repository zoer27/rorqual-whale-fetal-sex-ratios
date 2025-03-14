#Code to run adaptive sex ratio analysis for sperm whales
#Last updated 3/14/25

# Libraries ---------------------------------------------------------------
library(tidyverse)
#use brms version 2.21.0
#brmsurl<-"http://cran.r-project.org/src/contrib/Archive/brms/brms_2.21.0.tar.gz"
#install.packages(brmsurl, type = "source")
library(brms)
library(bayesplot)
library(patchwork)
library(tidybayes)
options(mc.cores = parallel::detectCores())
options(brms.backend = "cmdstanr")
library(ggthemes)


# Fetal sex correction ----------------------------------------------------
L99_res<-read_csv("Code/Results/L99_table.csv")

# Data --------------------------------------------------------------------

#from MotherLenghtFetalSexDataSummary.R
MLenFDat<-read_csv("Data/MotherLengthFetalsex_62923.csv")
table(MLenFDat$SpName)

#isolating sperm whales
SpermDat<-MLenFDat %>% 
  filter(SpName  == "Sperm")

#adding Male 1 0 column

SpermDat<-SpermDat %>% mutate(Male = ifelse(F_S == 1, 1, 0)) #males are coded as 1s and females as 2s
head(SpermDat)

# Removing small fetuses --------------------------------------------------
CutoffLen<-L99_res %>% filter(Species == "Sperm") %>% select(`L99[1]`) %>% unlist()
FilteredDat<-SpermDat %>% filter(SpName == "Sperm") %>% filter(F_dec_ft >= round(CutoffLen, 4))

# Binning Data ------------------------------------------------------------

c(1.5, 4.5, 7.5, 10.5)/12

#binning
BinVals<-FilteredDat  %>% filter(M_dec_ft != 0) %>% 
  group_by(SpName) %>% mutate(minL = floor(min(M_dec_ft)), maxL = ceiling(max(M_dec_ft))) %>% 
  select(SpName, minL, maxL) %>% unique() 

Bins<-mapply(seq,BinVals$minL - 0.125, BinVals$maxL + 0.125, by = 0.25)
Bins

#labels <- sapply(Bins, function(x){format(x[-which(x == tail(x, 1))], digits = 5)}) # Custom labels as lower value of each bin
labels<-format(Bins[-which(Bins == tail(Bins[,1],1)), 1], digits = 5)
#center scale function
scale_this <- function(x){ #need to write your own in order to use it for grouped data
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


Binned3Iinches<-FilteredDat %>% filter(M_dec_ft != 0) %>% 
  group_by(SpName) %>%
  mutate(binned = cut(M_dec_ft, breaks = Bins[,1], include.lowest = TRUE, labels = labels)) %>%
  mutate(binned_center_sc = scale_this(as.numeric(as.character(binned)))) %>%
  group_by(SpName, binned, binned_center_sc) %>%
  summarise(Total = n(), Nmale = sum(Male)) %>% #create binomial data
  filter(Total != 0) %>% #remove empty bins
  filter(binned_center_sc <= 4 && binned_center_sc >= -4)

ggplot(Binned3Iinches) + geom_point(aes(x = binned_center_sc, y = Total)) +
  facet_wrap(~SpName, scales = "free")




# Run Model ---------------------------------------------------------------

#no random effects because only one species
Mod1_Sperm<-brm(Nmale | trials(Total) ~ binned_center_sc,
                data = Binned3Iinches, family = "beta_binomial", 
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 10), class = b)),
                warmup = 5000, 
                iter   = 6000, 
                chains = 4, control = list(adapt_delta = 0.95), 
                file = "Code/Fetal Sex Mother Length/Results/Modelfit_Sperm.rds")
plot(Mod1_Sperm)
summary(Mod1_Sperm)

Mod1_Sperm$fit

ndraws<-4000
Sperm_pneg<-Mod1_Sperm %>% spread_draws(b_binned_center_sc) %>%  
  mutate(neg = ifelse(b_binned_center_sc < 0, 1, 0)) %>% 
  summarise(TotNeg = sum(neg), prob_neg = TotNeg/ndraws)

Sperm_ppos<-Mod1_Sperm %>% spread_draws(b_binned_center_sc) %>%  
  mutate(pos = ifelse(b_binned_center_sc > 0, 1, 0)) %>% 
  summarise(TotPos = sum(pos), prob_pos = TotPos/ndraws)

Sperm_intercept<-Mod1_Sperm %>% spread_draws(b_Intercept) %>%  
  mutate(real_intercept = plogis(b_Intercept))

quantile(Sperm_intercept$real_intercept, p = c(0.025, 0.5, 0.975))


#plot
col<-"#72DC9A"

#get predictions
par_draws<-Mod1_Sperm %>% spread_draws(b_binned_center_sc, b_Intercept)

Sp_lengths<-seq(-2, 2, by = 0.05)

out<-mapply(function(x,y) plogis(x + y*Sp_lengths), par_draws$b_Intercept, par_draws$b_binned_center_sc)

preds_quants<-apply(out, 1, mode_hdci) %>% bind_rows() %>% add_column(mat_length = Sp_lengths)


p1<-ggplot(preds_quants) +
  geom_hline(aes(yintercept = 0.5), color = 'black')+
  geom_line(aes(x = mat_length, y = y), alpha = 1, linewidth = 0.6, color = col) + 
  geom_ribbon(aes(x = mat_length, ymin = ymin, ymax = ymax), alpha = 0, linetype = "dashed", linewidth = 0.6, color = col) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.46, 0.54)) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Maternal length (scaled)", y = "Fetal sex ratio") + 
  theme_tufte(base_size = 12, base_family = "Arial") + 
  theme(legend.position = "none", 
        #text = element_text(family = "Arial", size = 12), 
        strip.text = element_text(color = "black", size = 11, vjust = -1), 
        panel.grid.major.y = element_blank(), 
        panel.spacing = unit(0.5, "lines"))


#ggsave("Figures/FigureS3SpermAdaptive.png", p1, width = 5, height = 5, units = 'in')
