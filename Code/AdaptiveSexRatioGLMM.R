#Bayesian GLMM Mother Length Fetal Sex
#Zoe Rand
#Last updated: 7/5/24


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(brms)
library(bayesplot)
library(patchwork)
library(tidybayes)
options(mc.cores = parallel::detectCores())
options(brms.backend = "cmdstanr")


# Fetal sex correction ----------------------------------------------------
#L99s estimated from mis-ID model
L99_res<-read_csv("Results/L99_table.csv")
#Pygmy blue whale is same as Antarctic
#common minke is same as Antarctic
#Bryde's same as Sei

#create cutoff dataframe
Cutoff<-L99_res %>% add_row(Species = "Pygmy blue", `L99[1]` = L99_res$`L99[1]`[L99_res$Species == "Antarctic blue"]) %>%
  add_row(Species = "Brydes", `L99[1]` = L99_res$`L99[1]`[L99_res$Species == "Sei"]) %>%
  add_row(Species = "Common minke", `L99[1]` = L99_res$`L99[1]`[L99_res$Species == "Antarctic minke"])

Cutoff
# Data --------------------------------------------------------------------

#read in Mother and fetal length data
MLenFDat<-read_csv("Data/MotherLengthFetalsex.csv")
table(MLenFDat$SpName)

#removing species with too small sample sizes:
MLenFDat2<-MLenFDat %>% 
  filter(SpName %in% c("Sperm", "Humpback","Fin","Sei","Pygmy blue","Brydes","Antarctic minke", 
                       "Common minke", "Antarctic blue"))
table(MLenFDat2$SpName)

#removing sperm whales because not roquals
Roq_MLen<-MLenFDat2 %>% filter(SpName != "Sperm")

#adding Male 1 0 column

Roq_MLen<-Roq_MLen %>% mutate(Male = ifelse(F_S == 1, 1, 0)) #males are coded as 1s and females as 2s
head(Roq_MLen)

SpList<-Cutoff$Species
SpList<-SpList[SpList != "Sperm"]#removing sperm whales
# Removing small fetuses --------------------------------------------------

SubsetFetuses<-function(Sp_dat){
  CutoffLen<-Cutoff %>% filter(Species == Sp_dat) %>% select(`L99[1]`) %>% unlist()
  FilteredDat<-Roq_MLen %>% filter(SpName == Sp_dat) %>% filter(F_dec_ft >= round(CutoffLen, 4))
  return(FilteredDat)
  #return(Cutoff)
}

#subsetting all species
FilteredDat<-lapply(SpList,SubsetFetuses) %>% bind_rows()

#test
FilteredDat %>% filter(SpName == "Antarctic minke") %>% summary(F_dec_ft)

# Binning Data ------------------------------------------------------------
#binning data and accounting for differences in length mesurement practices

#binning
BinVals<-FilteredDat  %>% filter(M_dec_ft != 0) %>% 
  group_by(SpName) %>% mutate(minL = floor(min(M_dec_ft)), maxL = ceiling(max(M_dec_ft))) %>% 
  select(SpName, minL, maxL) %>% unique() 

Bins<-mapply(seq,BinVals$minL - 0.125, BinVals$maxL + 0.125, by = 0.25)
Bins

labels <- lapply(Bins, function(x){format(x[-11], digits = 5)}) # Custom labels as lower value of each bin
labels
#center scale function
scale_this <- function(x){ #need to write your own in order to use it for grouped data
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


#bin data
Binned3Iinches<-FilteredDat %>% filter(M_dec_ft != 0) %>% 
  group_by(SpName) %>%
  mutate(binned = cut(M_dec_ft, breaks = Bins[[which(SpList == SpName[1])]], include.lowest = TRUE, labels = labels[[which(SpList == SpName[1])]])) %>%
  mutate(binned_center_sc = scale_this(as.numeric(as.character(binned)))) %>%
  group_by(SpName, binned, binned_center_sc) %>%
  summarise(Total = n(), Nmale = sum(Male)) %>% #create binomial data
  filter(Total != 0) %>% #remove empty bins
  filter(binned_center_sc <= 4 && binned_center_sc >= -4)

#plot
ggplot(Binned3Iinches) + geom_point(aes(x = binned_center_sc, y = Total)) +
  facet_wrap(~SpName, scales = "free")


# Model -------------------------------------------------------------------
#Random slope and intercept for species, no correlation parameter
#Full results from model used in paper are saved as "Results/GLMM_Modelfit_44.rds" so no need to run this
Mod1_ALLRE<-brm(Nmale | trials(Total) ~ binned_center_sc + (1+binned_center_sc||SpName), #double bar means no correlation
                       data = Binned3Iinches, family = "beta_binomial", 
                       prior = c(prior(normal(0, 10), class = Intercept),
                                 prior(normal(0, 10), class = b), 
                                 prior(exponential(1), class = sd)),
                       warmup = 5000, 
                       iter   = 6000, 
                       chains = 4, control = list(adapt_delta = 0.95), 
                      file = "Results/GLMM_Modelfit_44.rds")



plot(Mod1_ALLRE)
summary(Mod1_ALLRE)

# Numerical Results -------------------------------------------------------
#read in model fit if needed
Mod1<-readRDS("Results/GLMM_Modelfit_44.rds")

#number of saved draws
ndraws<-4*1000

#probability slope is negative
species_slopes<-Mod1 %>% spread_draws(b_binned_center_sc, r_SpName[c,t]) %>% mutate(Species = gsub("\\.", replacement = " ", c)) %>%
  filter(t == "binned_center_sc") %>% mutate(Sp_slope = b_binned_center_sc + r_SpName)

species_neg<-species_slopes %>% mutate(neg = ifelse(Sp_slope < 0, 1, 0)) %>% 
  group_by(Species) %>% summarise(TotNeg = sum(neg), prob_neg = TotNeg/ndraws)

pop_slope_neg<-spread_draws(Mod1, b_binned_center_sc) %>% mutate(Species = "Group") %>% 
  mutate(neg = ifelse(b_binned_center_sc < 0, 1, 0)) %>% summarise(TotNeg = sum(neg), prob_neg = TotNeg/ndraws)

sd_slope<-spread_draws(Mod1, sd_SpName__binned_center_sc) %>% reframe(mode_hdci(sd_SpName__binned_center_sc)) %>%
  add_column(Par = "sd_slope") %>% add_column(Rhat = summary(Mod1)$random$SpName$Rhat[2], ESS = summary(Mod1)$random$SpName$Tail_ESS[2])


#intercepts

species_intercepts <- Mod1 %>% spread_draws(b_Intercept, r_SpName[c,t]) %>% mutate(Species = gsub("\\.", replacement = " ", c)) %>%
  filter(t == "Intercept") %>% mutate(Sp_Intercept_logit = (b_Intercept + r_SpName), Sp_Intercept_real = plogis(Sp_Intercept_logit)) %>%
  ungroup() %>% group_by(c) %>% reframe(mode_hdci(Sp_Intercept_real))
  #summarise(med = median(Sp_Intercept_real), low = quantile(Sp_Intercept_real, 0.025), up = quantile(Sp_Intercept_real, 0.975))


pop_intercept_sum<-spread_draws(Mod1, b_Intercept) %>% 
  mutate(b_Intercept_real = plogis(b_Intercept), Species = "Group") %>%
  reframe(mode_hdci(b_Intercept_real)) %>% add_column(Par = "Intercept") %>%
  add_column(Rhat = summary(Mod1)$fixed$Rhat[1], ESS = summary(Mod1)$fixed$Tail_ESS[1])


#standard deviation
sd_intercept<-spread_draws(Mod1, sd_SpName__Intercept) %>% reframe(mode_hdci(sd_SpName__Intercept)) %>% add_column(Par = "sd_intercept") %>%
  add_column(Rhat = summary(Mod1)$random$SpName$Rhat[1], ESS = summary(Mod1)$random$SpName$Tail_ESS[1])
sd_intercept

#overdispersion
phi<-spread_draws(Mod1, phi) %>% reframe(mode_hdci(phi)) %>% add_column(Par = "phi") %>%
  add_column(Rhat = summary(Mod1)$spec_pars$Rhat[1], ESS = summary(Mod1)$spec_pars$Tail_ESS[1])
phi

#population slope for table
pop_slope<-spread_draws(Mod1, b_binned_center_sc) %>% reframe(mode_hdci(b_binned_center_sc)) %>% add_column(Par = "Slope") %>%
  add_column(Rhat = summary(Mod1)$fixed$Rhat[2], ESS = summary(Mod1)$fixed$Tail_ESS[2])

#number of data points
numFetus<-Binned3Iinches %>% group_by(SpName) %>% summarise(sum(Total))

#table for supplement
S4tab<-bind_rows(pop_intercept_sum, pop_slope, sd_intercept, sd_slope, phi) %>% select(Par, y, ymin, ymax, Rhat, ESS, .interval)
S4tab
#write_csv(S4tab, "Code/Results/GLMMTable_results.csv")

# Plots -------------------------------------------------------------------
library(ggthemes)
library(gghighlight)
library(png)


#colors

pal<-c("#3696ad", "#3494c5","#5e8ed5","#9482d6", "#c670c3", "#eb5e9f", "#fc576f", "#f56638", "#ffbd59")#add extra color for sperm whales

SpList<-c("Antarctic blue", "Sperm", "Pygmy blue","Fin","Humpback", "Sei",  "Brydes", "Antarctic minke", "Common minke")
col <-setNames(pal, SpList)


#get_variables(Mod1)

#Slope
species_slopes<-Mod1 %>% spread_draws(b_binned_center_sc, r_SpName[c,t]) %>% mutate(Species = gsub("\\.", replacement = " ", c)) %>%
  filter(t == "binned_center_sc") %>% mutate(Sp_slope = b_binned_center_sc + r_SpName)

pop_slope<-spread_draws(Mod1, b_binned_center_sc) %>% mutate(Species = "Group")

#annotating with probability of a negative slope
facet_names<-c("Antarctic blue", "Pygmy blue", "Fin",   
               "Humpback","Sei", "Brydes", "Antarctic minke",
               "Common minke")
neg_annots<-tibble(Species = facet_names) %>% left_join(species_neg, by = "Species") %>% select(prob_neg) %>% unlist()
neg_annots<-c(pop_slope_neg$prob_neg, neg_annots)
neg_annots<-paste("P(neg.) = ",  sprintf("%.3f", neg_annots))
names(neg_annots)<-c("Group", facet_names)


#plotting
GroupSlopes_allRE<-
  ggplot(species_slopes) +
  annotate("rect", xmin = -0.06, xmax = 0.042, ymin = "Group", ymax = "Antarctic blue", fill = "gray", alpha = 0.5)+
  stat_halfeye(aes(y = Species, x = Sp_slope, fill = Species), point_interval = mode_hdi, scale = 0.8, alpha = 0.8, normalize = "panels") + 
  stat_halfeye(data = pop_slope, aes(x = b_binned_center_sc, y = Species), point_interval = mode_hdi, scale= 0.9, fill = "gray20", interval_color = "red", point_color = "red") + 
  theme_classic() +  
  labs(x = "Slope", y = "Density") + 
  geom_vline(aes(xintercept = 0.0), linetype = "dashed") +
  scale_fill_manual(values = col, drop = TRUE) + 
  scale_y_discrete(expand = c(0.002, 0), limits = c("Group", "Antarctic blue", "Pygmy blue", "Fin",   
                                                    "Humpback","Sei", "Brydes", "Antarctic minke",
                                                    "Common minke"), 
                   labels = c("Group", "Antarctic blue", "Pygmy blue", "Fin",   
                              "Humpback","Sei", "Brydes", "Antarctic minke",
                              "Common minke")) + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), 
        axis.title.y = element_text(size = 14, vjust = 3),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 11), 
        axis.text.y = element_text(hjust = 0, vjust = -4.6, margin = margin(l = 0, r = -105), size = 14)) + 
  scale_x_continuous(limits = c(-0.06, 0.042), breaks = seq(-0.08, 0.08, by = 0.02), expand = c(0,0))

#adding silhouettes
#NOTE: pictures can be found at https://phylopic.org/
#Not included in repo, so this part will not run

#ut<-(-0.026--0.061)/24

#blue<-grid::rasterGrob(readPNG("Figures/Chris Huh Phylopic/blue.png", native = TRUE))
#fin<-grid::rasterGrob(readPNG("Figures/Chris Huh Phylopic/Fin.png", native = TRUE))
#hump<-grid::rasterGrob(readPNG("Figures/Chris Huh Phylopic/humpback.png", native = TRUE))
#sei<-grid::rasterGrob(readPNG("Figures/Chris Huh Phylopic/sei.png", native = TRUE))
#minke<-grid::rasterGrob(readPNG("Figures/Chris Huh Phylopic/minke.png", native = TRUE))
#brydes<-grid::rasterGrob(readPNG("Figures/Chris Huh Phylopic/Brydes.png", native = TRUE))

#GroupSlopes_allRE_annot<-GroupSlopes_allRE +  annotation_custom(blue, xmin = -0.061, xmax = -0.026,ymin = 0.5, ymax = 4) +
  #annotation_custom(blue, xmin = -0.061, xmax = -0.030,ymin = 1.5, ymax =5) +
  #annotation_custom(fin, xmin = -0.061, xmax = -0.032,ymin = 4, ymax = 6.5) +
  #annotation_custom(hump, xmin = -0.061, xmax = -0.041,ymin = 3.8, ymax = 7.9) + 
  #annotation_custom(sei, xmin = -0.061, xmax = -0.042,ymin = 6.1, ymax = 7.9) +
  #annotation_custom(minke, xmin = -0.060, xmax = -0.051,ymin = 8.5, ymax = 11) + 
  #annotation_custom(minke, xmin = -0.060, xmax = -0.050,ymin = 7.5, ymax = 10) + 
  #annotation_custom(brydes, xmin = -0.061, xmax = -0.042,ymin = 7.1, ymax =8.8) + 
  #annotate("text", x = rep(0.03, 9), y = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5), label = neg_annots, parse = FALSE, size = 4) + 
  #theme(axis.title = element_text(size = 14), 
        #axis.text.x = element_text(size = 14))

#GroupSlopes_allRE_annot
#ggsave("Figures/SlopeResults_Silh.png", dpi = 900, width = 7, height = 7, units = "in")



#Intercepts
species_intercepts<-Mod1 %>% spread_draws(b_Intercept, r_SpName[c,t]) %>% mutate(Species = gsub("\\.", replacement = " ", c)) %>%
  filter(t == "Intercept") %>% mutate(Sp_Intercept_logit = (b_Intercept + r_SpName), Sp_Intercept_real = plogis(Sp_Intercept_logit))


pop_intercept<-spread_draws(Mod1, b_Intercept) %>% mutate(b_Intercept_real = plogis(b_Intercept), Species = "Group")

#plotting intercepts
GroupIntercepts<-
  ggplot(species_intercepts) +
  annotate("rect", xmin = 0.475, xmax = 0.530, ymin = "Group", ymax = "Antarctic blue", fill = "gray", alpha = 0.5)+
  stat_halfeye(aes(y = Species, x = Sp_Intercept_real, fill = Species), point_interval = mode_hdi, scale = 0.8, alpha = 0.8, normalize = "panels") + 
  stat_halfeye(data = pop_intercept, aes(x = b_Intercept_real, y = Species), point_interval = mode_hdi, scale= 0.9, fill = "gray20", interval_color = "red", point_color = "red") + 
  theme_classic() +  
  labs(x = "Intercept", y = "Species") + 
  geom_vline(aes(xintercept = 0.5), linetype = "dashed") +
  scale_fill_manual(values = col, drop = TRUE) + 
  scale_y_discrete(expand = c(0.002, 0),  limits = c("Group", "Antarctic blue", "Pygmy blue", "Fin",   
                                                     "Humpback","Sei", "Brydes", "Antarctic minke",
                                                     "Common minke"), 
                   labels = c("Group", "Antarctic blue", "Pygmy blue", "Fin",   
                              "Humpback","Sei", "Brydes", "Antarctic minke",
                              "Common minke")) + 
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(0.475, 0.530), breaks = seq(0.47, 0.52, by = 0.01), expand = c(0,0)) 


#combining slopes and intercepts
GroupIntercepts
lay<-"
AAAAAABBBBB
AAAAAABBBBB
AAAAAABBBBB"

slop_int<-GroupSlopes_allRE_annot + GroupIntercepts + plot_layout(design = lay)
slop_int
#ggsave("Figures/Slope_Intercepts.png", dpi = 900, width = 10, height = 8, units = "in")

#linear plot
species_vals <- Mod1 %>% spread_draws(b_Intercept, b_binned_center_sc, r_SpName[c,t]) %>% pivot_wider(names_from = "t", values_from = "r_SpName") %>% mutate(Species = gsub("\\.", replacement = " ", c))

#function to get predictions
get_pred<-function(species_val){
  Sp_lengths<-seq(-2, 2, by = 0.25)
  sr_pred_logit<-species_val["b_Intercept"] + species_val["Intercept"] + (species_val["b_binned_center_sc"] + species_val["binned_center_sc"])*Sp_lengths
  #return(sr_pred_logit)
  return(plogis(sr_pred_logit))
}

#creating dataframe of predictions
species_preds<-apply(species_vals[,c("b_Intercept", "b_binned_center_sc", "Intercept", "binned_center_sc")], 1, get_pred)
id<-length(seq(-2,2, by = 0.25))

pred_tibble_all<-tibble(draw = rep(species_vals$`.draw`, each = id), sp = rep(species_vals$Species, each = id), 
                        mat_length = rep(seq(-2, 2, by = 0.25), nrow(species_vals)), ypred = as.vector(species_preds))


pred_tibble_summ<-pred_tibble_all %>% group_by(sp, mat_length) %>% mode_hdci(ypred) %>% rename("mode" = ypred, "2.5%" = `.lower`, "97.5%" = `.upper`)


#add sample size to facet labels
facet_names<-c("Common minke", "Antarctic minke", "Brydes", 
              "Sei", "Humpback", "Fin", "Pygmy blue", "Antarctic blue")
n_fetus<-tibble(SpName = facet_names) %>% left_join(numFetus, by = "SpName") %>% select(`sum(Total)`) %>% unlist()
facet_labs<-paste0(facet_names, "\n n=", formatC(n_fetus, big.mark = ","))
names(facet_labs)<-facet_names

#plotting
linesplot<-ggplot(pred_tibble_summ) +
  geom_hline(aes(yintercept = 0.5), color = 'black')+
  geom_line(aes(x = mat_length, y = mode, color = sp), alpha = 1, linewidth = 0.6) + 
  geom_ribbon(aes(x = mat_length, ymin = `2.5%`, ymax = `97.5%`, color = sp), alpha = 0, linetype = "dashed", linewidth = 0.6) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.46, 0.54)) + 
  scale_x_continuous(expand = c(0, 0)) +
  gghighlight(use_direct_label = FALSE, unhighlighted_params = list(linewidth = 0.1, color = "gray60")) +
  scale_color_manual(values = col, drop = TRUE) + 
  labs(x = "Maternal length (scaled)", y = "Fetal sex ratio") + 
  facet_wrap(~factor(sp, levels = c("Common minke", "Antarctic minke", "Brydes", 
                                    "Sei", "Humpback", "Fin", "Pygmy blue", "Antarctic blue")), 
             labeller = as_labeller(facet_labs), nrow = 4) + 
  theme_tufte(base_size = 12, base_family = "Arial") + 
  theme(legend.position = "none", 
        #text = element_text(family = "Arial", size = 12), 
        strip.text = element_text(color = "black", size = 11, vjust = -1), 
        panel.grid.major.y = element_blank(), 
        panel.spacing = unit(0.5, "lines"))


linesplot

#ggsave("Figures/predictions_species.png", linesplot, width = 7, height = 8, units = "in", dpi = 600)





