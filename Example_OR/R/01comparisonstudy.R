################################################################################
### Analyze OR data and comparison study  ######################################
################################################################################

library(epitools)
# small provides estimate when there are 0-cells, 
# Wald, Fisher, midp do not 
library(pairwiseCI)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(forcats)

source("./Example_OR/R/help_functions.R")
dir.create("./Example_OR/Figures")
dir.create("./Example_OR/Data")

# set simulation parameters as specified in manuscript
OR <- c(2, 3)
n.obs <- 50
p_x <- c(0.25, 0.5)


scenarios <- expand_grid(n.obs, p_x, OR) 

# generate XY data for each of the eight scenarios
datXY_sim <- apply(X = scenarios ,MARGIN = 1, function(scen){
    sim_OR.data(n.sim = 100000, n.obs = scen[1], p_x = scen[2], 
                OR = scen[3],  beta0 = 0, seed = 1)
  }) 


# save results 
# save(datXY_sim, file = "./Example_OR/Data/Simulated_XYData/datXY_sim_allscen.RData")

# load results 
#load("./Example_OR/Data/Simulated_XYData/datXY_sim_allscen.RData")


################################################################################
### Comparison study: discard for all. vs. failing methods only ################
################################################################################

# Note: abbreviations
# - "ACA" stands for discarding data sets with sampling zeros for the failing 
#   methods only
# - "CCA" stands for discarding data sets with sampling zeros for ALL methods 

######################
# get OR estimations #
######################

est.OR <- lapply(X = datXY_sim, FUN = est.OR_df, haldane.correct = FALSE)



# save files
#save(est.OR, file = "./Example_OR/Data/Simulated_XYData/est.OR.RData")
# load 
#load("./Example_OR/Data/Simulated_XYData/est.OR.RData")

##########
# Biases #
##########

# for bias at this point, remove Haldane correction for manual 
for(i in 1:4){
  
  est.OR[[i]]$estimates$manual <- ifelse(est.OR[[i]]$estimates$fisher == 0, 0, 
                                         ifelse(is.finite(est.OR[[i]]$estimates$fisher), est.OR[[i]]$estimates$manual, Inf))
  
  
}


# get bias (on log scale) 
bias <- lapply(X = est.OR, FUN = function(est.OR_df) compute.bias_logscale(est.OR_df$estimates, est.OR_df$true.OR)) %>%
                     do.call(what = rbind) %>% mutate(params = as.factor(params), deal.NA = as.factor(deal.NA))



##########################################
# rank by bias (on log scale); raw scale # 
##########################################

# make longer data frame 
bias.long <- pivot_longer(bias, cols= c("fisher", "midp", "small", "Woolf", "manual"), 
                                       names_to = "est.method", values_to = "bias") %>% as.data.frame()

# add rank (for each parameter setting and aca. vs. cca)
bias.long<- bias.long %>% group_by(params, deal.NA) %>% mutate(rank = rank(abs(bias)))

# add column indicating whether method can deal with sampling zeros or not 
bias.long$can.deal <- as.factor(ifelse(!bias.long$est.method %in% c("Woolf", "small"), TRUE, FALSE))


bias.long$est.method <- fct_recode(bias.long$est.method,
                                            "Fis" = "fisher",
                                            "Mid" = "midp",  
                                            "Sma" = "small", 
                                            "Woo" = "Woolf",  
                                            "Man" = "manual") 


##########################################
### Comparison studies using fallbacks ###
##########################################


est.OR_fallback <- lapply(X = datXY_sim, FUN = est.OR_df, haldane.correct = TRUE)


#save(est.OR_fallback, file = "./Example_OR/Data/Simulated_XYData/est.OR_fallback_neu.RData")
# load("./Example_OR/Data/Simulated_XYData/est.OR_fallback.RData")

# calculate bias
bias.fallback <- lapply(X = est.OR_fallback, FUN = function(est.OR_df) compute.bias_logscale(est.OR_df$estimates, est.OR_df$true.OR)) %>%
  do.call(what = rbind) %>% 
  mutate(params = as.factor(params), deal.NA = as.factor(deal.NA)) %>%
  filter(deal.NA == "ACA") %>% subset(select = -deal.NA)







