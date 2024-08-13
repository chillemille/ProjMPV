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

source("./Examples/Example_OR/R/sim_ORdata.R")

# Idee: ORs bei realistischeren Werten belassen, dafür aber n.obs und p_x variieren, 
# um MPVs zu erhalten

# n.obs 50 -> auch die Sample Size fix lassen, da die Auswirkungen auf die 
# Genauigkeit der Schätzung hat 

OR <- 2:5
n.obs <- 50
p_x <- c(0.25, 0.5)


scenarios <- expand_grid(n.obs, p_x, OR) 

# generate XY data for each of the eight scenarios
datXY_sim <- apply(X = scenarios ,MARGIN = 1, function(scen){
    sim_OR.data(n.sim = 100000, n.obs = scen[1], p_x = scen[2], 
                OR = scen[3],  beta0 = 0, seed = 1)
  }) 


# save results 
# save(datXY_sim, file = "~/Examples/Example_OR/Data/Simulated_XYData/datXY_sim_allscen.RData")

# load results 
load("~/MPVExamples/MPVExamples/Example_OR/Data/Simulated_XYData/datXY_sim_allscen.RData")
################################################################################
### Check proportions of contingency tables with zero cells in each scenario ###
################################################################################

lapply(X = datXY_sim, FUN = count.zerotables)

################################################################################
### Run methods on simulated data: Compare ACA vs. CCA #########################
################################################################################

##############################
# get OR estimations and CIs #
##############################

est.OR <- lapply(X = datXY_sim, FUN = est.OR_df, haldane.correct = FALSE)



# save files
save(est.OR, file = "~/Examples/Example_OR/Data/Simulated_XYData/est.OR.RData")
# load 
load("~/MPVExamples/MPVExamples/Example_OR/Data/Simulated_XYData/est.OR.RData")

##########
# Biases #
##########

# for bias at this point, remove Haldane correction for manual 
for(i in 1:8){
  
  est.OR[[i]]$estimates$manual <- ifelse(est.OR[[i]]$estimates$fisher == 0, 0, 
                                         ifelse(is.finite(est.OR[[i]]$estimates$fisher), est.OR[[i]]$estimates$manual, Inf))
  
  
}


# get bias (on log scale) 
bias <- lapply(X = est.OR, FUN = function(est.OR_df) compute.bias_logscale(est.OR_df$estimates, est.OR_df$true.OR)) %>%
                     do.call(what = rbind) %>% mutate(params = as.factor(params), deal.NA = as.factor(deal.NA))

MSE <- lapply(X = est.OR, FUN = function(est.OR_df) compute.bias_logscale(est.OR_df$estimates, est.OR_df$true.OR, absolute = TRUE)) %>%
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





######################
# make plot on ranks #
######################

add_labels_xaxis <- paste("Scenario ", 1:8)

plot.rank_bias <- ggplot(data = bias.long,
               aes(x = interaction(params, deal.NA, lex.order = TRUE),
                   y = rank, group = 1, label = est.method)) +
  geom_line(aes(group=as.factor(paste0(params, "-", est.method)), linetype = can.deal, col = est.method), size = 0.65)+
  geom_label(aes(group=as.factor(paste0(params, "-", est.method)), fill = est.method, fontface = "bold"), col = "white",size = 3.5) +
  # geom_text(aes(group=as.factor(paste0(params, "-", est.method, size = 12)), col = est.method)) + 
  scale_y_continuous(breaks = 5:1, transform = "reverse") +
  scale_x_discrete(labels= rep(c("Failing", "All"),
                               times = 8)) +
  # facet_wrap_paginate(~ page, ncol=1, nrow = 1, page = 1)
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 11),
        axis.text.y=element_text(size = 12),
        plot.margin = margin(t=0, b =1.5, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
        axis.title.y = element_text(size =14), 
        axis.title.x= element_text(size = 12, vjust = -14)) +
  # # Add the tool names to the plot:
  annotate(geom = "text",
           x = 1.5 + 2*(0:(length(unique(bias.long$params))-1)),
           y = 6,
           label = add_labels_xaxis, size = 4) +
  coord_cartesian(ylim=c(5,1),clip = "off") + # clip = "off" required to add scenario number below the plot
  xlab("Simulation scenario") +
 ylab("Rank (bias)") + 
  labs(colour='OR estimation method', linetype = "Fails with sampling zeros ")  + 
  # theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
guides(fill="none", col = "none") +
  # scale_colour_hue(l = 60) + 
  theme(legend.position = "top")


plot.rank_bias

##############################
# rank by mse (on log scale) # 
##############################

# make longer data frame 
mse.long <- pivot_longer(MSE, cols= c("fisher", "midp", "small", "Woolf", "manual"), 
                          names_to = "est.method", values_to = "bias") %>% as.data.frame()

# add rank (for each parameter setting and aca. vs. cca)
mse.long<- mse.long %>% group_by(params, deal.NA) %>% mutate(rank = rank(bias))

# add column indicating whether method can deal with sampling zeros or not 
mse.long$can.deal <- as.factor(ifelse(!mse.long$est.method %in% c("Woolf", "small"), TRUE, FALSE))


mse.long$est.method <- fct_recode(mse.long$est.method,
                                   "Fis" = "fisher",
                                   "Mid" = "midp",  
                                   "Sma" = "small", 
                                   "Woo" = "Woolf",  
                                   "Man" = "manual") 





######################
# make plot on ranks #
######################

add_labels_xaxis <- paste("Scenario ", 1:8)

plot.rank_mse <- ggplot(data = mse.long,
                        aes(x = interaction(params, deal.NA, lex.order = TRUE),
                            y = rank, group = 1, label = est.method)) +
  geom_line(aes(group=as.factor(paste0(params, "-", est.method)), linetype = can.deal, col = est.method), size = 0.65)+
  geom_label(aes(group=as.factor(paste0(params, "-", est.method)), fill = est.method, fontface = "bold"), col = "white",size = 3.5) +
  # geom_text(aes(group=as.factor(paste0(params, "-", est.method, size = 12)), col = est.method)) + 
  scale_y_continuous(breaks = 5:1, transform = "reverse") +
  scale_x_discrete(labels= rep(c("Failing", "All"),
                               times = 8)) +
  # facet_wrap_paginate(~ page, ncol=1, nrow = 1, page = 1)
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 11),
        axis.text.y=element_text(size = 12),
        plot.margin = margin(t=1, b =1.5, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
        axis.title.y = element_text(size =14), 
        axis.title.x= element_text(size = 12, vjust = -12)) +
  # # Add the tool names to the plot:
  annotate(geom = "text",
           x = 1.5 + 2*(0:(length(unique(mse.long$params))-1)),
           y = 6,
           label = add_labels_xaxis, size = 4) +
  coord_cartesian(ylim=c(5,1),clip = "off") + # clip = "off" required to add scenario number below the plot
   xlab("Simulation scenario") +
  ylab("Rank (MSE)") + 
  labs(colour='OR estimation method', linetype = "Fails with sampling zeros")  + 
  # theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  guides(fill="none", col = "none") +
  # scale_colour_hue(l = 60) + 
  theme(legend.position = "top")




plot.rank_mse




plot_grid(plot.rank_bias, plot.rank_mse, labels = c('A', 'B'), ncol = 1)
ggsave2("./Example_OR/Figures/Figure2.pdf",
        width = 11,
        height = 13)

################################################################################
### repeat analysis using fallback principle ###################################
################################################################################


est.OR_fallback <- lapply(X = datXY_sim, FUN = est.OR_df, haldane.correct = TRUE)


save(est.OR_fallback, file = "~/MPVExamples/MPVExamples/Example_OR/Data/Simulated_XYData/est.OR_fallback.RData")
# load("~/MPVExamples/MPVExamples/Example_OR/Data/Simulated_XYData/est.OR_fallback.RData")

# calculate bias
bias.fallback <- lapply(X = est.OR_fallback, FUN = function(est.OR_df) compute.bias_logscale(est.OR_df$estimates, est.OR_df$true.OR)) %>%
  do.call(what = rbind) %>% 
  mutate(params = as.factor(params), deal.NA = as.factor(deal.NA)) %>%
  filter(deal.NA == "ACA") %>% subset(select = -deal.NA)

mse.fallback <- lapply(X = est.OR_fallback, FUN = function(est.OR_df) compute.bias_logscale(est.OR_df$estimates, est.OR_df$true.OR, absolute = TRUE)) %>%
  do.call(what = rbind) %>% 
  mutate(params = as.factor(params), deal.NA = as.factor(deal.NA)) %>%
  filter(deal.NA == "ACA") %>% subset(select = -deal.NA)


#########################################
# rank by bias (on log scale); fallback # 
#########################################

cols.fallback <- names(bias.fallback)[-length(names(bias.fallback))]

# make longer data frame 
bias_fallback.long <- pivot_longer(bias.fallback, cols= cols.fallback, 
                          names_to = "est.method", values_to = "bias") %>% as.data.frame()

# add rank (for each parameter setting and aca. vs. cca)
bias_fallback.long<- bias_fallback.long %>% group_by(params) %>% mutate(rank = rank(abs(bias)))


bias_fallback.long$est.method <- fct_recode(bias_fallback.long$est.method,
                                 "Fis/+0.5" = "fisher",
                                 "Mid/+0.5" = "midp",  
                                 "Sma" = "small", 
                                 "Woo" = "Woolf",  
                                 "Man/+0.5" = "manual",
                                 "Mid/Sma" = "midp_fallback.small", 
                                 "Mid/Woo" = "midp_fallback.Woolf", 
                                 "Fis/Sma" = "fisher_fallback.small", 
                                 "Fis/Woo" = "fisher_fallback.Woolf") 

plot.bias_fallback <- ggplot(data = bias_fallback.long,
                    aes(x = params, y = rank, label = est.method)) +
  # geom_point(aes(x = params, y = rank, col = est.method)) +
  geom_label(aes(x = params, y = rank, fill = est.method, fontface = "bold"), col = "white",size = 3.5) +
  scale_y_continuous(breaks = 9:1, transform = "reverse") +
  scale_x_discrete(labels= paste0("Scen. ", 1:8)) +
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 12),
        axis.text.y=element_text(size = 12),
        plot.margin = margin(t=1, b =2, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
        axis.title.y = element_text(size =14), 
        axis.title.x= element_text(size = 13, vjust = -11)) +
  coord_cartesian(ylim=c(9,1),clip = "off") + # clip = "off" required to add scenario number below the plot
  xlab("Simulation scenario") +
  ylab("Rank (bias)") + 
  theme(legend.position="none", panel.grid.minor = element_blank()) +
  scale_colour_hue(l = 90) 
  
## uncomment to save

plot.bias_fallback



# make longer data frame 
mse_fallback.long <- pivot_longer(mse.fallback, cols= cols.fallback, 
                                   names_to = "est.method", values_to = "mse") %>% as.data.frame()

# add rank (for each parameter setting and aca. vs. cca)
mse_fallback.long<- mse_fallback.long %>% group_by(params) %>% mutate(rank = rank(abs(mse)))


mse_fallback.long$est.method <- fct_recode(mse_fallback.long$est.method,
                                            "Fis/+0.5" = "fisher",
                                            "Mid/+0.5" = "midp",  
                                            "Sma" = "small", 
                                            "Woo" = "Woolf",  
                                            "Man/+0.5" = "manual",
                                            "Mid/Sma" = "midp_fallback.small", 
                                            "Mid/Woo" = "midp_fallback.Woolf", 
                                            "Fis/Sma" = "fisher_fallback.small", 
                                            "Fis/Woo" = "fisher_fallback.Woolf") 

plot.mse_fallback <- ggplot(data = mse_fallback.long,
                            aes(x = params, y = rank, label = est.method)) +
  # geom_point(aes(x = params, y = rank, col = est.method)) +
  geom_label(aes(x = params, y = rank, fill = est.method, fontface = "bold"), col = "white",size = 3.5) +
  scale_y_continuous(breaks = 9:1, transform = "reverse") +
  scale_x_discrete(labels= paste0("Scen. ", 1:8)) +
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 12),#
        axis.text.y=element_text(size = 12),
        plot.margin = margin(t=1, b =2, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
        axis.title.y = element_text(size =14), 
        axis.title.x= element_text(size = 13, vjust = -11)) +
  coord_cartesian(ylim=c(9,1),clip = "off") + # clip = "off" required to add scenario number below the plot
  xlab("Simulation scenario") +
  ylab("Rank (MSE)") + 
  theme(legend.position="none", panel.grid.minor = element_blank()) +
  scale_colour_hue(l = 90) 
## uncomment to save

plot.mse_fallback


plot_grid(plot.bias_fallback, plot.mse_fallback, labels = c('A', 'B'), ncol = 1)
ggsave2("./Example_OR/Figures/Figure3.pdf",
        width = 11,
        height = 13)

# thoughts on fallback strategy 
# - might require additional considerations in practice: which combinations of 
# original and fallback make sense or are, which don't?
# -> for instance, users using a method from the epitools package (e.g. midp, fisher)
# are more likely to use small when sampling zeros occur

# -> some fallbacks might be quite difficult to apply in the first place
# -> less likely to be applied by users 