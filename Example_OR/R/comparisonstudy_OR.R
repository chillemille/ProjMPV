################################################################################
### Analyze OR data and comparison study  ######################################
################################################################################

library(epitools)
library(PropCIs) # for oscoreci
# small provides estimate when there are 0-cells, 
# Wald, Fisher, midp do not 
library(pairwiseCI)
library(tidyr)
library(dplyr)
library(ggplot2)


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
### Check propprtion of contingency table with zero cells in each scenario #####
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

# get bias (on log scale) 
bias <- lapply(X = est.OR, FUN = function(est.OR_df) compute.bias_logscale(est.OR_df$estimates, est.OR_df$true.OR)) %>%
                     do.call(what = rbind) %>% mutate(params = as.factor(params), deal.NA = as.factor(deal.NA))

bias.absolute <- lapply(X = est.OR, FUN = function(est.OR_df) compute.bias_logscale(est.OR_df$estimates, est.OR_df$true.OR, absolute = TRUE)) %>%
                    do.call(what = rbind) %>% mutate(params = as.factor(params), deal.NA = as.factor(deal.NA))



###############################
# rank by bias (on log scale) # 
###############################

# make longer data frame 
bias.long <- pivot_longer(bias, cols= c("fisher", "midp", "small", "Woolf", "manual"), 
                                       names_to = "est.method", values_to = "bias") %>% as.data.frame()

# add rank (for each parameter setting and aca. vs. cca)
bias.long<- bias.long %>% group_by(params, deal.NA) %>% mutate(rank = rank(abs(bias)))

# add column indicating whether method can deal with sampling zeros or not 
bias.long$can.deal <- as.factor(ifelse(bias.long$est.method %in% c("Woolf", "small", "manual"), TRUE, FALSE))

# add unique ID to each column 



######################
# make plot on ranks #
######################

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

add_labels_xaxis <- paste("Scen. ", 1:8)

plot.rank <- ggplot(data = bias.long,
               aes(x = interaction(params, deal.NA, lex.order = TRUE),
                   y = rank, group = 1)) +
  geom_line(aes(group=as.factor(paste0(params, "-", est.method)),  col = est.method))+
  geom_point(aes(group=as.factor(paste0(params, "-", est.method, size = 12)), col = est.method, shape = est.method)) +
  scale_y_continuous(breaks = 5:1, transform = "reverse") +
  scale_x_discrete(labels= rep(c("Discard 1", "Discard 2"),
                               times = 8)) +
  # facet_wrap_paginate(~ page, ncol=1, nrow = 1, page = 1)
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 11),
        plot.margin = margin(t=1, b =2, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
        axis.title.y = element_text(size =14), 
        axis.title.x= element_text(size = 13, vjust = -14)) +
  # # Add the tool names to the plot:
  annotate(geom = "text",
           x = 1.5 + 2*(0:(length(unique(bias.long$params))-1)),
           y = 5.7,
           label = add_labels_xaxis, size = 4) +
  coord_cartesian(ylim=c(5,1),clip = "off") + # clip = "off" required to add scenario number below the plot
  xlab("Simulation scenarios") +
 ylab("Rank") + 
  labs(color='OR estimation method')  + 
  #theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  scale_fill_manual(values=cbPalette) + 
  scale_shape_manual(values = c(15, 17, 18, 15, 11))


plot.rank



# make plot on raw (but absolute) "bias values" 
# plot.raw <- ggplot(data = bias.long,
#                     aes(x = interaction(params, deal.NA, lex.order = TRUE),
#                         y = abs(bias), group = 1)) +
#   geom_line(aes(group=as.factor(paste0(params, "-", est.method)),  col = est.method))+
#   geom_point(aes(group=as.factor(paste0(params, "-", est.method)), col = est.method, shape = est.method)) +
#   scale_x_discrete(labels= rep(c("Discard for failing", "Discard for all"),
#                                times = 8)) +
#   # facet_wrap_paginate(~ page, ncol=1, nrow = 1, page = 1)
#   theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 11),
#         plot.margin = margin(t=1, b =2, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
#         axis.title.y = element_text(size =14), 
#         axis.title.x= element_text(size = 13, vjust = -14)) +
#   # # Add the tool names to the plot:
#   annotate(geom = "text",
#            x = 1.5 + 2*(0:(length(unique(bias.long$params))-1)),
#            y = -0.13,
#            label = add_labels_xaxis, size = 4) +
#   coord_cartesian(ylim=c(0, 0.6),clip = "off") + # clip = "off" required to add scenario number below the plot
#   xlab("Simulation scenarios") +
#   ylab("Absolute bias") + 
#   labs(color='OR estimation \n method') 


################################################################################
### repeat analysis using fallback principle ###################################
################################################################################


est.OR_fallback <- lapply(X = datXY_sim, FUN = est.OR_df, haldane.correct = TRUE)

save(est.OR_fallback, file = "~/MPVExamples/MPVExamples/Example_OR/Data/Simulated_XYData/est.OR_fallback.RData")
# load("~/MPVExamples/MPVExamples/Example_OR/Data/Simulated_XYData/est.OR_fallback.RData")

# calculate bias
bias.fallback <- lapply(X = est.OR_fallback, FUN = function(est.OR_df) compute.bias_logscale(est.OR_df$estimates, est.OR_df$true.OR)) %>%
  do.call(what = rbind) %>% mutate(params = as.factor(params), deal.NA = as.factor(deal.NA))



plot.rank_fallback <- ggplot(data = bias.long,
                    aes(x = interaction(params, deal.NA, lex.order = TRUE),
                        y = rank, group = 1)) +
  geom_line(aes(group=as.factor(paste0(params, "-", est.method)),  col = est.method))+
  geom_point(aes(group=as.factor(paste0(params, "-", est.method)), col = est.method)) +
  scale_y_continuous(breaks = 1:6) +
  scale_x_discrete(labels= rep(c("ACA", "CCA"),
                               times = 8)) +
  # facet_wrap_paginate(~ page, ncol=1, nrow = 1, page = 1)
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 11),
        plot.margin = margin(t=1, b =2, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
        axis.title.y = element_text(size =14), 
        axis.title.x= element_text(size = 13, vjust = -14)) +
  # # Add the tool names to the plot:
  annotate(geom = "text",
           x = 1.5 + 2*(0:(length(unique(bias.long$params))-1)),
           y = 0,
           label = add_labels_xaxis, size = 4) +
  coord_cartesian(ylim=c(1,5),clip = "off") + # clip = "off" required to add scenario number below the plot
  xlab("Simulation scenarios") +
  ylab("Rank") + 
  labs(color='OR estimation method')  
## uncomment to save

# make plot on raw (but absolute) "bias values" 
plot.raw <- ggplot(data = bias.long,
                   aes(x = interaction(params, deal.NA, lex.order = TRUE),
                       y = abs(bias), group = 1)) +
  geom_line(aes(group=as.factor(paste0(params, "-", est.method)),  col = est.method))+
  geom_point(aes(group=as.factor(paste0(params, "-", est.method)), col = est.method)) +
  scale_x_discrete(labels= rep(c("ACA", "CCA"),
                               times = 8)) +
  # facet_wrap_paginate(~ page, ncol=1, nrow = 1, page = 1)
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 11),
        plot.margin = margin(t=1, b =2, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
        axis.title.y = element_text(size =14), 
        axis.title.x= element_text(size = 13, vjust = -14)) +
  # # Add the tool names to the plot:
  annotate(geom = "text",
           x = 1.5 + 2*(0:(length(unique(bias.long$params))-1)),
           y = -0.13,
           label = add_labels_xaxis, size = 4) +
  coord_cartesian(ylim=c(0, 0.6),clip = "off") + # clip = "off" required to add scenario number below the plot
  xlab("Simulation scenarios") +
  ylab("Absolute bias") + 
  labs(color='OR estimation \n method') 


