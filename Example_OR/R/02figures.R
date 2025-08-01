###############################
# Reproduce Figures 1 and 2 ###
###############################

source("./Example_OR/R/01comparisonstudy.R")

library(viridis)
library(cowplot)



###############
### Table 5 ###
###############

# check proportions of contingency tables with zero cells
lapply(X = datXY_sim, FUN = count.zerotables)


######################
# Figure 1 (A and B) #
######################

##############
# Figure 1 A # 
##############


add_labels_xaxis <- paste("Scenario ", 1:4)

Figure1_A <- ggplot(data = bias.long,
                    aes(x = interaction(params, deal.NA, lex.order = TRUE),
                        y = rank, group = 1, label = est.method)) +
  geom_line(aes(group=as.factor(paste0(params, "-", est.method)), linetype = can.deal, col = est.method), linewidth = 0.65)+
  geom_label(aes(group=as.factor(paste0(params, "-", est.method)), fill = est.method, fontface = "bold"), col = "white",size = 4.5) +
  # geom_text(aes(group=as.factor(paste0(params, "-", est.method, size = 12)), col = est.method)) + 
  scale_y_continuous(breaks = 5:1, transform = "reverse") +
  scale_x_discrete(labels= rep(c("Single", "All"),
                               times = 4)) +
  # facet_wrap_paginate(~ page, ncol=1, nrow = 1, page = 1)
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 14),
        axis.text.y=element_text(size = 14),
        plot.margin = margin(t=0, b =1.5, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
        axis.title.y = element_text(size =16), 
        axis.title.x= element_text(size = 16, vjust = -6)) +
  # # Add the tool names to the plot:
  annotate(geom = "text",
           x = 1.5 + 2*(0:(length(unique(bias.long$params))-1)),
           y = 6,
           label = add_labels_xaxis, size = 5) +
  coord_cartesian(ylim=c(5,1),clip = "off") + # clip = "off" required to add scenario number below the plot
  xlab("Simulation scenario") +
  ylab("Rank (bias)") + 
  labs(colour='OR estimation method', linetype = "Fails with sampling zeros ")  + 
  # theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  guides(fill="none", col = "none") +
  # scale_colour_hue(l = 60) + 
  theme(legend.position = "top",
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_color_viridis(discrete = TRUE) +
scale_fill_viridis(discrete = TRUE)

plot_grid(Figure1_A, Figure1_B, labels = c('A', 'B'), ncol = 2, align ="h")


##############
# Figure 1 B # 
##############

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

Figure1_B <- ggplot(data = bias_fallback.long,
                    aes(x = params, y = rank, label = est.method)) +
  # geom_point(aes(x = params, y = rank, col = est.method)) +
  geom_label(aes(x = params, y = rank, fill = est.method, fontface = "bold"), col = "white",size = 4.5) +
  scale_y_continuous(breaks = 9:1, transform = "reverse") +
  scale_x_discrete(labels= paste0("Scenario ", 1:4)) +
  theme(axis.text.x=element_text(angle = 50, vjust = 1, hjust = 1, size = 14),
        axis.text.y=element_text(size = 14),
        plot.margin = margin(t=1, b =2, l=1, r=1, unit="cm"),  ## add space below the actual plot (needed for the GSA tool names)
        axis.title.y = element_text(size =16), 
        axis.title.x= element_text(size = 16, vjust = -6)) +
  coord_cartesian(ylim=c(9,1),clip = "off") + # clip = "off" required to add scenario number below the plot
  xlab("Simulation scenario") +
  ylab("Rank (bias)") + 
  theme(legend.position="none", panel.grid.minor = element_blank()) +
  # scale_colour_hue(l = 90)  +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)


# combine to obtain Figure 1 
plot_grid(Figure1_A, Figure1_B, labels = c('A', 'B'), ncol = 2)
# vertical alignment 
plot_grid(Figure1_A, Figure1_B, labels = c('A', 'B'), ncol = 2, align ="h")


# uncomment to save the figure: 

ggsave2("./Example_OR/Figures/Figure1.pdf",
        width = 13,
        height = 7)




