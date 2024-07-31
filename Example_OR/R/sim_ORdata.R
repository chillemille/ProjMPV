# note: the higher beta0, the higher the number of zeros in the contingency table 

library(dplyr) # for function between()



################################################################################
### Function to generate OR data ###############################################
################################################################################

# generate X and Y data with a pre-specified OR 
sim_OR.data <- function(n.sim, n.obs, p_x, OR, beta0 , seed){
  
  # n.sim = 1000
  # n.obs=200
  # p_x=0.7
  # OR = 4 
  # beta0 = 0.3
  # seed = 4 
  
  set.seed(seed)
  
  # empty list to be filled with the XY-data from the simulation
  datXY_list <- list()
  
  # store random states 
  statesdat_list <- list()
  
  for(i in 1:n.sim){
    
  # store random number generator state
  statesdat_list[[i]] = .Random.seed
  
  # beta1 (effect of X) can be directly derived from log(OR)
  beta1 <- log(OR)
  
  # initialise X
  X <- c()

  # sample n.obs binary values with the pre-specified p_xability p_x
  X <- sample(c(0,1), n.obs, replace=TRUE, 
                       prob= c(1-p_x, p_x))
  
  # for each X-value, calculate corresponding linear predictor 
  linpred<- beta0 + X * beta1
  # simulate Y-values using logistic regression
  Y <- rbinom(n.obs,size = 1, prob = plogis(linpred))
  
  # convert X and Y into factors 
  X <- factor(X, levels = c(0,1))
  Y <- factor(Y, levels = c(0,1))
  
  datXY_list[[i]] <- data.frame(X = X, Y = Y)
  
  }
  
  return(list(datXY_list = datXY_list, 
              true.OR = OR, 
              params = paste0(n.obs, "_",p_x, "_", OR)))

  
}

################################################################################
### Further help functions #####################################################
################################################################################

# generate contingency table FOR A SINGLE XY data set 
# perform haldane correction if needed (and specify value to add the cells, default is 0.5)
contingencytable <- function(datXY, haldane.correct = FALSE, add.bc = 0.5, add.ad = 0.5){
  
  # datXY <- datXY_sim1$datXY_list[[798]]
  
  # generate contingency table 
  contingencytable <- table(datXY)
  
  # Check if any of the cells in the contingency table contain 0, 
  # if so, add 0.5 to each of  the cells 
  # note: we make an exception for small (see below)
  if((haldane.correct == TRUE) && ((contingencytable[2] == 0) ||  (contingencytable[3] == 0))){
    
    contingencytable <- contingencytable + add.bc
    
    
    # note to myself: if the first condition is true (i.e. b = 0 or c = 0), then 0.5 is added 
    # to all entries in the contingency table 
    # The second, i.e. else-if, condition, is then evaluated on the "corrected" contingencytable in 
    # which all entries as > 0 
    # this means that the 0.5 correction is never performed twice 
  }
  
  if((haldane.correct == TRUE) && ((contingencytable[1] == 0) ||  (contingencytable[4] == 0))){
    
    contingencytable <- contingencytable + add.ad
    
  }

  return(contingencytable)
}

################################################################################
### indicate whether a contingency table contains sampling zeros ###############
################################################################################

ind.zerotables <- function(sim.data){
  
  ind.zerotables <- which(unlist(lapply(X = sim.data$datXY_list, FUN = function(datXY.list) !all(contingencytable(datXY.list) != 0))))
  
  # return
  ind.zerotables
  
}


################################################################################
### count proportion of zero tables from simulation run ########################
################################################################################

count.zerotables <- function(sim.data){
  
  prop.zerotables <- mean(sapply(X=sim.data$datXY_list, FUN = function(datXY.list) !all(contingencytable(datXY.list) !=0)))
  
  return(prop.zerotables)
  
}

  

################################################################################
### Obtain the OR estimations and CIs for all simulation iterations ############
################################################################################

est.OR_df <- function(simXY.data_list, haldane.correct = FALSE){
  
  # simXY.data_list <- datXY_sim[[1]]
  # haldane.correct = TRUE

  # specify which methods are to be included in the comparison  
  methods <- c("fisher", "midp", "small", "Woolf", "manual") 
  
  # empty data frame to be filled with the OR estimations
  OR.est <- data.frame(matrix(NA, ncol = length(methods), nrow =length(simXY.data_list$datXY_list)))
  names(OR.est) <- methods
  
  # empty data frame to be filled with lower and upper bounds of the CIs 
  #CI.est <- data.frame(matrix(NA, ncol = 2*length(methods), nrow=length(simXY.data_list$datXY_list))) 
  #names(CI.est) <- paste0(rep(methods, each = 2), ".", rep(c("lower", "upper"), times = length(methods)))
  
  
  for(i in which(!1:length(methods) %in% grep("manual", methods))){
    
    # Wald does not round off so that we can use 0.5 as offset 
    # Small can handle b=0 and c=0 internally such that we do not need an offset
    # however: when a =0 or d=0, it outputs 0 as estimation so that in this case, we offset by 0.5 
    # for the remaining methods, we use 0.5 as offset 
    add.bc <- ifelse(methods[i] == "wald", 0.5, ifelse(methods[i] == "small", 0, 1))
    add.ad <- ifelse(methods[i] %in% c("wald","small"), 0.5,1)
    
    # get contingency table for each simulation iteration (haldane correction performed
    # if specified in function arguments)
    conttable <- lapply(X = simXY.data_list$datXY_list, FUN = contingencytable, haldane.correct = haldane.correct, add.bc = add.bc, add.ad = add.ad)
    
    
    est.results <- lapply(X = conttable, FUN = function(x) tryCatch(oddsratio(x, method = methods[i])$measure,
                                                                    error = function(e) NA))
    
    
    OR.est[,i]  <- unlist(lapply(X = est.results, FUN = function(x) tryCatch(x["1", "estimate"],
                                                                             error = function(e) NA)))
    # get contingency table with Haldane Correction 
    # CI.est[, grep(methods[i], names(CI.est))] <- do.call(what = rbind, lapply(X = est.results, FUN = function(x) tryCatch(x["1", c("lower", "upper")],
    #                                                                                                                      error = function(e) c(NA, NA))))
    
  }
  
  # estimate OR manually  
  for(i in grep("manual", methods)){
    
    # Perform Haldane Correction in case of zerotables 
    OR.est[,i] <- lapply(X = simXY.data_list$datXY_list, FUN = contingencytable, haldane.correct = TRUE, add.bc = 0.5, add.ad=0.5) %>% 
      lapply(FUN = estOR.manual) %>% 
      unlist()
    
    # we do not compute the confidence interval manually as it is most likely too complicated for a real user 
    
  }
  
  
  # for remaining methods/pipelines: get uncorrected contingencytables #
  # Woolf corrects automatically for zeros in the contingency table 
  conttable.uncorrect <- lapply(X = simXY.data_list$datXY_list, FUN = contingencytable, haldane.correct = FALSE)  
    
    for(i in grep("Woolf", methods)){
      
      
      est.results <- lapply(X = conttable.uncorrect, FUN = function(z) tryCatch(Prop.or(z[1,], z[2,], CImethod = "Woolf"),
                                                                      error = function(e) NA))
      
      # point estimate
      OR.est[,i]  <- unlist(lapply(X = est.results, FUN = function(z) tryCatch(z$estimate,
                                                                               error = function(e) NA)))
      
      # confidence interval
      # CI.est[, grep(methods[i], names(CI.est))] <- do.call(what = rbind, lapply(X = est.results, FUN = function(x) tryCatch(as.numeric(x$conf.int),
      #                                                                                                                      error = function(e) c(NA, NA))))
      
    }
  
  
    # in the following, we perform the preceding calculations again for midp and fisher, 
    # now, however, using small or Woolf for those contingency tables with sampling zeros (fallback)
    if(haldane.correct == TRUE){
      
      # check whether contingencytable contains (i) no sampling zeros, (ii) sampling
      # zeros in the off-diagonals, or (iii) sampling zeros in the main diagonals
      cat.zero <- check.zeros(conttable.uncorrect)
      
      # define pipelines 
      OR.est$midp_fallback.small <- ifelse(cat.zero == 0 | cat.zero == 2, OR.est$midp, OR.est$small)
      OR.est$midp_fallback.Woolf <- ifelse(cat.zero == 0 , OR.est$midp, OR.est$Woolf)
      OR.est$fisher_fallback.small <- ifelse(cat.zero == 0 | cat.zero == 2, OR.est$fisher, OR.est$small)
      OR.est$fisher_fallback.Woolf <- ifelse(cat.zero == 0 , OR.est$fisher, OR.est$Woolf)
      
      
      # OR.est$midp_fallback.small <- ifelse(cat.zero == 0 | cat.zero == 2, OR.est$midp, OR.est$small)
      # unlist(lapply(FUN = midp_fisher.fallback, X = simXY.data_list$datXY_list, orig.method = "midp", fallback.method = "small"))
      # OR.est$midp_fallback.Woolf <- unlist(lapply(FUN = midp_fisher.fallback, X = simXY.data_list$datXY_list, orig.method = "midp", fallback.method = "Woolf"))
      # OR.est$fisher_fallback.small <- unlist(lapply(FUN = midp_fisher.fallback, X = simXY.data_list$datXY_list, orig.method = "fisher", fallback.method = "small"))
      # OR.est$fisher_fallback.Woolf <- unlist(lapply(FUN = midp_fisher.fallback, X = simXY.data_list$datXY_list, orig.method = "fisher", fallback.method = "Woolf"))
      # 
      
      
    }  
  
  # für die manuelle Schätzung des ORs erstmal das KI nicht berechnen, da das ein Anwender 
  # nicht unbedingt tun würde 
  # CI.est <- CI.est[, !grepl("manual", names(CI.est))]
  
  # add variable that uniquely defines parameter combination of data generating mechanism
  OR.est$params <- simXY.data_list$params
  # CI.est$params <- simXY.data_list$params
  
  
  
  return(list(estimates = OR.est,
              # confidence.intervals = CI.est, 
              true.OR = simXY.data_list$true.OR))
  
}




################################################################################
### Performance Measures #######################################################
################################################################################

####################
### Compute bias ###
####################

# Here, we compute the to the true OR 

# deal_Inf: method to deal with missing values in the OR estimations 
# (i) deal_Inf <- "available_case": compute bias based on available estimations 
#                                   for each INDIVIDUAL method
# (ii) deal_Inf <- "complete_case": for the computation of the bias include only the estimations from 
#                                   those iterations 1:R for which valid estimations are made for ALL methods 

# Note that we compute the empirical bias on the log scale 
compute.bias_logscale <- function(OR.est.df, true.OR, absolute = FALSE){
  
  # OR.est.df <- est.OR_fallback[[1]]$estimates
  # true.OR <- est.OR_fallback[[1]]$true.OR
  
  if(!all(names(OR.est.df) != "wald")){
  
  OR.est.df <- subset(OR.est.df, select = - wald)
  
  }

  
  # extract the names of the methods 
  methods <- names(OR.est.df)[names(OR.est.df) != "params"]
  # for each method, indicate whether the OR estimate is finite and non-zero
  ind_notInf <- list()
  for(i in 1:length(methods)){
    
    ind_notInf[[i]] <- which(is.finite(OR.est.df[, methods[i]])
                                      & (OR.est.df[, methods[i]] !=0))
    
  }
  
  # data frame to be filled with the coverage for method 
  bias.df <- as.data.frame(matrix(NA, ncol = length(methods), nrow =2))
  
  # for complete case analysis: get intersection of successful iterations across all methods 
  ind_completecase <- Reduce(intersect,  ind_notInf)


    
    for(i in 1:(ncol(OR.est.df)-1)){
      
      
      
      ###########################
      # available case analysis # 
      ###########################
      
      bias.df_aca <- OR.est.df[ind_notInf[[i]], i]
      
      # get coverage 
      if(absolute == FALSE){
      
      bias.df[1,i] <- mean(log(bias.df_aca) - log(true.OR))
      
      }else if(absolute == TRUE){
        
        bias.df[1,i] <- mean((log(bias.df_aca) - log(true.OR))^2)
        
      }
    
      ##########################
      # complete case analysis # 
      ##########################
    
      # reduce estimates df to iterations successful for ALL methods 
      est.OR_clean <- OR.est.df[ind_completecase,i]
      
      # calculate mean (absolute) bias of available OR estimations
      
      if(absolute == FALSE){
        
      bias.df[2,i] <- mean(log(est.OR_clean) - log(true.OR))
      
      }else if(absolute == TRUE){
        
      bias.df[2,i] <- mean((log(est.OR_clean) - log(true.OR))^2)
        
        
      }
  
      
    }
    
  # add info
  names(bias.df) <- methods
  # add parameters from underlying data generating mechanism 
  bias.df$params <- unique(OR.est.df$params)
  bias.df$deal.NA <- c("ACA", "CCA")
  
  return(bias.df)
  
}


########################
### Compute coverage ###
########################

# CI.df: a list with each element containing the CIs resulting from the R simulations
# two columns: "lower" and "upper" which correspond to the lower and upper bound 
# of the CI 
# deal_Inf: how to deal with cases in which one bound of the CI contains an Inf 
# - "deal_Inf = "compl_case_analysis": ignore these cases in the  computation of the coverage 
# - "deal_Inf = "Set_nonCovering": set the corresponding CI to "not covering the 
# true OR: true OR
# nom_coverage: nominal coverage based on which the true OR is built  
compute.coverage <- function(CI.df, true_OR, nominal.coverage = 0.95){
  
  # CI.df <- est.OR_df(datXY_sim1, FALSE)$confidence.interval
  # true_OR = 30
  
  # extract the names of the methods 
  methods <- unique(gsub("\\..*","",names(CI.df)))
  # for each method, indicate whether lower AND upper bound of the CI exist  
  ind_notInf <- lapply(X = methods, function(method) which((is.finite(CI.df[, min(grep(method, names(CI.df)))])) & (is.finite(CI.df[, max(grep(method, names(CI.df)))]))))
  
  # data frame to be filled with the coverage for method 
  coverage.df <- as.data.frame(matrix(NA, ncol = length(methods), nrow =2))
  
  
  ind_completecase <- Reduce(intersect,  ind_notInf)
  
  
  
  # proceed according to the method to handling the NAs 

    for(i in 1:length(methods)){
      
      # indicate columns of CI.df that correspond to current method
      ind_method <- grep(methods[i], names(CI.df))
      
      ###########################
      # available case analysis # 
      ###########################
      
      CI.df_aca <- CI.df[ind_notInf[[i]], ind_method]
      
      # get coverage 
      coverage.df[1,i] <- mean(apply(CI.df_aca, MARGIN = 1, 
                                     function(CI.df) between(true_OR, CI.df[1], CI.df[2])))
      
      ##########################
      # complete case analysis #
      ##########################
    
      CI.df_cca <- CI.df[ind_completecase, ind_method]
      
      # get coverage 
      coverage.df[2,i] <- mean(apply(CI.df_cca, MARGIN = 1, 
                                     function(CI.df) between(true_OR, CI.df[1], CI.df[2])))
      
    }
  

  # add info
  names(coverage.df) <- methods
  rownames(coverage.df) <- c("available case analysis", "complete case analysis")
  coverage.df$nominal.level <- rep(nominal.coverage,2)
  # add info on parameters of underlying data generating mechanism 
  coverage.df$params <- CI.df$params
  
  return(coverage.df)
  
}

#CI.df <- est.OR_df(datXY_sim1)$confidence.interval
#testtt <- compute.coverage(CI.df, 6.71, nominal.coverage = 0.95)


# compute empirical odds ratio from a contingency table manually as ad/bc
estOR.manual <- function(conttable_corrected){
  
  OR <- (conttable_corrected[1,1]*conttable_corrected[2,2]) / (conttable_corrected[1,2] * conttable_corrected[2,1])
  
  return(OR)
  
}


midp_fisher.fallback <- function(datXY, orig.method, fallback.method){
  
  # datXY <- datXY_sim[[1]]$datXY_list[[6815]]
  # orig.method <- "midp"
  # fallback.method <- "Woolf"
  
  # generate contingencytable from simulated XY data (NO Haldane correction)
  conttable <- contingencytable(datXY, haldane.correct = FALSE)
  
  # if none of the entries in the contingencytable is 0, run the chosen "original" method
  if(all(conttable != 0)){
    
    OR.estimate <- oddsratio(conttable, method = orig.method)$measure["1", "estimate"]
    
    
    # use estimator small as fallback
  }else if( (!all(conttable != 0)) && (fallback.method == "small")){
    
    
    # small returns estimate of 0 if diagonal elements are 0 
    # this leads to issues when computing bias on log scale 
    # we therefore resort to Haldance correction in this case (but: with original method)
    if(conttable[1,1] == 0 | conttable[2,2] == 0){
      
      conttable.correct <- contingencytable(datXY, haldane.correct = TRUE, add.ad = 1)
      
      OR.estimate <- oddsratio(conttable.correct, method = orig.method)$measure["1", "estimate"]
      
    }else{# if sampling zeros are in off-diagonal, use small as 
    
    
    OR.estimate <- oddsratio(conttable, method = "small")$measure["1", "estimate"]
    
    }
  
  # use estimator Woolf as fallback(works with zeros in main and off diagonals!) 
  }else if((!all(conttable != 0)) & (fallback.method == "Woolf")){
    
    conttable.correct <- contingencytable(datXY, haldane.correct = TRUE, add.ad = 1)
    
    OR.estimate <- Prop.or(conttable.correct[1,], conttable.correct[2,], CImethod = "Woolf")$estimate %>% 
                    unname()
    
  }
  
  return(OR.estimate)
}


check.zeros <- function(contingencytable.list){
  
  # contingencytable.list <- conttables
  
  cat.zero <- c()
  for(i in 1:length(contingencytable.list)){
    
    if(all(contingencytable.list[[i]] != 0)){
      
      cat.zero[i] <- 0
      
    }else if((contingencytable.list[[i]][1,2] == 0) | (contingencytable.list[[i]][2,1] == 0)){
      
      
      cat.zero[i] <- 1
      
      
    }else if(contingencytable.list[[i]][1,1] == 0 | contingencytable.list[[i]][2,2] == 0){
    
    
    cat.zero[i] <- 2
    
    }
    
  }
  
  
  return(as.factor(cat.zero))
  
}











