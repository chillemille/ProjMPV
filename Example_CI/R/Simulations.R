################################################################################
### Playaround #################################################################
################################################################################

library(mlr3)
library(rpart)
library(dplyr)

# get SAPA data set 
load("./Example_CI/Data/sapa_noltenius.RData")
# load help functions 
source("./Example_CI/R/help_functions.R")




################################################################################
### Run ########################################################################
################################################################################

# run simulation with 1000 iterations 
simulation = experiment(sapa[,c("smoking", "neuro", "consci", "agree", "open", "extra")], 
                        "smoking", 1000, 15, 0.8, 1)

# save(simulation, file = "./Example_CI/Data/Results/simulation.RData")
# load("./Example_CI/Data/Results/simulation.RData")


# get number of iterations in which method N failed
sum(is.na(simulation$conf.int_naive))


#####################################################
# first three values in Table 2 in main manuscript ##
#####################################################
# get coverage for both methods, for available and complete case analysis 
handle.missCI(simulation$check.cover_naive, simulation$check.cover_correct)

################################################
# fourth column in Table 2 in main manuscript ##
################################################

# get coverage resulting from setting undefined "naive" CIs to point estimator 
check.cover_naive.singleton <- c()

for(i in 1:length(simulation$check.cover_naive)){
  
  # if the CI is NOT undefined, take the corresponding covering indicator
  if(!is.na(simulation$check.cover_naive[[i]])){
    
    check.cover_naive.singleton[i] <- simulation$check.cover_naive[[i]]
  
  }else{# else: set undefined CI as corresponding point estimator and check 
    # whether it is equal to the approx. true AUC
      
    check.cover_naive.singleton[i] <- simulation$true.AUC[[i]] == simulation$subsampl.AUC_aggr[[i]]
    
    
    
    }
  
}
# resulting coverage 
mean(check.cover_naive.singleton)

################################################################################
### Investigate source of missingness ##########################################
################################################################################

# indicate which iteration leads to failure of method N
ind.fail <- which(is.na(simulation$conf.int_naive))

# inspect true and estimated AUCs in these failed iterations 

# true AUC: all equal to 0.5
all(simulation$true.AUC[ind.fail] == 0.5)

# check whether all of the 15 estimated AUCs are = 0.5 in all "problematic" simulation iterations
table(unlist(lapply(X = simulation$subsampl.AUC_all[ind.fail], function(AUC_all) all(AUC_all == 0.5))))


# count number of splits for the decision trees in the "problematic" iterations 
n.node <- unlist(lapply(simulation$model[ind.fail], FUN = function(dec.tree) dec.tree$parms$split))
# get table of frequencies 
table(n.node) #->> in all of these iterations, the decision tree consists of a single node, i.e. no splitting is performed 

  
  
