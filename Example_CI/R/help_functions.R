################################################################################
### function to obtain CI based on corrected resampled t-test ##################
################################################################################

# t_test_result: results from function t.test()
# ratio_testtotrain: ratio of the test set to the training set (n.test/n.train)
# ->> as proposed by Nadeau and Bengio (1999)

corrected_resampled_t.test <- function(t_test_result, ratio_testtotrain) {
  
  stopifnot(t_test_result$alternative == "two.sided")
  alpha <- 1 - attr(t_test_result$conf.int, "conf.level")
  J <- unname(t_test_result$parameter + 1L)
  # get standard error from regular t-test
  stderr <- t_test_result$stderr
  # correct standard error
  stderr_corrected <- sqrt(((1 / J) + (ratio_testtotrain)) * (J * (stderr^2)))
  # get corrected resampled t-statistic
  statistic <- (t_test_result$estimate - t_test_result$null.value) / stderr_corrected
  t_quantile_stderr_corrected <- qt(1 - (alpha / 2), df = J - 1L) * stderr_corrected
  conf.int <- c(t_test_result$estimate - t_quantile_stderr_corrected, t_test_result$estimate + t_quantile_stderr_corrected)
  attr(conf.int, which = "conf.level") <- 1 - alpha
  
  t_test_result$statistic <- setNames(statistic, nm = "t_corrected")
  t_test_result$p.value <- unname(2 * pt(-abs(statistic), df = J - 1L))
  #' Addition (otherwise the bounds are both called "mean of x") ----
  t_test_result$conf.int <- setNames(conf.int, NULL)
  #'   --------------------------------------------------------------
  t_test_result$stderr <- stderr_corrected
  t_test_result$method <- paste0("Corrected ", t_test_result$method)
  
  t_test_result
}


CI.correct <- function(est.AUC_vec, ratio_test.train, alpha = 0.05){
  
  # correction term 
  c <- (1/length(est.AUC_vec)) + ratio_test.train
  # sample variance 
  sample.var <- sd(est.AUC_vec)^2
  # lower bound of corrected confidence interval 
  CI_lower <- mean(est.AUC_vec) - qt(1-(alpha/2), length(est.AUC_vec) - 1) * sqrt((c*sample.var))
  # upper bound of corrected confidence interval 
  CI_upper <- mean(est.AUC_vec) + qt(1-(alpha/2), length(est.AUC_vec) - 1) * sqrt((c*sample.var))
  
  # corrected confidence interval 
  conf.int_corrected <- c(CI_lower, CI_upper)
  
  return(conf.int_corrected)
  
}


################################################################################
experiment <- function(data, target_var, n.rep, n_it.subsampl = 15, ratio.subsampl = 0.8, seed){
  
  
  # data = sapa[,c("smoking", "neuro", "consci", "agree", "open", "extra")]
  # target_var = "smoking"
  # n_it.subsampl = 15
  # ratio.subsampl = 0.8
  # n.rep = 20
  # seed = 1
  
  
  set.seed(seed)
  
  
  #n.test <- round(nrow(data)*0.8810161)
  n.train <- round(nrow(data)*0.20)
  
  # 2. get ratio of test set to training set (needed for corrected CI)
  ratio_testtotrain <- (nrow(data) - n.train) / n.train
  
  
  # initiate list to store state of random number generator 
  statesdat_list <- list()
  
  
  true.auc_vec <- c()
  model_list <- list()
  
  subsampl.AUC_list <- list()
  subsampl.AUC_aggr <- list()
  
  # vector to check whether the corrected confidence interval contains the true AUC 
  check.cover_naive <- c()
  check.cover_correct <- c()
  
  conf.int_naive <- list()
  conf.int_correct <- list()
  
  
  # step 1: Partition data data set into training (20%) and test set (80%) 
  for(j in 1:n.rep){
    
    # store random number generator state
    statesdat_list[[j]] = .Random.seed
    
    #########
    # step 1: Partition data data into training and test set 
    #########
    
    
    # sample indices of observations to be placed in the test set
    ind.train = sample(1:nrow(data), n.train, replace = FALSE)
    # assign remaining indices to training set
    ind.test = (1:nrow(data))[! (1:nrow(data)) %in% ind.train]
    
    # reduce the training data set once more 
    #ind.train2 <- sample(ind.train, round(length(ind.train)*0.5), replace = FALSE)
    
    # subset data set to test set
    data.test = data[ind.test, ]
    data.train = data[ind.train, ]
    
    
    
    # generate task from the training set 
    task = TaskClassif$new(
      id = "SAPA",
      backend = data.train,
      target = target_var,
    )
    
    # set learner as classification tree 
    learner = lrn("classif.rpart")
    
    
    
    #########
    # step 2: Approximate "true" generalization error 
    #########
    
    
    # train on training data set 
    model = learner$train(task)
    # store model 
    model_list[[j]] = model$model
    
    # predict on test set 
    model$predict_type = "prob"
    pred = model$predict_newdata(as.data.table(data.test))
    
    # obtain approximated true AUC 
    true.auc_vec[j]<- pred$score(msr("classif.auc"))
    
    ########
    # step 3: obtain true AUC
    ########
    
    # set subsampling (15 iterations, ratio 4:1) as resampling strategy 
    #subsampl = rsmp("cv", folds = 5)
    subsampl = rsmp("subsampling", repeats = n_it.subsampl, ratio = ratio.subsampl)
    
    # perform repeated subsampling
    subsampl.run = resample(task, learner, subsampl)
    
    # store AUCs from 15 repetitions
    subsampl.AUC_list[[j]] = subsampl.run$score(msr("classif.auc"))$classif.auc
    # aggregate AUCs from the 15 repetitions
    subsampl.AUC_aggr[[j]] = subsampl.run$aggregate(msr("classif.auc"))
    
    #########
    # step 4: apply corrected resampled t-test to estimated AUCs from resampling
    ########
    
    # 1. obtain naive CI from function t.test() 
    conf.int_naive[[j]] = tryCatch(as.numeric(t.test(subsampl.AUC_list[[j]])$conf.int), 
                                   error = function(e) NA)
    
    # 2. obtain confidence interval stemming from the corrected resampled t-test
    conf.int_correct[[j]] = CI.correct(subsampl.AUC_list[[j]], ratio_testtotrain)
    
    # 3. check whether the confidence interval contains the true AUC 
    check.cover_naive[j] <- tryCatch(between(true.auc_vec[j], conf.int_naive[[j]][1], conf.int_naive[[j]][2]), 
                                     error = function(e) NA)
    
    # 3. check whether the confidence interval contains the true AUC 
    check.cover_correct[j] <- tryCatch(between(true.auc_vec[j], conf.int_correct[[j]][1], conf.int_correct[[j]][2]), 
                                       error = function(e) NA)
    
  }
  
  
  
  
  return(list(true.AUC = true.auc_vec, # vector of approx. true AUCs 
              model = model_list, # vector of desicion tree models 
              subsampl.AUC_all = subsampl.AUC_list, # 15 estimated AUCs from repeated subsampling
              subsampl.AUC_aggr = subsampl.AUC_aggr, # aggregated AUC estimation 
              conf.int_naive = conf.int_naive, 
              conf.int_correct = conf.int_correct, 
              check.cover_naive = check.cover_naive, 
              check.cover_correct = check.cover_correct, # indication whether confidence interval contains true AUC 
              # emp.coverage = emp.coverage), #empirical coverage (single value or error message)
              states = statesdat_list)) # random states 
}




handle.missCI <- function(check.cover_method1, check.cover_method2){
  
  # check.cover_method1 <- simulation$check.cover_naive
  # check.cover_method2 <- simulation$check.cover_correct
  
  
  coverage.df <- data.frame()
  
  
  ###########################
  # Available case analysis #
  ###########################
  
  coverage.df["Naive", "avail.case"] <- mean(check.cover_method1, na.rm = TRUE)  
  coverage.df["Corrected", "avail.case"] <- mean(check.cover_method2, na.rm = TRUE)  
  
  ##########################
  # Complete case analysis #
  ##########################
  
  # indicate iterations in which both methods produce a CI     
  success.iter <- !is.na(check.cover_method1) & !is.na(check.cover_method2)
  
  # compute empirical coverage based only those iterations where the
  # confidence interval could be computed
  coverage.df["Naive", "compl.case"] <- mean(check.cover_method1[success.iter])
  coverage.df["Corrected", "compl.case"] <- mean(check.cover_method2[success.iter])
  
  #############################################
  # Set missing CIs to "not covering true AUC # 
  #############################################
  
  # in computation of the coverage, set missing confidence intervals
  # to "not covering the true value"
  coverage.df["Naive", "set.FALSE"] = mean(ifelse(is.na(check.cover_method1),
                                                  FALSE,
                                                  check.cover_method1))
  
  coverage.df["Corrected", "set.FALSE"] = mean(ifelse(is.na(check.cover_method2),
                                                      FALSE,
                                                      check.cover_method2))
  
  return(coverage.df)
  
  
}