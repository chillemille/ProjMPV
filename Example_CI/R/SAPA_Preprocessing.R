################################################################################
### Preprocess raw SAPA data set ###############################################
################################################################################

library(mlr3)

# # load raw SAPA data set 
# load("./Example_CI/Data/SAPA_raw/SAPAdata01jan2022thru31dec2022.rdata")
# # remove item info 
# rm(ItemInfo)
# 
# # rename SAPA data set 
# SAPA <- SAPAdata01jan2022thru31dec2022
# rm(SAPAdata01jan2022thru31dec2022)

# load SAPA data set 
load("./Example_CI/Data/SAPA_raw/sapa_favorite_new.RData")


################################################################################
### only BIG5 variables as predictors ##########################################
################################################################################

SAPA.BIG5 <- sapa_favorite



# subset data set: we only want to keep the targe (education) and the BIG5 variables 
SAPA.BIG5 = SAPA.BIG5[, c("smoking", "extra", "neuro", "consci", "agree", "open")]

# remove all observations with missing values 
SAPA.BIG5 = SAPA.BIG5[complete.cases(SAPA.BIG5),]
dim(SAPA.BIG5)

# transform target (education) to factor
SAPA.BIG5$smoking <- factor(SAPA.BIG5$smoking, levels = c(FALSE,TRUE), labels = c(0,1))


################################################################################
### Use full set of predictor variables ########################################
################################################################################

SAPA.Full <- sapa_favorite

# remove unnecessary features and remaining 4 targets 
SAPA.Full <- subset(SAPA.Full,select = c(-continent, - BMI, -education, -relstatus, -phys_act))

# only consider complete cases 
SAPA.Full <- SAPA.Full[complete.cases(SAPA.Full), ]

# set categorical variables as factors 
SAPA.Full$gender <- as.factor(SAPA.Full$gender)
SAPA.Full$age <- as.factor(SAPA.Full$age)
SAPA.Full$jobstatus <- as.factor(SAPA.Full$jobstatus)
SAPA.Full$smoking <- as.factor(SAPA.Full$smoking)





################################################################################
### remove initial data set ####################################################
################################################################################

rm(sapa_favorite)







