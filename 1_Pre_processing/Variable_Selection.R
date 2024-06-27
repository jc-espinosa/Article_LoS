
# Roots

root_input <- "./Datasets/input"
root_output <- "./Datasets/output"

# Libraries with version in comments

library(dplyr) # 1.0.10
library(BeSS) # 2.0.3
library(fastDummies) # 1.6.3
library(glmnet) # 4.1-7

###########################
##     Read datasets     ##
###########################

# Original datasets
dat_resume_sex_age <- read.csv2(paste0(root_input,"/train_outer_dataset.csv"))
dat_resume_sex_age <- read.csv2(paste0(root_input,"/test_dataset.csv"))

# Imputed data set by MICE (saved in rds)
data_sum_imputed_by_mice<- readRDS(paste0(root_output,"/data_sum_imputed_by_mice.rds"))
data_ts_sum_imputed_by_mice <- readRDS(paste0(root_output,"/data_ts_sum_imputed_by_mice.rds"))

#############################
##     LASSO Cox model     ##
#############################

# Deterioro = deterioration
# `Episodio (único)` = Unique identifier (Id)
# los = length of stay
# Sexo = Sex

# Database preparation

data_for_lasso <- data_sum_imputed_by_mice %>% dplyr::select(-`Episodio (único)`)
data_for_lasso$los <- as.numeric(data_for_lasso$los)
data_for_lasso$Deterioro <- as.numeric(data_for_lasso$Deterioro)
y = cbind(time = data_for_lasso$los, status = data_for_lasso$Deterioro)

est_neur_bin <- fastDummies::dummy_cols(data_for_lasso$last_est_neur)
est_neur_bin <- est_neur_bin[,-c(1,5)]
names(est_neur_bin) <- c("neur_1","neur_2","neur_3")
x <- data_for_lasso %>% 
  dplyr::select(-Deterioro,-los,-last_est_neur)
x <- cbind(x,est_neur_bin)
x$Fem <- as.numeric(x$Sexo) - 1 # Cathegory 0 for Male
x$Sexo <- NULL
x <- as.matrix(x)

# 1. Deterioration as main event

set.seed(1)
surv_lasso_cv <- cv.glmnet(x, y, family = "cox", type.measure = "C")
plot(surv_lasso_cv)
vars_lasso <- coef(surv_lasso_cv)

# 2. Favourable discharge as main event

y2 <- y
y2[,2] <- -1 * ( y2[,2] - 1 )
set.seed(1)
surv_lasso_cv_2 <- cv.glmnet(x, y2, family = "cox", type.measure = "C")
plot(surv_lasso_cv_2)
vars_lasso_2 <- coef(surv_lasso_cv_2)

# Selection of final dataset (s_Det and s_Fav)

var_def_lasso <- data.frame(Variables=names(data_for_lasso))
s_Det_lasso <- data.frame(as.matrix(vars_lasso)) %>% dplyr::filter(X1 !=0)
s_Det_lasso <- data.frame(Variables=row.names(s_Det_lasso),Var_Det=s_Det_lasso$X1)
var_def_lasso <- var_def_lasso %>% 
  left_join(s_Det_lasso)

s_Fav_lasso <- data.frame(as.matrix(vars_lasso_2)) %>% dplyr::filter(X1 !=0)
s_Fav_lasso <- data.frame(Variables=row.names(s_Fav_lasso),Var_Alta=s_Fav_lasso$X1)
var_def_lasso <- var_def_lasso %>% 
  left_join(s_Fav_lasso)

# Final set of variables under LASSO (set s)

var_def_lasso <- var_def_lasso %>% mutate(
  Var_Det=ifelse(Variables %in% c("last_est_neur","Sexo","Deterioro","los"),1,Var_Det),
  Var_Alta=ifelse(Variables %in% c("last_est_neur","Sexo","Deterioro","los"),1,Var_Alta)
)

var_def_lasso

############################
##     BeSS cox model     ##
############################

# 1. Deterioration as main event

surv_bess <- bess(x, y, method = "sequential",
                  family = "cox")
summary(surv_bess)
surv_bess_one <- bess.one(x, y,s=40,
                          family = "cox")
s_Det_bess <- names(surv_bess_one$bestmodel$coefficients)
s_Det_bess <- gsub("xbest","",s_Det_bess)

# 2. Favourable discharge as main event

surv_bess2 <- bess(x, y2, method = "sequential",
                   family = "cox")
summary(surv_bess2)
surv_bess_one2 <- bess.one(x, y2,s=40,
                           family = "cox")
s_Fav_bess <- names(surv_bess_one2$bestmodel$coefficients)
s_Fav_bess <- gsub("xbest","",s_Fav_bess)

# Final set of variables selected by BeSS (set s)

var_def_bess <- data.frame(Variables=names(data_for_lasso))
var_def_bess <- var_def_bess %>% left_join(data.frame(Variables=s_Det_bess,Var_Det=1))
var_def_bess <- var_def_bess %>% left_join(data.frame(Variables=s_Fav_bess,Var_Alta=1))
var_def_bess <- var_def_bess %>% mutate(
  Var_Det=ifelse(Variables %in% c("last_est_neur","Sexo","Deterioro","los"),1,Var_Det),
  Var_Alta=ifelse(Variables %in% c("last_est_neur","Sexo","Deterioro","los"),1,Var_Alta)
)

var_def_bess


####################
##     Output     ##
####################

### Final Output: two database with variable selection from LASSO and BeSS

# Save data-sets
saveRDS(var_def_lasso,
        paste0(root_output,"/var_def_lasso.rds"))

saveRDS(var_def_bess,
        paste0(root_output,"/var_def_bess.rds"))

