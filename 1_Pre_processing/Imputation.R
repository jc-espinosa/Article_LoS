
# Roots

root_input <- "./Datasets/input"
root_output <- "./Datasets/output"

# Libraries with version in comments

library(mice) # 3.15.0
library(caret) # 6.0-93
library(data.table) # 1.14.6
library(dplyr) # 1.0.10
library(missMethods) # 0.4.0
library(pcaMethods) # 1.90.0

# Functions

#' percNA calculates the missing data rate in a vector.
#' 
#' @param x A numeric vector.
#' @returns A numeric value.
#' @examples
#' x <- c(1,2,3,NA,5,6,NA)
#' percNA(x)
percNA <- function(x){
  return(length(which(is.na(x)))/length(x))
}

##########################
##     Read datasets    ##
##########################

# Train outer data-set
dat_resume_sex_age <- read.csv2(paste0(root_input,"/train_outer_dataset.csv"))

########################
##     Data split     ##
########################

# tr sufix here is assigned to train-inner
# ts sufix here is asigned to validation
# Deterioro = deterioration
# `Episodio (único)` = Unique identifier (Id)

dat_resume_sex_age$Deterioro <- as.factor(dat_resume_sex_age$Deterioro)
n_t <- dim(dat_resume_sex_age)[1]
per_tr <- 0.7 # 70%
n_tr <- floor(n_t*per_tr)
set.seed(12)
train.index <- createDataPartition(dat_resume_sex_age$Deterioro, p = per_tr, list = FALSE)
train.index <- dat_resume_sex_age[train.index,] %>% dplyr::select(`Episodio (único)`)
dat_resume_sex_age$Deterioro <- as.numeric(dat_resume_sex_age$Deterioro)

# Train-inner data-set
dat_tr <- dat_resume_sex_age %>% inner_join(train.index)
# Validation data-set
dat_ts <- dat_resume_sex_age %>% anti_join(train.index)

#########################
##     Preliminars     ##
#########################

# Create a copy of databases

dat_tr_no_na <- copy(dat_tr)
summary(dat_tr_no_na)
dat_ts_no_na <- copy(dat_ts)
summary(dat_ts_no_na)

# Filter only numeric variables in databases
# and delete Id and outcomes (deterioration and los)
# los = length of stay

dat_tr_no_na_pca1 <- dat_tr_no_na %>% 
  select_if(is.numeric) %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los)

dat_ts_no_na_pca1 <- dat_ts_no_na %>% 
  select_if(is.numeric) %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los)

# To ensure the same order of variables in testing (w.r.t training)

dat_ts_no_na_pca1 <- dat_ts_no_na_pca1 %>%
  dplyr::select(names(dat_tr_no_na_pca1))

# Calculate missing data rate

perc_of_na <- apply(dat_ts_no_na_pca1,2,percNA)

# Delete missing data to work with complete cases

base_comp <- na.omit(dat_ts_no_na_pca1)

# Void function to create data-sets with artificial missing under MCAR approach

base_na_artif <- function(){
  set.seed(NULL)
  delete_MCAR(base_comp,perc_of_na)
}

# Iterations 
R = 500

################################
##     Imputation by mean     ##
################################

# Calculate global mean

means_train <- apply(dat_tr_no_na_pca1,2,mean,na.rm=TRUE)

# Arrays to allocate results at each iteration

rmse_mean_red_data <- array(0,dim=c(R,dim(base_comp)[2]))
mae_mean_red_data <- array(0,dim=c(R,dim(base_comp)[2]))

# Iterations to calculate NRMSE and NMAE

for(i in 1:R){
  base_na_artif_imput1 <- base_na_artif()
  for(j in colnames(base_na_artif_imput1)){
    base_na_artif_imput1[,j][is.na(base_na_artif_imput1[,j])] <- means_train[j]
  }
  temp_num_rmse <- apply((base_na_artif_imput1 - base_comp)^2,2,function(x) sqrt(mean((x))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  rmse_mean_red_data[i,] <- temp_num_rmse/temp_den_rmse 
  temp_num_rmse <- apply((base_na_artif_imput1 - base_comp),2,function(x) mean(abs((x))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  mae_mean_red_data[i,] <- temp_num_rmse/temp_den_rmse
  print(i)
}

# NRMSE data-set
rmse_mean_red_data <- data.frame(rmse_mean_red_data)
names(rmse_mean_red_data) <- names(base_comp)

# NMAE data-set
mae_mean_red_data <- data.frame(mae_mean_red_data)
names(mae_mean_red_data) <- names(base_comp)
mae_mean_red_data

# Mean NRMSE
temp_imp_mean_rmse <- apply(rmse_mean_red_data,2,mean)
mean(temp_imp_mean_nrmse[temp_imp_mean_nrmse>0])
sd(temp_imp_mean_nrmse[temp_imp_mean_nrmse>0])

# Mean NMAE
temp_imp_mean_mae <- apply(rmse_mean_red_data,2,mean)
mean(temp_imp_mean_mae[temp_imp_mean_mae>0])
sd(temp_imp_mean_mae[temp_imp_mean_mae>0])

#################################
##     Imputation by median    ##
#################################

# Calculate global median

medians_train <- apply(dat_tr_no_na_pca1,2,median,na.rm=TRUE)

# Arrays to allocate results at each iteration

rmse_median_red_data <- array(0,dim=c(R,dim(base_comp)[2]))
mae_median_red_data <- array(0,dim=c(R,dim(base_comp)[2]))

# Iterations to calculate NRMSE and NMAE

for(i in 1:R){
  base_na_artif_imput1 <- base_na_artif()
  for(j in colnames(base_na_artif_imput1)){
    base_na_artif_imput1[,j][is.na(base_na_artif_imput1[,j])] <- medians_train[j]
  }
  temp_num_rmse <- apply((base_na_artif_imput1 - base_comp)^2,2,function(x) sqrt(mean((x))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  rmse_median_red_data[i,] <- temp_num_rmse/temp_den_rmse
  temp_num_rmse <- apply((base_na_artif_imput1 - base_comp),2,function(x) mean(abs((x))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  mae_median_red_data[i,] <- temp_num_rmse/temp_den_rmse # 
  print(i)
}

# NRMSE data-set
rmse_median_red_data <- data.frame(rmse_median_red_data)
names(rmse_median_red_data) <- names(base_comp)

# NMAE data-set
mae_median_red_data <- data.frame(mae_median_red_data)
names(mae_median_red_data) <- names(base_comp)

# Mean NRMSE
temp_imp_median_rmse <- apply(rmse_median_red_data,2,mean)
mean(temp_imp_median_rmse[temp_imp_median_rmse>0])
sd(temp_imp_median_rmse[temp_imp_median_rmse>0])

# Mean NMAE
temp_imp_median_mae <- apply(rmse_median_red_data,2,mean)
mean(temp_imp_median_mae[temp_imp_median_mae>0])
sd(temp_imp_median_mae[temp_imp_median_mae>0])

###############################
##     Imputation by MICE    ##
###############################
# Number of iterations 

iter_mice <- 5

# Arrays to allocate results at each iteration

rmse_mice_red_data <- array(0,dim=c(R,dim(base_comp)[2]))
mae_mice_red_data <- array(0,dim=c(R,dim(base_comp)[2]))

# Separe databases in two: numerical and categorical
# and delete Id and outcomes (deterioration and los)

data_sum_imputed_by_mice_num <- dat_tr %>%  
  dplyr::select(-`Episodio (único)`,-Deterioro,-los) %>% 
  select_if(is.numeric)
data_sum_imputed_by_mice_cat <- dat_tr %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los) %>% 
  select_if(is.factor)
data_ts_sum_imputed_by_mice_cat <- dat_ts %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los) %>% 
  select_if(is.factor)

# Iterations to calculate NRMSE and NMAE

# For parallel processing, indicate the number of cores for MICE
ncores_mice <- 4

for(i in 1:R){
  data_ts_sum_imputed_by_mice_num <- base_na_artif()
  ignore <- c(rep(TRUE,dim(data_ts_sum_imputed_by_mice_num)[1]),
              rep(FALSE,dim(data_sum_imputed_by_mice_num)[1]))
  data_ts_sum_imputed_by_mice_prev <- futuremice(rbind(data_ts_sum_imputed_by_mice_num,
                                                       data_sum_imputed_by_mice_num),
                                                 remove.collinear = FALSE,
                                                 m=iter_mice,maxit=10,meth='pmm',
                                                 ignore=ignore,
                                                 parallelseed=1234,
                                                 n.core	= ncores_mice)
  data_ts_sum_imputed_by_mice <- complete(data_ts_sum_imputed_by_mice_prev,1)
  data_ts_sum_imputed_by_mice <- data_ts_sum_imputed_by_mice/iter_mice
  for(j in 2:iter_mice){
    data_ts_sum_imputed_by_mice <- data_ts_sum_imputed_by_mice + 
      (complete(data_ts_sum_imputed_by_mice_prev,j))/iter_mice
  }
  data_ts_sum_imputed_by_mice <- 
    data_ts_sum_imputed_by_mice[1:(dim(data_ts_sum_imputed_by_mice_num)[1]),]
  temp_num_rmse <- as.vector(apply((data_ts_sum_imputed_by_mice - base_comp)^2,2,
                                   function(x) sqrt(mean((x)))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  rmse_mice_red_data[i,] <- temp_num_rmse/temp_den_rmse
  temp_num_rmse <- as.vector(apply((data_ts_sum_imputed_by_mice - base_comp),2,
                                   function(x) mean(abs((x)))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  mae_mice_red_data[i,] <- temp_num_rmse/temp_den_rmse
}

# NRMSE data-set
rmse_mice_red_data <- data.frame(rmse_mice_red_data)
names(rmse_mice_red_data) <- names(base_comp)

# NMAE data-set
mae_mice_red_data <- data.frame(mae_mice_red_data)
names(mae_mice_red_data) <- names(base_comp)

# Mean NRMSE
temp_imp_mice_rmse <- apply(rmse_mice_red_data,2,mean)
mean(temp_imp_mice_rmse[temp_imp_mice_rmse>0])
sd(temp_imp_mice_rmse[temp_imp_mice_rmse>0])

# Mean NMAE
temp_imp_mice_mae <- apply(rmse_mice_red_data,2,mean)
mean(temp_imp_mice_mae[temp_imp_mice_mae>0])
sd(temp_imp_mice_mae[temp_imp_mice_mae>0])

#######################################
##     Imputation by Bayesian PCA    ##
#######################################

# Filter numerical variables
# and delete Id and outcomes (deterioration and los)

dat_tr_for_nlpca <- dat_tr %>% 
  dplyr::select_if(is.numeric) %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los)
dat_ts_for_nlpca <- dat_ts %>% 
  dplyr::select_if(is.numeric) %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los)

# Train a BPCA with 20 principal components and a maxSteps of 100. 
# Data is scaled using unit variance (uv)

bpca1 <- pca(dat_tr_for_nlpca, nPcs = 20,maxSteps=100, method="bpca",
             scale="uv",completeObs = TRUE)

# Calculate missing data rates

perc_of_na_pca <- apply(dat_ts_for_nlpca,2,percNA)
base_comp_pca <- na.omit(dat_ts_for_nlpca)

# Void function to create missing data under MCAR

base_na_artif2 <- function(){
  delete_MCAR(base_comp_pca,perc_of_na_pca)
}

# Arrays to allocate results at each iteration

rmse_bpca_red_data <- array(0,dim=c(R,dim(base_comp)[2]))
mae_bpca_red_data <- array(0,dim=c(R,dim(base_comp)[2]))

# Iterations to calculate NRMSE and NMAE

for(i in 1:R){
  dat_ts_resume_sex_age_for_bpca <- base_na_artif2()
  imputed_dat_ts_bpca <- predict(bpca1, 
                                 newdata = dat_ts_resume_sex_age_for_bpca,
                                 pre = TRUE,
                                 pcs = nP(bpca1))
  temp_num_rmse <- apply((imputed_dat_ts_bpca[[2]] - base_comp)^2,2,
                         function(x) sqrt(mean((x))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  rmse_bpca_red_data[i,] <- temp_num_rmse/temp_den_rmse
  temp_num_rmse <- apply((imputed_dat_ts_bpca[[2]] - base_comp),2,
                         function(x) mean(abs((x))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  mae_bpca_red_data[i,] <- temp_num_rmse/temp_den_rmse
  print(i)
}

# NRMSE data-set
rmse_bpca_red_data <- data.frame(rmse_bpca_red_data)
names(rmse_bpca_red_data) <- names(base_comp)

# NMAE data-set
mae_bpca_red_data <- data.frame(mae_bpca_red_data)
names(mae_bpca_red_data) <- names(base_comp)

# Mean NRMSE
temp_imp_bpca_rmse <- apply(rmse_bpca_red_data,2,mean)
mean(temp_imp_bpca_rmse[temp_imp_bpca_rmse>0])
sd(temp_imp_bpca_rmse[temp_imp_bpca_rmse>0])

# Mean NMAE
temp_imp_bpca_rmse <- apply(rmse_bpca_red_data,2,mean)
mean(temp_imp_bpca_rmse[temp_imp_bpca_rmse>0])
sd(temp_imp_bpca_rmse[temp_imp_bpca_rmse>0])

#################################
##     Imputation by NIPALS    ##
#################################

# Filter numerical variables
# and delete Id and outcomes (deterioration and los)

dat_tr_for_nlpca <- dat_tr %>% 
  dplyr::select_if(is.numeric) %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los)
dat_ts_for_nlpca <- dat_ts %>% 
  dplyr::select_if(is.numeric) %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los)

# Train a non-linear PCA with 20 principal components employing NIPALS. 
# Data is scaled using unit variance (uv)

nipals1 <- pca(dat_tr_for_nlpca, nPcs = 20, 
               method="nipals",scale="uv",completeObs = TRUE)

# Arrays to allocate results at each iteration

rmse_nipals_red_data <- array(0,dim=c(R,dim(base_comp)[2]))
mae_nipals_red_data <- array(0,dim=c(R,dim(base_comp)[2]))

# Iterations to calculate NRMSE and NMAE

for(i in 1:R){
  dat_ts_resume_sex_age_for_nipals <- base_na_artif2()
  imputed_dat_ts_nipals <- predict(nipals1, 
                                   newdata = dat_ts_resume_sex_age_for_nipals,
                                   pre = TRUE,
                                   pcs = nP(nipals1))
  temp_num_rmse <- apply((imputed_dat_ts_nipals[[2]] - base_comp)^2,2,
                         function(x) sqrt(mean((x))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  rmse_nipals_red_data[i,] <- temp_num_rmse/temp_den_rmse
  temp_num_rmse <- apply((imputed_dat_ts_nipals[[2]] - base_comp),2,
                         function(x) mean(abs((x))))
  temp_den_rmse <- apply((base_comp),2,function(x) sd((x)))
  mae_nipals_red_data[i,] <- temp_num_rmse/temp_den_rmse
  print(i)
}

# NRMSE data-set
rmse_nipals_red_data <- data.frame(rmse_nipals_red_data)
names(rmse_nipals_red_data) <- names(base_comp)

# NMAE data-set
mae_nipals_red_data <- data.frame(mae_nipals_red_data)
names(mae_nipals_red_data) <- names(base_comp)

# Mean NRMSE
temp_imp_nipals_rmse <- apply(rmse_nipals_red_data,2,mean)
mean(temp_imp_nipals_rmse[temp_imp_nipals_rmse>0])
sd(temp_imp_nipals_rmse[temp_imp_nipals_rmse>0])

# Mean NMAE
temp_imp_nipals_mae <- apply(rmse_nipals_red_data,2,mean)
mean(temp_imp_nipals_mae[temp_imp_nipals_mae>0])
sd(temp_imp_nipals_mae[temp_imp_nipals_mae>0])


## Average NRMSE and NMAE over validation data-set

mean_imp <- array(0,dim=c(5,dim(rmse_mean_red_data)[2]))
mean_imp[1,] <- apply(rmse_mean_red_data,2,mean)
mean_imp[2,] <- apply(rmse_median_red_data,2,mean)
mean_imp[3,] <- apply(rmse_bpca_red_data,2,mean)
mean_imp[4,] <- apply(rmse_nipals_red_data,2,mean)
mean_imp[5,] <- apply(rmse_mice_red_data,2,mean)

apply(mean_imp,1,mean)
apply(mean_imp,1,sd)

mean_imp_mae <- array(0,dim=c(5,dim(rmse_mean_red_data)[2]))
mean_imp_mae[1,] <- apply(mae_mean_red_data,2,mean)
mean_imp_mae[2,] <- apply(mae_median_red_data,2,mean)
mean_imp_mae[3,] <- apply(mae_bpca_red_data,2,mean)
mean_imp_mae[4,] <- apply(mae_nipals_red_data,2,mean)
mean_imp_mae[5,] <- apply(mae_mice_red_data,2,mean)

apply(mean_imp_mae,1,mean)
apply(mean_imp_mae,1,sd)

##########################
##     Final Output     ##
##########################

### Final Output: one database imputed by MICE for training and one for testing

# Read train-outer and test data-sets
# tr sufix here is assigned to train-outer
# ts sufix here is asigned to testing

# Train-outer
dat_resume_sex_age <- readRDS(paste0(root_input,"/dat_resume_sex_age2.rds"))
# Test
dat_resume_ts_sex_age <- readRDS(paste0(root_input,"/dat_resume_ts_sex_age2.rds"))

# Define MICE parameters

iter_mice <- 5

# Subset numeric and categorical variables from train

data_sum_imputed_by_mice_num <- dat_resume_sex_age %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los) %>% 
  select_if(is.numeric)
data_sum_imputed_by_mice_cat <- dat_resume_sex_age %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los) %>% 
  select_if(is.factor)

# MICE in parallel over train data-set

data_sum_imputed_by_mice_prev <- futuremice(data_sum_imputed_by_mice_num,
                                            m=iter_mice,maxit=10,
                                            remove.collinear = FALSE,
                                            meth='pmm',
                                            n.core = ncores_mice, 
                                            parallelseed=1234)
data_sum_imputed_by_mice <- complete(data_sum_imputed_by_mice_prev,1)
data_sum_imputed_by_mice <- data_sum_imputed_by_mice/iter_mice
for(i in 2:iter_mice){
  data_sum_imputed_by_mice <- data_sum_imputed_by_mice + 
    (complete(data_sum_imputed_by_mice_prev,i))/iter_mice
}

# Prepare data-set
data_sum_imputed_by_mice <- cbind(data_sum_imputed_by_mice,data_sum_imputed_by_mice_cat)
data_sum_imputed_by_mice$`Episodio (único)` <- dat_resume_sex_age$`Episodio (único)`
data_sum_imputed_by_mice$Deterioro <- dat_resume_sex_age$Deterioro
data_sum_imputed_by_mice$los <- dat_resume_sex_age$los

# Save data-set
saveRDS(data_sum_imputed_by_mice,
        paste0(root_output,"/data_sum_imputed_by_mice.rds"))

### Testing ###

# MICE in parallel over test data-set

# Subset numeric and categorical variables from test

data_ts_sum_imputed_by_mice_num <- dat_resume_ts_sex_age %>%  
  dplyr::select(-`Episodio (único)`,-Deterioro,-los) %>% 
  select_if(is.numeric)
data_ts_sum_imputed_by_mice_cat <- dat_resume_ts_sex_age %>% 
  dplyr::select(-`Episodio (único)`,-Deterioro,-los) %>% 
  select_if(is.factor)

# Identify which elements are in train and test

ignore <- c(rep(TRUE,dim(data_ts_sum_imputed_by_mice_num)[1]),
            rep(FALSE,dim(data_sum_imputed_by_mice_num)[1]))

# MICE in parallel over test data-set

data_ts_sum_imputed_by_mice_prev <- futuremice(rbind(data_ts_sum_imputed_by_mice_num,
                                                     data_sum_imputed_by_mice_num),
                                               remove.collinear = FALSE,
                                               m=iter_mice,maxit=10,meth='pmm',
                                               ignore=ignore,
                                               n.core = 4, 
                                               parallelseed=1234)

data_ts_sum_imputed_by_mice <- complete(data_ts_sum_imputed_by_mice_prev,1)
data_ts_sum_imputed_by_mice <- data_ts_sum_imputed_by_mice/iter_mice
for(i in 2:iter_mice){
  data_ts_sum_imputed_by_mice <- data_ts_sum_imputed_by_mice + 
    (complete(data_ts_sum_imputed_by_mice_prev,i))/iter_mice
}
data_ts_sum_imputed_by_mice <- 
  data_ts_sum_imputed_by_mice[1:dim(data_ts_sum_imputed_by_mice_num)[1],]
summary(data_ts_sum_imputed_by_mice)

# Prepare data-set
data_ts_sum_imputed_by_mice <- cbind(data_ts_sum_imputed_by_mice,
                                     data_ts_sum_imputed_by_mice_cat)
data_ts_sum_imputed_by_mice$`Episodio (único)` <- dat_resume_ts_sex_age$`Episodio (único)`
data_ts_sum_imputed_by_mice$Deterioro <- dat_resume_ts_sex_age$Deterioro
data_ts_sum_imputed_by_mice$los <- dat_resume_ts_sex_age$los

# Save data-set
saveRDS(data_ts_sum_imputed_by_mice,
        paste0(root_output,"/data_ts_sum_imputed_by_mice.rds"))

