
# Roots

root_input <- "./Datasets/input"
root_output <- "./Datasets/output"

# Libraries with version in comments

library(dplyr) # 1.0.10
library(pec) # 2022.05.04
library(randomForestSRC) # 3.2.3
library(riskRegression) # 2022.11.28
library(survival) # 3.5-0

# Functions

#' compRiskModelsTr trains CSC and FG models from a dataset.
#' 
#' @param data_set A dataframe with at least this variables: los (length of stay);
#' Deterioro (a binary 1/0 variable, in this case means deterioration); 
#' one or more covariates. 
#' If dataframe contains some ID, it must be called "Episodio (único)".
#' @param cause Cause of interest.
#' @returns A list with two elements: model_csc, the CSC trained model; 
#' model_fgr, the FG model trained.
#' @examples
#' data_set <- data.frame(los=rexp(100,1/50),
#' Deterioro=rbinom(100,1,0.2),
#' cov1=rnorm(100,10,2),
#' cov2=rnorm(100,100,10))
#' compRiskModelsTr(data_set,1)
compRiskModelsTr <- function(data_set,cause){
  if("Episodio (único)" %in% names(data_set)){
    data_set <- data_set %>% dplyr::select(-`Episodio (único)`)
  }
  data_set$los <- as.numeric(data_set$los)
  data_set$Deterioro <- as.factor(data_set$Deterioro)
  data_set$Deterioro <- as.numeric(data_set$Deterioro)
  form_surv <- formula(paste0("Hist(time=los, event=Deterioro) ~ ",paste(names(data_set)[
    !grepl("Deterioro|los",names(data_set))
  ],collapse = " + ")))
  model_csc <- CSC(form_surv,
                   data=data_set)
  model_fgr <- FGR(form_surv,cause=cause,
                   data=data_set)
  return(list(model_csc=model_csc,model_fgr=model_fgr))
}

###########################
##     Read datasets     ##
###########################

data_sum_imputed_by_mice <- readRDS(paste0(root_output,"/data_sum_imputed_by_mice.rds"))
data_ts_sum_imputed_by_mice <- readRDS(paste0(root_output,"/data_ts_sum_imputed_by_mice.rds"))

var_def_lasso <- readRDS(paste0(root_output,"/var_def_lasso.rds"))
var_def_bess <- readRDS(paste0(root_output,"/var_def_bess.rds"))

#####################################
##     Competing risk training     ##
#####################################

#### 1. Full data-set ####

## 1.1 Event: Favourable discharge

## Linear models (Both CSC and FG)

compRiskModels_imp_mice_disc <- compRiskModelsTr(data_sum_imputed_by_mice,
                                                 cause=1)
saveRDS(compRiskModels_imp_mice_disc,
        paste0(root_output,"/compRiskModels_imp_mice_disc.rds"))

## Ensemble learning

# Prepare data-set

data_for_rf_full <- data_sum_imputed_by_mice %>% 
  dplyr::select(-`Episodio (único)`) %>% 
  mutate_if(is.character,as.factor)
data_for_rf_full$Deterioro <- as.numeric(data_for_rf_full$Deterioro)

# Options for parallel computing

rf.cores_rsf = 8
mc.cores_rsf = 8

options(rf.cores=rf.cores_rsf, mc.cores=mc.cores_rsf)

# LR-RSF 

data_for_rf_full2 <- data_for_rf_full
data_for_rf_full2$Deterioro <- as.factor(as.numeric(data_for_rf_full$Deterioro) + 1)

# Tuning hyper-parameters

rf_tune2 <- tune(Surv(los, Deterioro) ~ ., cause = c(1,0),
                 splitrule = "logrank", 
                 data = data_for_rf_full2, 
                 trace = TRUE,
                 ntreeTry = 100, seed = 1234)

# Optimal
# nodesize     mtry 
#    1          17

# Traing the Best RF based on hyper-param obtained

arfsurv.obj_disch_full <- rfsrc(Surv(los, Deterioro) ~ ., 
                                data = data_for_rf_full2,
                                cause = c(1,0),
                                splitrule = "logrank",
                                mtry = 17,
                                nodesize = 1,
                                importance = TRUE,
                                ntree = 100,
                                seed = 1234,
                                ntime=24*(2:15))  

saveRDS(arfsurv.obj_disch_full,
        paste0(root_output,"/arfsurv.obj_disch_full.rds"))

# GT-RSF
# For both fav. discharge and deterioration is just one training.

# Tuning

rf_tune_fg <- tune(Surv(los, Deterioro) ~ ., cause = c(1,1),
                   splitrule = "logrankCR", 
                   data = data_for_rf_full2, 
                   trace = TRUE,
                   ntreeTry = 100,
                   seed=1234)
# Optimal
# nodesize     mtry 
#    1          21 

# Training the Best RF based on hyper-param obtained

arfsurv.obj_disch_full_fg <- rfsrc(Surv(los, Deterioro) ~ ., 
                                   data = data_for_rf_full2,
                                   cause = c(1,1),
                                   #splitrule = "logrank",
                                   splitrule = "logrankCR",
                                   mtry = 21,
                                   nodesize = 1,
                                   importance = TRUE,
                                   ntree = 100,
                                   seed = 1234,
                                   ntime=24*(2:15))  

saveRDS(arfsurv.obj_disch_full_fg,
        paste0(root_output,"/arfsurv.obj_disch_full_fg.rds"))

# Null model: LR-RSF
# Prepare data-set

data_for_rf_full2_null <- data_for_rf_full2 %>% 
  dplyr::select(los,Deterioro)
data_for_rf_full2_null$Base <- 1

arfsurv.obj_disch_full_null <- rfsrc(Surv(los, Deterioro) ~ ., 
                                     data = data_for_rf_full2_null,
                                     cause = c(1,0),
                                     splitrule = "logrank",
                                     mtry = 1,
                                     nodesize = 4,
                                     importance = TRUE,
                                     ntree = 100,
                                     seed = 1234,
                                     ntime=24*(2:15))

saveRDS(arfsurv.obj_disch_full_null,
        paste0(root_output,"/arfsurv.obj_disch_full_null.rds"))

# Null model: GT-RSF
# For both fav. discharge and deterioration is just one training.

arfsurv.obj_disch_full_null_fg <- rfsrc(Surv(los, Deterioro) ~ ., 
                                        data = data_for_rf_full2_null,
                                        cause = c(1,1),
                                        splitrule = "logrankCR",
                                        mtry = 1,
                                        nodesize = 4,
                                        importance = TRUE,
                                        ntree = 100,
                                        seed = 1234,
                                        ntime=24*(2:15))

saveRDS(arfsurv.obj_disch_full_null_fg,
        paste0(root_output,"/arfsurv.obj_disch_full_null_fg.rds"))

## 1.2. Event: Deterioration

# Linear models

compRiskModels_imp_mice_deter <- compRiskModelsTr(data_sum_imputed_by_mice,
                                                  cause=2)
saveRDS(compRiskModels_imp_mice_deter,
        paste0(root_output,"/compRiskModels_imp_mice_deter.rds"))

# Ensemble learning

# LR-RSF tuning

rf_tune <- tune(Surv(los, Deterioro) ~ ., cause = c(0,1),splitrule = "logrank", 
                data = data_for_rf_full2, trace = TRUE,
                ntreeTry = 100,seed=1234)

# Optimal
# nodesize     mtry 
#    3          50
# Best LR-RSF based on tuning hyper-param

arfsurv.obj_deter_full <- rfsrc(Surv(los, Deterioro) ~ ., 
                                cause = c(0,1),
                                splitrule = "logrank", 
                                data = data_for_rf_full2,
                                mtry = 50,
                                nodesize = 3,
                                importance = TRUE,
                                ntree = 100,
                                seed = 1234,
                                ntime=24*(2:15))

saveRDS(arfsurv.obj_deter_full,
        paste0(root_output,"/arfsurv.obj_deter_full.rds"))

# Null model: LR-RSF

arfsurv.obj_deter_full_null <- rfsrc(Surv(los, Deterioro) ~ ., 
                                     data = data_for_rf_full2_null,
                                     cause = c(0,1),
                                     splitrule = "logrank",
                                     mtry = 50,
                                     nodesize = 3,
                                     importance = TRUE,
                                     ntree = 10,
                                     seed = 1234,
                                     ntime=24*(2:15))

saveRDS(arfsurv.obj_deter_full_null,
        paste0(root_output,"/arfsurv.obj_deter_full_null.rds"))

#### 2. Reduced LASSO data-set ####

# Select the variables selected by LASSO
vars_comp_lasso <- var_def_lasso %>% dplyr::filter(!is.na(Var_Det) | !is.na(Var_Alta)) %>% 
  dplyr::select(Variables)

### 2.1. Event: Favourable discharge

## Linear models (Both CSC and FG)

compRiskModels_imp_mice_lasso_disc <- compRiskModelsTr(data_sum_imputed_by_mice %>% 
                                                         dplyr::select(as.vector(t(vars_comp_lasso))),
                                                       cause=1)
saveRDS(compRiskModels_imp_mice_lasso_disc,
        paste0(root_output,"/compRiskModels_imp_mice_lasso_disc.rds"))

## Ensemble learning

# LR-RSF 

# Prepare data-set

data_for_rf_lasso <- data_sum_imputed_by_mice %>% 
  dplyr::select(-`Episodio (único)`) %>% 
  dplyr::select(as.vector(t(vars_comp_lasso))) %>% 
  mutate_if(is.character,as.factor)
data_for_rf_lasso$Deterioro <- as.numeric(data_for_rf_lasso$Deterioro)
data_for_rf_lasso2 <- data_for_rf_lasso
data_for_rf_lasso2$Deterioro <- as.factor(as.numeric(data_for_rf_lasso2$Deterioro)+1)

# Parallel options

options(rf.cores=rf.cores_rsf, mc.cores=mc.cores_rsf)

# Tuning hyper-parameters for LR-RSF

rf_tune_lasso2 <- tune(Surv(los, Deterioro) ~ ., 
                       cause = c(1,0),
                       splitrule = "logrank", 
                       data = data_for_rf_lasso2, trace = TRUE,
                       ntreeTry = 100,seed = 1234)
# Optimal
# nodesize     mtry 
#    10         12

# Best LR-RSF based on tuning

arfsurv.obj_lasso_disc <- rfsrc(Surv(los, Deterioro) ~ ., 
                                cause = c(1,0),
                                splitrule = "logrank", 
                                data = data_for_rf_lasso2,
                                mtry = 12,
                                nodesize = 10,
                                importance = TRUE,
                                ntree = 100,
                                seed = 1234,
                                ntime=24*(2:15))  

saveRDS(arfsurv.obj_lasso_disc,
        paste0(root_output,"/arfsurv.obj_lasso_disc.rds"))

# Tuning hyper-parameters for GT-RSF

rf_tune_lasso2_fg <- tune(Surv(los, Deterioro) ~ ., 
                          cause = c(1,1),
                          splitrule = "logrankCR", 
                          data = data_for_rf_lasso2, 
                          trace = TRUE,
                          ntreeTry = 100,
                          seed = 1234)
# Optimal
# nodesize     mtry 
#    6          12 

# Best GT-RSF based on tuning

arfsurv.obj_lasso_disc_fg <- rfsrc(Surv(los, Deterioro) ~ ., 
                                   splitrule = "logrankCR", 
                                   data = data_for_rf_lasso2,
                                   mtry = 12,
                                   nodesize =6,
                                   importance = TRUE,
                                   seed = 1234,
                                   ntree = 100,
                                   ntime=24*(2:15))  
saveRDS(arfsurv.obj_lasso_disc_fg,
        paste0(root_output,"/arfsurv.obj_lasso_disc_fg.rds"))

### 2.2. Event: Deterioration

## Linear models

compRiskModels_imp_mice_lasso_deter <- compRiskModelsTr(data_sum_imputed_by_mice %>% 
                                                          dplyr::select(as.vector(t(vars_comp_lasso))),
                                                        cause=2)

saveRDS(compRiskModels_imp_mice_lasso_deter,
        paste0(root_output,"/compRiskModels_imp_mice_lasso_deter.rds"))

## Ensemble learning
# LR-RSF 

# Tuning hyper-parameters for GT-RSF

rf_tune_lasso <- tune(Surv(los, Deterioro) ~ ., cause = c(0,1),splitrule = "logrank", 
                      data = data_for_rf_lasso2, trace = TRUE,
                      ntreeTry = 100,seed = 1234)

# Optimal 
# nodesize     mtry 
# 1       32

# Best LR-RSF based on tuning

arfsurv.obj_lasso_deter <- rfsrc(Surv(los, Deterioro) ~ ., cause = c(0,1),
                                 splitrule = "logrank", 
                                 data = data_for_rf_lasso2,
                                 mtry = 32,
                                 nodesize = 1,
                                 importance = TRUE,
                                 ntree = 100,
                                 seed = 1234,
                                 ntime=24*(2:15))  

saveRDS(arfsurv.obj_lasso_deter,
        paste0(root_output,"/arfsurv.obj_lasso_deter.rds"))

#### 3. Reduced BeSS data-set ####

# Select the variables selected by LASSO

vars_comp_bess <- var_def_bess %>% dplyr::filter(!is.na(Var_Det) | !is.na(Var_Alta)) %>% 
  dplyr::select(Variables)

### 3.1. Event: Deterioration

## Linear models (Both CSC and FG)

compRiskModels_imp_mice_bess_disc <- compRiskModelsTr(data_sum_imputed_by_mice %>% 
                                                        dplyr::select(as.vector(t(vars_comp_bess))),
                                                      cause=1)
saveRDS(compRiskModels_imp_mice_bess_disc,
        paste0(root_output,"/compRiskModels_imp_mice_bess_disc.rds"))

## Ensemble learning

# Prepare data-set

data_for_rf_bess <- data_sum_imputed_by_mice %>% 
  dplyr::select(-`Episodio (único)`) %>% 
  dplyr::select(as.vector(t(vars_comp_bess))) %>% 
  mutate_if(is.character,as.factor)
data_for_rf_bess$Deterioro <- as.numeric(data_for_rf_bess$Deterioro)
data_for_rf_bess2 <- data_for_rf_bess
data_for_rf_bess2$Deterioro <- as.factor(as.numeric(data_for_rf_bess2$Deterioro)+1)

# Tuning hyper-parameters for LR-RSF

rf_tune_bess2 <- tune(Surv(los, Deterioro) ~ ., cause = c(1,0),
                      splitrule = "logrank", 
                      data = data_for_rf_bess2, trace = TRUE,
                      ntreeTry = 100,seed = 1234)

# Optimal
# nodesize     mtry 
#    1          21

# Best LR-RSF based on tuning

arfsurv.obj_bess2_disc <- rfsrc(Surv(los, Deterioro) ~ ., 
                                cause = c(1,0),
                                splitrule = "logrank", 
                                data = data_for_rf_bess2,
                                mtry = 21,
                                nodesize = 1,
                                importance = TRUE,
                                ntree = 100,
                                seed = 1234,
                                ntime=24*(2:15))  

saveRDS(arfsurv.obj_bess2_disc,
        paste0(root_output,"/arfsurv.obj_bess2_disc.rds"))

# Tuning hyper-parameters for GT-RSF

rf_tune_bess2_fg <- tune(Surv(los, Deterioro) ~ ., cause = c(1,1),
                         splitrule = "logrankCR", 
                         data = data_for_rf_bess2, trace = TRUE,
                         ntreeTry = 100,seed = 1234)

# Optimal
# nodesize     mtry 
#    3          14

# Best GT-RSF based on tuning

arfsurv.obj_bess2_disc_fg <- rfsrc(Surv(los, Deterioro) ~ ., 
                                   splitrule = "logrankCR", 
                                   data = data_for_rf_bess2,
                                   mtry = 14,
                                   nodesize = 3,
                                   importance = TRUE,
                                   ntree = 100,
                                   seed = 1234,
                                   ntime=24*(2:15)) 

saveRDS(arfsurv.obj_bess2_disc_fg,
        paste0(root_output,"/arfsurv.obj_bess2_disc_fg.rds"))

### 3.1. Event: Deterioration

## Linear models (Both CSC and FG)

compRiskModels_imp_mice_bess_deter <- compRiskModelsTr(data_sum_imputed_by_mice %>% 
                                                         dplyr::select(as.vector(t(vars_comp_bess))),
                                                       cause=2)
saveRDS(compRiskModels_imp_mice_bess_deter,
        paste0(root_output,"/compRiskModels_imp_mice_bess_deter.rds"))

## Ensemble learning

# Tuning hyper-parameters for LR-RSF

rf_tune_bess <-  tune(Surv(los, Deterioro) ~ ., cause = c(0,1),
                      splitrule = "logrank", 
                      data = data_for_rf_bess2, trace = TRUE,
                      ntreeTry = 100,seed=1234)

# Optimal
# nodesize     mtry 
#    1          7 

# Best LR-RSF based on tuning

arfsurv.obj_bess_deter <- rfsrc(Surv(los, Deterioro) ~ ., 
                                cause = c(0,1),
                                splitrule = "logrank", 
                                data = data_for_rf_bess2,
                                mtry = 7,
                                nodesize = 1,
                                importance = TRUE,
                                ntree = 100,
                                seed = 1234,
                                ntime=24*(2:15))  

saveRDS(arfsurv.obj_bess_deter,
        paste0(root_output,"/arfsurv.obj_bess_deter.rds"))

