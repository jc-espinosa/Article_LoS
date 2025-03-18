
# Roots

root_input <- "./Datasets/input"
root_output <- "./Datasets/output"

# Libraries with version in comments

library(dplyr) # 1.0.10
library(pec) # 2022.05.04
library(randomForestSRC) # 3.2.3
library(riskRegression) # 2022.11.28
library(survival) # 3.5-0
library(timeROC) # 0.4

# Functions

#' compRiskModelsTest test linear models CSC and FG, previously trained.
#' It calculates predictions from models, IBS, BS and C-Index.
#' 
#' @param data_set A test dataframe with at least this variables: los (length of stay);
#' Deterioro (a binary 1/0 variable, in this case means deterioration); 
#' one or more covariates. 
#' If dataframe contains some ID, it must be called "Episodio (único)".
#' @param train_mod a list obtained from compRiskModelsTr function
#' (see Competing_Risk_Training.R).
#' @param cause Cause of interest.
#' @param times Time window to evaluate IBS, BS and C-Index.
#' @returns A list with: predictions in pred_csc amd pred_fgr;
#' BS in IPA_csc and IPA_fgr; C-Index in c_index; 
#' IBS in ibs_csc, ibs_fg and ibs_null.
#' @examples
#' source("Competing_Risk_Training.R")
#' set.seed(123)
#' data_set_train <- data.frame(los=rexp(5000,1/50),
#'                              Deterioro=rbinom(5000,1,0.2),
#'                              cov1=rnorm(5000,10,2),
#'                              cov2=runif(5000,100,1000))
#' cr_train <- compRiskModelsTr(data_set_train,2)
#' times <- rexp(1000,1/50)
#' data_set_test <- data.frame(los=times,
#'                          Deterioro=rbinom(1000,1,0.2),
#'                          cov1=rnorm(1000,10,2),
#'                          cov2=runif(1000,100,1000))
#' compRiskModelsTest(data_set_test,cr_train,cause=2,times=quantile(times))
compRiskModelsTest <- function(data_set,train_mod,cause=2,times=24*(2:15)){
  #browser()
  if("Episodio (único)" %in% names(data_set)){
    data_set <- data_set %>% dplyr::select(-`Episodio (único)`)
  }
  data_set$los <- as.numeric(data_set$los)
  data_set$Deterioro <- as.factor(data_set$Deterioro)
  data_set$Deterioro <- as.numeric(data_set$Deterioro)
  form_surv <- formula(paste0("Hist(time=los, event=Deterioro) ~ ",paste(names(data_set)[
    !grepl("Deterioro|los",names(data_set))
  ],collapse = " + ")))
  IPA_csc <- IPA(train_mod$model_csc,form_surv,
                 # null.model = T,
                 newdata = data_set,times=times,
                 cause=cause)
  IPA_fgr <- IPA(train_mod$model_fgr,form_surv,newdata = data_set,
                 times=times,cause=cause)
  c_index = pec::cindex(list(train_mod$model_csc,train_mod$model_fgr),
                        form_surv,data=data_set,eval.times=times,
                        cause=cause)
  ibs <-  crps(pec::pec(list(train_mod$model_csc,train_mod$model_fgr),
                        formula=form_surv,
                        data=data_set,
                        cause=cause,
                        start=24*2,
                        maxtime=24*15,
                        exactness=13,
                        exact=F),
               times=c(24*15),
               start=24*2)
  ibs_csc <- ibs[2,1]
  ibs_fg <- ibs[3,1]
  ibs_null <- ibs[1,1]
  pred_csc <- predictRisk(train_mod$model_csc,cause=cause,
                          newdata=data_set,times=times)
  pred_fgr <- predictRisk(train_mod$model_fgr,cause=cause,
                          newdata=data_set,times=times)
  return(list(IPA_csc=data.frame(IPA_csc) %>% dplyr::select(-IPA,-IPA.drop),
              IPA_fgr=data.frame(IPA_fgr) %>% dplyr::select(-IPA),
              c_index=c_index$AppCindex,
              pred_csc=pred_csc,
              pred_fgr=pred_fgr,
              ibs_csc=ibs_csc,
              ibs_fg=ibs_fg,
              ibs_null=ibs_null))
}

#' brierCR calculates BS for Competing risk models without censored observations
#' 
#' @param time_samp The times until event observed in data-set.
#' @param t The time of interest to evaluate BS.
#' @param event_samp The endpoint for each subject (Competing risk).
#' @param event_int The event (cause) of interest.
#' @param cif The estimated CIF.
#' @returns A numeric value with the BS.
#' @examples 
#' set.seed(1234)
#' time_samp <- rexp(1000,1/50)
#' t <- median(time_samp)
#' event_samp <- rbinom(1000,1,0.2)
#' event_int <- 1
#' cif <- rbeta(1000,1,0.5)
#' brierCR(time_samp,t,event_samp,event_int,cif)
brierCR <- function(time_samp,t,event_samp,event_int,cif){
  time_samp <- as.double(time_samp)
  Ind <- ifelse(time_samp <= t & event_samp == event_int,1,0) 
  BS <- mean((Ind-cif)^2)
  return(BS)
}

#################################
##     Read trained models     ##
#################################

### Models over Full data-set ###

# CSC-Full and FG-Full favourable discharge
compRiskModels_imp_mice_disc <- readRDS(paste0(root_output,"/compRiskModels_imp_mice_disc.rds"))
# CSC-Full and FG-Full deterioration
compRiskModels_imp_mice_deter <- readRDS(paste0(root_output,"/compRiskModels_imp_mice_deter.rds"))
# LR-RSF-Full favourable discharge
arfsurv.obj_disch_full <- readRDS(paste0(root_output,"/arfsurv.obj_disch_full.rds"))
# LR-RSF-Full deterioration
arfsurv.obj_deter_full <- readRDS(paste0(root_output,"/arfsurv.obj_deter_full.rds"))
# GT-RSF-Full favourable discharge and deterioration
arfsurv.obj_disch_full_fg <- readRDS(paste0(root_output,"/arfsurv.obj_disch_full_fg.rds"))

### Models over LASSO data-set ###

# CSC-LASSO and FG-LASSO favourable discharge
compRiskModels_imp_mice_lasso_disc <- readRDS(paste0(root_output,"/compRiskModels_imp_mice_lasso_disc.rds"))
# CSC-LASSO and FG-LASSO deterioration
compRiskModels_imp_mice_lasso_deter <- readRDS(paste0(root_output,"/compRiskModels_imp_mice_lasso_deter.rds"))
# LR-RSF-LASSO favourable discharge
arfsurv.obj_lasso_disc <- readRDS(paste0(root_output,"/arfsurv.obj_lasso_disc.rds"))
# LR-RSF-LASSO deterioration
arfsurv.obj_lasso_deter <- readRDS(paste0(root_output,"/arfsurv.obj_lasso_deter.rds"))
# GT-RSF-LASSO favourable discharge and deterioration
arfsurv.obj_lasso_disc_fg <- readRDS(paste0(root_output,"/arfsurv.obj_lasso_disc_fg.rds"))

### Models over BeSS data-set ###

# CSC-BeSS and FG-BeSS favourable discharge
compRiskModels_imp_mice_bess_disc <- readRDS(paste0(root_output,"/compRiskModels_imp_mice_bess_disc.rds"))
# CSC-BeSS and FG-BeSS deterioration
compRiskModels_imp_mice_bess_deter <- readRDS(paste0(root_output,"/compRiskModels_imp_mice_bess_deter.rds"))
# LR-RSF-BeSS favourable discharge
arfsurv.obj_bess2_disc <- readRDS(paste0(root_output,"/arfsurv.obj_bess2_disc.rds"))
# LR-RSF-BeSS deterioration
arfsurv.obj_bess_deter <- readRDS(paste0(root_output,"/arfsurv.obj_bess_deter.rds"))
# GT-RSF-BeSS favourable discharge and deterioration
arfsurv.obj_bess2_disc_fg <- readRDS(paste0(root_output,"/arfsurv.obj_bess2_disc_fg.rds"))

### Models over Null dataset  ###

# LR-RSF-Null favourable discharge
arfsurv.obj_disch_full_null <- readRDS(paste0(root_output,"/arfsurv.obj_disch_full_null.rds"))
# LR-RSF-Null deterioration
arfsurv.obj_deter_full_null <- readRDS(paste0(root_output,"/arfsurv.obj_deter_full_null.rds"))
# GT-RSF-Null favourable discharge and deterioration
arfsurv.obj_disch_full_null_fg <- readRDS(paste0(root_output,"/arfsurv.obj_disch_full_null_fg.rds"))

#####################################
##     Competing risk training     ##
#####################################

######### 1. Full data-set #########

## 1.1 Event: Favourable discharge

## Linear models (Both CSC and FG)

compRiskModels_mice_test_full_disc <- compRiskModelsTest(data_ts_sum_imputed_by_mice,
                                                         compRiskModels_imp_mice_disc,
                                                         cause = 1)
# IBS CSC
compRiskModels_mice_test_full_disc$ibs_csc
# IBS FG
compRiskModels_mice_test_full_disc$ibs_fg
# IBS Null
compRiskModels_mice_test_full_disc$ibs_null

# Brier Score CSC
brier_csc_ts_disc <- compRiskModels_mice_test_full_disc$IPA_csc %>% filter(Variable=="Full model")
# Brier Score FG
brier_fgr_ts_disc <- compRiskModels_mice_test_full_disc$IPA_fgr %>% filter(model=="FGR")
# Brier Score CSC Null model
brier_csc_ts_disc_null <- compRiskModels_mice_test_full_disc$IPA_csc %>% filter(Variable=="Null model")
# Brier Score FG Null model
brier_fgr_ts_disc_null <- compRiskModels_mice_test_full_disc$IPA_fgr %>% filter(model=="Null model")

# C-Index CSC
cindex_csc_ts_disc <- compRiskModels_mice_test_full_disc$c_index$CauseSpecificCox
# C-Index FG
cindex_fgr_ts_disc <- compRiskModels_mice_test_full_disc$c_index$FGR
# Cumulative C-Index CSC
sum(cindex_csc_ts_disc)
# Cumulative C-Index FG
sum(cindex_fgr_ts_disc)

# Time-dependent AUROC

# Ensure states starting with 1 (since we do not have censored data)
cause_cr <- data_ts_sum_imputed_by_mice$Deterioro + 1
# Recover time covariate
time_cr <- data_ts_sum_imputed_by_mice$los
auroc_t_FG_Full_disc <- rep(0,14)

for(i in 1:14){
  auroc_temp <- timeROC(T=time_cr,
                        delta=cause_cr,
                        weighting="marginal",
                        marker=compRiskModels_mice_test_full_disc$pred_fgr[,i],
                        cause=1,
                        times=c(i+1)*24)
  auroc_t_FG_Full_disc[i] <- as.numeric(auroc_temp$AUC_2[2])
}

## Ensemble learning

# Prepare data-set 

data_ts_for_rf_full <- data_ts_sum_imputed_by_mice %>% 
  dplyr::select(-`Episodio (único)`) %>% 
  mutate_if(is.character,as.factor)
data_ts_for_rf_full2 <- data_ts_for_rf_full
data_ts_for_rf_full2$Deterioro <- as.factor(as.numeric(data_ts_for_rf_full2$Deterioro)+1)

# Prediction for LR-RSF

pred_rf_full_alt <- predict.rfsrc(arfsurv.obj_disch_full,data_ts_for_rf_full2,
                                  outcome="test")

# Time-dependent AUROC

auroc_t_LR_RSF_Full_disc <- rep(0,14)

for(i in 1:14){
  auroc_temp <- timeROC(T=time_cr,
                        delta=cause_cr,
                        weighting="marginal",
                        marker=pred_rf_full_alt$cif[,i,1],
                        cause=1,
                        times=c(i+1)*24)
  auroc_t_LR_RSF_Full_disc[i] <- as.numeric(auroc_temp$AUC_1[2])
}

# Generic formula for RSF

form_rsf_full_alt <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                    paste(names(data_ts_for_rf_full2)[!grepl("Deterioro|los",
                                                                             names(data_ts_for_rf_full2))
                                    ],collapse = " + ")))
# IBS LR-RSF

ibs_full_alt <-  crps(pec::pec(as.matrix(pred_rf_full_alt$cif[,,1]),
                               formula=form_rsf_full_alt,
                               data=data_ts_for_rf_full2,
                               cause=1,
                               start=24*2,
                               maxtime=24*15,
                               exactness=13,
                               exact=F),
                      times=c(24*2,24*15),
                      start=24*2)
ibs_full_alt <- ibs_full_alt[2,2]

# Prediction for GT-RSF

pred_rf_full_alt_fg <- predict.rfsrc(arfsurv.obj_disch_full_fg,data_ts_for_rf_full2,
                                     outcome="test")

# Time-dependent AUROC

auroc_t_GT_RSF_Full_disc <- rep(0,14)

for(i in 1:14){
  auroc_temp <- timeROC(T=time_cr,
                        delta=cause_cr,
                        weighting="marginal",
                        marker=pred_rf_full_alt_fg$cif[,i,1],
                        cause=1,
                        times=c(i+1)*24)
  auroc_t_GT_RSF_Full_disc[i] <- as.numeric(auroc_temp$AUC_2[2])
}

# IBS GT-RSF

ibs_full_alt_fg <-  crps(pec::pec(as.matrix(pred_rf_full_alt_fg$cif[,,1]),
                                  formula=form_rsf_full_alt,
                                  data=data_ts_for_rf_full2,
                                  cause=1,
                                  start=24*2,
                                  maxtime=24*15,
                                  exactness=13,
                                  exact=F),
                         times=c(24*2,24*15),
                         start=24*2)
ibs_full_alt_fg <- ibs_full_alt_fg[2,2]

# LR-RSF Null

data_ts_for_rf_full2_null <- data_ts_for_rf_full2 %>% 
  dplyr::select(Deterioro,los)
data_ts_for_rf_full2_null$Base <- 1

pred_rf_full_null <- predict.rfsrc(arfsurv.obj_disch_full_null,data_ts_for_rf_full2_null,
                                   outcome="test")

# IBS LR-RSF-Null

form_rsf_full_null <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                     paste(names(data_ts_for_rf_full2_null)[!grepl("Deterioro|los",
                                                                                   names(data_ts_for_rf_full2_null))
                                     ],collapse = " + ")))
ibs_full_alt_null <-  crps(pec::pec(as.matrix(pred_rf_full_null$cif[,,1]),
                                    formula=form_rsf_full_null,
                                    data=data_ts_for_rf_full2_null,
                                    cause=1,
                                    start=24*2,
                                    maxtime=24*15,
                                    exactness=13,
                                    exact=F),
                           times=c(24*2,24*15),
                           start=24*2)
ibs_full_alt_null <- ibs_full_alt_null[2,2]

# IBS GT-RSF-Null

pred_rf_full_null_fg <- predict.rfsrc(arfsurv.obj_disch_full_null_fg,data_ts_for_rf_full2_null,
                                      outcome="test")
form_rsf_null_alt_fg <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                       paste(names(data_ts_for_rf_full2_null)[!grepl("Deterioro|los",
                                                                                     names(data_ts_for_rf_full2_null))
                                       ],collapse = " + ")))
ibs_null_alt_fg <-  crps(pec::pec(as.matrix(pred_rf_full_null_fg$cif[,,1]),
                                  formula=form_rsf_null_alt_fg,
                                  data=data_ts_for_rf_full2_null,
                                  cause=1,
                                  start=24*2,
                                  maxtime=24*15,
                                  exactness=13,
                                  exact=F),
                         times=c(24*2,24*15),
                         start=24*2)

ibs_null_alt_fg <- ibs_null_alt_fg[2,2]

# Brier Score and C-Index 

# Vectors to allocate info

times <- 24*c(2:15)
c_ind_disc <- rep(0,14)
bs.rsf_full_alt <- rep(0,14)
c_ind_disc_fg <- rep(0,14)
bs.rsf_full_alt_fg <- rep(0,14)
bs.rsf_full_alt_null <- rep(0,14)
bs.rsf_full_alt_null_fg <- rep(0,14)

# Calculus of BS and C-Index

for(i in 1:14){
  # BS for LR-RSF
  bs.rsf_full_alt[i] <- brierCR(time_samp=data_ts_for_rf_full2$los,
                                t=times[i],
                                event_samp=data_ts_for_rf_full2$Deterioro,
                                event_int=1,   
                                cif=pred_rf_full_alt$cif[,i,1])
  # BS for GT-RSF
  bs.rsf_full_alt_fg[i] <- brierCR(time_samp=data_ts_for_rf_full2$los,
                                   t=times[i],
                                   event_samp=data_ts_for_rf_full2$Deterioro,
                                   event_int=1,   
                                   cif=pred_rf_full_alt_fg$cif[,i,1])
  # BS for LR-RSF-Null
  bs.rsf_full_alt_null[i] <- brierCR(time_samp=data_ts_for_rf_full2$los,
                                     t=times[i],
                                     event_samp=data_ts_for_rf_full2$Deterioro,
                                     event_int=1,   
                                     cif=pred_rf_full_null$cif[,i,1])
  # BS for GT-RSF-Null
  bs.rsf_full_alt_null_fg[i] <- brierCR(time_samp=data_ts_for_rf_full2$los,
                                        t=times[i],
                                        event_samp=data_ts_for_rf_full2$Deterioro,
                                        event_int=1,   
                                        cif=pred_rf_full_null_fg$cif[,i,1])
  # C-index of LR-RSF
  c_ind_disc[i] <- as.double(pec::cindex(matrix(pred_rf_full_alt$cif[,i,1]),
                                         Surv(los, Deterioro) ~ .,
                                         data_ts_for_rf_full2,
                                         cause=1,eval.times=24*(i+1))$AppCindex)
  # C-index of GT-RSF
  c_ind_disc_fg[i] <- as.double(pec::cindex(matrix(pred_rf_full_alt_fg$cif[,i,1]),
                                            Surv(los, Deterioro) ~ .,
                                            data_ts_for_rf_full2,
                                            cause=1,eval.times=24*(i+1))$AppCindex)
  print(i)
}

# CC-Index of LR-RSF
sum(c_ind_disc)
# CC-Index of GT-RSF
sum(c_ind_disc)

# Table with BS and C-Index for favourable discharge
# with Full data-set

full_data_discharge_cr_c_br <- data.frame(times=brier_csc_ts_disc$times,
                                          csc_brier=brier_csc_ts_disc$Brier,
                                          fgr_brier=brier_fgr_ts_disc$Brier,
                                          rf_brier=bs.rsf_full_alt,
                                          csc_c_index=cindex_csc_ts_disc,
                                          fgr_c_index=cindex_fgr_ts_disc,
                                          rf_c_index=c_ind_disc,
                                          csc_brier_null=brier_csc_ts_disc_null$Brier,
                                          fgr_brier_null=brier_fgr_ts_disc_null$Brier,
                                          rf_c_index_fg=c_ind_disc_fg,
                                          rf_brier_fg=bs.rsf_full_alt_fg,
                                          rsf_brier_null=bs.rsf_full_alt_null,
                                          rsf_brier_null_fg=bs.rsf_full_alt_null_fg
)


## 1.2 Event: Deterioration

## Linear models (Both CSC and FG)

compRiskModels_mice_test_full_deter <- compRiskModelsTest(data_ts_sum_imputed_by_mice,
                                                          compRiskModels_imp_mice_deter,
                                                          cause = 2)
# IBS CSC
compRiskModels_mice_test_full_deter$ibs_csc
# IBS FG
compRiskModels_mice_test_full_deter$ibs_fg
# IBS Null
compRiskModels_mice_test_full_deter$ibs_null

# BS CSC
brier_csc_ts_deter <- compRiskModels_mice_test_full_deter$IPA_csc %>% filter(Variable=="Full model")
# BS FG
brier_fgr_ts_deter <- compRiskModels_mice_test_full_deter$IPA_fgr %>% filter(model=="FGR")
# BS CSC-Null
brier_csc_ts_deter_null <- compRiskModels_mice_test_full_deter$IPA_csc %>% filter(Variable=="Null model")
# BS FG-Null
brier_fgr_ts_deter_null <- compRiskModels_mice_test_full_deter$IPA_fgr %>% filter(model=="Null model")

# C-Index CSC 
cindex_csc_ts_deter <- compRiskModels_mice_test_full_deter$c_index$CauseSpecificCox
# C-Index FG
cindex_fgr_ts_deter <- compRiskModels_mice_test_full_deter$c_index$FGR

# CC-Index CSC
sum(cindex_csc_ts_deter)

# CC-Index FG
sum(cindex_fgr_ts_deter)

## Ensemble learning

# Prediction

pred_rf_full_det <- predict.rfsrc(arfsurv.obj_deter_full,data_ts_for_rf_full2)

# Time-dependent AUROC LR-RSF

auroc_t_LR_RSF_Full_deter <- rep(0,14)
for(i in 1:14){
  auroc_temp <- timeROC(T=time_cr,
                        delta=cause_cr,
                        weighting="marginal",
                        marker=pred_rf_full_det$cif[,i,2],
                        cause=2,
                        times=c(i+1)*24)
  auroc_t_LR_RSF_Full_deter[i] <- as.numeric(auroc_temp$AUC_1[2])
}

# Time-dependent AUROC GT-RSF

auroc_t_GT_RSF_Full_deter <- rep(0,14)
for(i in 1:14){
  auroc_temp <- timeROC(T=time_cr,
                        delta=cause_cr,
                        weighting="marginal",
                        marker=pred_rf_full_alt_fg$cif[,i,2],
                        cause=2,
                        times=c(i+1)*24)
  auroc_t_GT_RSF_Full_deter[i] <- as.numeric(auroc_temp$AUC_2[2])
}

# Generic formula for RSF
form_rsf_full_det <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                    paste(names(data_ts_for_rf_full2)[!grepl("Deterioro|los",
                                                                             names(data_ts_for_rf_full2))
                                    ],collapse = " + ")))

# IBS LR-RSF

ibs_full_alt_det <-  crps(pec::pec(as.matrix(pred_rf_full_det$cif[,,2]),
                                   formula=form_rsf_full_det,
                                   data=data_ts_for_rf_full2,
                                   cause=2,
                                   start=24*2,
                                   maxtime=24*15,
                                   exactness=13,
                                   exact=F),
                          times=c(24*2,24*15),
                          start=24*2)
ibs_full_alt_det <- ibs_full_alt_det[2,2]

# # IBS GT-RSF

ibs_full_deter_fg <-  crps(pec::pec(as.matrix(pred_rf_full_alt_fg$cif[,,2]),
                                    formula=form_rsf_full_alt,
                                    data=data_ts_for_rf_full2,
                                    cause=2,
                                    start=24*2,
                                    maxtime=24*15,
                                    exactness=13,
                                    exact=F),
                           times=c(24*2,24*15),
                           start=24*2)
ibs_full_deter_fg <- ibs_full_deter_fg[2,2]

# Prediction for Null

pred_rf_full_det_null <- predict.rfsrc(arfsurv.obj_deter_full_null,data_ts_for_rf_full2_null)

# Generic formula for RSF Null
form_rsf_full_det_null <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                         paste(names(data_ts_for_rf_full2_null)[!grepl("Deterioro|los",
                                                                                       names(data_ts_for_rf_full2_null))
                                         ],collapse = " + ")))

# IBS for LR-RSF-Null

ibs_full_alt_det_null <-  crps(pec::pec(as.matrix(pred_rf_full_det_null$cif[,,2]),
                                        formula=form_rsf_full_det_null,
                                        data=data_ts_for_rf_full2_null,
                                        cause=2,
                                        start=24*2,
                                        maxtime=24*15,
                                        exactness=13,
                                        exact=F),
                               times=c(24*2,24*15),
                               start=24*2)
ibs_full_alt_det_null <- ibs_full_alt_det_null[2,2]

# IBS for GT-RSF-Null

ibs_null_det_fg <-  crps(pec::pec(as.matrix(pred_rf_full_null_fg$cif[,,2]),
                                  formula=form_rsf_null_alt_fg,
                                  data=data_ts_for_rf_full2_null,
                                  cause=2,
                                  start=24*2,
                                  maxtime=24*15,
                                  exactness=13,
                                  exact=F),
                         times=c(24*2,24*15),
                         start=24*2)
ibs_null_det_fg <- ibs_null_det_fg[2,2]

# Brier Score and C-Index 

# Vectors to allocate info

times <- 24*c(2:15)
c_ind_deter <- rep(0,14)
bs.rsf_full_det <- rep(0,14)
c_ind_deter_fg <- rep(0,14)
bs.rsf_full_det_fg <- rep(0,14)
c_ind_deter_null <- rep(0,14)
bs.rsf_full_det_null <- rep(0,14)
bs.rsf_full_alt_null_fg <- rep(0,14)

for(i in 1:14){
  # BS LR-RSF
  bs.rsf_full_det[i] <- brierCR(time_samp=data_ts_for_rf_full2$los,
                                t=times[i],
                                event_samp=data_ts_for_rf_full2$Deterioro,
                                event_int=2,   
                                cif=pred_rf_full_det$cif[,i,2])
  # BS GT-RSF
  bs.rsf_full_det_fg[i] <- brierCR(time_samp=data_ts_for_rf_full2$los,
                                   t=times[i],
                                   event_samp=data_ts_for_rf_full2$Deterioro,
                                   event_int=2,   
                                   cif=pred_rf_full_alt_fg$cif[,i,2])
  # BS LR-RSF Null
  bs.rsf_full_det_null[i] <- brierCR(time_samp=data_ts_for_rf_full2$los,
                                     t=times[i],
                                     event_samp=data_ts_for_rf_full2$Deterioro,
                                     event_int=2,   
                                     cif=pred_rf_full_det_null$cif[,i,2])
  # BS GT-RSF-Null
  bs.rsf_full_alt_null_fg[i] <- brierCR(time_samp=data_ts_for_rf_full2$los,
                                        t=times[i],
                                        event_samp=data_ts_for_rf_full2$Deterioro,
                                        event_int=2,   
                                        cif=pred_rf_full_null_fg$cif[,i,2])
  # C-Index LR-RSF
  c_ind_deter[i] <- as.double(pec::cindex(matrix(pred_rf_full_det$cif[,i,2]),
                                          Surv(los, Deterioro) ~ .,
                                          data_ts_for_rf_full2,
                                          cause=2,eval.times=24*(i+1))$AppCindex)
  # C-Index GT-RSF
  c_ind_deter_fg[i] <- as.double(pec::cindex(matrix(pred_rf_full_alt_fg$cif[,i,2]),
                                             Surv(los, Deterioro) ~ .,
                                             data_ts_for_rf_full2,
                                             cause=2,eval.times=24*(i+1))$AppCindex)
  print(i)
}

# CC-Index LR-RSF
sum(c_ind_deter)

# CC-Index LR-RSF
sum(c_ind_deter_fg)

# Table with BS and C-Index for deterioration
# with Full data-set

full_data_deter_cr_c_br <- data.frame(times=brier_csc_ts_deter$times,
                                      csc_brier=brier_csc_ts_deter$Brier,
                                      fgr_brier=brier_fgr_ts_deter$Brier,
                                      rf_brier=bs.rsf_full_det,
                                      csc_c_index=cindex_csc_ts_deter,
                                      fgr_c_index=cindex_fgr_ts_deter,
                                      rf_c_index=c_ind_deter,
                                      csc_brier_null=brier_csc_ts_deter_null$Brier,
                                      fgr_brier_null=brier_fgr_ts_deter_null$Brier,
                                      rf_c_index_fg=c_ind_deter_fg,
                                      rf_brier_fg=bs.rsf_full_det_fg,
                                      null_c_index=c_ind_deter_null,
                                      rsf_brier_null=bs.rsf_full_det_null,
                                      rsf_brier_null_fg=bs.rsf_full_alt_null_fg
)


######### 2. LASSO data-set #########

## 2.1 Event: Favourable discharge

## Linear models

compRiskModels_imp_mice_lasso_disc_ts <- compRiskModelsTest(data_ts_sum_imputed_by_mice %>% 
                                                              dplyr::select(as.vector(t(vars_comp_lasso))),
                                                            compRiskModels_imp_mice_lasso_disc,
                                                            cause = 1)
# IBS CSC
compRiskModels_imp_mice_lasso_disc_ts$ibs_csc
# IBS FG
compRiskModels_imp_mice_lasso_disc_ts$ibs_fg
# IBS Null
compRiskModels_imp_mice_lasso_disc_ts$ibs_null

# BS CSC
brier_lasso_csc_disc <- compRiskModels_imp_mice_lasso_disc_ts$IPA_csc %>% filter(Variable=="Full model")
# BS FG
brier_lasso_fgr_disc <- compRiskModels_imp_mice_lasso_disc_ts$IPA_fgr %>% filter(model=="FGR")

# C-Index CSC
cindex_lasso_csc_disc <- compRiskModels_imp_mice_lasso_disc_ts$c_index$CauseSpecificCox
# C-Index FG
cindex_lasso_fgr_disc <- compRiskModels_imp_mice_lasso_disc_ts$c_index$FGR

## Ensemble learning

# Prepare data-set

data_ts_for_rf_lasso <- data_ts_sum_imputed_by_mice %>% 
  dplyr::select(-`Episodio (único)`) %>% 
  dplyr::select(as.vector(t(vars_comp_lasso))) %>% 
  mutate_if(is.character,as.factor)
data_ts_for_rf_lasso$Deterioro <- as.factor(as.numeric(data_ts_for_rf_lasso$Deterioro)+1)
data_ts_for_rf_lasso2 <- data_ts_for_rf_lasso

# Prediction LR-RSF
pred_rf_lasso_alt <- predict.rfsrc(arfsurv.obj_lasso_disc,data_ts_for_rf_lasso2)

# General formula
form_rsf_lasso_alt <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                     paste(names(data_ts_for_rf_lasso2)[!grepl("Deterioro|los",
                                                                               names(data_ts_for_rf_lasso2))
                                     ],collapse = " + ")))

# IBS LR-RSF
ibs_full_alt_lasso <-  crps(pec::pec(as.matrix(pred_rf_lasso_alt$cif[,,1]),
                                     formula=form_rsf_lasso_alt,
                                     data=data_ts_for_rf_lasso2,
                                     cause=1,
                                     start=24*2,
                                     maxtime=24*15,
                                     exactness=13,
                                     exact=F),
                            times=c(24*2,24*15),
                            start=24*2)
ibs_full_alt_lasso <- ibs_full_alt_lasso[2,2]

# Prediction GT-RSF
pred_rf_lasso_alt_fg <- predict.rfsrc(arfsurv.obj_lasso_disc_fg,data_ts_for_rf_lasso2)

# IBS GT-RSF
ibs_full_alt_lasso_fg <-  crps(pec::pec(as.matrix(pred_rf_lasso_alt_fg$cif[,,1]),
                                        formula=form_rsf_lasso_alt,
                                        data=data_ts_for_rf_lasso2,
                                        cause=1,
                                        start=24*2,
                                        maxtime=24*15,
                                        exactness=13,
                                        exact=F),
                               times=c(24*2,24*15),
                               start=24*2)
ibs_full_alt_lasso_fg <- ibs_full_alt_lasso_fg[2,2]

# Brier Score and C-Index 

# Vectors to allocate info

bs.rsf_lasso_alt <- rep(0,14)
bs.rsf_lasso_alt_fg <- rep(0,14)
c_ind_disc_lasso_alt <- rep(0,14)
c_ind_disc_fg_lasso_alt <- rep(0,14)

for(i in 1:14){
  # BS LR-RSF
  bs.rsf_full_alt[i] <- brierCR(time_samp=data_ts_for_rf_lasso2$los,
                                t=times[i],
                                event_samp=data_ts_for_rf_lasso2$Deterioro,
                                event_int=1,   
                                cif=pred_rf_lasso_alt$cif[,i,1])
  # BS GT-RSF
  bs.rsf_lasso_alt_fg[i] <- brierCR(time_samp=data_ts_for_rf_lasso2$los,
                                    t=times[i],
                                    event_samp=data_ts_for_rf_lasso2$Deterioro,
                                    event_int=1,   
                                    cif=pred_rf_lasso_alt_fg$cif[,i,1])
  # C-Index LR-RSF
  c_ind_disc_lasso_alt[i] <- as.double(pec::cindex(matrix(pred_rf_lasso_alt$cif[,i,1]),
                                         Surv(los, Deterioro) ~ .,
                                         data_ts_for_rf_lasso2,
                                         cause=1,eval.times=24*(i+1))$AppCindex)
  # C-Index GT-RSF
  c_ind_disc_fg_lasso_alt[i] <- as.double(pec::cindex(matrix(pred_rf_lasso_alt_fg$cif[,i,1]),
                                            Surv(los, Deterioro) ~ .,
                                            data_ts_for_rf_lasso2,
                                            cause=1,eval.times=24*(i+1))$AppCindex)
  print(i)
}

# CC-Index LR-RSF
sum(c_ind_disc_lasso_alt)
# CC-Index GT-RSF
sum(c_ind_disc_fg_lasso_alt)

# Table with BS and C-Index for favourable discharge
# with LASSO data-set

lasso_data_disc_cr_c_br <- data.frame(times=brier_lasso_csc_disc$times,
                                      csc_brier=brier_lasso_csc_disc$Brier,
                                      fgr_brier=brier_lasso_fgr_disc$Brier,
                                      rf_brier=bs.rsf_full_alt,
                                      csc_c_index=cindex_lasso_csc_disc,
                                      fgr_c_index=cindex_lasso_fgr_disc,
                                      rf_c_index=c_ind_disc_lasso_alt,
                                      rf_c_index_fg=c_ind_disc_fg_lasso_alt,
                                      rf_brier_fg=bs.rsf_lasso_alt_fg
)

## 2.2 Event: Deterioration

## Linear models

compRiskModels_imp_mice_lasso_deter_ts <- compRiskModelsTest(data_ts_sum_imputed_by_mice %>% 
                                                               dplyr::select(as.vector(t(vars_comp_lasso))),
                                                             compRiskModels_imp_mice_lasso_deter,
                                                             cause = 2)
# IBS CSC
compRiskModels_imp_mice_lasso_deter_ts$ibs_csc
# IBS FG
compRiskModels_imp_mice_lasso_deter_ts$ibs_fg
# IBS Null
compRiskModels_imp_mice_lasso_deter_ts$ibs_null

# BS CSC
brier_lasso_csc_deter <- compRiskModels_imp_mice_lasso_deter_ts$IPA_csc %>% filter(Variable=="Full model")
# BS FG
brier_lasso_fgr_deter <- compRiskModels_imp_mice_lasso_deter_ts$IPA_fgr %>% filter(model=="FGR")

# C-Index CSC
cindex_lasso_csc_deter <- compRiskModels_imp_mice_lasso_deter_ts$c_index$CauseSpecificCox
# C-Index FG
cindex_lasso_fgr_deter <- compRiskModels_imp_mice_lasso_deter_ts$c_index$FGR

# Time-dependent AUROC CSC-LASSO

auroc_t_CSC_LASSO_deter <- rep(0,14)
for(i in 1:14){
  auroc_temp <- timeROC(T=time_cr,
                        delta=cause_cr,
                        weighting="marginal",
                        marker=compRiskModels_imp_mice_lasso_deter_ts$pred_csc[,i],
                        cause=2,
                        times=c(i+1)*24)
  auroc_t_CSC_LASSO_deter[i] <- as.numeric(auroc_temp$AUC_1[2])
}

## Ensemble learning

# Prediction LR-RSF
pred_rf_lasso_det <- predict.rfsrc(arfsurv.obj_lasso_deter,data_ts_for_rf_lasso)

# General formula
form_rsf_full_lasso_det <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                          paste(names(data_ts_for_rf_lasso)[!grepl("Deterioro|los",
                                                                                   names(data_ts_for_rf_lasso))
                                          ],collapse = " + ")))

# IBS LR-RSF
ibs_full_lasso_det<-  crps(pec::pec(as.matrix(pred_rf_lasso_det$cif[,,2]),
                                    formula=form_rsf_full_lasso_det,
                                    data=data_ts_for_rf_lasso,
                                    cause=2,
                                    start=24*2,
                                    maxtime=24*15,
                                    exactness=13,
                                    exact=F),
                           times=c(24*2,24*15),
                           start=24*2)
ibs_full_lasso_det <- ibs_full_lasso_det[2,2]

# IBS GT-RSF
ibs_full_alt_lasso_deter_fg <-  crps(pec::pec(as.matrix(pred_rf_lasso_alt_fg$cif[,,2]),
                                              formula=form_rsf_full_lasso_det,
                                              data=data_ts_for_rf_lasso2,
                                              cause=2,
                                              start=24*2,
                                              maxtime=24*15,
                                              exactness=13,
                                              exact=F),
                                     times=c(24*2,24*15),
                                     start=24*2)
ibs_full_alt_lasso_deter_fg <- ibs_full_alt_lasso_deter_fg[2,2]

# Brier Score and C-Index 

# Vectors to allocate info

bs.rsf_lasso_det <- rep(0,14)
bs.rsf_lasso_det_fg <- rep(0,14)
c_ind_deter_lasso_det <- rep(0,14)
c_ind_deter_fg_lasso_det <- rep(0,14)

for(i in 1:14){
  # BS LR_RSF
  bs.rsf_lasso_det[i] <- brierCR(time_samp=data_ts_for_rf_lasso2$los,
                                 t=times[i],
                                 event_samp=data_ts_for_rf_lasso2$Deterioro,
                                 event_int=2,   
                                 cif=pred_rf_lasso_det$cif[,i,2])
  # BS GT-RSF
  bs.rsf_lasso_det_fg[i] <- brierCR(time_samp=data_ts_for_rf_lasso2$los,
                                    t=times[i],
                                    event_samp=data_ts_for_rf_lasso2$Deterioro,
                                    event_int=2,   
                                    cif=pred_rf_lasso_alt_fg$cif[,i,2])
  # C-Index LR-RSF
  c_ind_deter_lasso_det[i] <- as.double(pec::cindex(matrix(pred_rf_lasso_det$cif[,i,2]),
                                          Surv(los, Deterioro) ~ .,
                                          data_ts_for_rf_lasso2,
                                          cause=2,eval.times=24*(i+1))$AppCindex)
  # C-Index GT-RSF
  c_ind_deter_fg_lasso_det[i] <- as.double(pec::cindex(matrix(pred_rf_lasso_alt_fg$cif[,i,2]),
                                             Surv(los, Deterioro) ~ .,
                                             data_ts_for_rf_lasso2,
                                             cause=2,eval.times=24*(i+1))$AppCindex)
  print(i)
}

# CC-Index LR-RSF
sum(c_ind_deter_lasso_det)

# CC-Index LR-RSF
sum(c_ind_deter_fg_lasso_det)

# Table with BS and C-Index for deterioration
# with LASSO data-set

lasso_data_deter_cr_c_br <- data.frame(times=brier_lasso_csc_deter$times,
                                       csc_brier=brier_lasso_csc_deter$Brier,
                                       fgr_brier=brier_lasso_fgr_deter$Brier,
                                       rf_brier=bs.rsf_lasso_det,
                                       csc_c_index=cindex_lasso_csc_deter,
                                       fgr_c_index=cindex_lasso_fgr_deter,
                                       rf_c_index=c_ind_deter_lasso_det,
                                       rf_c_index_fg=c_ind_deter_fg_lasso_det,
                                       rf_brier_fg=bs.rsf_lasso_det_fg
)

######### 3. BeSS data-set #########

## 3.1 Event: Favourable discharge

## Linear models

compRiskModels_imp_mice_bess_disc_ts <- compRiskModelsTest(data_ts_sum_imputed_by_mice %>% 
                                                             dplyr::select(as.vector(t(vars_comp_bess))),
                                                           compRiskModels_imp_mice_bess_disc,
                                                           cause = 1)

# IBS CSC
compRiskModels_imp_mice_bess_disc_ts$ibs_csc
# IBS FG
compRiskModels_imp_mice_bess_disc_ts$ibs_fg
# IBS Null
compRiskModels_imp_mice_bess_disc_ts$ibs_null

# BS CSC
brier_bess_csc_disc <- compRiskModels_imp_mice_bess_disc_ts$IPA_csc %>% filter(Variable=="Full model")
# BS FG
brier_bess_fgr_disc <- compRiskModels_imp_mice_bess_disc_ts$IPA_fgr %>% filter(model=="FGR")

# C-Index CSC
cindex_bess_csc_disc <- compRiskModels_imp_mice_bess_disc_ts$c_index$CauseSpecificCox
# C-Index FG
cindex_bess_fgr_disc <- compRiskModels_imp_mice_bess_disc_ts$c_index$FGR

# Time-dependent AUROC CSC-BeSS

auroc_t_CSC_BeSS_disc <- rep(0,14)
for(i in 1:14){
  auroc_temp <- timeROC(T=time_cr,
                        delta=cause_cr,
                        weighting="marginal",
                        marker=compRiskModels_imp_mice_bess_disc_ts$pred_csc[,i],
                        cause=1,
                        times=c(i+1)*24)
  auroc_t_CSC_BeSS_disc[i] <- as.numeric(auroc_temp$AUC_1[2])
}

## Ensemble learning

# Prepare data-set

data_ts_for_rf_bess <- data_ts_sum_imputed_by_mice %>% 
  dplyr::select(-`Episodio (único)`) %>% 
  dplyr::select(as.vector(t(vars_comp_bess))) %>% 
  mutate_if(is.character,as.factor)
data_ts_for_rf_bess$Deterioro <- as.factor(as.numeric(data_ts_for_rf_bess$Deterioro)+1)
data_ts_for_rf_bess2 <- data_ts_for_rf_bess

# Prediction LR-RSF
pred_rf_bess_alt <- predict.rfsrc(arfsurv.obj_bess2_disc,data_ts_for_rf_bess2)

# General formula
form_rsf_bess_alt <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                    paste(names(data_ts_for_rf_bess2)[!grepl("Deterioro|los",
                                                                             names(data_ts_for_rf_bess2))
                                    ],collapse = " + ")))

# IBS LR-RSF
ibs_full_bess_alt <-  crps(pec::pec(as.matrix(pred_rf_bess_alt$cif[,,1]),
                                    formula=form_rsf_bess_alt,
                                    data=data_ts_for_rf_bess2,
                                    cause=1,
                                    start=24*2,
                                    maxtime=24*15,
                                    exactness=13,
                                    exact=F),
                           times=c(24*2,24*15),
                           start=24*2)
ibs_full_bess_alt <- ibs_full_bess_alt[2,2]

# Prediction GT-RSF
pred_rf_bess_alt_fg <- predict.rfsrc(arfsurv.obj_bess2_disc_fg,data_ts_for_rf_bess2)

# IBS GT-RSF
ibs_full_bess_alt_fg <-  crps(pec::pec(as.matrix(pred_rf_bess_alt_fg$cif[,,1]),
                                       formula=form_rsf_bess_alt,
                                       data=data_ts_for_rf_bess2,
                                       cause=1,
                                       start=24*2,
                                       maxtime=24*15,
                                       exactness=13,
                                       exact=F),
                              times=c(24*2,24*15),
                              start=24*2)
ibs_full_bess_alt_fg <- ibs_full_bess_alt_fg[2,2]

# Brier Score and C-Index 

# Vectors to allocate info

bs.rsf_bess_alt <- rep(0,14)
bs.rsf_bess_alt_fg <- rep(0,14)
c_ind_disc_alt_bess <- rep(0,14)
c_ind_disc_fg_alt_bess <- rep(0,14)

for(i in 1:14){
  # BS LR-RSF
  bs.rsf_bess_alt[i] <- brierCR(time_samp=data_ts_for_rf_bess2$los,
                                t=times[i],
                                event_samp=data_ts_for_rf_bess2$Deterioro,
                                event_int=1,   
                                cif=pred_rf_bess_alt$cif[,i,1])
  # BS GT-RSF
  bs.rsf_bess_alt_fg[i] <- brierCR(time_samp=data_ts_for_rf_bess2$los,
                                   t=times[i],
                                   event_samp=data_ts_for_rf_bess2$Deterioro,
                                   event_int=1,   
                                   cif=pred_rf_bess_alt_fg$cif[,i,1])
  # C-Index LR
  c_ind_disc_alt_bess[i] <- as.double(pec::cindex(matrix(pred_rf_bess_alt$cif[,i,1]),
                                         Surv(los, Deterioro) ~ .,
                                         data_ts_for_rf_bess2,
                                         cause=1,eval.times=24*(i+1))$AppCindex)
  # C-Index GT-RSF
  c_ind_disc_fg_alt_bess[i] <- as.double(pec::cindex(matrix(pred_rf_bess_alt_fg$cif[,i,1]),
                                            Surv(los, Deterioro) ~ .,
                                            data_ts_for_rf_bess2,
                                            cause=1,eval.times=24*(i+1))$AppCindex)
  
  print(i)
}

# CC-Index LR-RSF
sum(c_ind_disc_alt_bess)

# CC-Index LR-RSF
sum(c_ind_disc_fg_alt_bess)

# Table with BS and C-Index for deterioration
# with BeSS data-set

bess_data_disc_cr_c_br <- data.frame(times=brier_bess_csc_disc$times,
                                     csc_brier=brier_bess_csc_disc$Brier,
                                     fgr_brier=brier_bess_fgr_disc$Brier,
                                     rf_brier=bs.rsf_bess_alt,
                                     csc_c_index=cindex_bess_csc_disc,
                                     fgr_c_index=cindex_bess_fgr_disc,
                                     rf_c_index=c_ind_disc_alt_bess,
                                     rf_c_index_fg=c_ind_disc_fg_alt_bess,
                                     rf_brier_fg=bs.rsf_bess_alt_fg
)

## 3.2 Event: Deterioration

## Linear models

compRiskModels_imp_mice_bess_deter_ts <- compRiskModelsTest(data_ts_sum_imputed_by_mice %>% 
                                                              dplyr::select(as.vector(t(vars_comp_bess))),
                                                            compRiskModels_imp_mice_bess_deter,
                                                            cause = 2)

# IBS CSC
compRiskModels_imp_mice_bess_deter_ts$ibs_csc
# IBS FG
compRiskModels_imp_mice_bess_deter_ts$ibs_fg
# IBS Null
compRiskModels_imp_mice_bess_deter_ts$ibs_null

# BS CSC
brier_bess_csc_deter <- compRiskModels_imp_mice_bess_deter_ts$IPA_csc %>% filter(Variable=="Full model")
# BS FG
brier_bess_fgr_deter <- compRiskModels_imp_mice_bess_deter_ts$IPA_fgr %>% filter(model=="FGR")

# C-Index CSC
cindex_bess_csc_deter <- compRiskModels_imp_mice_bess_deter_ts$c_index$CauseSpecificCox
# C-Index FG
cindex_bess_fgr_deter <- compRiskModels_imp_mice_bess_deter_ts$c_index$FGR

# Time-dependent AUROC FG-BeSS

auroc_t_FG_BeSS_deter <- rep(0,14)
for(i in 1:14){
  auroc_temp <- timeROC(T=time_cr,
                        delta=cause_cr,
                        weighting="marginal",
                        marker=compRiskModels_imp_mice_bess_deter_ts$pred_fgr[,i],
                        cause=2,
                        times=c(i+1)*24)
  auroc_t_FG_BeSS_deter[i] <- as.numeric(auroc_temp$AUC_2[2])
}

## Ensemble learning

# Prepare data-set
data_ts_for_rf_bess$Deterioro <- as.numeric(data_ts_for_rf_bess$Deterioro)

# Prediction LR-RSF
pred_rf_bess_det <- predict.rfsrc(arfsurv.obj_bess_deter,data_ts_for_rf_bess)

# General formula
form_rsf_bess_det <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                    paste(names(data_ts_for_rf_bess)[!grepl("Deterioro|los",
                                                                            names(data_ts_for_rf_bess))
                                    ],collapse = " + ")))

# IBS LR-RSF
ibs_full_bess_det <-  crps(pec::pec(as.matrix(pred_rf_bess_alt_fg$cif[,,2]),
                                    formula=form_rsf_bess_det,
                                    data=data_ts_for_rf_bess,
                                    cause=2,
                                    start=24*2,
                                    maxtime=24*15,
                                    exactness=13,
                                    exact=F),
                           times=c(24*2,24*15),
                           start=24*2)
ibs_full_bess_det <- ibs_full_bess_det[2,2]

# IBS GT-RSF

ibs_full_bess_deter_fg <-  crps(pec::pec(as.matrix(pred_rf_bess_alt_fg$cif[,,2]),
                                         formula=form_rsf_bess_det,
                                         data=data_ts_for_rf_bess2,
                                         cause=2,
                                         start=24*2,
                                         maxtime=24*15,
                                         exactness=13,
                                         exact=F),
                                times=c(24*2,24*15),
                                start=24*2)
ibs_full_bess_deter_fg <- ibs_full_bess_deter_fg[2,2]

# Brier Score and C-Index 

# Vectors to allocate info

bs.rsf_bess_det <- rep(0,14)
bs.rsf_bess_det_fg <- rep(0,14)
c_ind_deter_det_bess <- rep(0,14)
c_ind_deter_fg_det_bess  <- rep(0,14)

for(i in 1:14){
  # BS LR-RSF
  bs.rsf_bess_det[i] <- brierCR(time_samp=data_ts_for_rf_bess$los,
                                t=times[i],
                                event_samp=data_ts_for_rf_bess$Deterioro,
                                event_int=2,   
                                cif=pred_rf_bess_det$cif[,i,2])
  # BS GT-RSF
  bs.rsf_bess_det_fg[i] <- brierCR(time_samp=data_ts_for_rf_bess2$los,
                                   t=times[i],
                                   event_samp=data_ts_for_rf_bess2$Deterioro,
                                   event_int=2,   
                                   cif=pred_rf_bess_alt_fg$cif[,i,2])
  # C-Index LR-RSF
  c_ind_deter_det_bess[i] <- as.double(pec::cindex(matrix(pred_rf_bess_det$cif[,i,2]),
                                          Surv(los, Deterioro) ~ .,
                                          data_ts_for_rf_bess2,
                                          cause=2,eval.times=24*(i+1))$AppCindex)
  # C-index GT-RSF
  c_ind_deter_fg_det_bess[i] <- as.double(pec::cindex(matrix(pred_rf_bess_alt_fg$cif[,i,2]),
                                             Surv(los, Deterioro) ~ .,
                                             data_ts_for_rf_bess2,
                                             cause=2,eval.times=24*(i+1))$AppCindex)
  print(i)
}

# CC-Index LR-RSF
sum(c_ind_deter_det_bess)

# CC-Index LR-RSF
sum(c_ind_deter_fg_det_bess)

# Table with BS and C-Index for deterioration
# with BeSS data-set

bess_data_deter_cr_c_br <- data.frame(times=brier_bess_csc_deter$times,
                                      csc_brier=brier_bess_csc_deter$Brier,
                                      fgr_brier=brier_bess_fgr_deter$Brier,
                                      rf_brier=bs.rsf_bess_det,
                                      csc_c_index=cindex_bess_csc_deter,
                                      fgr_c_index=cindex_bess_fgr_deter,
                                      rf_c_index=c_ind_deter_det_bess,
                                      rf_c_index_fg=c_ind_deter_fg_det_bess,
                                      rf_brier_fg=bs.rsf_bess_det_fg
)
