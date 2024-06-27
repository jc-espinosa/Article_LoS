
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

source("Competing_Risk_Training.R")
source("Competing_Risk_Testing.R")

###################################
######## Data preparation  ########
###################################

# Read training database with last observations for NEWS

dat_resume_sat <- read.csv2(paste0(root_input,"/dataset_last_observations.csv"))

# Calculate NEWS

news_score <- dat_resume_sat %>% mutate(resp_news = cut(last_freq_resp,c(-Inf,8,11,20,24,Inf),labels = c(3,1,0,2,3),right = TRUE),
                                        o2sat_news = cut(last_satur,c(-Inf,91,93,95,Inf),labels = c(3,2,1,0),right = TRUE),
                                        supp_o2_news = ifelse(last_satur_aux==0,0,ifelse(last_satur_aux>0,2,NA)),
                                        temp_news = cut(last_temp,c(-Inf,35,36,38,39,Inf),labels = c(3,1,0,1,2),right = TRUE), # celsius
                                        sbp_news = cut(last_pres_sisto,c(-Inf,90,100,110,219,Inf),labels = c(3,2,1,0,3),right = TRUE),
                                        hr_news = cut(last_freq_card,c(-Inf,40,50,90,110,130,Inf),labels = c(3,1,0,1,2,3),right = TRUE),
                                        avpu_news = ifelse(last_est_neur %in% c('1','nm'),0,3)
)
news_score <- news_score %>% select(resp_news,o2sat_news,supp_o2_news,temp_news,
                                    sbp_news,hr_news,avpu_news)
news_score <- apply(news_score,2,as.numeric)
news_score[is.na(news_score)] <- 0
news_score <- apply(news_score,1,sum,na.rm=T)

# Definite data-set

dat_resume_sat <- dat_resume_sat %>% 
  dplyr::select(`Episodio (único)`,los,Deterioro,Sexo,Edad)
dat_resume_sat$NEWS <- news_score

###########################
######## Training  ########
###########################

### 1. Event: Favourable discharge

## Linear models

compRiskModels_disc_ews <- compRiskModelsTr(dat_resume_sat,
                                            cause=1)

## Ensemble learning

# Prepare data-set
data_for_rf_ews <- dat_resume_sat %>% 
  dplyr::select(-`Episodio (único)`) %>% 
  mutate_if(is.character,as.factor)
data_for_rf_ews$Deterioro <- as.factor(as.numeric(data_for_rf_ews$Deterioro) + 1)
data_for_rf_ews <- as.data.frame(data_for_rf_ews)

# Parallel options
rf.cores_rsf = 8
mc.cores_rsf = 8
options(rf.cores=rf.cores_rsf, mc.cores=mc.cores_rsf)

# Tuning of hyper-parameters of LR-RSF
rf_tune_ews_disc <- tune(Surv(los, Deterioro) ~ ., cause = c(1,0),
                         splitrule = "logrank", 
                         data = data_for_rf_ews, trace = TRUE,
                         ntreeTry = 100, seed = 1234)
# Optimal
# nodesize     mtry 
#    80          1

# LR-RSF based on hyper-parameters tuning
arfsurv.obj_disch_ews <- rfsrc(Surv(los, Deterioro) ~ ., 
                               data = data_for_rf_ews,
                               cause = c(1,0),
                               splitrule = "logrank",
                               mtry = 1,
                               nodesize = 80,
                               importance = TRUE,
                               ntree = 100,
                               seed = 1234,
                               ntime=24*(2:15))  

# Tuning of hyper-parameters of GT-RSF
rf_tune_fg_ews <- tune(Surv(los, Deterioro) ~ ., cause = c(1,1),
                       splitrule = "logrankCR", 
                       data = data_for_rf_ews, trace = TRUE,
                       ntreeTry = 100,seed=1234)

# Optimal
#nodesize     mtry 
#  2           1 

# GT-RSF based on hyper-parameters tuning
arfsurv.obj_disch_fg_ews <- rfsrc(Surv(los, Deterioro) ~ ., 
                                  data = data_for_rf_ews,
                                  cause = c(1,1),
                                  splitrule = "logrankCR",
                                  mtry = 1,
                                  nodesize = 2,
                                  importance = TRUE,
                                  ntree = 100,
                                  seed = 1234,
                                  ntime=24*(2:15))  

### 2. Event: Deterioration

## Linear models
compRiskModels_det_ews <- compRiskModelsTr(dat_resume_sat,
                                           cause=2)

# Tuning of hyper-parameters of LR-RSF
rf_tune_ews_det <- tune(Surv(los, Deterioro) ~ ., cause = c(0,1),
                        splitrule = "logrank", 
                        data = data_for_rf_ews, trace = TRUE,
                        ntreeTry = 100, seed = 1234)

# Optimal
# nodesize     mtry 
#   25           1

# LR-RSF based on hyper-parameters tuning
arfsurv.obj_det_ews <- rfsrc(Surv(los, Deterioro) ~ ., 
                             data = data_for_rf_ews,
                             cause = c(0,1),
                             splitrule = "logrank",
                             mtry = 1,
                             nodesize = 25,
                             importance = TRUE,
                             ntree = 100,
                             seed = 1234,
                             ntime=24*(2:15))  

##########################
######## Testing  ########
##########################

# Read testing database with last observations for NEWS

dat_resume_sat_ts <- readRDS(paste0(root_input,"/dat_resume_sat_ts.rds"))

# Calculate NEWS

news_score_ts <- dat_resume_sat_ts %>% mutate(resp_news = cut(last_freq_resp,c(-Inf,8,11,20,24,Inf),labels = c(3,1,0,2,3),right = TRUE),
                                              o2sat_news = cut(last_satur,c(-Inf,91,93,95,Inf),labels = c(3,2,1,0),right = TRUE),
                                              supp_o2_news = ifelse(last_satur_aux==0,0,2),
                                              temp_news = cut(last_temp,c(-Inf,35,36,38,39,Inf),labels = c(3,1,0,1,2),right = TRUE), # celsius
                                              sbp_news = cut(last_pres_sisto,c(-Inf,90,100,110,219,Inf),labels = c(3,2,1,0,3),right = TRUE),
                                              hr_news = cut(last_freq_card,c(-Inf,40,50,90,110,130,Inf),labels = c(3,1,0,1,2,3),right = TRUE),
                                              avpu_news = ifelse(last_est_neur %in% c('1','nm'),0,3)
)
news_score_ts <- news_score_ts %>% select(resp_news,o2sat_news,supp_o2_news,temp_news,
                                          sbp_news,hr_news,avpu_news)
news_score_ts <- apply(news_score_ts,2,as.numeric)
news_score_ts[is.na(news_score_ts)] <- 0
news_score_ts <- apply(news_score_ts,1,sum,na.rm=T)

# Definite data-set

dat_resume_sat_ts <- dat_resume_sat_ts %>% 
  dplyr::select(`Episodio (único)`,los,Deterioro,Sexo,Edad)
dat_resume_sat_ts <- as.data.frame(dat_resume_sat_ts)
dat_resume_sat_ts$NEWS <- news_score_ts

### 1. Event: Favourable discharge

## Linear models

compRiskModels_mice_test_disc_ews <- compRiskModelsTest(dat_resume_sat_ts,
                                                        compRiskModels_disc_ews,
                                                        state = 1)
# IBS CSC
compRiskModels_mice_test_disc_ews$ibs_csc
# IBS FG
compRiskModels_mice_test_disc_ews$ibs_fg

# BS CSC
brier_csc_ews_ts_disc <- compRiskModels_mice_test_disc_ews$IPA_csc %>% filter(Variable=="Full model")
# BS FG
brier_fg_ews_ts_disc <- compRiskModels_mice_test_disc_ews$IPA_fgr %>% filter(model=="FGR")

# C-Index CSC
cindex_csc_ews_ts_disc <- compRiskModels_mice_test_disc_ews$c_index$CauseSpecificCox
# C-Index FG
cindex_fg_ews_ts_disc <- compRiskModels_mice_test_disc_ews$c_index$FGR

# CC-Index CSC
sum(cindex_csc_ews_ts_disc)
# CC-Index BS
sum(cindex_fg_ews_ts_disc)

## Ensemble learning

# Prepare data-set

data_ts_for_rf_ews <- dat_resume_sat_ts %>% 
  dplyr::select(-`Episodio (único)`) %>% 
  mutate_if(is.character,as.factor)
data_ts_for_rf_ews$Deterioro <- as.factor(as.numeric(data_ts_for_rf_ews$Deterioro)+1)

# Prediction LR-RSF
pred_rf_ews_alt <- predict.rfsrc(arfsurv.obj_disch_ews,data_ts_for_rf_ews)

# Prediction GT-RSF
pred_rf_ews_alt_fg <- predict.rfsrc(arfsurv.obj_disch_fg_ews,data_ts_for_rf_ews)

# General formula
form_rsf_ews_disc <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                    paste(names(data_ts_for_rf_ews)[!grepl("Deterioro|los",
                                                                           names(data_ts_for_rf_ews))
                                    ],collapse = " + ")))

# IBS LR-RSF
ibs_ews_disc <-  crps(pec::pec(as.matrix(pred_rf_ews_alt$cif[,,1]),
                               formula=form_rsf_ews_disc,
                               data=data_ts_for_rf_ews,
                               cause=1,
                               start=24*2,
                               maxtime=24*15,
                               exactness=13,
                               exact=F),
                      times=c(24*2,24*15),
                      start=24*2)
ibs_ews_disc <- ibs_ews_disc[2,2]

# IBS GT-RSF
ibs_ews_disc_fg <-  crps(pec::pec(as.matrix(pred_rf_ews_alt_fg$cif[,,1]),
                                  formula=form_rsf_ews_disc,
                                  data=data_ts_for_rf_ews,
                                  cause=1,
                                  start=24*2,
                                  maxtime=24*15,
                                  exactness=13,
                                  exact=F),
                         times=c(24*2,24*15),
                         start=24*2)
ibs_ews_disc_fg <- ibs_ews_disc_fg[2,2]

# Brier Score and C-Index 
# Vectors to allocate info

bs_disc_ews  <- rep(0,14)
bs_disc_ews_fg  <- rep(0,14)
c_ind_disc_ews <- rep(0,14)
c_ind_disc_ews_fg <- rep(0,14)
times_ews <- 24*c(2:15)

for(i in 1:14){
  # BS LR-RSF
  bs_disc_ews[i] <- brierCR(time_samp=data_ts_for_rf_ews$los,
                            t=times_ews[i],
                            event_samp=data_ts_for_rf_ews$Deterioro,
                            event_int=1,   
                            cif=pred_rf_ews_alt$cif[,i,1])
  # BS GT-RSF
  bs_disc_ews_fg[i] <- brierCR(time_samp=data_ts_for_rf_ews$los,
                               t=times_ews[i],
                               event_samp=data_ts_for_rf_ews$Deterioro,
                               event_int=1,   
                               cif=pred_rf_ews_alt_fg$cif[,i,1])
  # C-Index LR-RSF
  c_ind_disc_ews[i] <- as.double(pec::cindex(matrix(pred_rf_ews_alt$cif[,i,1]),
                                             Surv(los, Deterioro) ~ .,
                                             data_ts_for_rf_ews,
                                             cause=1,
                                             eval.times=24*(i+1))$AppCindex)
  # C-Index GT-RSF
  c_ind_disc_ews_fg[i] <- as.double(pec::cindex(matrix(pred_rf_ews_alt_fg$cif[,i,1]),
                                                Surv(los, Deterioro) ~ .,
                                                data_ts_for_rf_ews,
                                                cause=1,
                                                eval.times=24*(i+1))$AppCindex)
}

# BS LR-RSF
bs_disc_ews
# BS GT-RSF
bs_disc_ews_fg
# C-Index LR-RSF
c_ind_disc_ews
# C-Index GT-RSF
c_ind_disc_ews_fg
# CC-Index LR-RSF
sum(c_ind_disc_ews)
# CC-Index GT-RSF
sum(c_ind_disc_ews_fg)

### 2. Event: deterioration

## Linear models

compRiskModels_mice_test_det_ews <- compRiskModelsTest(dat_resume_sat_ts,
                                                       compRiskModels_det_ews,
                                                       state = 2)
# IBS CSC
compRiskModels_mice_test_det_ews$ibs_csc
# IBS FG
compRiskModels_mice_test_det_ews$ibs_fg

# BS CSC
brier_csc_ews_ts_det <- compRiskModels_mice_test_det_ews$IPA_csc %>% filter(Variable=="Full model")
# BS FG
brier_fg_ews_ts_det <- compRiskModels_mice_test_det_ews$IPA_fgr %>% filter(model=="FGR")

# C-Index CSC
cindex_csc_ews_ts_det <- compRiskModels_mice_test_det_ews$c_index$CauseSpecificCox
# C-Index FG
cindex_fg_ews_ts_det <- compRiskModels_mice_test_det_ews$c_index$FGR

# CC-Index CSC
sum(cindex_csc_ews_ts_det)
# CC-Index FG
sum(cindex_fg_ews_ts_det)

## Ensemble learning

# Prediction LR-RSF
pred_rf_ews_det <- predict.rfsrc(arfsurv.obj_det_ews,data_ts_for_rf_ews)

# Prediction GT-RSF
pred_rf_ews_det_fg <- predict.rfsrc(arfsurv.obj_disch_fg_ews,data_ts_for_rf_ews)

# General formula
form_rsf_ews_det <- formula(paste0("Surv(time=los, event=Deterioro) ~ ",
                                   paste(names(data_ts_for_rf_ews)[!grepl("Deterioro|los",
                                                                          names(data_ts_for_rf_ews))
                                   ],collapse = " + ")))
# IBS LR-RSF
ibs_ews_det <-  crps(pec::pec(as.matrix(pred_rf_ews_det$cif[,,2]),
                              formula=form_rsf_ews_det,
                              data=data_ts_for_rf_ews,
                              cause=2,
                              start=24*2,
                              maxtime=24*15,
                              exactness=13,
                              exact=F),
                     times=c(24*2,24*15),
                     start=24*2)
ibs_ews_det <- ibs_ews_det[2,2]

# IBS GT-RSF
ibs_ews_det_fg <-  crps(pec::pec(as.matrix(pred_rf_ews_det_fg$cif[,,2]),
                                 formula=form_rsf_ews_det,
                                 data=data_ts_for_rf_ews,
                                 cause=2,
                                 start=24*2,
                                 maxtime=24*15,
                                 exactness=13,
                                 exact=F),
                        times=c(24*2,24*15),
                        start=24*2)
ibs_ews_det_fg <- ibs_ews_det_fg[2,2]

# Brier Score and C-Index 
# Vectors to allocate info

bs_det_ews  <- rep(0,14)
bs_det_ews_fg  <- rep(0,14)
c_ind_det_ews <- rep(0,14)
c_ind_det_ews_fg <- rep(0,14)
times_ews <- 24*c(2:15)

for(i in 1:14){
  # BS LR-RSF
  bs_det_ews[i] <- brierCR(time_samp=data_ts_for_rf_ews$los,
                           t=times_ews[i],
                           event_samp=data_ts_for_rf_ews$Deterioro,
                           event_int=2,   
                           cif=pred_rf_ews_det$cif[,i,2])
  # BS GT-RSF
  bs_det_ews_fg[i] <- brierCR(time_samp=data_ts_for_rf_ews$los,
                              t=times_ews[i],
                              event_samp=data_ts_for_rf_ews$Deterioro,
                              event_int=2,   
                              cif=pred_rf_ews_det_fg$cif[,i,2])
  # C-Index LR-RSF
  c_ind_det_ews[i] <- as.double(pec::cindex(matrix(pred_rf_ews_det$cif[,i,2]),
                                            Surv(los, Deterioro) ~ .,
                                            data_ts_for_rf_ews,
                                            cause=2,
                                            eval.times=24*(i+1))$AppCindex)
  # C-Index GT-RSF
  c_ind_det_ews_fg[i] <- as.double(pec::cindex(matrix(pred_rf_ews_det_fg$cif[,i,2]),
                                               Surv(los, Deterioro) ~ .,
                                               data_ts_for_rf_ews,
                                               cause=2,
                                               eval.times=24*(i+1))$AppCindex)
}

# BS LR-RSF
bs_det_ews
# BS GT-RSF
bs_det_ews_fg
# C-Index LR-RSF
c_ind_det_ews
# C-Index GT-RSF
c_ind_det_ews_fg
# CC-Index LR-RSF
sum(c_ind_det_ews)
# CC-Index GT-RSF
sum(c_ind_det_ews_fg)

