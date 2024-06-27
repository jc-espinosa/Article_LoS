# LoS_CR_Modelling
Scripts developed to obtain the main results of article "Modelling in-hospital length of stay: A comparison of linear and ensemble models for competing risk analysis".

# Structure

1. Folder "1_Pre_processing" contains two files: <br />

   **1.1. Imputation.R**: Imputation methods mean, median, MICE, BPCA and NIPALS. It needs the covariates ID, LoS, endpoint, sex, age, and transversal statistics of vital signs, allocated in the csv "train_outer_dataset.csv" in the "Datasets" folder.  <br />
   **1.2. Variable_Selection.R**: Variable selection with BeSS and LASSO relugarized Cox regression. It needs previous compilation of "Imputation.R" script.  <br />
2. Folder "2_Modelling" contains two files:  <br />

   **2.1. Competing_Risk_Training.R**: Cause-Specific Cox, Subdistribution Hazard and Random Survival Forest with Generalized Log-Rank test and Gray's test as splitting rules, are trained in this file. It needs previous compilation of both "Imputation.R" and "Variable_Selection.R" script, and uses "train_outer_dataset.csv" and "test_dataset.csv", in the "Datasets" folder.  <br />
   **2.2. Competing_Risk_Testing.R**: Testing of models trained in "Competing_Risk_Training.R", evaluated via Brier Score, Integrated Brier Score, C-Index and a custom version of Cumulative C-Index .  <br />

3. Folder "3_Baseline" contains two files:  <br />

   **3.1. MEWS.R**: Calculation of Modified Early Warning Score. It includes training and testing of CR models. It only needs the dataset of last observations of vital signs, as well as ID, LoS and endpoint, named "dataset_last_observations.csv" in folder "Datasets". Instead of MICE, use 0 for missing data. <br />
   **3.2. NEWS.R**: Calculation of National Early Warning Score. It includes training and testing of CR models. It only needs the dataset of last observations, described in 3.1.  <br />

Note: Null models are estimated in scipts from folder "2_Modelling", not in "3_Baseline".
