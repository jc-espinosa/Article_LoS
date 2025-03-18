# LoS_CR_Modelling
Scripts developed to obtain the main results of article "Modelling in-hospital length of stay: A comparison of linear and ensemble models for competing risk analysis".

# Structure

1. Folder "1_Pre_processing". This folder contains two files: <br />

   **1.1. Imputation.R**: Implements imputation methods, including mean, median, MICE, BPCA, and NIPALS. This script requires the following covariates: ID, LoS, endpoint, sex, age, and transversal statistics of vital signs. These covariates should be stored in the CSV file train_outer_dataset.csv within the "Datasets" folder.  <br />
   **1.2. Variable_Selection.R**: Performs variable selection using BeSS and LASSO-regularized Cox regression. This script requires prior execution of "Imputation.R".  <br />
2. Folder "2_Modelling". This folder contains two files:  <br />

   **2.1. Competing_Risk_Training.R**: Trains competing risk models, including the Cause-Specific Cox model, Subdistribution Hazard model, and Random Survival Forest with Generalized Log-Rank and Gray’s test as splitting rules. This script requires prior execution of both "Imputation.R" and "Variable_Selection.R", and uses the datasets train_outer_dataset.csv and test_dataset.csv from the "Datasets" folder.  <br />
   **2.2. Competing_Risk_Testing.R**: Evaluates the models trained in "Competing_Risk_Training.R" using metrics such as Brier Score, Integrated Brier Score, C-Index, and a custom version of the Cumulative C-Index.  <br />

3. Folder "3_Baseline". This folder contains two files:  <br />

   **3.1. MEWS.R**: Computes the Modified Early Warning Score (MEWS) and includes training and testing of competing risk models. It requires only the dataset containing the last recorded observations of vital signs, along with ID, LoS, and endpoint. This dataset is stored as dataset_last_observations.csv in the "Datasets" folder. Unlike "Imputation.R", this script replaces missing values with 0 instead of using MICE. <br />
   **3.2. NEWS.R**: Computes the National Early Warning Score (NEWS) and includes training and testing of competing risk models. It requires the same dataset described in 3.1.  <br />

Note: Null models are estimated in the scripts within the "2_Modelling" folder, not in "3_Baseline".
