# QRT-data-Challenge
My participation to the QRT challenge : The goal is to predict overall survival (OS) in patients diagnosed with acute myeloid leukemia (AML).

## Stucture

- My participation file (.csv) is located in the data folder.
- The notebook contains the code used to construct my predictive model.
- The scripts folder stores the reusable functions.

## Model

- I cleaned the data following a simple rule: if too much information was missing, I removed the entire column (e.g., MONOCYTES).
- If less than 10% of the values were missing, I applied imputation, after performing a statistical test to confirm that the imputed distribution matched the original one.
- I engineered features from both clinical and molecular data, supported by a review of the scientific literature to identify which variables are informative and which are not.
- I compared several survival models, including Cox proportional hazards, Random Survival Forests (RSF), and deep learning–based approaches. RSF performed best, as it captures non-linear dependencies, while deep learning suffered from limited data availability.
- I used a conservative methodology with K-fold cross-validation combined with Optuna for hyperparameter optimization to select the best model configuration.

## Result

IPCW-C-index : 0.7440 with only one submission, expected at least top 10% as the first one has 0.7208 in private leader bord, and i should not suffer from leader bord fitting as i submitted only once. 
