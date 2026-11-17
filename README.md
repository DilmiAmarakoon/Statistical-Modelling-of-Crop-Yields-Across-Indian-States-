# Statistical Modelling of Crop Yield Across Indian States

An R project for IT3081 that compares statistical models for predicting crop yield and identifying above-median yields across Indian states. The analysis covers rice, wheat, maize, sugarcane, cotton (lint), and jute.

## Dataset

- **Source:** *Crop Yield in Indian States* by Akshat Gupta, published on Kaggle.
- **Scope:** 3,973 records across six crops from 1997 to 2019.
- **Predictors:** Crop, state, season, year, log-transformed cultivated area, and log-transformed annual rainfall.
- **Regression target:** Log-relative yield, calculated using each crop's median yield in the training data as its reference.
- **Classification target:** Whether yield exceeds the training median for that crop.

## Methodology

1. Cleaned and audited the dataset, then filtered it to the six crops.
2. Used a chronological split, with 1997–2016 for training, 2017 for validation, 2018 for calibration, and 2019 for final testing.
3. Compared historical baselines with regression models including ordinary least squares, ridge, LASSO, and interaction ridge. Also evaluated logistic regression, ridge logistic regression, and naïve Bayes for classification.
4. Checked prediction interval coverage and compared performance across crops.

## Key Results

| Task | Best approach on the chosen metric | 2019 test result |
| --- | --- | --- |
| Regression | Interaction Ridge | Macro RMSE of **0.1661**, a **23.05%** reduction from the crop–state historical baseline (0.2158) |
| Classification | Crop–state probability baseline | Macro balanced accuracy of **80.03%** |

The 2019 test set contained 196 records. Regression errors are measured on the log-relative yield scale.

## Repository Files

| File | Description |
| --- | --- |
| `Six_Crops.R` | R script for preprocessing, modelling, evaluation, and plots |
| `crop_yield.csv` | Input dataset |
| [Statistical Modelling of Crop Yield Across Indian States.pdf](Statistical%20Modelling%20of%20Crop%20Yield%20Across%20Indian%20States.pdf) | Project report covering the methodology and findings |

## How to Run

1. Install **R 4.3.0 or newer**.
2. Download the repository and open R or RStudio with the repository folder as your working directory. Keep `Six_Crops.R` and `crop_yield.csv` together.
3. Run these commands in the R console.

```r
install.packages(c("dplyr", "ggplot2", "car", "glmnet", "MASS", "e1071", "rstatix", "lmtest"))
source("Six_Crops.R")
```

The script prints result tables in the console and displays plots during the analysis. 
