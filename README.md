## SurvPancreas

> [!NOTE]
> We provide a synthetic example data in the `data/` directory for reproducibility checks.
> Any analytical results based on it do not represent the real data.
> For the actual model performance and prognostic associations, please refer to the original research article.

All analyses were performed using R version 4.5.1.

### Required R Packages

- [`{tidyverse}`](https://tidyverse.org/packages/):
  A collection of packages for tidy data manipulation and functional programming.
- [`{survival}`](https://github.com/therneau/survival) (v3.8.3):
  Create censored time-to-event object and fit the Cox models.
- [`{grpreg}`](https://github.com/pbreheny/grpreg) (v3.5.0):
  Apply group Lasso to Cox models.
- [`{aorsf}`](https://github.com/ropensci/aorsf/) (v0.1.5):
  Fit the oblique random survival forests.
- [`{xgboost}`](https://xgboost.readthedocs.io/en/latest/R-package/) (v3.1.1):
  Fit XGBoost models.
- [`{survivalmodels}`](https://github.com/RaphaelS1/survivalmodels/) (v0.1.19):
  Fit deep-learning survival models.
- [`{riskRegression}`](https://github.com/tagteam/riskRegression) (v2025.9.17):
  Calculate time-dependent performance metrics.
- [`{mice}`](https://github.com/amices/mice) (v3.18.0):
  Perform multiple imputation by chained equations.
- [`{plotshap}`](https://github.com/darentsai/plotshap):
  Visualize SHAP values; available only on Github for now.

Install the key packages above, along with other auxiliary packages:

```r
pkgs <- c("tidyverse", "survival", "grpreg", "aorsf", "xgboost", "survivalmodels",
          "riskRegression", "furrr", "mice", "gbm", "furrr", "progressr", "patchwork", "cli")

install.packages(pkgs)
remotes::install_github("darentsai/plotshap")
```

### Source code

The source code for this project is available in the `R/` directory:

- `R/` directory
  - `model_fit.R`: Wrapper functions for each model that reparameterize the original models, tune their hyperparameters, and output fitted models.
  - `model_pred.R`: Wrapper functions for each model that predict survival probabilities for each new instance at specified times using the fitted models.
  - `bench.R`: Benchmark function that performs the repeated and nested cross-validation, and outputs multiple time-dependent performance metrics.
  - `utils.R`: Other auxiliary functions.

For analysis, load these scripts first to import all the functions.
The vignettes below demonstrate their usage.

### Vignettes

- `analysis1.qmd`: [SHAP Analysis and Visualization](https://darentsai.github.io/SurvPancreas/analysis1)
- `analysis2.qmd`: [Model Performance Evaluation](https://darentsai.github.io/SurvPancreas/analysis2)
