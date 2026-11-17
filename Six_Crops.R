# ============================================================
# IT3081 MASTER PROJECT: RELATIVE YIELD IN SIX CROPS
# Rice, Wheat, Maize, Sugarcane, Cotton(lint), Jute
# Run Modules 0-8 in order, or source the whole file in a fresh R session.
# Each module prints tables/plots immediately. No automatic exports.
# Requires R >= 4.3.0. Never describe unverified units as tonnes/ha.
# ============================================================

# ============================================================
## Module 0: Helper Functions & Setup
# ============================================================

# ------------------------------------------------------------
# 0.1 Input and display settings
# ------------------------------------------------------------

csv_path <- "crop_yield.csv"
make_plots <- TRUE
apply_provisional_fibre_conversion <- FALSE
# FALSE is the primary analysis. TRUE is an explicitly unverified scenario.
# Relative yields are invariant to constant within-crop conversion factors.

# ------------------------------------------------------------
# 0.2 Package and file checks
# ------------------------------------------------------------

required <- c("dplyr", "ggplot2", "car", "glmnet", "MASS", "e1071", "rstatix", "lmtest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: install.packages(c(", paste(shQuote(missing),
  collapse = ", "), "))")
if (getRversion() < "4.3.0") stop("R 4.3.0 or newer is required.")
library(dplyr)
library(ggplot2)
if (!file.exists(csv_path) && interactive()) csv_path <- file.choose()
if (!file.exists(csv_path)) stop("Set csv_path to the full path of crop_yield.csv.")

# ------------------------------------------------------------
# 0.3 Reproducibility and model settings
# ------------------------------------------------------------

set.seed(123)
options(width = 130)
tables <- list()
crops <- c("Rice", "Wheat", "Maize", "Sugarcane", "Cotton", "Jute")
lambda_grid <- 10^seq(2, -6, length.out = 25)
shrinkage_k <- 10
alpha_interval <- 0.05
numeric_features <- c("Year", "Log_Area", "Log_Rain")
main_formula <- Z ~ Crop + State + Season + Year + Log_Area + Log_Rain
interaction_formula <- ~ Crop * State + Season + Year + Log_Area + Log_Rain + Crop:Year + Crop:Log_Rain
main_xformula <- delete.response(terms(main_formula))
class_formula <- High ~ Crop + State + Season + Year + Log_Area + Log_Rain
project <- list(tables = tables)

# Utilities: reporting; preprocessing; model fitting; metrics; intervals.

# ------------------------------------------------------------
# 0.4 Reporting helpers
# ------------------------------------------------------------

print_table <- function(name, value, n = 24) {
  tables[[name]] <<- value
  project$tables[[name]] <<- value
  cat("\n=====", name, "=====\n")
  if (is.data.frame(value)) print(as.data.frame(utils::head(value, n)), row.names = FALSE) else print(value)
  if (is.data.frame(value) && nrow(value) > n) cat("Full table retained in project$tables$", name,
    "\n", sep = "")
  flush.console()
  invisible(value)
}

print_plot <- function(p) {
  if (make_plots) {
    print(p);
    flush.console()
  }
  invisible(p)
}

# ------------------------------------------------------------
# 0.5 Target construction and preprocessing
# ------------------------------------------------------------

add_targets <- function(d, ref = reference) {
  d <- left_join(d, ref, by = "Crop")
  if (anyNA(d$Crop_Median)) stop("Crop absent from the training reference.")
  d$Relative_Yield <- d$Yield / d$Crop_Median
  d$Log_Yield_Reported <- log1p(d$Yield) # Within-crop description only; not a common physical scale.
  d$Log_Relative <- log1p(d$Relative_Yield)
  d$Z <- d$Log_Relative
  d$High <- as.integer(d$Yield > d$Crop_Median)
  d$Yield_Class <- factor(ifelse(d$High == 1, "High", "Low"), levels = c("Low", "High"))
  d$Year <- d$Crop_Year
  d$Log_Area <- log1p(d$Area)
  d$Log_Rain <- log1p(d$Annual_Rainfall)
  d
}

learn_preprocessing <- function(d) {
  list(levels = lapply(d[c("Crop", "State", "Season")], function(x) sort(unique(as.character(x)))),
    centre = vapply(d[numeric_features], mean, numeric(1)), spread = vapply(d[numeric_features], sd,
    numeric(1)))
}

apply_preprocessing <- function(d, preprocessing) {
  for (v in names(preprocessing$levels)) {
    unseen <- setdiff(as.character(d[[v]]), preprocessing$levels[[v]])
    if (length(unseen)) stop("Unsupported new category in ", v, ": ", paste(unique(unseen), collapse = ", "))
    d[[v]] <- factor(d[[v]], levels = preprocessing$levels[[v]])
  }
  if (any(!is.finite(preprocessing$spread) | preprocessing$spread <= 0)) stop("A numeric feature is constant.")
  for (v in numeric_features) d[[v]] <- (d[[v]] - preprocessing$centre[v]) / preprocessing$spread[v]
  d
}

make_design_matrix <- function(d, formula) model.matrix(formula, data = d)[, -1, drop = FALSE]

# ------------------------------------------------------------
# 0.6 Safe fitting and shared model helpers
# ------------------------------------------------------------

fail <- function(e) list(type = "failed", reason = conditionMessage(e))

safe_fit <- function(expression) tryCatch(expression, error = fail)

historical_prediction <- function(a, b, response, state = FALSE) {
  crop <- a %>% group_by(Crop) %>% summarise(Value = mean(.data[[response]]), .groups = "drop")
  fallback <- crop$Value[match(b$Crop, crop$Crop)]
  if (!state) return(fallback)
  grouped <- a %>% group_by(Crop, State) %>% summarise(N = n(), Total = sum(.data[[response]]), .groups = "drop") %>%
    left_join(crop, by = "Crop") %>% mutate(Estimate = (Total + shrinkage_k * Value) / (N + shrinkage_k))
  i <- match(paste(b$Crop, b$State, sep = "|||"), paste(grouped$Crop, grouped$State, sep = "|||"))
  p <- grouped$Estimate[i]
  p[is.na(p)] <- fallback[is.na(p)]
  p
}

glmnet_fit <- function(a, response, family, alpha, formula) {
  x <- make_design_matrix(a, formula)
  keep <- apply(x, 2, function(v) length(unique(v)) > 1)
  x <- x[, keep, drop = FALSE]
  if (ncol(x) < 2) stop("Insufficient variable predictors.")
  list(type = "glmnet", formula = formula, columns = colnames(x),
    fit = glmnet::glmnet(x, a[[response]], family = family, alpha = alpha, lambda = lambda_grid, standardize = TRUE))
}

# ------------------------------------------------------------
# 0.7 Regression models and evaluation
# ------------------------------------------------------------

fit_regression <- function(a) {
  models <- list(Crop_median = list(type = "median"), Crop_history = list(type = "history", a = a,
    state = FALSE),
    Crop_state_history = list(type = "history", a = a, state = TRUE))
  models$OLS <- safe_fit(list(type = "lm", fit = lm(main_formula, data = a)))
  if (models$OLS$type != "failed") {
    models$Backward <- safe_fit(list(type = "lm", fit = MASS::stepAIC(models$OLS$fit, direction = "backward",
      trace = FALSE)))
    models$Forward <- safe_fit(list(type = "lm", fit = MASS::stepAIC(lm(Z ~ 1, data = a), scope = list(lower = ~1,
      upper = main_formula), direction = "forward", trace = FALSE)))
  }
  models$Ridge <- safe_fit(glmnet_fit(a, "Z", "gaussian", 0, main_xformula))
  models$LASSO <- safe_fit(glmnet_fit(a, "Z", "gaussian", 1, main_xformula))
  models$Interaction_Ridge <- safe_fit(glmnet_fit(a, "Z", "gaussian", 0, interaction_formula))
  models
}

predict_regression <- function(object, b, lambda = NA_real_) {
  if (object$type == "failed") stop(object$reason)
  if (object$type == "median") p <- rep(log(2), nrow(b))
  if (object$type == "history") p <- historical_prediction(object$a, b, "Z", object$state)
  if (object$type == "lm") p <- predict(object$fit, newdata = b, rankdeficient = "NA")
  if (object$type == "glmnet") p <- predict(object$fit, newx = make_design_matrix(b, object$formula)[,
    object$columns, drop = FALSE], s = lambda)
  p <- pmax(0, as.numeric(p))
  if (length(p) != nrow(b) || any(!is.finite(p))) stop("Non-finite/non-estimable regression prediction.")
  p
}

mean_by_crop <- function(v, crop) mean(vapply(split(v, as.character(crop)), mean, numeric(1)))

regression_metrics <- function(b, p) {
  e <- b$Z - p
  relative_error <- abs(b$Relative_Yield - expm1(p))
  if (any(!is.finite(relative_error))) stop("Unstable inverse transformation.")
  r2 <- if (sum((b$Z - mean(b$Z))^2) > 0) 1 - sum(e^2) / sum((b$Z - mean(b$Z))^2) else NA_real_
  c(Macro_RMSE_Z = mean(sqrt(vapply(split(e^2, as.character(b$Crop)), mean, numeric(1)))),
    RMSE_Z = sqrt(mean(e^2)), MAE_Z = mean(abs(e)), R2_Z = r2, Macro_MAE_Relative = mean_by_crop(relative_error,
    b$Crop))
}

specifications <- function(models) bind_rows(lapply(names(models), function(m) {
    data.frame(Model = m, Lambda = if (models[[m]]$type == "glmnet") lambda_grid else NA_real_)
  }))

score_regression <- function(models, specs, b) bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    s <- specs[i, ]
    tryCatch(cbind(s, as.data.frame(as.list(regression_metrics(b, predict_regression(models[[s$Model]],
    b, s$Lambda)))), Failure = ""),
    error = function(e) cbind(s, Macro_RMSE_Z = NA_real_, RMSE_Z = NA_real_, MAE_Z = NA_real_, R2_Z = NA_real_,
    Macro_MAE_Relative = NA_real_, Failure = conditionMessage(e)))
  }))

# ------------------------------------------------------------
# 0.8 Classification models and evaluation
# ------------------------------------------------------------

fit_classification <- function(a) {
  out <- list(Majority = list(type = "majority", probability = mean(a$High)),
    Crop_probability = list(type = "history", a = a, state = FALSE),
    Crop_state_probability = list(type = "history", a = a, state = TRUE))
  out$Logistic <- safe_fit({
      fit <- glm(class_formula, data = a, family = binomial(), control = glm.control(maxit = 100))
      if (!fit$converged || anyNA(coef(fit)) || any(fitted(fit) < 1e-8 | fitted(fit) > 1 - 1e-8)) stop("Unpenalized logistic fit unstable/separation suspected; use penalized comparison.")
      list(type = "glm", fit = fit)
    })
  out$Ridge_logistic <- safe_fit(glmnet_fit(a, "High", "binomial", 0, main_xformula))
  out$Naive_Bayes <- safe_fit(list(type = "nb", fit = e1071::naiveBayes(Yield_Class ~ Crop + State +
    Season + Year + Log_Area + Log_Rain, data = a, laplace = 1)))
  out
}

predict_classification <- function(object, b, lambda = NA_real_) {
  if (object$type == "failed") stop(object$reason)
  if (object$type == "majority") p <- rep(object$probability, nrow(b))
  if (object$type == "history") p <- historical_prediction(object$a, b, "High", object$state)
  if (object$type == "glm") p <- predict(object$fit, newdata = b, type = "response")
  if (object$type == "nb") p <- predict(object$fit, newdata = b, type = "raw")[, "High"]
  if (object$type == "glmnet") p <- predict(object$fit, newx = make_design_matrix(b, object$formula)[,
    object$columns, drop = FALSE], s = lambda, type = "response")
  p <- as.numeric(p)
  if (length(p) != nrow(b) || any(!is.finite(p) | p < 0 | p > 1)) stop("Invalid class probability.")
  p
}

roc_auc <- function(y, p) {
  n1 <- sum(y == 1);
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

binary_metrics <- function(y, p) {
  predicted <- as.integer(p >= 0.5)
  sensitivity <- if (any(y == 1)) mean(predicted[y == 1] == 1) else NA_real_
  specificity <- if (any(y == 0)) mean(predicted[y == 0] == 0) else NA_real_
  c(Accuracy = mean(predicted == y), Sensitivity = sensitivity, Specificity = specificity,
    Balanced_Accuracy = (sensitivity + specificity) / 2, AUC = roc_auc(y, p), Brier = mean((y - p)^2))
}

class_metrics <- function(b, p) {
  ba <- vapply(split(seq_len(nrow(b)), as.character(b$Crop)), function(i) binary_metrics(b$High[i],
    p[i])["Balanced_Accuracy"], numeric(1))
  c(Macro_BA = if (any(is.finite(ba))) mean(ba[is.finite(ba)]) else NA_real_, binary_metrics(b$High,
    p), Crops_Both_Classes = sum(is.finite(ba)))
}

score_classification <- function(models, specs, b) bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    s <- specs[i, ]
    tryCatch(cbind(s, as.data.frame(as.list(class_metrics(b, predict_classification(models[[s$Model]],
    b, s$Lambda)))), Failure = ""),
    error = function(e) cbind(s, Macro_BA = NA_real_, Accuracy = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
    Balanced_Accuracy = NA_real_, AUC = NA_real_, Brier = NA_real_, Crops_Both_Classes = NA_integer_,
    Failure = conditionMessage(e)))
  }))

# ------------------------------------------------------------
# 0.9 Calibration, ROC and interval helpers
# ------------------------------------------------------------

conformal_quantile <- function(scores) {
  if (any(!is.finite(scores))) stop("Invalid calibration score.")
  n <- length(scores);
  k <- ceiling((n + 1) * (1 - alpha_interval))
  if (k > n) Inf else sort(scores)[k]
}

roc_points <- function(y, p, model) {
  if (length(unique(y)) < 2) return(data.frame())
  order_p <- order(p, decreasing = TRUE);
  y <- y[order_p];
  p <- p[order_p]
  last <- !duplicated(p, fromLast = TRUE)
  data.frame(Model = model, FPR = c(0, (cumsum(y == 0) / sum(y == 0))[last]), TPR = c(0, (cumsum(y == 1) / sum(y == 1))[last]))
}

make_intervals <- function(q, method) {
  low_z <- pmax(0, test_p - q);
  high_z <- test_p + q
  data.frame(Row_ID = test$Row_ID, Crop = as.character(test$Crop), Method = method,
    Lower_Z = low_z, Upper_Z = high_z, Covered = test$Z >= low_z & test$Z <= high_z,
    Lower_Reported = test$Crop_Median * expm1(low_z), Upper_Reported = test$Crop_Median * expm1(high_z),
    Width_Relative = expm1(high_z) - expm1(low_z))
}

# ============================================================
## Module 1: Task 3 - Data Audit & Preprocessing
# ============================================================
cat("\n>>> Module 1: Task 3 - Data Audit & Preprocessing\n")

# ------------------------------------------------------------
# 1.1 Import, clean and audit the six-crop scope
# ------------------------------------------------------------

raw <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
raw <- raw %>% mutate(across(where(is.character), trimws))
required_columns <- c("Crop", "State", "Season", "Crop_Year", "Area", "Annual_Rainfall", "Yield")
if (!all(required_columns %in% names(raw))) stop("Missing required CSV columns.")
numeric_columns <- c("Crop_Year", "Area", "Annual_Rainfall", "Yield")
if (!all(vapply(raw[numeric_columns], is.numeric, logical(1)))) stop("A numeric field was read as text; inspect the CSV.")
raw$Row_ID <- seq_len(nrow(raw))
raw$Crop_Source <- raw$Crop
raw$Crop[raw$Crop == "Cotton(lint)"] <- "Cotton"
scoped <- raw[raw$Crop %in% crops, ]
print_table("Missing_values", data.frame(Variable = required_columns, Missing = colSums(is.na(scoped[required_columns]))))
valid <- complete.cases(scoped[required_columns]) & Reduce(`&`, lapply(scoped[numeric_columns],
  is.finite))
valid <- valid & scoped$Area > 0 & scoped$Annual_Rainfall >= 0 & scoped$Yield >= 0
valid <- valid & scoped$Crop_Year == floor(scoped$Crop_Year)
valid <- valid & nzchar(scoped$State) & nzchar(scoped$Season)
valid[is.na(valid)] <- FALSE
invalid_records <- scoped[!valid, ]
data <- scoped[valid & scoped$Crop_Year >= 1997 & scoped$Crop_Year <= 2019, ]
duplicate_count <- sum(duplicated(data[required_columns]))
if (duplicate_count > 0) stop("Exact analytical duplicates found: adjudicate them before modelling.")
print_table("Data_audit", data.frame(Item = c("Source rows", "Six-crop rows before year restriction",
  "Invalid six-crop rows", "1997-2019 study rows", "Exact duplicates"),
  Count = c(nrow(raw), nrow(scoped), nrow(invalid_records), nrow(data), duplicate_count)))
if (!setequal(unique(data$Crop), crops)) stop("The selected six crops are not all present.")

# ------------------------------------------------------------
# 1.2 Unit status and provisional conversion
# ------------------------------------------------------------

data$Raw_Yield <- data$Yield
data$Conversion_Factor <- ifelse(data$Crop == "Cotton", 0.17, ifelse(data$Crop == "Jute", 0.18,
  1))
if (apply_provisional_fibre_conversion) data$Yield <- data$Raw_Yield * data$Conversion_Factor
data$Unit_Status <- ifelse(data$Crop %in% c("Cotton", "Jute"), "Fibre unit unverified", "Reported unit; metadata verification needed")
unit_audit <- unique(data[c("Crop", "Conversion_Factor", "Unit_Status")])
unit_audit$Conversion_Applied <- apply_provisional_fibre_conversion & unit_audit$Crop %in% c("Cotton",
  "Jute")
print_table("Unit_audit", unit_audit)

# ------------------------------------------------------------
# 1.3 Chronological data splits
# ------------------------------------------------------------

train_raw <- data[data$Crop_Year <= 2016, ]
val_raw <- data[data$Crop_Year == 2017, ]
cal_raw <- data[data$Crop_Year == 2018, ]
test_raw <- data[data$Crop_Year == 2019, ]
print_table("Splits", data.frame(Set = c("Training", "Validation", "Calibration", "Retrospective test"),
  Years = c("1997-2016", "2017", "2018", "2019"), N = c(nrow(train_raw), nrow(val_raw), nrow(cal_raw),
  nrow(test_raw))))
if (any(c(nrow(train_raw), nrow(val_raw), nrow(cal_raw), nrow(test_raw)) == 0)) stop("An evaluation split is empty.")
# Production is outcome-side information; no claim that Yield equals Production/Area.
# Fertilizer/pesticide are not in the primary model: no verified local input resolution.

# ------------------------------------------------------------
# 1.4 Training references and preprocessing
# ------------------------------------------------------------

reference <- train_raw %>% group_by(Crop) %>% summarise(Crop_Median = median(Yield), .groups = "drop")
if (any(reference$Crop_Median <= 0)) stop("A crop has a nonpositive training median.")
print_table("Training_crop_references", reference)
training_unscaled <- add_targets(train_raw)
preprocessing <- learn_preprocessing(training_unscaled)
train <- apply_preprocessing(training_unscaled, preprocessing)
validation <- apply_preprocessing(add_targets(val_raw), preprocessing)
# Calibration/test targets and feature matrices are constructed only after selection.

# ------------------------------------------------------------
# 1.5 Descriptive statistics and outlier review
# ------------------------------------------------------------

print_table("Training_descriptives", training_unscaled %>% group_by(Crop) %>% summarise(N = n(),
  Zeros = sum(Yield == 0),
  Mean_Reported = mean(Yield), Median_Reported = median(Yield), SD_Reported = sd(Yield), Max_Reported = max(Yield),
  Mean_Relative = mean(Relative_Yield), Median_Relative = median(Relative_Yield), .groups = "drop"))
limits <- train %>% group_by(Crop) %>% summarise(Q1 = quantile(Z, 0.25), Q3 = quantile(Z, 0.75),
  Low = quantile(Z, 0.01), High_Cap = quantile(Z, 0.99), .groups = "drop")
review <- left_join(train, limits, by = "Crop")
review$Flag <- review$Z < review$Q1 - 3 * (review$Q3 - review$Q1) | review$Z > review$Q3 + 3 * (review$Q3 - review$Q1)
print_table("Outlier_review", review %>% group_by(Crop) %>% summarise(N = n(), Flagged = sum(Flag),
  .groups = "drop"))
print_table("Flagged_records", review %>% filter(Flag) %>% dplyr::select(Row_ID, Crop, State, Crop_Year,
  Season, Raw_Yield, Z))

# ------------------------------------------------------------
# 1.6 Training-data visualizations
# ------------------------------------------------------------

print_plot(ggplot(training_unscaled, aes(Yield)) + geom_histogram(bins = 25, fill = "steelblue",
  colour = "white") +
  facet_wrap(~Crop, scales = "free") + labs(title = "Training yields: separate crop scales", x = "Reported/scenario yield",
  y = "Records") + theme_bw())
print_plot(ggplot(train, aes(Crop, Z)) + geom_boxplot(fill = "lightblue") + labs(title = "Training crop-relative yields",
  y = "log(1 + yield / training crop median)") + theme_bw())
print_plot(ggplot(training_unscaled, aes(Annual_Rainfall, Z)) + geom_point(alpha = 0.25) + facet_wrap(~Crop) +
  labs(title = "Rainfall and relative yield: training data", y = "Log-relative yield") + theme_bw())

# ============================================================
## Module 2: Task 4 - Exploratory Inference & Plots
# ============================================================
cat("\n>>> Module 2: Task 4 - Exploratory Inference & Plots\n")

# H0 Welch: all group means of Z are equal. H1: at least one differs.
# H0 Levene: group variances of Z are equal. Median-centred version used.
# Repeated crop/state/year records and estimated reference medians mean these
# ordinary p-values are exploratory, not cluster-corrected causal evidence.

# ------------------------------------------------------------
# 2.1 Welch, Levene and Games-Howell comparisons
# ------------------------------------------------------------

inference_rows <- list()
posthoc <- list()
for (group in c("Crop", "Season", "State")) {
  counts <- table(train[[group]])
  accepted <- names(counts[counts >= 6])
  d <- droplevels(train[as.character(train[[group]]) %in% accepted, ])
  variances <- tapply(d$Z, d[[group]], var)
  accepted <- names(variances)[is.finite(variances) & variances > 0]
  d <- droplevels(d[as.character(d[[group]]) %in% accepted, ])
  if (length(unique(d[[group]])) < 2) next
  f <- reformulate(group, response = "Z")
  welch <- oneway.test(f, data = d, var.equal = FALSE)
  levene <- car::leveneTest(f, data = d, center = median)
  ordinary <- lm(f, data = d)
  eta <- 1 - sum(residuals(ordinary)^2) / sum((d$Z - mean(d$Z))^2)
  inference_rows[[group]] <- data.frame(Group = group, N = nrow(d), Excluded_Small_Constant = nrow(train) - nrow(d),
    Welch_F = unname(welch$statistic), DF1 = unname(welch$parameter[1]), DF2 = unname(welch$parameter[2]),
    Welch_p = welch$p.value, Levene_F = levene[1, "F value"], Levene_p = levene[1, "Pr(>F)"], Descriptive_Eta2 = eta)
  # GH already adjusts pairwise p-values within each grouping family.
  posthoc[[group]] <- rstatix::games_howell_test(d, f, conf.level = 0.95, detailed = TRUE)
}
inference <- bind_rows(inference_rows)
inference$Welch_Holm <- p.adjust(inference$Welch_p, "holm")
inference$Levene_Holm <- p.adjust(inference$Levene_p, "holm")
print_table("Inference", inference)
for (g in names(posthoc)) print_table(paste0("Games_Howell_", g), posthoc[[g]])
cat("Eta2 is ordinary one-way descriptive explained variation, not Welch's effect size or a causal effect.\n")
cat("Crop tests concern normalized distribution means, not absolute physical crop productivity.\n")

# Optional proportion comparison, with explicit independence limitations.

# ------------------------------------------------------------
# 2.2 Exploratory comparison of proportions
# ------------------------------------------------------------

prop_counts <- train %>% group_by(Season) %>% summarise(High = sum(High), N = n(), .groups = "drop")
prop_counts$Expected_High <- prop_counts$N * sum(prop_counts$High) / sum(prop_counts$N)
prop_counts$Expected_Low <- prop_counts$N - prop_counts$Expected_High
print_table("Proportion_counts", prop_counts)
proportion_test <- NULL
if (all(prop_counts$Expected_High >= 5 & prop_counts$Expected_Low >= 5)) {
  proportion_test <- prop.test(prop_counts$High, prop_counts$N, correct = FALSE)
  print_table("Exploratory_proportions_test", proportion_test)
} else cat("Proportion chi-square approximation skipped: sparse expected counts.\n")
cat("Proportions H0: all season-level High proportions are equal. Confounding/dependence remain.\n")

# Display the tested response distribution immediately after the inference tables.
print_plot(ggplot(train, aes(Season, Z)) + geom_boxplot(fill = "lightblue") + theme_bw() +
  labs(title = "Training log-relative yield by season", y = "Log-relative yield") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)))
print_plot(ggplot(inference, aes(Group, Descriptive_Eta2)) + geom_col(fill = "steelblue") + theme_bw() +
  labs(title = "Descriptive one-way eta-squared", y = "Eta-squared (not causal)"))

# ============================================================
## Module 3: Task 5 - Regression Modeling & Validation
# ============================================================
cat("\n>>> Module 3: Task 5 - Regression Modeling & Validation\n")

# ------------------------------------------------------------
# 3.1 Fit candidates and select using validation
# ------------------------------------------------------------

models <- fit_regression(train)
reg_validation <- score_regression(models, specifications(models), validation)
print_table("Regression_failures", reg_validation %>% filter(Failure != ""))
reg_best <- reg_validation %>% filter(is.finite(Macro_RMSE_Z)) %>% group_by(Model) %>%
  slice_min(Macro_RMSE_Z, n = 1, with_ties = FALSE) %>% ungroup() %>% arrange(Macro_RMSE_Z, Model)
if (!nrow(reg_best)) stop("No valid regression candidates.")
winner <- reg_best[1, ]
print_table("Regression_validation_best_per_family", reg_best)
print_table("Regression_lambda_boundaries", reg_best %>% filter(!is.na(Lambda)) %>% mutate(At_Boundary = Lambda %in% range(lambda_grid)))
cat("Regression selected on 2017 crop-balanced RMSE:", winner$Model, "\n")

# ------------------------------------------------------------
# 3.2 OLS collinearity and variance diagnostics
# ------------------------------------------------------------

if (models$OLS$type != "failed") {
  vif_result <- tryCatch(car::vif(models$OLS$fit), error = function(e) conditionMessage(e))
  if (is.matrix(vif_result)) {
    vif_result <- data.frame(Term = rownames(vif_result), GVIF = vif_result[, 1], DF = vif_result[,
      2], Adjusted_GVIF = vif_result[, 3], row.names = NULL)
    vif_result$Ordinary_VIF_if_DF1 <- ifelse(vif_result$DF == 1, vif_result$GVIF, NA_real_)
  }
  print_table("VIF_GVIF", vif_result)
  print_table("OLS_Breusch_Pagan_exploratory", tryCatch(lmtest::bptest(models$OLS$fit), error = function(e) conditionMessage(e)))
}
cat("Adjusted GVIF is not ordinary VIF; no automatic 5/10 threshold is applied to multi-DF factors.\n")

# ------------------------------------------------------------
# 3.3 Inspect selected LASSO coefficients
# ------------------------------------------------------------

if (models$LASSO$type == "glmnet") {
  s <- reg_best %>% filter(Model == "LASSO")
  if (nrow(s)) {
    co <- as.matrix(coef(models$LASSO$fit, s = s$Lambda))
    print_table("LASSO_coefficients", data.frame(Term = rownames(co), Coefficient = as.numeric(co)))
    cat("LASSO nonzero predictors:", sum(co[-1, 1] != 0), "of", nrow(co) - 1, "\n")
  }
}

print_plot(ggplot(reg_best, aes(reorder(Model, Macro_RMSE_Z), Macro_RMSE_Z)) + geom_col(fill = "steelblue") +
  coord_flip() + theme_bw() + labs(title = "2017 regression validation", x = "Model", y = "Crop-balanced RMSE on Z"))

# ============================================================
## Module 4: Task 7 - PCA Analysis
# ============================================================
cat("\n>>> Module 4: Task 7 - PCA Analysis\n")

# ------------------------------------------------------------
# 4.1 Fit PCA and summarize variance
# ------------------------------------------------------------

pca_input <- training_unscaled[numeric_features]
pca <- prcomp(pca_input, center = TRUE, scale. = TRUE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
print_table("PCA_variance", data.frame(PC = paste0("PC", seq_along(variance)), Percent = 100 * variance,
  Cumulative = 100 * cumsum(variance)))
print_table("PCA_loadings", data.frame(Variable = rownames(pca$rotation), pca$rotation, row.names = NULL))
cat("PCA components for at least 90% variance:", which(cumsum(variance) >= 0.9)[1], "of", ncol(pca_input),
  "\n")
sampled <- sample(seq_len(nrow(pca$x)), min(300, nrow(pca$x)))
print_plot(ggplot(data.frame(PC = seq_along(variance), Percent = 100 * variance), aes(PC, Percent)) +
  geom_line() + geom_point() + theme_bw() + labs(title = "PCA scree: three numeric predictors"))
print_plot(ggplot(as.data.frame(pca$x[sampled, 1:2, drop = FALSE]), aes(PC1, PC2)) + geom_point(alpha = 0.5) +
  theme_bw() +
  labs(title = "PCA scores: up to 300 training records", x = sprintf("PC1 %.1f%%", 100 * variance[1]),
  y = sprintf("PC2 %.1f%%", 100 * variance[2])))
cat("PCA excludes Yield and categorical predictors. Compression is limited with only three inputs; predictive benefit is untested.\n")

# ============================================================
## Module 5: Task 8 - Classification Modeling
# ============================================================
cat("\n>>> Module 5: Task 8 - Classification Modeling\n")

# ------------------------------------------------------------
# 5.1 Fit and compare classifiers
# ------------------------------------------------------------

class_models <- fit_classification(train)
class_validation <- score_classification(class_models, specifications(class_models), validation)
print_table("Classification_failures", class_validation %>% filter(Failure != ""))
class_best <- class_validation %>% filter(is.finite(Macro_BA)) %>% group_by(Model) %>%
  slice_max(Macro_BA, n = 1, with_ties = FALSE) %>% ungroup() %>% arrange(desc(Macro_BA), Model)
if (!nrow(class_best)) stop("No valid classifiers.")
class_winner <- class_best[1, ]
print_table("Classification_validation_best_per_family", class_best)
print_table("Classification_lambda_boundaries", class_best %>% filter(!is.na(Lambda)) %>% mutate(At_Boundary = Lambda %in% range(lambda_grid)))
for (label in levels(train$Yield_Class)) print_table(paste0("NB_correlations_", label), cor(train[train$Yield_Class == label,
  numeric_features]))
cat("Low correlations do not prove conditional independence; categorical dependence also matters.\n")

print_plot(ggplot(class_best, aes(reorder(Model, Macro_BA), Macro_BA)) + geom_col(fill = "darkseagreen") +
  coord_flip() + theme_bw() + labs(title = "2017 classification validation", x = "Model", y = "Crop-balanced accuracy"))

# ------------------------------------------------------------
# 5.2 Validation ROC curves
# ------------------------------------------------------------

validation_roc <- list()
for (i in seq_len(nrow(class_best))) {
  specification <- class_best[i, ]
  probability <- predict_classification(class_models[[specification$Model]], validation, specification$Lambda)
  validation_roc[[specification$Model]] <- roc_points(validation$High, probability, specification$Model)
}
print_plot(ggplot(bind_rows(validation_roc), aes(FPR, TPR, colour = Model)) + geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + coord_equal() + theme_bw() + labs(title = "2017 validation ROC curves"))

# ============================================================
## Module 6: Task 5 Support - Outliers & Conformal Calibration
# ============================================================
cat("\n>>> Module 6: Task 5 Support - Outliers & Conformal Calibration\n")

# Fixed 2017 evaluation; no calibration/test outcomes used. Main retains all rows.

# ------------------------------------------------------------
# 6.1 Training-only outlier sensitivity
# ------------------------------------------------------------

capped <- review
capped$Z <- pmax(capped$Low, pmin(capped$High_Cap, capped$Z))
scenarios <- list(Keep_all = train, Cap_training_Z_1_99 = capped, Remove_training_3IQR = review[!review$Flag,
  ])
best_fitted <- reg_best %>% filter(!Model %in% c("Crop_median", "Crop_history", "Crop_state_history")) %>%
  slice_head(n = 1)
sensitivity_specs <- unique(bind_rows(winner[c("Model", "Lambda")], best_fitted[c("Model", "Lambda")]))
sensitivity <- list()
for (scenario in names(scenarios)) {
  alternative <- fit_regression(scenarios[[scenario]])
  s <- score_regression(alternative, sensitivity_specs, validation)
  sensitivity[[scenario]] <- cbind(Scenario = scenario, Training_N = nrow(scenarios[[scenario]]),
    s)
}
print_table("Outlier_sensitivity", bind_rows(sensitivity))

# ------------------------------------------------------------
# 6.2 Constant-unit conversion sensitivity
# ------------------------------------------------------------

converted_raw <- train_raw
converted_raw$Yield <- train_raw$Raw_Yield * train_raw$Conversion_Factor
converted_ref <- converted_raw %>% group_by(Crop) %>% summarise(Crop_Median = median(Yield), .groups = "drop")
raw_ref <- train_raw %>% group_by(Crop) %>% summarise(Crop_Median = median(Raw_Yield), .groups = "drop")
raw_scenario <- train_raw;
raw_scenario$Yield <- raw_scenario$Raw_Yield
converted_z <- add_targets(converted_raw, converted_ref)$Z
raw_z <- add_targets(raw_scenario, raw_ref)$Z
unit_difference <- max(abs(converted_z - raw_z))
print_table("Unit_invariance", data.frame(Max_Absolute_Z_Difference = unit_difference, Pass_1e_10 = unit_difference < 1e-10))
if (unit_difference >= 1e-10) stop("Unexpected failure of constant-scale invariance.")
cat("Unit check proves algebraic invariance, NOT that fibre conversions are correct.\n")

print_plot(ggplot(bind_rows(sensitivity) %>% filter(is.finite(Macro_RMSE_Z)), aes(Scenario, Macro_RMSE_Z,
  fill = Model)) +
  geom_col(position = "dodge") + theme_bw() + labs(title = "2017 outlier sensitivity", y = "Crop-balanced RMSE on Z") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)))
# Freeze the original 1997-2016 crop medians for targets AND class labels.

# ------------------------------------------------------------
# 6.3 Refit selected models through 2017
# ------------------------------------------------------------

final_raw <- bind_rows(train_raw, val_raw)
final_unscaled <- add_targets(final_raw)
final_preprocessing <- learn_preprocessing(final_unscaled)
final_train <- apply_preprocessing(final_unscaled, final_preprocessing)
calibration <- apply_preprocessing(add_targets(cal_raw), final_preprocessing)
final_models <- fit_regression(final_train)
final_class_models <- fit_classification(final_train)
selected_model <- final_models[[winner$Model]]
selected_classifier <- final_class_models[[class_winner$Model]]
if (selected_model$type == "failed" || selected_classifier$type == "failed") stop("Selected method failed on refit; do not choose a substitute using test results.")

# ------------------------------------------------------------
# 6.4 Selected-model residual diagnostics
# ------------------------------------------------------------

fitted_z <- predict_regression(selected_model, final_train, winner$Lambda)
diagnostics <- data.frame(Fitted = fitted_z, Residual = final_train$Z - fitted_z)
print_plot(ggplot(diagnostics, aes(Fitted, Residual)) + geom_point(alpha = 0.25) + geom_hline(yintercept = 0,
  colour = "red") + theme_bw() + labs(title = paste("Training diagnostics:", winner$Model)))
print_plot(ggplot(diagnostics, aes(sample = Residual)) + stat_qq(alpha = 0.3) + stat_qq_line(colour = "red") +
  theme_bw())
print_plot(ggplot(diagnostics, aes(Fitted, sqrt(abs(Residual)))) + geom_point(alpha = 0.25) + theme_bw() +
  labs(y = "Sqrt absolute residual", title = "Residual spread; not standardized scale-location"))

# ------------------------------------------------------------
# 6.5 Calibrate prediction intervals using 2018
# ------------------------------------------------------------

calibration_prediction <- predict_regression(selected_model, calibration, winner$Lambda)
calibration_scores <- abs(calibration$Z - calibration_prediction)
pooled_q <- conformal_quantile(calibration_scores)
crop_q <- bind_rows(lapply(crops, function(crop) {
    scores <- calibration_scores[as.character(calibration$Crop) == crop]
    data.frame(Crop = crop, Calibration_N = length(scores), Rank = ceiling((length(scores) + 1) * 0.95),
    Q = conformal_quantile(scores))
  }))
print_table("Conformal_calibration", data.frame(N = nrow(calibration), Rank = ceiling((nrow(calibration) +
  1) * 0.95), Nominal = 0.95, Pooled_Margin_Z = pooled_q))
print_table("Crop_conformal_calibration", crop_q)
cat("Infinite crop-specific margins are reported honestly when too few calibration observations exist.\n")
cat("95% is nominal: chronological drift and repeated records mean exchangeability is unverified.\n")

# ============================================================
## Module 7: Task 5/8 - 2019 Test Evaluation
# ============================================================
cat("\n>>> Module 7: Task 5/8 - 2019 Test Evaluation\n")

# ------------------------------------------------------------
# 7.1 Regression evaluation on 2019
# ------------------------------------------------------------

test <- apply_preprocessing(add_targets(test_raw), final_preprocessing)
reg_test <- score_regression(final_models, reg_best[c("Model", "Lambda")], test)
print_table("Regression_2019", reg_test)
baseline_error <- reg_test$Macro_RMSE_Z[reg_test$Model == "Crop_state_history"]
if (length(baseline_error) == 1 && is.finite(baseline_error) && baseline_error > 0) {
  print_table("Regression_gains_vs_crop_state", reg_test %>% mutate(Gain_Percent = 100 * (baseline_error - Macro_RMSE_Z) / baseline_error))
}

# ------------------------------------------------------------
# 7.2 Crop-specific errors in reported units
# ------------------------------------------------------------

reported_metrics <- list()
for (i in seq_len(nrow(reg_best))) {
  s <- reg_best[i, ]
  p <- tryCatch(predict_regression(final_models[[s$Model]], test, s$Lambda), error = function(e) NULL)
  if (!is.null(p)) {
    reported_metrics[[s$Model]] <- data.frame(Crop = test$Crop, Actual = test$Yield, Predicted = test$Crop_Median * expm1(p)) %>%
      group_by(Crop) %>% summarise(N = n(), RMSE_Reported = sqrt(mean((Actual - Predicted)^2)),
      MAE_Reported = mean(abs(Actual - Predicted)), .groups = "drop") %>% mutate(Model = s$Model)
  }
}
print_table("Reported_unit_errors_by_crop", bind_rows(reported_metrics))
test_p <- predict_regression(selected_model, test, winner$Lambda)
print_plot(ggplot(data.frame(Predicted = test_p, Actual = test$Z), aes(Predicted, Actual)) + geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "red") + theme_bw() + labs(title = "2019 actual versus predicted log-relative yield"))

# ------------------------------------------------------------
# 7.3 Classification evaluation and confusion matrices
# ------------------------------------------------------------

class_test <- score_classification(final_class_models, class_best[c("Model", "Lambda")], test)
print_table("Classification_2019", class_test)
test_prob <- predict_classification(selected_classifier, test, class_winner$Lambda)
confusion <- table(Actual = test$Yield_Class, Predicted = factor(ifelse(test_prob >= 0.5, "High",
  "Low"), levels = c("Low", "High")))
print_table("Selected_confusion_matrix", confusion)
all_confusions <- list()
roc_rows <- list()
for (i in seq_len(nrow(class_best))) {
  s <- class_best[i, ]
  p <- tryCatch(predict_classification(final_class_models[[s$Model]], test, s$Lambda), error = function(e) NULL)
  if (!is.null(p)) {
    all_confusions[[s$Model]] <- table(Actual = test$Yield_Class, Predicted = factor(ifelse(p >= 0.5,
      "High", "Low"), levels = c("Low", "High")))
    roc_rows[[s$Model]] <- roc_points(test$High, p, s$Model)
  }
}
print_table("All_confusion_matrices", all_confusions)
print_plot(ggplot(bind_rows(roc_rows), aes(FPR, TPR, colour = Model)) + geom_line() + geom_abline(slope = 1,
  intercept = 0, linetype = 2) + coord_equal() + theme_bw() + labs(title = "2019 ROC curves: validation-fixed models"))

# ------------------------------------------------------------
# 7.4 Performance by crop
# ------------------------------------------------------------

crop_metrics <- bind_rows(lapply(split(seq_len(nrow(test)), as.character(test$Crop)), function(i) {
    cbind(data.frame(Crop = as.character(test$Crop[i[1]]), N = length(i)),
    as.data.frame(as.list(regression_metrics(test[i, ], test_p[i]))), as.data.frame(as.list(binary_metrics(test$High[i],
    test_prob[i]))))
  }))
  print_table("Selected_methods_by_crop", crop_metrics)

  pooled_intervals <- make_intervals(rep(pooled_q, nrow(test)), "Pooled_Z")
  crop_intervals <- make_intervals(crop_q$Q[match(as.character(test$Crop), crop_q$Crop)], "Crop_specific_Z")
  intervals <- bind_rows(pooled_intervals, crop_intervals)
  print_table("Interval_overall", intervals %>% group_by(Method) %>% summarise(N = n(), Coverage = mean(Covered),
  Mean_Relative_Width = mean(Width_Relative), Infinite_Intervals = sum(!is.finite(Upper_Z)), .groups = "drop"))
  coverage <- intervals %>% group_by(Method, Crop) %>% summarise(N = n(), Coverage = mean(Covered),
  Mean_Relative_Width = mean(Width_Relative), Mean_Reported_Width = mean(Upper_Reported - Lower_Reported),
  Infinite_Intervals = sum(!is.finite(Upper_Z)), .groups = "drop")
  print_table("Interval_by_crop", coverage)
  print_plot(ggplot(coverage, aes(Crop, Coverage, fill = Method)) + geom_col(position = "dodge") +
  geom_hline(yintercept = 0.95, linetype = 2) +
  scale_y_continuous(limits = c(0, 1)) + theme_bw() + labs(title = "2019 interval coverage: small crop samples are uncertain"))

# ============================================================
## Module 8: Task 9-11 - Dashboard & Framework Templates
# ============================================================
  cat("\n>>> Module 8: Task 9-11 - Dashboard & Framework Templates\n")

  annual <- training_unscaled %>% group_by(Crop, Crop_Year) %>% summarise(N = n(), States = n_distinct(State),
  Median_Relative = median(Relative_Yield), .groups = "drop")
  print_table("Annual_training_summary", annual)
  print_plot(ggplot(annual, aes(Crop_Year, Median_Relative)) + geom_line() + facet_wrap(~Crop) + theme_bw() +
  labs(title = "Descriptive annual medians; changing state composition may affect trends"))
  time_outline <- data.frame(Topic = c("Trend", "Seasonality", "ARIMA", "Stochastic modelling", "Decision use"),
  Discussion = c("Choose a consistent crop-state series and audit changing coverage.",
  "Annual observations cannot identify within-year seasonal cycles; Season labels are not equally spaced timestamps.",
  "Require a regular series, stationarity/differencing checks, residual checks and rolling forecast evaluation; no ARIMA fitted here.",
  "Consider uncertainty in rainfall, yield and costs; document assumed distributions and dependence.",
  "Use validated forecasts for planning with uncertainty, not automatic fertilizer prescriptions."))
  print_table("Task9_outline", time_outline)

  dashboard <- data.frame(Row_ID = test$Row_ID, Crop = as.character(test$Crop), State = as.character(test$State),
  Year = test$Crop_Year,
  Predicted_Relative = expm1(test_p), Predicted_Reported = test$Crop_Median * expm1(test_p), Probability_High = test_prob,
  Lower_Relative = expm1(pooled_intervals$Lower_Z), Upper_Relative = expm1(pooled_intervals$Upper_Z),
  Unit_Status = test$Unit_Status, Model = winner$Model, Classifier = class_winner$Model)
  dashboard <- left_join(dashboard, coverage %>% filter(Method == "Pooled_Z") %>% dplyr::select(Crop,
  Retrospective_Crop_Coverage = Coverage), by = "Crop")
  print_table("Dashboard_preview", dashboard, n = 12)
  cat("Framework: source audit -> historical references -> validation -> fixed model -> calibration -> dashboard -> human review -> outcome monitoring.\n")
  cat("Dashboard values are retrospective demonstrations; no operational alert threshold or intervention benefit has been validated.\n")

  expert_template <- data.frame(Question = c("Is the six-crop problem relevant?", "Are yield units and aggregation interpretable?",
  "When are the proposed inputs actually available?", "Are the baseline and crop-relative target useful?",
  "Are uncertainty intervals useful for any concrete decision?", "What are the costs of missed versus false alerts?",
  "What implementation barriers exist?", "What recommendation should change?"),
  Expert_Response = NA_character_, Evidence_Reference = NA_character_, Team_Action = NA_character_,
  Status = "NOT CONSULTED")
  print_table("Expert_feedback_template", expert_template)
  cat("Sensitivity analysis is model validation, not industry expert validation. No expert feedback is fabricated.\n")

  literature_matrix <- data.frame(Paper_ID = sprintf("P%02d", 1:15), Citation_DOI = NA_character_,
  Year = NA_integer_,
  Journal = NA_character_, Peer_Reviewed = NA, Dataset_Context = NA_character_, Method = NA_character_,
  Validation_Design = NA_character_, Findings = NA_character_, Limitations = NA_character_, Relevance_To_Our_Study = NA_character_)
  print_table("Literature_matrix_template", literature_matrix)
  design_outline <- data.frame(Design = c("CRD", "RCBD"),
  Proposal = c("Randomly assign approved treatments to comparable plots within one crop and site.",
  "Block plots by a measured fertility/terrain gradient, randomize all treatments within each block."),
  Strength = c("Simple if experimental units are reasonably homogeneous.", "Controls a known nuisance source of variation."),
  Limitation = c("Field heterogeneity can reduce precision.", "Needs suitable blocks, replication, space, and careful logistics."))
  print_table("Experimental_design_outline", design_outline)
  bayesian_outline <- data.frame(Method = c("Naive Bayes", "Bayesian regression", "Bayesian decision analysis"),
  Status = c("Implemented and evaluated", "Proposed; not fitted", "Proposed; not implemented"),
  Discussion = c("Check conditional dependence, probability quality and temporal transfer.",
  "A hierarchical crop-state model could partially pool sparse groups; specify priors, convergence and predictive checks.",
  "Combine validated probabilities with expert-supplied consequences and costs; no invented utility values."))
  print_table("Bayesian_methods_outline", bayesian_outline)
  cat("\nTASK 1: Describe a proposed agricultural planning client, stakeholders and decisions; do not invent an engagement.\n")
  cat("TASK 12: Link each recommendation to results, literature, risks and actual expert feedback when obtained.\n")
  cat("Task 2 needs 15 peer-reviewed papers, at least 10 journal papers and at least 5 within five years of submission.\n")
  cat("2019 was previously inspected: this is retrospective, not a new untouched confirmatory test.\n")
  cat("The new six-crop scope differs from previous work: all numerical results must be regenerated.\n")
  cat("Do not claim a causal fertilizer effect, verified tonnes/ha, normal residuals, or guaranteed time-series coverage.\n")
  cat("\n===== ANALYSIS COMPLETE: EXPERT VALIDATION AND LITERATURE VERIFICATION STILL REQUIRED =====\n")
  print(sessionInfo())
  project <- list(tables = tables, data = data, invalid_records = invalid_records, reference = reference,
  regression_validation_all = reg_validation, classification_validation_all = class_validation,
  regression_selected = winner, classification_selected = class_winner, fitted_regression_models = final_models,
  fitted_classification_models = final_class_models, final_preprocessing = final_preprocessing,
  test_predictions_Z = test_p, test_probabilities = test_prob, intervals = intervals, dashboard = dashboard,
  pca = pca, expert_template = expert_template, literature_matrix = literature_matrix,
  source_md5 = unname(tools::md5sum(csv_path)), session = sessionInfo())

# Full results are in project; tables were also available immediately after printing.
# Examples: project$tables$Inference; project$tables$Regression_2019
