#' Cox Regression + Group Regularization
#' 
#' @param alpha elastic-net mixing parameter (1: lasso / 0: ridge)
#' @param lambda.opt choice of the penalty parameter lambda
#' \itemize{
#'   \item 1: value of lambda that gives the minimum mean cross-validated error
#'   \item 2: largest value of lambda such that error is within 1 standard error of the minimum;
#'   set as default to follow the convention of the {glmnet} package
#' }
#' @param interaction generate all possible two-way interaction terms
#' @param hierarchy hierarchy/marginality principle
#' 
cox_fit <- function(formula, data, alpha = 1, n_fold = 5, lambda.opt = 2,
                    interaction = FALSE, hierarchy = FALSE, control = coxph.control(iter.max = 50), ...) {
  
  start_time <- Sys.time()
  lambda <- paste0("lambda.", c("min", "1se"))[lambda.opt]
  
  if(interaction) {
    formula <- update(terms(formula, data = data), ~ .^2)
  }
  
  selector <- grpreg_cox_fit(formula, data, alpha = alpha, n_fold = n_fold, ...)$fit
  xvar <- as.character(predict(selector, type = "groups", lambda = selector[[lambda]]))
  
  if(hierarchy) {
    xvar <- union(xvar, unlist(strsplit(xvar, ':', fixed = TRUE)))
  }
  
  formula.new <- reformulate(xvar, formula[[2]])
  fit <- coxph(formula.new, data, control = control)
  
  end_time <- Sys.time()
  list(fit = fit,
       time_fit = end_time - start_time)
}

grpreg_cox_fit <- function(formula, data, alpha = 1, n_fold = 5,
                           # var_trans = NULL, class_trans = NULL, trans = ~ .x,
                           ...) {
  
  start_time <- Sys.time()
  mf <- model.frame(formula, data)
  xvar <- labels(terms(mf))
  
  # Not used:
  #
  # if(!is.null(c(var_trans, class_trans))) {
  #   cl <- attr(terms(mf), "dataClasses")
  #   var_trans <- unique(c(var_trans, names(cl)[cl %in% class_trans]))
  #   
  #   fmt <- gsub("\\.x?", "%s", update(trans, NULL ~ .))[-1]
  #   xvar <- sapply(strsplit(xvar, ':', fixed = TRUE), \(x) {
  #     ind_trans <- x %in% var_trans
  #     x[ind_trans] <- sprintf(fmt, x[ind_trans])
  #     paste(x, collapse = ':')
  #   })
  #   
  #   formula <- reformulate(xvar, formula[[2]])
  # }
  
  designMat <- model.matrix(formula, data)
  cvfit <- grpreg::cv.grpsurv(X = designMat[, -1], # exclude intercept
                              y = model.response(mf),
                              group = xvar[attr(designMat, "assign")],
                              alpha = alpha,
                              nfolds = n_fold,
                              returnY = TRUE, # to calculate AUC
                              ...)
  
  cvfit$lambda.1se <- with(cvfit, max(lambda[cve < (cve + cvse)[which.min(cve)]]))
  # cvfit$lambda.auc <- cvfit$lambda[which.max(grpreg::AUC(cvfit))]
  
  end_time <- Sys.time()
  list(fit = cvfit,
       formula = formula,
       time_fit = end_time - start_time)
}

#' Oblique Random Survival Forests
#' 
#' @param alpha significance level of a log-rank test required to split a node
#' 
orsf_fit <- function(formula, data, param_expand = FALSE,
                     n_tree = 500, colsample = NULL, alpha = 0.01,
                     n_split = 10, n_retry = 5,
                     doCV = TRUE, n_fold = 5, n_repeat = 1, eval_horizon = NULL,
                     n_core = availableCores(), progress = FALSE) {
  
  start_time <- Sys.time()
  
  nvar <- length(labels(terms(formula, data = data)))
  if(is.null(colsample)) {
    colsample <- 1 / sqrt(nvar)
  }
  
  parm_lst <- list(n_tree = as.integer(n_tree),
                   colsample = colsample,
                   alpha = alpha)
  
  param_grid <- if(param_expand) {
    expand.grid(parm_lst, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  } else {
    as.data.frame(parm_lst)
  }
  
  param_grid <- param_grid |>
    transform(mtry = ceiling(colsample * nvar),
              split_min_stat = qchisq(1 - alpha, df = 1)) |>
    unique()
  
  best_ind <- 1L
  if(doCV) {
    formula.null <- update(formula, . ~ 1)
    n_grid <- nrow(param_grid)
    fold_id <- unname(kfold_cv_ind(nrow(data), k = n_fold, repeats = n_repeat))
    
    # Not used:
    # future_pmap() does not speed up in this case.
    #
    # if(n_grid > 1) {
    #   plan(multisession, workers = min(n_grid, n_core))
    #   on.exit(plan(sequential), add = TRUE)
    # }
    
    progressr::with_progress({
      p <- progressr::progressor(steps = n_grid)
      cv_tuner <- purrr::pmap_dfr(param_grid, ~ { # future_pmap_dfr()
        params <- list(...)
        cv <- try({
          score <- purrr::map_dbl(fold_id, \(id) {
            fit <- aorsf::orsf(data = data[-id, ],
                               formula = formula,
                               n_tree = params$n_tree,
                               mtry = params$mtry,
                               n_split = n_split,
                               n_retry = n_retry,
                               split_rule = "logrank",
                               split_min_stat = params$split_min_stat,
                               oobag_pred_type = "none",
                               importance = "none")
    
            if(is.null(eval_horizon)) {
              eval_horizon <- fit$pred_horizon
            }
            
            test <- data[id, ]
            pred <- predict(fit, test, pred_type = "risk", pred_horizon = eval_horizon)
            score <- riskRegression::Score(
              object = list(pred),
              formula = formula.null,
              data = test[all.vars(formula.null)],
              metrics = c("auc"),
              times = eval_horizon
            )$AUC$score$AUC |> mean()
            
            return(score)
          })
          
          data.frame(c_index = mean(score))
        })
        
        p()
        
        if(inherits(cv, "try-error"))
          list(c_index = NA)
        else cv
      })
    }, enable = progress)
    
    best_ind <- which.max(cv_tuner$c_index)
  }
  
  best_param <- c(param_grid[best_ind, ])
  fit <- aorsf::orsf(data = data,
                     formula = formula,
                     n_tree = best_param$n_tree,
                     mtry = best_param$mtry,
                     n_split = n_split,
                     n_retry = n_retry,
                     split_rule = "logrank",
                     split_min_stat = best_param$split_min_stat,
                     oobag_pred_type = "none",
                     importance = "none")
  
  end_time <- Sys.time()
  list(fit = fit,
       tuner = if(doCV) cbind(param_grid, cv_tuner) else NULL,
       time_fit = end_time - start_time)
}

#' Extreme Gradient Boosting (XGBoost)
#' 
#' @param penalty amount of regularization
#' @param mixture elastic-net mixing parameter
#' 
xgbst_cox_fit <- function(formula, data, x = TRUE, param_expand = FALSE,
                          learn_rate = 0.01, tree_depth = 5, min_n = 10,
                          colsample = NULL, subsample = 2/3,
                          penalty = 1, mixture = 0, n_round = 5000, device = "cpu",
                          doCV = TRUE, n_fold = 5, n_repeat = 1,
                          n_core = availableCores(), progress = FALSE) {
  
  start_time <- Sys.time()
  mf <- model.frame(formula, data, na.action = NULL)
  resp <- model.response(mf)
  time <- resp[, "time"]; status <- resp[, "status"]
  
  xdat <- mf[labels(terms(mf))]
  ydat <- time * ifelse(status == 1, 1, -1)
  xy.DMat <- xgb.DMatrix(data = xdat, label = ydat)

  if(is.null(colsample)) {
    colsample <- 1 / sqrt(ncol(xdat))
  }
  
  parm_lst <- list(eta = learn_rate,
                   max_depth = as.integer(tree_depth),
                   min_child_weight = as.integer(min_n),
                   colsample_bynode = colsample,
                   subsample = subsample,
                   penalty = penalty,
                   mixture = mixture)
  
  param_grid <- if(param_expand) {
    expand.grid(parm_lst, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  } else {
    as.data.frame(parm_lst)
  }
  
  param_grid <- param_grid |>
    transform(alpha = penalty * mixture,        # L1 regularization term on weights
              lambda = penalty * (1-mixture),   # L2 regularization term on weights
              objective = 'survival:cox',
              eval_metric = 'cox-nloglik',
              tree_method = "hist",
              device = device) |> 
    subset(select = -c(penalty, mixture)) |>
    unique()
  
  best_ind <- 1L
  if(doCV) {
    n_grid <- nrow(param_grid)
    if(n_grid > 1) {
      plan(multisession, workers = min(n_grid, n_core))
      on.exit(plan(sequential), add = TRUE)
    }
    
    cv_tuner <- replicate(n_repeat, {
      fold_id <- unname(kfold_cv_ind(nrow(xdat), k = n_fold))
      progressr::with_progress({
        p <- progressr::progressor(steps = n_grid)
        future_pmap_dfr(param_grid, ~ {
          cv <- try(
            xgb.cv(
              params = list(...),
              data = xgb.DMatrix(data = xdat, label = ydat),
              folds = fold_id,
              nrounds = n_round,
              early_stopping_rounds = round(sqrt(n_round)),
              maximize = FALSE,
              verbose = FALSE
            )$early_stop
          )
          
          p()
          
          if(inherits(cv, "try-error"))
            list(best_iteration = NA, best_score = NA)
          else
            cv[c("best_iteration", "best_score")]
          
        }, .options = furrr_options(seed = TRUE))
      }, enable = progress)
    }, simplify = FALSE)
    
    cv_tuner <- Reduce(`+`, cv_tuner) / n_repeat
    cv_tuner$best_iteration <- as.integer(cv_tuner$best_iteration)
    best_ind <- which.min(cv_tuner$best_score)
  }
  
  best_param <- c(param_grid[best_ind, ])
  fit <- xgb.train(params = best_param,
                   data = xy.DMat,
                   nrounds = if(doCV) cv_tuner$best_iteration[best_ind] else n_round)
  
  # baseline cumulative hazard H0(t) estimated using Breslow method
  lin_preds_fit <- predict(fit, newdata = xy.DMat, outputmargin = TRUE)
  H0 <- function(tt, smooth = TRUE) {
    gbm::basehaz.gbm(t = time,
                     delta = status,
                     f.x = lin_preds_fit,
                     t.eval = tt,
                     smooth = smooth,
                     cumulative = TRUE)
  }
  
  end_time <- Sys.time()
  list(fit = fit,
       H0 = H0,
       x = if(x) xy.DMat else NULL,
       tuner = if(doCV) cbind(param_grid, cv_tuner) else NULL,
       time_fit = end_time - start_time)
}

#' Neural Networks
#' 
#' @param n_layer number of hidden layers in the network
#' @param frac_layer1 ratio of the width in the first hidden layer to that in the input layer
#' @param width_ratio common ratio between layer widths
#' 
nn_surv_fit <- function(formula, data, nn_fun, param_expand = FALSE,
                        n_layer = 1, frac_layer1 = 2/3, width_ratio = 0.8,
                        epochs = 1000, learn_rate = 0.01, dropout = 0.1, weight_decay = 0,
                        batch_size = 64L, batch_norm = TRUE, activation = "relu", device = "cpu",
                        val_frac = 0.2, n_core = availableCores(), progress = FALSE) {
  
  # Prevent multi-threading conflicts between PyTorch and riskRegression
  Sys.setenv(OMP_NUM_THREADS = 1)
  on.exit(Sys.unsetenv("OMP_NUM_THREADS"), add = TRUE)
  
  start_time <- Sys.time()
  mf <- model.frame(formula, data)
  xvar <- labels(terms(mf))
  yvar <- all.vars(formula)[1:2]
  
  # n0 <- length(xvar)
  n0 <- ncol(model.matrix(formula, data))
  
  parm_lst <- list(n_layer = as.integer(n_layer),
                   frac_layer1 = frac_layer1,
                   width_ratio = width_ratio,
                   learning_rate = learn_rate,
                   dropout = dropout,
                   weight_decay = weight_decay,
                   batch_size = batch_size)
  
  param_grid <- if(param_expand) {
    expand.grid(parm_lst, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  } else {
    as.data.frame(parm_lst)
  }
  param_grid <- unique(param_grid)
  
  n_grid <- nrow(param_grid)
  seed <- sample(.Machine$integer.max, 1)
  
  if(n_grid > 1) {
    plan(multisession, workers = min(n_grid, n_core))
    on.exit(plan(sequential), add = TRUE)
  }
  
  progressr::with_progress({
    p <- progressr::progressor(steps = n_grid)
    
    val_tuner <- future_pmap_dfr(param_grid, ~ {
      
      # ensure the same validation set is used when evaluating each set of parameters
      survivalmodels::set_seed(seed)
      params <- list(...)
      
      val <- try(
        nn_fun(
          data = data[c(xvar, yvar)],
          time_variable = yvar[1],
          status_variable = yvar[2],
          num_nodes = pmax(2, round(
            (n0 * params$frac_layer1) * params$width_ratio^(seq_len(params$n_layer) - 1)
          )),
          epochs = epochs,
          learning_rate = params$learning_rate,
          dropout = params$dropout,
          weight_decay = params$weight_decay,
          batch_size = params$batch_size,
          batch_norm = batch_norm,
          activation = activation,
          frac = val_frac,
          early_stopping = TRUE,
          patience = round(sqrt(epochs)),
          device = device
        )$model$val_metrics$scores$loss |> as.data.frame()
      )
      
      p()
      
      if(inherits(val, "try-error"))
        list(epoch = NA, score = NA)
      else
        val[which.min(val$score), ]
      
    }, .options = furrr_options(seed = TRUE))
  }, enable = progress)
  
  best_ind <- which.min(val_tuner$score)
  best_param <- c(param_grid[best_ind, ])
  fit <- nn_fun(
    data = data[c(xvar, yvar)],
    time_variable = yvar[1],
    status_variable = yvar[2],
    num_nodes = pmax(2, round(
      (n0 * best_param$frac_layer1) * best_param$width_ratio^(seq_len(best_param$n_layer) - 1)
    )),
    epochs = val_tuner$epoch[best_ind],
    learning_rate = best_param$learning_rate,
    dropout = best_param$dropout,
    weight_decay = best_param$weight_decay,
    batch_size = best_param$batch_size,
    batch_norm = batch_norm,
    activation = activation,
    device = device
  )
  
  end_time <- Sys.time()
  list(fit = fit,
       tuner = cbind(param_grid, val_tuner),
       time_fit = end_time - start_time)
}

nn_deepsurv_fit <- function(formula, data, ...) {
  nn_surv_fit(formula, data, nn_fun = survivalmodels::deepsurv, ...)
}

nn_coxtime_fit <- function(formula, data, ...) {
  nn_surv_fit(formula, data, nn_fun = survivalmodels::coxtime, ...)
}

nn_pchazard_fit <- function(formula, data, ...) {
  nn_surv_fit(formula, data, nn_fun = survivalmodels::pchazard, ...)
}
