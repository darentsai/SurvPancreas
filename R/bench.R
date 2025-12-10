#' Cutpoint searcher for a ROC curve
#' 
roc_cutpoint <- function(tpr, tnr, cutpoint.opt) {
  switch(cutpoint.opt,
         which.max(tpr + tnr),
         which.min((1-tpr)^2 + (1-tnr)^2),
         which.max(tpr * tnr))
}

#' Benchmark for survival models
#' 
#' @param test_ind a vector of test set's indices. When test_ind is supplied, test_prop is ignored.
#' @param test_prop fraction of data to use for test set.
#' @param cutpoint.opt cut-point selection criteria
#' \itemize{
#'   \item 1: maximize Youden's index (sensitivity + specificity - 1)
#'   \item 2: minimize the distance between the coordinate (0, 1) and the ROC curve
#'   \item 3: maximum product of sensitivity and specificity
#' }
#' @param roc if TRUE, return the average ROC curve
#' @param return_pred if TRUE, return actual outcomes and predicted mortality risks (1 - survival)
#'   of the validation set without evaluating performance metrics; used for calibration assessment.
#' 
bench <- function(formula, data, model_fn, pred_horiz,
                  test_ind = NULL, test_prop, cutpoint.opt = 2,
                  roc = FALSE, return_pred = FALSE) {
  
  if(is.null(test_ind)) {
    nr <- nrow(data)
    test_ind <- sample(nr, nr * test_prop)
  }
  
  train <- data[-test_ind, ]
  test <- data[test_ind, ]
  formula.null <- update(formula, . ~ 1)
  
  res <- purrr::map2(model_fn, format(names(model_fn)), \(fn, nms) {
    cat("- Model:", nms,
        cli::col_blue(cli::symbol$arrow_right), "training... ")
    model_fit <- fn[["fit"]](formula, train)
    time_fit <- model_fit$time_fit
    cat(prettyunits::pretty_dt(time_fit),
        cli::col_blue(cli::symbol$arrow_right), "testing... ")
    model_pred <- try(
      fn[["pred"]](model_fit, test, pred_horiz)
    )
    
    if(inherits(model_pred, "try-error")) {
      model_pred <- NA
      time_fit <- NA
    }
    cat(if(anyNA(model_pred)) cli::col_yellow(cli::symbol$warning)
        else cli::col_green(cli::symbol$tick), "\n")
    
    return(list(
      time_fit = time_fit,
      risk_preds = 1 - model_pred
    ))
  }) |> purrr::list_transpose(simplify = FALSE)
  
  risk_preds <- res$risk_preds
  time_fit <- res$time_fit
  
  fail <- sapply(risk_preds, anyNA)
  risk_preds <- risk_preds[!fail]
  time_fit <- time_fit[!fail]
  time_fit <- setNames(stack(lapply(time_fit, as.double, units = "secs"))[2:1],
                       c("model", "time"))
  
  if(return_pred) {
    return(list(
      actual_outcome = test[all.vars(formula.null)],
      risk_preds = risk_preds,
      time_fit = time_fit
    ))
  }
  
  perf1 <- lapply(risk_preds, \(pred) {
    lapply(pred_horiz, \(t) {
      
      pred_t <- list(pred[, match(t, pred_horiz)])
      roc_t <- riskRegression::Score(
        object = pred_t, 
        formula = formula.null, 
        data = test[all.vars(formula.null)],
        plots = "ROC",
        times = t
      )$ROC$plotframe
      
      cp <- roc_t[roc_cutpoint(TPR, 1-FPR, cutpoint.opt), ][["risk"]]
      riskRegression::Score(
        object = pred_t, 
        formula = formula.null, 
        data = test[all.vars(formula.null)],
        times = t,
        cutpoints = cp
      )$AUC$cutpoints[, model := NULL]
      
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist(idcol = "model")
  
  sc <- riskRegression::Score(
    object = risk_preds, 
    formula = formula.null, 
    data = test[all.vars(formula.null)],
    metrics = c("auc", "brier"),
    summary = c('IPA', 'ibs'),
    plots = if(roc) "ROC" else NULL,
    times = pred_horiz
  )
  
  perf2 <- cbind(sc$AUC$score[, .(model, times, AUC, se.AUC = se)],
                 sc$IPA$score[model != "Null model", .(IPA)])
  
  # if(roc) {
  #   roc.grid <- seq(0, 1, by = 0.01)
  #   roc.dt <- sc$ROC$plotframe[, {
  #     roc.interp <- approx(x = FPR,
  #                          y = TPR,
  #                          xout = roc.grid,
  #                          ties = mean,
  #                          yleft = 0,
  #                          yright = 1)
  #     list(FPR = roc.interp$x,
  #          TPR = roc.interp$y)
  #   }, by = .(model, times)]
  # }
  
  return(list(
    perf_metric = as.data.frame(perf1[perf2, on = c("model", "times"), ]),
    time_fit = time_fit
    # roc = if(roc) as.data.frame(roc.dt) else NULL
  ))
}
