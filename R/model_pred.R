cox_pred <- function(object, newdata, pred_horiz) {
  
  # Not used:
  # pred <- t(
  #   summary(survfit(object$fit, newdata = newdata),
  #           times = pred_horiz)$surv
  # )
  #
  # 1. Slow!
  # 2. survfit() is not able to create a curve for models that
  #    contain an interaction without the lower order effect
  
  coef <- object$fit$coefficients
  object$fit$coefficients <- replace(coef, is.na(coef), 0)
  pred <- sapply(pred_horiz, function(t) {
    newdata[[all.vars(formula(object$fit)[[2]])[1]]] <- t
    predict(object$fit, newdata, type = "survival")
  })
  
  return(pred)
}

grpreg_cox_pred <- function(object, newdata, pred_horiz, lambda.opt = 2) {
  
  lambda <- paste0("lambda.", c("min", "1se"))[lambda.opt]
  S.lst <- predict(object$fit$fit,
                   model.matrix(object$formula, newdata)[, -1], # exclude intercept
                   type = "survival",
                   lambda = object$fit[[lambda]])
  
  pred <- t(sapply(S.lst, \(S) S(pred_horiz)))
  return(pred)
}

orsf_pred <- function(object, newdata, pred_horiz, type = "surv") {
  
  pred <- predict(object$fit,
                  new_data = newdata,
                  pred_type = type,
                  pred_horizon = pred_horiz)
  return(pred)
}

xgbst_cox_pred <- function(object, newdata, pred_horiz, smooth = TRUE) {
  
  # The order of the columns must match with that of the data from which the model was fitted
  # (i.e. columns will not be referenced by their names, just by their order in the data)
  xdat.new <- newdata[getinfo(object$fit, "feature_name")]
  lin_preds_pred <- predict(object$fit, newdata = xdat.new, outputmargin = TRUE)
  
  # baseline hazard H0(t) estimate at pred_horiz
  base_haz <- object$H0(pred_horiz, smooth = smooth)
  pred <- outer(lin_preds_pred, base_haz,
                FUN = \(x, H0) exp(exp(x) * -H0))
  return(pred)
}

nn_surv_pred <- function(object, newdata, pred_horiz, lin_interp = TRUE) {
  
  # {survivalmodels} gives predictions only at the observed event time
  # The predicted probabilities are linearly interpolated to any prediction horizons
  pred_alltime <- predict(object$fit, newdata, type = "survival")
  pmat_time <- as.numeric(colnames(pred_alltime))
  pmat <- unname(pred_alltime)
  pred <- mat_interp(pmat, pmat_time, pred_horiz, lin_interp = lin_interp)
  return(pred)
}
