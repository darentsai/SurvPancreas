#' ID generator for repeated k-fold cross-validation
#' 
#' @param k number of folds
#' @param repeats number of repeats
#' 
kfold_cv_ind <- function(nrow, k = 5, repeats = 1) {
  replicate(repeats, simplify = FALSE, {
    split(sample(nrow), cut(seq_len(nrow), k, labels = paste0("Fold.", 1:k)))
  }) |>
    setNames(paste0("Repeat.", 1:repeats)) |>
    purrr::list_flatten()
}

#' Interpolate a matrix with timestamps to any time horizons
#' 
#' @examples
#' x <- matrix(1 - sort(runif(15)), 3, 5)
#' t <- sort(runif(5, 0, 5))
#' t_interp <- seq(min(t)-1, max(t)+1, len = 200)
#' x_interp1 <- mat_interp(x, t, t_interp, lin_interp = FALSE)
#' x_interp2 <- mat_interp(x, t, t_interp, lin_interp = TRUE)
#' matplot(t, t(x), pch = 1)
#' matplot(t_interp, t(x_interp1), type = 'l', add = TRUE)
#' matplot(t_interp, t(x_interp2), type = 'l', add = TRUE)
#' 
mat_interp <- function(mat, time, interp_time, lin_interp = TRUE) {
  
  if(lin_interp) {
    pos <- approx(x = time, y = seq_along(time), xout = interp_time, rule = 2)$y
    w <- pos %% 1
    interp <- t(t(mat[, floor(pos)]) * (1-w) + t(mat[, ceiling(pos)]) * w)
  } else {
    pos <- outer(interp_time, time, FUN = ">=")
    pos[rowSums(pos) == 0, 1] <- TRUE
    mat_col <- max.col(pos, ties.method = "last")
    interp <- mat[, mat_col]
  }
  return(interp)
}
