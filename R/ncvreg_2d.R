ncvreg_2d = function(x, y, penalty = c("SCAD", "MCP"),
                     ngamma = 9, max.gamma = 150, min.gamma = NULL,
                     nlambda = 50, intercept = TRUE, normalize = TRUE) {
  # Check for ncvreg package
  if (!require("ncvreg",quietly=TRUE)) {
    stop("Package ncvreg not installed (required here)!")
  }
  
  penalty = match.arg(penalty)
  
  # Set min.gamma based on penalty constraints
  if (is.null(min.gamma)) {
    if (penalty == "SCAD") {
      min.gamma = 2.001
    } else {
      min.gamma = 1.001
    }
  }
  
  # Build gamma grid like sparsenet
  gamma_grid = exp(seq(from = log(max.gamma),
                       to = log(min.gamma),
                       length.out = ngamma))
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  # Check input data
  # check.xy(x=x,y=y)
  
  # Save original x and y
  x0 = x
  y0 = y
  
  # # Center and scale, etc.
  obj = standardize(x, y, intercept, normalize)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx
  
  family = "gaussian"

  
  # Storage: p rows, ngamma * nlambda columns
  beta = matrix(0, p, ngamma * nlambda)
  beta0 = numeric(ngamma * nlambda)
  gamma_vec = numeric(ngamma * nlambda)
  lambda_vec = numeric(ngamma * nlambda)
  
  col_idx = 1
  
  # Fit for each gamma
  for (i in seq_along(gamma_grid)) {
    g = gamma_grid[i]
    
    fit = ncvreg(x0, y0, family = family, penalty = penalty,
                 gamma = g, nlambda = nlambda)
    
    nlam_actual = length(fit$lambda)
    
    # Fill in columns for this gamma
    for (j in 1:nlam_actual) {
      beta0[col_idx] = fit$beta[1, j]
      beta[, col_idx] = fit$beta[-1, j]
      gamma_vec[col_idx] = g
      lambda_vec[col_idx] = fit$lambda[j]
      col_idx = col_idx + 1
    }
  }
  
  # Trim if fewer lambdas were used
  if (col_idx - 1 < ngamma * nlambda) {
    beta = beta[, 1:(col_idx - 1), drop = FALSE]
    beta0 = beta0[1:(col_idx - 1)]
    gamma_vec = gamma_vec[1:(col_idx - 1)]
    lambda_vec = lambda_vec[1:(col_idx - 1)]
  }
  
  # Column names
  colnames(beta) = paste0("g", rep(1:ngamma, each = nlambda)[1:ncol(beta)],
                          "_l", rep(1:nlambda, ngamma)[1:ncol(beta)])
  
  out = list(
    beta = beta,
    beta0 = beta0,
    gamma = gamma_vec,
    lambda = lambda_vec,
    gamma_grid = gamma_grid,
    x = x0,
    y = y0,
    bx = bx,
    by = by,
    sx = sx,
    intercept = intercept,
    normalize = normalize,
    family = family,
    penalty = penalty
  )
  class(out) = "ncvreg_2d"
  return(out)
}

coef.ncvreg_2d = function(object, s, ...) {
  beta = object$beta
  if (object$intercept) {
    beta = rbind("(Intercept)" = object$beta0, beta)
  }
  return(beta)
}

predict.ncvreg_2d = function(object, newx, s, ...) {
  beta = coef.ncvreg_2d(object, s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx, ncol = ncol(object$x))
  
  if (object$intercept) newx = cbind(rep(1, nrow(newx)), newx)
  
  eta = newx %*% beta
  return(eta)
}
