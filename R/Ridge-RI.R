# GD
library(relaimpo)
ridgeGD = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000), nlambda = 10,
                   lambda = NULL, intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
  # Check for glmnet package
  if (!require("glmnet",quietly=TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  # Check input data
  check.xy(x=x,y=y)
  
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
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps, 500)
  action = numeric(buf)                      # Actions taken
  df = c() # Degrees of freedom
  beta = matrix(0, p, 1)                     # FS estimates
  
  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  GD = unname(calc.relimp(cbind(y,x))@lmg)
  GD_order = order(GD, decreasing = T)
  
  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  lambda_storage = NULL
  
  while (k <= maxsteps) {
    r = r + 1
    A = c(A, GD_order[r])
    action[k] = GD_order[r]
    
    if(k == 1){
      beta_temp = matrix(0, p, 1)
      if (intercept) {
        beta_temp[A, 1] = lsfit(x0[, A], y0)$coef 
      } else {
        beta_temp[A, 1] = lsfit(x0[, A], y0, int = FALSE)$coef
      }
      df[k + 1] = 1
    } else {
      beta_temp = matrix(0, p, nlambda)
      glmnet.obj = glmnet(x0[, A], y0, alpha = 0, nlambda = nlambda, dfmax = p,
                          lambda.min.ratio = ifelse(nrow(x0[, A]) < ncol(x0[, A]), 0.01, 0.0001), lambda = lambda,
                          intercept = intercept, standardize = TRUE)
      beta_temp[A, 1:nlambda] = as.matrix(glmnet.obj$beta)
      df = c(df, rep(k, nlambda))
      if( is.null(lambda_storage) ) lambda_storage = glmnet.obj$lambda
    }
    
    beta = cbind(beta, beta_temp)
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }
    
    # Update counter
    k = k + 1
  }
  
  # we do not include the OLS
  completepath = FALSE
  bls = NULL
  
  if (verbose) cat("\n")
  
  # Assign column names
  colnames(beta) = as.character(c(0, 1, rep(Seq(2, k-1), each = nlambda)))
  
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,lambda=lambda_storage,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "ridgeGD"
  return(out)
}

coef.ridgeGD = function(object, s, ...) {
  return( object$beta )
}

predict.ridgeGD = function(object, newx, s, ...) {
  beta = coef.ridgeGD(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}




# CRI
ridgecri = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000), nlambda = 10,
                    lambda = NULL, intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
  # Check for glmnet package
  if (!require("glmnet",quietly=TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)

  # Check input data
  check.xy(x=x,y=y)

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

  # Things to keep track of, and return at the end
  buf = min(maxsteps, 500)
  action = numeric(buf)                      # Actions taken
  df = c() # Degrees of freedom
  beta = matrix(0, p, 1)                     # FS estimates

  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  cri_order = order(cri(x, y), decreasing = T)

  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  lambda_storage = NULL
  
  while (k <= maxsteps) {
    r = r + 1
    A = c(A, cri_order[r])
    action[k] = cri_order[r]
    
    if(k == 1){
      beta_temp = matrix(0, p, 1)
      if (intercept) {
        beta_temp[A, 1] = lsfit(x0[, A], y0)$coef 
      } else {
        beta_temp[A, 1] = lsfit(x0[, A], y0, int = FALSE)$coef
      }
      df[k + 1] = 1
    } else {
      beta_temp = matrix(0, p, nlambda)
      glmnet.obj = glmnet(x0[, A], y0, alpha = 0, nlambda = nlambda, dfmax = p,
                          lambda.min.ratio = ifelse(nrow(x0[, A]) < ncol(x0[, A]), 0.01, 0.0001), lambda = lambda,
                          intercept = intercept, standardize = TRUE)
      beta_temp[A, 1:nlambda] = as.matrix(glmnet.obj$beta)
      df = c(df, rep(k, nlambda))
      if( is.null(lambda_storage) ) lambda_storage = glmnet.obj$lambda
    }
    
    beta = cbind(beta, beta_temp)
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }

    # Update counter
    k = k + 1
  }
  
  # we do not include the OLS
  completepath = FALSE
  bls = NULL

  if (verbose) cat("\n")

  # Assign column names
  colnames(beta) = as.character(c(0, 1, rep(Seq(2, k-1), each = nlambda)))

  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,lambda=lambda_storage,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "ridgecri"
  return(out)
}

coef.ridgecri = function(object, s, ...) {
  return( object$beta )
}

predict.ridgecri = function(object, newx, s, ...) {
  beta = coef.ridgecri(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))

  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}



# CRI-Z
ridgecriz = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000), nlambda = 10,
                    lambda = NULL, intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
  # Check for glmnet package
  if (!require("glmnet",quietly=TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  # Check input data
  check.xy(x=x,y=y)
  
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
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps, 500)
  action = numeric(buf)                      # Actions taken
  df = c() # Degrees of freedom
  beta = matrix(0, p, 1)                     # FS estimates
  
  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  criz_order = order(criz(x, y), decreasing = T)
  
  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  lambda_storage = NULL
  
  while (k <= maxsteps) {
    r = r + 1
    A = c(A, criz_order[r])
    action[k] = criz_order[r]
    
    if(k == 1){
      beta_temp = matrix(0, p, 1)
      if (intercept) {
        beta_temp[A, 1] = lsfit(x0[, A], y0)$coef 
      } else {
        beta_temp[A, 1] = lsfit(x0[, A], y0, int = FALSE)$coef
      }
      df[k + 1] = 1
    } else {
      beta_temp = matrix(0, p, nlambda)
      glmnet.obj = glmnet(x0[, A], y0, alpha = 0, nlambda = nlambda, dfmax = p,
                          lambda.min.ratio = ifelse(nrow(x0[, A]) < ncol(x0[, A]), 0.01, 0.0001), lambda = lambda,
                          intercept = intercept, standardize = TRUE)
      beta_temp[A, 1:nlambda] = as.matrix(glmnet.obj$beta)
      df = c(df, rep(k, nlambda))
      if( is.null(lambda_storage) ) lambda_storage = glmnet.obj$lambda
    }
    
    beta = cbind(beta, beta_temp)
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }
    
    # Update counter
    k = k + 1
  }
  
  # we do not include the OLS
  completepath = FALSE
  bls = NULL
  
  if (verbose) cat("\n")
  
  # Assign column names
  colnames(beta) = as.character(c(0, 1, rep(Seq(2, k-1), each = nlambda)))
  
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,lambda=lambda_storage,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "ridgecriz"
  return(out)
}

coef.ridgecriz = function(object, s, ...) {
  return( object$beta )
}

predict.ridgecriz = function(object, newx, s, ...) {
  beta = coef.ridgecriz(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}



# CAR
ridgecar = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000), nlambda = 10,
                     lambda = NULL, intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
  # Check for glmnet package
  if (!require("glmnet",quietly=TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  # Check input data
  check.xy(x=x,y=y)
  
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
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps, 500)
  action = numeric(buf)                      # Actions taken
  df = c() # Degrees of freedom
  beta = matrix(0, p, 1)                     # FS estimates
  
  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  car_order = order(abs(carscore(x, y, verbose = FALSE)), decreasing = T)
  
  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  lambda_storage = NULL
  
  while (k <= maxsteps) {
    r = r + 1
    A = c(A, car_order[r])
    action[k] = car_order[r]
    
    if(k == 1){
      beta_temp = matrix(0, p, 1)
      if (intercept) {
        beta_temp[A, 1] = lsfit(x0[, A], y0)$coef 
      } else {
        beta_temp[A, 1] = lsfit(x0[, A], y0, int = FALSE)$coef
      }
      df[k + 1] = 1
    } else {
      beta_temp = matrix(0, p, nlambda)
      glmnet.obj = glmnet(x0[, A], y0, alpha = 0, nlambda = nlambda, dfmax = p,
                          lambda.min.ratio = ifelse(nrow(x0[, A]) < ncol(x0[, A]), 0.01, 0.0001), lambda = lambda,
                          intercept = intercept, standardize = TRUE)
      beta_temp[A, 1:nlambda] = as.matrix(glmnet.obj$beta)
      df = c(df, rep(k, nlambda))
      if( is.null(lambda_storage) ) lambda_storage = glmnet.obj$lambda
    }
    
    beta = cbind(beta, beta_temp)
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }
    
    # Update counter
    k = k + 1
  }
  
  # we do not include the OLS
  completepath = FALSE
  bls = NULL
  
  if (verbose) cat("\n")
  
  # Assign column names
  colnames(beta) = as.character(c(0, 1, rep(Seq(2, k-1), each = nlambda)))
  
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,lambda=lambda_storage,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "ridgecar"
  return(out)
}

coef.ridgecar = function(object, s, ...) {
  return( object$beta )
}

predict.ridgecar = function(object, newx, s, ...) {
  beta = coef.ridgecar(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}



# SIS
ridgeSIS = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000), nlambda = 10,
                    lambda = NULL, intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
  # Check for glmnet package
  if (!require("glmnet",quietly=TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }
  
  # Set up data
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  # Check input data
  check.xy(x=x,y=y)
  
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
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps, 500)
  action = numeric(buf)                      # Actions taken
  df = c() # Degrees of freedom
  beta = matrix(0, p, 1)                     # FS estimates
  
  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  SIS_order = order(abs(cor(x, y)), decreasing = T)
  
  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  lambda_storage = NULL
  
  while (k <= maxsteps) {
    r = r + 1
    A = c(A, SIS_order[r])
    action[k] = SIS_order[r]
    
    if(k == 1){
      beta_temp = matrix(0, p, 1)
      if (intercept) {
        beta_temp[A, 1] = lsfit(x0[, A], y0)$coef 
      } else {
        beta_temp[A, 1] = lsfit(x0[, A], y0, int = FALSE)$coef
      }
      df[k + 1] = 1
    } else {
      beta_temp = matrix(0, p, nlambda)
      glmnet.obj = glmnet(x0[, A], y0, alpha = 0, nlambda = nlambda, dfmax = p,
                          lambda.min.ratio = ifelse(nrow(x0[, A]) < ncol(x0[, A]), 0.01, 0.0001), lambda = lambda,
                          intercept = intercept, standardize = TRUE)
      beta_temp[A, 1:nlambda] = as.matrix(glmnet.obj$beta)
      df = c(df, rep(k, nlambda))
      if( is.null(lambda_storage) ) lambda_storage = glmnet.obj$lambda
    }
    
    beta = cbind(beta, beta_temp)
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }
    
    # Update counter
    k = k + 1
  }
  
  # we do not include the OLS
  completepath = FALSE
  bls = NULL
  
  if (verbose) cat("\n")
  
  # Assign column names
  colnames(beta) = as.character(c(0, 1, rep(Seq(2, k-1), each = nlambda)))
  
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,lambda=lambda_storage,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "ridgeSIS"
  return(out)
}

coef.ridgeSIS = function(object, s, ...) {
  return( object$beta )
}

predict.ridgeSIS = function(object, newx, s, ...) {
  beta = coef.ridgeSIS(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}



