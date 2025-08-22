# GD
library(relaimpo)
lsGD = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000),
                intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
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
  
  # Center and scale, etc.
  obj = standardize(x,y,intercept,normalize)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps + 1, 500)
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0, p, buf)     # FS estimates
  
  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  GD = unname(calc.relimp(cbind(y,x))@lmg)
  GD_order = order(GD, decreasing = T)
  
  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  
  while (k <= maxsteps) {
    df[k] = r
    r = r + 1
    A = c(A, GD_order[r])
    action[k] = GD_order[r]
    if (intercept) beta[A, k + 1] = lsfit(x[, A], y)$coef
    else beta[A, k + 1] = lsfit(x[, A], y, int = FALSE)$coef
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }
    
    # Update counter
    k = k + 1
  }
  
  df[k] = k-1
  # Trim
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k)]
  
  # If we stopped short of the complete path, then note this
  if (k-1 < min(n-intercept,p)) {
    completepath = FALSE
    bls = NULL
  } else { # Else we computed the complete path, so record LS solution
    completepath = TRUE
    bls = beta[,k]
  }
  
  if (verbose) cat("\n")
  
  # Adjust for the effect of centering and scaling
  if (intercept) df = df+1
  if (normalize) beta = beta/sx
  if (normalize && completepath) bls = bls/sx
  
  # Assign column names
  colnames(beta) = as.character(Seq(0,k-1))
  
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "lsGD"
  return(out)
}



coef.lsGD = function(object, s, ...) {
  return( object$beta )
}


predict.lsGD = function(object, newx, s, ...) {
  beta = coef.lsGD(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}



# CRI
lscri = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000),
              intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
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
  
  # Center and scale, etc.
  obj = standardize(x,y,intercept,normalize)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps + 1, 500)
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0, p, buf)     # FS estimates
  
  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  cri_order = order(cri(x, y), decreasing = T)
  
  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  
  while (k <= maxsteps) {
    df[k] = r
    r = r + 1
    A = c(A, cri_order[r])
    action[k] = cri_order[r]
    if (intercept) beta[A, k + 1] = lsfit(x[, A], y)$coef
    else beta[A, k + 1] = lsfit(x[, A], y, int = FALSE)$coef
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }
    
    # Update counter
    k = k + 1
  }
  
  df[k] = k-1
  # Trim
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k)]

  # If we stopped short of the complete path, then note this
  if (k-1 < min(n-intercept,p)) {
    completepath = FALSE
    bls = NULL
  } else { # Else we computed the complete path, so record LS solution
    completepath = TRUE
    bls = beta[,k]
  }
  
  if (verbose) cat("\n")
  
  # Adjust for the effect of centering and scaling
  if (intercept) df = df+1
  if (normalize) beta = beta/sx
  if (normalize && completepath) bls = bls/sx
  
  # Assign column names
  colnames(beta) = as.character(Seq(0,k-1))
  
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "lscri"
  return(out)
}

coef.lscri = function(object, s, ...) {
  return( object$beta )
}

predict.lscri = function(object, newx, s, ...) {
  beta = coef.lscri(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}

Standardize = function(x){
  x = x - mean(x)
  x = x / sqrt(sum(x^2))
  return(x)
}
cri = function(x, y){
  x = scale(x, center=TRUE, scale=FALSE)
  norms = sqrt(colSums(x^2))
  x = sweep(x, 2, norms, "/")  
  y = Standardize(y)
  r = Matrix::rankMatrix(x)
  SVD.x = svd(x, nu = r, nv = r)
  cri = (SVD.x$v %*% (SVD.x$d[1:r] * t(SVD.x$v)))^2 %*% (SVD.x$v %*% t(SVD.x$u) %*% y)^2
  cri = as.vector(cri)
  names(cri) = colnames(x)
  return( cri )
}



# CRI-Z
lscriz = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000),
                 intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
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
  
  # Center and scale, etc.
  obj = standardize(x,y,intercept,normalize)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps + 1, 500)
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0, p, buf)     # FS estimates
  
  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  criz_order = order(criz(x, y), decreasing = T)
  
  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  
  while (k <= maxsteps) {
    df[k] = r
    r = r + 1
    A = c(A, criz_order[r])
    action[k] = criz_order[r]
    if (intercept) beta[A, k + 1] = lsfit(x[, A], y)$coef
    else beta[A, k + 1] = lsfit(x[, A], y, int = FALSE)$coef
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }
    
    # Update counter
    k = k + 1
  }
  
  df[k] = k-1
  # Trim
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k)]
  
  # If we stopped short of the complete path, then note this
  if (k-1 < min(n-intercept,p)) {
    completepath = FALSE
    bls = NULL
  } else { # Else we computed the complete path, so record LS solution
    completepath = TRUE
    bls = beta[,k]
  }
  
  if (verbose) cat("\n")
  
  # Adjust for the effect of centering and scaling
  if (intercept) df = df+1
  if (normalize) beta = beta/sx
  if (normalize && completepath) bls = bls/sx
  
  # Assign column names
  colnames(beta) = as.character(Seq(0,k-1))
  
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "lscriz"
  return(out)
}

coef.lscriz = function(object, s, ...) {
  return( object$beta )
}

predict.lscriz = function(object, newx, s, ...) {
  beta = coef.lscriz(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}

criz = function(x, y){
  x = scale(x, center=TRUE, scale=FALSE)
  norms = sqrt(colSums(x^2))
  x = sweep(x, 2, norms, "/")  
  y = Standardize(y)
  r = Matrix::rankMatrix(x)
  SVD.x = svd(x, nu = r, nv = r)
  criz = (SVD.x$v %*% t(SVD.x$u) %*% y)^2
  criz = as.vector(criz)
  names(criz) = colnames(x)
  return( criz )
}



# CAR
library(care)
lscar = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000),
                  intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
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
  
  # Center and scale, etc.
  obj = standardize(x,y,intercept,normalize)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps + 1, 500)
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0, p, buf)     # FS estimates
  
  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  car_order = order(abs(carscore(x, y, verbose = FALSE)), decreasing = T)
  
  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  
  while (k <= maxsteps) {
    df[k] = r
    r = r + 1
    A = c(A, car_order[r])
    action[k] = car_order[r]
    if (intercept) beta[A, k + 1] = lsfit(x[, A], y)$coef
    else beta[A, k + 1] = lsfit(x[, A], y, int = FALSE)$coef
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }
    
    # Update counter
    k = k + 1
  }
  
  df[k] = k-1
  # Trim
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k)]
  
  # If we stopped short of the complete path, then note this
  if (k-1 < min(n-intercept,p)) {
    completepath = FALSE
    bls = NULL
  } else { # Else we computed the complete path, so record LS solution
    completepath = TRUE
    bls = beta[,k]
  }
  
  if (verbose) cat("\n")
  
  # Adjust for the effect of centering and scaling
  if (intercept) df = df+1
  if (normalize) beta = beta/sx
  if (normalize && completepath) bls = bls/sx
  
  # Assign column names
  colnames(beta) = as.character(Seq(0,k-1))
  
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "lscriz"
  return(out)
}

coef.lscar = function(object, s, ...) {
  return( object$beta )
}

predict.lscar = function(object, newx, s, ...) {
  beta = coef.lscar(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}



# SIS 
lsSIS = function(x, y, maxsteps = min(nrow(x) - intercept, ncol(x), 2000),
                 intercept = TRUE, normalize = TRUE, verbose = FALSE) {
  
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
  
  # Center and scale, etc.
  obj = standardize(x,y,intercept,normalize)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps + 1, 500)
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0, p, buf)     # FS estimates
  
  # Record action, df, solution (df and solution are here
  df[1] = 0
  beta[,1] = 0
  SIS_order = order(abs(cor(x, y)), decreasing = T)
  
  # Other things to keep track of, but not return
  r = 0                       # Size of active set
  A = c()                     # Active set
  k = 1                       # Step counter
  
  while (k <= maxsteps) {
    df[k] = r
    r = r + 1
    A = c(A, SIS_order[r])
    action[k] = SIS_order[r]
    if (intercept) beta[A, k + 1] = lsfit(x[, A], y)$coef
    else beta[A, k + 1] = lsfit(x[, A], y, int = FALSE)$coef
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...", k, A[r], r))
    }
    
    # Update counter
    k = k + 1
  }
  
  df[k] = k-1
  # Trim
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k)]
  
  # If we stopped short of the complete path, then note this
  if (k-1 < min(n-intercept,p)) {
    completepath = FALSE
    bls = NULL
  } else { # Else we computed the complete path, so record LS solution
    completepath = TRUE
    bls = beta[,k]
  }
  
  if (verbose) cat("\n")
  
  # Adjust for the effect of centering and scaling
  if (intercept) df = df+1
  if (normalize) beta = beta/sx
  if (normalize && completepath) bls = bls/sx
  
  # Assign column names
  colnames(beta) = as.character(Seq(0,k-1))
  
  out = list(action=action,df=df,beta=beta,completepath=completepath,bls=bls,
             x=x0,y=y0,bx=bx,by=by,intercept=intercept,normalize=normalize)
  class(out) = "lsSIS"
  return(out)
}



coef.lsSIS = function(object, s, ...) {
  return( object$beta )
}


predict.lsSIS = function(object, newx, s, ...) {
  beta = coef.lsSIS(object,s)
  if (missing(newx)) newx = object$x
  else newx = matrix(newx,ncol=ncol(object$x))
  
  newx = scale(newx,object$bx,FALSE)
  if (object$intercept) newx = cbind(rep(1,nrow(newx)),newx)
  return(newx %*% beta)
}
