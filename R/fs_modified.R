# This file contains a modified version of the `fs` function and its helper
# methods from the public repository of Hastie et al. (2020).
#
# MODIFICATION RATIONALE:
# The original code relies on calls to pre-compiled C functions within a
# private package ("bestsubset") that is not publicly available. This modified
# version replaces the C calls in the `updateQR` function with a pure-R
# equivalent to ensure full reproducibility.
#
# All other logic is preserved from the original implementation.

fs.mod = function(x, y, maxsteps=min(nrow(x)-intercept,ncol(x),2000),
              intercept=TRUE, normalize=TRUE, verbose=FALSE) {
  
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
  
  #####
  # Find the first variable to enter and its sign
  z = scale(x,center=F,scale=sqrt(colSums(x^2)))
  u = t(z) %*% y
  j.hit = which.max(abs(u))   # Hitting coordinate
  sign.hit = Sign(u[j.hit])   # Hitting sign
  
  if (verbose) {
    cat(sprintf("1. Added variable %i, |A|=%i...",j.hit,1))
  }
  
  # Now iterate to find the sequence of FS estimates
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps+1,500)
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0,p,buf)     # FS estimates
  
  # Record action, df, solution (df and solution are here
  # correspond to step 0; always a step behind)
  action[1] = j.hit
  df[1] = 0
  beta[,1] = 0
  
  # Other things to keep track of, but not return
  r = 1                       # Size of active set
  A = j.hit                   # Active set
  I = Seq(1,p)[-j.hit]        # Inactive set
  sign = sign.hit             # Active signs
  X1 = x[,j.hit,drop=FALSE]   # Matrix X[,A]
  X2 = x[,-j.hit,drop=FALSE]  # Matrix X[,I]
  k = 2                       # Step counter
  
  # Compute a skinny QR decomposition of X1
  qr.obj = qr(X1)
  Q = qr.Q(qr.obj,complete=TRUE)
  Q1 = Q[,1,drop=FALSE];
  Q2 = Q[,-1,drop=FALSE]
  R = qr.R(qr.obj)
  
  # Throughout the algorithm, we will maintain
  # the decomposition X1 = Q1*R. Dimensions:
  # X1: n x r
  # Q1: n x r
  # Q2: n x (n-r)
  # R:  r x r
  
  while (k <= maxsteps) {
    ##########
    # Check if we've reached the end of the buffer
    if (k > length(action)) {
      buf = length(action)
      action = c(action,numeric(buf))
      df = c(df,numeric(buf))
      beta = cbind(beta,matrix(0,p,buf))
    }
    
    # Key quantities for the next entry
    a = backsolve(R,t(Q1) %*% y)
    b = backsolve(R,t(Q1) %*% X2)
    X2.resid = X2 - X1 %*% b
    z = scale(X2.resid,center=F,scale=sqrt(colSums(X2.resid^2)))
    u = as.numeric(t(z) %*% y)
    
    # Otherwise find the next hitting time
    sign.u = Sign(u)
    abs.u = sign.u * u
    j.hit = which.max(abs.u)
    sign.hit = sign.u[j.hit]
    
    # Record action, df, solution
    action[k] = I[j.hit]
    df[k] = r
    beta[A,k] = a
    
    # Update rest of the variables
    r = r+1
    A = c(A,I[j.hit])
    I = I[-j.hit]
    sign = c(sign,sign.hit)
    X1 = cbind(X1,X2[,j.hit])
    X2 = X2[,-j.hit,drop=FALSE]
    
    # Update the QR decomposition
    updated.qr = updateQR.mod(Q1,Q2,R,X1[,r])
    Q1 = updated.qr$Q1
    Q2 = updated.qr$Q2
    R = updated.qr$R
    
    if (verbose) {
      cat(sprintf("\n%i. Added variable %i, |A|=%i...",k,A[r],r))
    }
    
    # Update counter
    k = k+1
  }
  
  # Record df and solution at last step
  df[k] = k-1
  beta[A,k] = backsolve(R,t(Q1) %*% y)
  
  # Trim
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k)]
  beta = beta[,Seq(1,k),drop=FALSE]
  
  # If we stopped short of the complete path, then note this
  if (k-1 < min(n-intercept,p)) {
    completepath = FALSE
    bls = NULL
  }
  
  # Else we computed the complete path, so record LS solution
  else {
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
  class(out) = "fs"
  return(out)
}

##############################

updateQR.mod = function(Q1,Q2,R,col) {
  m = nrow(Q1)
  n = ncol(Q1)
  k = ncol(Q2)
  
  a = .C("update1",
         Q2=as.double(Q2),
         w=as.double(t(Q2) %*% col),
         m=as.integer(m),
         k=as.integer(k),
         dup=FALSE)
  
  Q2 = matrix(a$Q2,nrow=m)
  w = c(t(Q1) %*% col,a$w)
  
  # Re-structure: delete a column from Q2, add one to
  # Q1, and expand R
  Q1 = cbind(Q1,Q2[,1])
  Q2 = Q2[,-1,drop=FALSE]
  R = rbind(R,rep(0,n))
  R = cbind(R,w[Seq(1,n+1)])
  
  return(list(Q1=Q1,Q2=Q2,R=R))
}
