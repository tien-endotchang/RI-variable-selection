# This `sim.xy.ext` function extends the original `sim.xy` in `sim.R` function from 
# Hastie et al. (2020) to also include simulation types from Fan & Lv (2008).
#
# Our modifications are:
# 1. Added a 'sim.type' argument to select between HTT and FL data generation.
# 2. Changed the 'Sigma.half' computation from svd to chol to reduce computation time.

sim.xy.ext = function(n, p, nval, rho=0, s=5, beta.type=1, snr=1, sim.type = c("HTT","FL")) {
  # Generate predictors
  x = matrix(rnorm(n*p),n,p)
  xval = matrix(rnorm(nval*p),nval,p)
  
  if (sim.type == "HTT") {
    # Introduce autocorrelation, if needed
    if (rho != 0) {
      inds = 1:p
      Sigma = stats::toeplitz(rho^(0:(p - 1)))
      Sigma.half = chol(Sigma)
      x = x %*% Sigma.half
      xval = xval %*% Sigma.half
    }
    else Sigma = diag(1,p)
    
    # Generate underlying coefficients
    s = min(s,p)
    beta = rep(0,p)
    if (beta.type==1) {
      beta[round(seq(1,p,length=s))] = 1
    } else if (beta.type==2) {
      beta[1:s] = 1
    } else if (beta.type==3) {
      beta[1:s] = seq(10,0.5,length=s)
    } else if (beta.type==4) {
      beta[1:6] = c(-10,-6,-2,2,6,10)
    } else {
      beta[1:s] = 1
      beta[(s+1):p] = 0.5^(1:(p-s))
    }
    
  } else if (sim.type == "FL") {
    if (beta.type==1) {
      Sigma = diag(p)
      Sigma[upper.tri(Sigma)|lower.tri(Sigma)] = rho
      Sigma.half = chol(Sigma)
      x = x %*% Sigma.half
      xval = xval %*% Sigma.half
      beta = rep(0,p)
      beta[1:3] = 5
    } else if (beta.type==2) {
      Sigma = diag(p)
      Sigma[upper.tri(Sigma)|lower.tri(Sigma)] = rho
      Sigma[4, ] = Sigma[, 4] = rho^(1/2)
      Sigma[4, 4] = 1
      Sigma.half = chol(Sigma)
      x = x %*% Sigma.half
      xval = xval %*% Sigma.half
      beta = rep(0,p)
      beta[1:3] = 5
      beta[4] = -15 * rho^(1/2)
    } else if (beta.type==3) {
      Sigma = diag(p)
      Sigma[upper.tri(Sigma)|lower.tri(Sigma)] = rho
      Sigma[4, ] = Sigma[, 4] = rho^(1/2)
      Sigma[4, 4] = 1
      Sigma[5, ] = Sigma[, 5] = 0
      Sigma[5, 5] = 1
      Sigma.half = chol(Sigma)
      x = x %*% Sigma.half
      xval = xval %*% Sigma.half
      beta = rep(0,p)
      beta[1:3] = 5
      beta[4] = -15 * rho^(1/2)
      beta[5] = 1
    }
  }
  
  # Set snr based on sample variance on infinitely large test set
  vmu = as.numeric(t(beta) %*% Sigma %*% beta)
  sigma = sqrt(vmu/snr)
  
  # Generate responses
  y = as.numeric(x %*% beta + rnorm(n)*sigma)
  yval = as.numeric(xval %*% beta + rnorm(nval)*sigma)
  
  enlist(x,y,xval,yval,Sigma,beta,sigma)
}


# Conduct simulation in part I.
#
# The `sim.master.sel` function is a modified version of the `sim.master` 
# function in `sim.R` from the public repository of Hastie et al. (2020). 
#
# Our modifications are:
# 1. Added a 'sim.type' argument to select between HTT and FL data generation.
# 2. Changed the internal call from `sim.xy` to our extended `sim.xy.ext` function.
# 3. The internal evaluation loop is simplified to only compute ranking metrics
#   (e.g., Minimum model size, Proportion of true predictors included).
# 4. Removed codes related to prediction error, relative risk, etc.

sim.master.sel = function(n, p, nval, select.funs, nrep=50, seed=NULL, verbose=FALSE, d=ifelse(n>p, p, 100),
                          file=NULL, file.rep=5, rho=0, s=NULL, beta.type=1, snr=1, sim.type=c("HTT","FL")) {
  
  this.call = match.call()
  if (!is.null(seed)) set.seed(seed)
  
  N = length(select.funs)
  select.names = names(select.funs)
  if (is.null(select.names)) select.names = paste("Method",1:N)
  
  minsize = prop = runtime = vector(mode="list",length=N)
  names(minsize) = names(prop) = names(runtime) = select.names
  for (j in 1:N) {
    minsize[[j]] = runtime[[j]] = matrix(NA,nrep,1)
    prop[[j]] = matrix(NA,nrep,d)
  }
  filled = rep(FALSE,N)
  
  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("  Generating data ...\n")
    }
    
    # Generate x, y, xval, yval
    xy.obj = sim.xy.ext(n,p,nval,rho,s,beta.type,snr,sim.type)
    true.ind = which(xy.obj$beta!=0)
    
    # Loop through the regression methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("  Applying regression method %i (of %i) ...\n",
                    j,N))
      }
      
      tryCatch({
        # Apply the regression method in hand
        runtime[[j]][i] = system.time({
          select.obj = select.funs[[j]](xy.obj$x,xy.obj$y)
        })[1]
        
        # Record all of our metrics
        minsize[[j]][i] = select_k(select.obj, true.ind)
        prop[[j]][i, ] = coverage_ratio(select.obj, true.ind, d)
        
      }, error = function(err) {
        if (verbose) {
          cat(paste("    Oops! Something went wrong, see error message",
                    "below; recording all metrics here as NAs ...\n"))
          cat("    ***** Error message *****\n")
          cat(sprintf("    %s\n",err$message))
          cat("    *** End error message ***\n")
        }
        # N.B. No need to do anything, the metrics are already filled with NAs
      })
    }
    
    # Save intermediate results?
    if (!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
      saveRDS(enlist(minsize,prop,runtime),file=file)
    }
  }
  
  # Save results now (in case of an error that might occur below)
  out = enlist(minsize,prop,runtime)
  if (!is.null(file)) saveRDS(out, file)
  
  # Tune according to validation error, and according to test error
  # Save final results
  out = c(out,list(rho=rho,s=s,beta.type=beta.type,snr=snr,call=this.call))
  class(out) = "sim.select"
  if (!is.null(file)) { saveRDS(out, file); invisible(out) }
  else return(out)
}


# Conduct simulation in part II.
#
# The `sim.master.ext` function is a modified version of the `sim.master` 
# function in `sim.R` from the public repository of Hastie et al. (2020). 
#
# Our modifications are:
# 1. Added a 'sim.type' argument to select between HTT and FL data generation.
# 2. Changed the internal call from `sim.xy` to our extended `sim.xy.ext` function.

sim.master.ext = function(n, p, nval, reg.funs, nrep=50, seed=NULL, verbose=FALSE,
                          file=NULL, file.rep=5, rho=0, s=5, beta.type=1, snr=1, sim.type=c("HTT","FL")) {
  
  this.call = match.call()
  if (!is.null(seed)) set.seed(seed)
  
  N = length(reg.funs)
  reg.names = names(reg.funs)
  if (is.null(reg.names)) reg.names = paste("Method",1:N)
  
  err.train = err.val = err.test = prop = risk = nzs = fpos = fneg = F1 = 
    opt = runtime = vector(mode="list",length=N)
  names(err.train) = names(err.val) = names(err.test) = names(prop) =
    names(risk) = names(nzs) = names(fpos) = names(fneg) = names(F1) = 
    names(opt) = names(runtime) = reg.names
  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] = risk[[j]] =
      nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] = runtime[[j]] =
      matrix(NA,nrep,1)
  }
  filled = rep(FALSE,N)
  err.null = risk.null = sigma = rep(NA,nrep)
  
  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("  Generating data ...\n")
    }
    
    # Generate x, y, xval, yval
    xy.obj = sim.xy.ext(n,p,nval,rho,s,beta.type,snr,sim.type) # Changed sim.xy to sim.xy.ext
    risk.null[i] = diag(t(xy.obj$beta) %*% xy.obj$Sigma %*% xy.obj$beta)
    err.null[i] = risk.null[i] + xy.obj$sigma^2
    sigma[i] = xy.obj$sigma
    
    # Loop through the regression methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("  Applying regression method %i (of %i) ...\n",
                    j,N))
      }
      
      tryCatch({
        # Apply the regression method in hand
        runtime[[j]][i] = system.time({
          reg.obj = reg.funs[[j]](xy.obj$x,xy.obj$y)
        })[1]
        
        # Grab the estimated coefficients, and the predicted values on the
        # training and validation sets
        betahat = as.matrix(coef(reg.obj))
        m = ncol(betahat); nc = nrow(betahat)
        
        # Check for intercept
        if (nc == p+1) {
          intercept = TRUE
          betahat0 = betahat[1,]
          betahat = betahat[-1,]
        }
        else intercept = FALSE
        
        muhat.train = as.matrix(predict(reg.obj,xy.obj$x))
        muhat.val = as.matrix(predict(reg.obj,xy.obj$xval))
        
        # Populate empty matrices for our metrics, of appropriate dimension
        if (!filled[j]) {
          err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] =
            risk[[j]] = nzs[[j]] = fpos[[j]] = fneg[[j]] = F1[[j]] = opt[[j]] =
            matrix(NA,nrep,m)
          filled[j] = TRUE
          # N.B. Filling with NAs is important, because the filled flag could
          # be false for two reasons: i) we are at the first iteration, or ii)
          # we've failed in all previous iters to run the regression method
        }
        
        # Record all of our metrics
        err.train[[j]][i,] = colMeans((muhat.train - xy.obj$y)^2)
        err.val[[j]][i,] = colMeans((muhat.val - xy.obj$yval)^2)
        delta = betahat - xy.obj$beta
        risk[[j]][i,] = diag(t(delta) %*% xy.obj$Sigma %*% delta)
        if (intercept) risk[[j]][i,] = risk[[j]][i,] + betahat0^2
        err.test[[j]][i,] = risk[[j]][i,] + xy.obj$sigma^2
        prop[[j]][i,] = 1 - err.test[[j]][i,] / err.null[i]
        nzs[[j]][i,] = colSums(betahat!=0)
        tpos = colSums((betahat!=0)*(xy.obj$beta!=0))
        fpos[[j]][i,] = nzs[[j]][i,]-tpos
        fneg[[j]][i,] = colSums((betahat==0)*(xy.obj$beta!=0))
        F1[[j]][i,] = 2*tpos/(2*tpos+fpos[[j]][i,]+fneg[[j]][i,])
        opt[[j]][i,] = (err.test[[j]][i,] - err.train[[j]][i,]) /
          err.train[[j]][i,]
      }, error = function(err) {
        if (verbose) {
          cat(paste("    Oops! Something went wrong, see error message",
                    "below; recording all metrics here as NAs ...\n"))
          cat("    ***** Error message *****\n")
          cat(sprintf("    %s\n",err$message))
          cat("    *** End error message ***\n")
        }
        # N.B. No need to do anything, the metrics are already filled with NAs
      })
    }
    
    # Save intermediate results?
    if (!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
      saveRDS(enlist(err.train,err.val,err.test,err.null,prop,risk,risk.null,
                     nzs,fpos,fneg,F1,opt,sigma,runtime),file=file)
    }
  }
  
  # Save results now (in case of an error that might occur below)
  out = enlist(err.train,err.val,err.test,err.null,prop,risk,risk.null,nzs,fpos,
               fneg,F1,opt,sigma,runtime)
  if (!is.null(file)) saveRDS(out, file)
  
  # Tune according to validation error, and according to test error
  out = choose.tuning.params(out)
  
  # Save final results
  out = c(out,list(rho=rho,s=s,beta.type=beta.type,snr=snr,call=this.call))
  class(out) = "sim"
  if (!is.null(file)) { saveRDS(out, file); invisible(out) }
  else return(out)
}