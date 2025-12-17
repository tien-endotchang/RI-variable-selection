## Real-world dataset Examples
# --------------------------------------------------------------------------
# SCRIPT SETUP: LOADING FUNCTIONS AND LIBRARIES
#
# This script should be run with the R working directory set to the root
# of the project folder (e.g., '~/RI-variable-selection/').
# setwd("G:/其他電腦/我的筆記型電腦/PhD/Journal Paper/CRI-feature-selection/RI-variable-selection")
# --------------------------------------------------------------------------
rm(list = ls())
# --- 1. Source R files from the external Hastie et al. code ---
# These are the original, unmodified benchmark functions.
message("Loading benchmark functions from external/hastie_code/R/...")
external_files <- list.files("external/hastie_code/R", 
                             pattern = "\\.R$", 
                             full.names = TRUE, 
                             ignore.case = TRUE)
sapply(external_files, source)

# --- 2. Source our own custom and modified R functions ---
# These include our implementations of RI methods and modified simulation runners.
message("Loading custom functions from R/...")
custom_files <- list.files("R", 
                           pattern = "\\.R$", 
                           full.names = TRUE, 
                           ignore.case = TRUE)
sapply(custom_files, source)

# --- 3. Load the compiled C/Fortran code (if necessary) ---
# This step is needed for the original benchmark code that relies on compiled functions.
# Note: This may require the user to have compilation tools installed.
dll_path <- file.path("src", "matrixcomps.dll")
if (file.exists(dll_path)) {
  message(paste("Loading dynamic library:", dll_path))
  dyn.load(dll_path)
} else {
  warning(paste("Dynamic library not found at:", dll_path, 
                "\nSome benchmark functions may not work without compilation."))
}
glmnet::glmnet.control(fdev=0)
message("Setup complete. Starting simulations...")

#------------------------------------------------------------------------------#

run.realdata = function(x, y, reg.funs, nrep=50, seed=NULL, verbose=FALSE,
                        file=NULL, file.rep=5) {
  
  this.call = match.call()
  if (!is.null(seed)) set.seed(seed)
  
  N = length(reg.funs)
  reg.names = names(reg.funs)
  if (is.null(reg.names)) reg.names = paste("Method",1:N)
  
  err.train = err.val = err.test = nzs = acc = bacc = F1 = runtime = vector(mode="list",length=N)
  names(err.train) = names(err.val) = names(err.test) = names(nzs) = 
    names(acc) = names(bacc) = names(F1) = names(runtime) = reg.names
  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = nzs[[j]] = acc[[j]] = 
      bacc[[j]] = F1[[j]] = runtime[[j]] =
      matrix(NA,nrep,1)
  }
  filled = rep(FALSE,N)

  x = as.matrix(x)
  y = as.matrix(y)
  n = nrow(x)
  p = ncol(x)
  n.tv = round(0.8 * n)
  n.tv0 = round(n.tv * sum(y == 0) / length(y))
  n.tv1 = n.tv - n.tv0
  
  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Real-world example %i (of %i) ...\n", i, nrep))
      cat("  Splitting data ...\n")
    }
    
    training.ind = c(sample(which(y == 1), n.tv1), sample(which(y == 0), n.tv0)) 
    xtv = x[training.ind, ]
    ytv = y[training.ind]
    xtest = x[-training.ind,]
    ytest = y[-training.ind]
    
    n.t = round(0.8 * n.tv)
    n.t0 = round(n.t * sum(ytv == 0) / length(ytv))
    n.t1 = n.t - n.t0
    
    train.ind = c(sample(which(ytv == 1), n.t1), sample(which(ytv == 0), n.t0)) 
    xtrain = xtv[train.ind, ]
    ytrain =  ytv[train.ind]
    xval = xtv[-train.ind,]
    yval = ytv[-train.ind]
    
    # Loop through the regression methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("  Applying regression method %i (of %i) ...\n",
                    j,N))
      }
      
      tryCatch({
        # Apply the regression method in hand
        runtime[[j]][i] = system.time({
          reg.obj = reg.funs[[j]](xtrain, ytrain)
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
        
        muhat.train = as.matrix(predict(reg.obj, xtrain))
        muhat.val = as.matrix(predict(reg.obj, xval))
        muhat.test = as.matrix(predict(reg.obj, xtest))
        
        # Populate empty matrices for our metrics, of appropriate dimension
        if (!filled[j]) {
          err.train[[j]] = err.val[[j]] = err.test[[j]] = nzs[[j]] = acc[[j]] = 
            bacc[[j]] = F1[[j]] = matrix(NA,nrep,m)
          filled[j] = TRUE
          # N.B. Filling with NAs is important, because the filled flag could
          # be false for two reasons: i) we are at the first iteration, or ii)
          # we've failed in all previous iters to run the regression method
        }
        
        # Record all of our metrics
        err.train[[j]][i,] = colMeans((muhat.train - ytrain)^2)
        err.val[[j]][i,] = colMeans((muhat.val - yval)^2)
        err.test[[j]][i,] = colMeans((muhat.test - ytest)^2)
        nzs[[j]][i,] = colSums(betahat!=0)
        
        yhat = ifelse(muhat.test > 0.5, 1, 0)
        acc[[j]][i, ] = colMeans(yhat == ytest)
        tpos = colSums(ytest * yhat)
        fpos = colSums(yhat) - tpos
        tneg = colSums((1 - ytest) * (1 - yhat))
        fneg = colSums(1 - yhat) - tneg
        sens = tpos / (tpos + fneg)
        spec = tneg / (fpos + tneg)
        bacc[[j]][i, ] = (sens + spec) / 2
        F1[[j]][i, ] = 2 * tpos / (2 * tpos + fpos + fneg)
        
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
      saveRDS(enlist(err.train,err.val,err.test,nzs,acc,bacc,F1,runtime),file=file)
    }
  }
  
  # Save results now (in case of an error that might occur below)
  out = enlist(err.train,err.val,err.test,nzs,acc,bacc,F1,runtime)
  if (!is.null(file)) saveRDS(out, file)
  
  # Tune according to validation error, and according to test error
  out = choose.tuning.params(out)
  
  # Save final results
  out = c(out,list(call=this.call))
  class(out) = "real"
  if (!is.null(file)) { saveRDS(out, file); invisible(out) }
  else return(out)
}

base.dir = "data/"
datasets = c("leukemia", "gli_85")
foldername = paste0("results/tab/")
dir.create(file.path(foldername), showWarnings = TRUE, recursive=TRUE)
nrep = 100
seed = 42
for(i in 1:length(datasets)){
  cat(sprintf("Real-world dataset example: %s (%i of %i) ...\n", datasets[i], i, length(datasets)))

  x = read.csv (paste0(base.dir, datasets[i], "/", datasets[i], "x.csv"), header = FALSE) 
  y = read.csv (paste0(base.dir, datasets[i], "/", datasets[i], "y.csv"), header = FALSE) 
  
  reg.funs = list()
  reg.funs[["Lasso"]] = function(x,y) lasso(x, y, intercept = T, nlam = 2*nrow(x))
  reg.funs[["Forward stepwise"]] = function(x,y) fs.mod(x, y, intercept = T)
  
  reg.funs[["Relaxed lasso"]] = function(x,y) lasso(x, y, intercept = T, nrelax = 10, nlam = 2*nrow(x))
  
  reg.funs[["LS-CRI"]] = function(x,y) lscri(x, y, intercept = T)
  reg.funs[["LS-CRI.Z"]] = function(x,y) lscriz(x, y, intercept = T)
  reg.funs[["LS-CAR"]] = function(x,y) lscar(x, y, intercept = T)
  reg.funs[["LS-SIS"]] = function(x,y) lsSIS(x, y, intercept = T)
  
  reg.funs[["Ridge-CRI"]] = function(x,y) ridgecri(x, y, intercept = T, nlam = 20)
  reg.funs[["Ridge-CRI.Z"]] = function(x,y) ridgecriz(x, y, intercept = T, nlam = 20)
  reg.funs[["Ridge-CAR"]] = function(x,y) ridgecar(x, y, intercept = T, nlam = 20)
  reg.funs[["Ridge-SIS"]] = function(x,y) ridgeSIS(x, y, intercept = T, nlam = 20)
  
  file = paste0(foldername, "res.", datasets[i], ".rds")
  run.obj = run.realdata(x, y, reg.funs=reg.funs, nrep=nrep, seed=seed, verbose=TRUE, file=file)
  
  tab = matrix(NA, nrow=length(reg.funs), ncol=6)
  row.names(tab) = names(reg.funs)
  metrics = c("err.test", "acc", "bacc", "F1", "nzs", "runtime")
  colnames(tab) = metrics
  
  for(j in 1:length(metrics)){
    if(j < length(metrics)){
      z = run.obj[[metrics[j]]]
      res = tune.and.aggregate(run.obj, z)
      tab[, j] = paste0(round(res$z.val.ave, 2), " (", round(res$z.val.std, 2), ")")
    }else{
      z = run.obj[[metrics[j]]]
      res = matrix(unlist(z), ncol=length(reg.funs))
      tab[, j] = paste0(round(colMeans(res*1000), 2), " (", round(apply(res*1000, 2, sd) / sqrt(nrep), 2), ")")
    }
  }
  write.csv(tab, paste0(foldername, "res.", datasets[i], ".csv")) 
}
