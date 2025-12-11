## High-dimensional simulation, n=50 and p=1000
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
# --------------------------------------------------------------------------

# Set some overall simulation parameters
n = 50; p = 1000 # Size of training set, and number of predictors
nval = n # Size of validation set
nrep = 30 # Number of repetitions for a given setting
seed = 42 # Random number generator seed
s = 5 # Number of nonzero coefficients
type.vec = c(1:3) # Simulation settings to consider
rho.vec = c(0,0.35,0.7,0.9) # Pairwise predictor correlations
snr.vec = exp(seq(log(0.05),log(6),length=10)) # Signal-to-noise ratios 
stem = paste0("sim.n",n,".p",p)
SIM.TYPE = "HTT"

# Regression functions: lasso, forward stepwise, and best subset selection
reg.funs = list()

reg.funs[["Lasso"]] = function(x,y) lasso(x, y, intercept = FALSE, nlam = 100)
reg.funs[["Forward stepwise"]] = function(x,y) fs.mod(x, y, intercept = FALSE, max = 50)

reg.funs[["Relaxed lasso"]] = function(x,y) lasso(x, y, intercept = FALSE, nrelax = 10, nlam = 100)

reg.funs[["LS-CRI"]] = function(x,y) lscri(x, y, intercept = FALSE, max = 50)
reg.funs[["LS-CRI.Z"]] = function(x,y) lscriz(x, y, intercept = FALSE, max = 50)
reg.funs[["LS-CAR"]] = function(x,y) lscar(x, y, intercept = FALSE, max = 50)
reg.funs[["LS-SIS"]] = function(x,y) lsSIS(x, y, intercept = FALSE, max = 50)

reg.funs[["Ridge-CRI"]] = function(x,y) ridgecri(x, y, intercept = FALSE, max = 50, nlam = 20)
reg.funs[["Ridge-CRI.Z"]] = function(x,y) ridgecriz(x, y, intercept = FALSE, max = 50, nlam = 20)
reg.funs[["Ridge-CAR"]] = function(x,y) ridgecar(x, y, intercept = FALSE, max = 50, nlam = 20)
reg.funs[["Ridge-SIS"]] = function(x,y) ridgeSIS(x, y, intercept = FALSE, max = 50, nlam = 20)

## NOTE: the loop below was not run in serial, it was in fact was split up
## and run on a Linux cluster

foldername = paste0("results/rds/part2/hi50/")
dir.create(file.path(foldername), showWarnings = TRUE, recursive=TRUE)
file.list = c() # Vector of files for the saved rds files
for (beta.type in type.vec) {
  for (rho in rho.vec) {
    name = paste0(stem, ".beta", beta.type, sprintf(".rho%0.2f", rho))
    for (snr in snr.vec) {
      file = paste0(foldername, name, ".snr", round(snr,2), ".rds")
      cat("..... NEW SIMULATION .....\n")
      cat("--------------------------\n")
      cat(paste0("File: ", file, "\n\n"))
      
      sim.master.ext(n, p, nval, reg.funs=reg.funs, nrep=nrep, seed=seed, s=s,
                 verbose=TRUE, file=file, rho=rho, beta.type=beta.type, snr=snr, sim.type=SIM.TYPE)
      
      file.list = c(file.list, file)
      cat("\n")
    }
  }
}

##############################
# Run the code below to reproduce the figures without rerunning the sims
n = 50; p = 1000
fig_dir = sprintf("results/fig/part2/")
dir.create(file.path(fig_dir), showWarnings = TRUE, recursive=TRUE)
file.list = paste0("results/rds/part2/hi50/", list.files("results/rds/part2/hi50/"))

method.nums = c(1:11)
method.names = c("Lasso","Forward stepwise",
                 "Relaxed lasso",
                 "LS-CRI","LS-CRI.Z","LS-CAR","LS-SIS",
                 "Ridge-CRI","Ridge-CRI.Z","Ridge-CAR","Ridge-SIS")
pch = c(rep(19, 3), rep(15, 3), 19, rep(15, 3), 19)

plot.from.file.mod(file.list, what="error", rel.to=NULL, tuning="val",
                   method.nums=method.nums, method.names=method.names,
                   main=paste0("n=",n,", p=",p,", s=",s), pch=pch, make.pdf=TRUE,
                   fig.dir=fig_dir, w=10, h=10, lwd=0.5, subset=FALSE,
                   file.name=paste0("sim.n",n,".p",p,".val.err.rel"))

plot.from.file.mod(file.list, what="F", tuning="val",
                   method.nums=method.nums, method.names=method.names,
                   main=paste0("n=",n,", p=",p,", s=",s), pch=pch, make.pdf=TRUE,
                   fig.dir=fig_dir, w=10, h=10, lwd=0.5, subset=FALSE,
                   file.name=paste0("sim.n",n,".p",p,".val.F"))