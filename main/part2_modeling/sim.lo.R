## Low-dimensional simulation, n=100 and p=10
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
n = 100; p = 10 # Size of training set, and number of predictors
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
reg.funs[["Best subset"]] = function(x,y) bs(x,y,intercept=FALSE,
                                             time.limit=1800,
                                             params=list(Threads=4))
reg.funs[["Lasso"]] = function(x,y) lasso(x, y, intercept = FALSE, nlam = 50)
reg.funs[["Forward stepwise"]] = function(x,y) fs.mod(x, y, intercept = FALSE)

reg.funs[["Relaxed lasso"]] = function(x,y) lasso(x, y, intercept = FALSE, nrelax = 10, nlam = 50)

reg.funs[["LS-GD"]] = function(x,y) lsGD(x, y, intercept = FALSE)
reg.funs[["LS-CRI"]] = function(x,y) lscri(x, y, intercept = FALSE)
reg.funs[["LS-CRI.Z"]] = function(x,y) lscriz(x, y, intercept = FALSE)
reg.funs[["LS-SIS"]] = function(x,y) lsSIS(x, y, intercept = FALSE)

reg.funs[["Ridge-GD"]] = function(x,y) ridgeGD(x, y, intercept = FALSE, nlam = 20)
reg.funs[["Ridge-CRI"]] = function(x,y) ridgecri(x, y, intercept = FALSE, nlam = 20)
reg.funs[["Ridge-CRI.Z"]] = function(x,y) ridgecriz(x, y, intercept = FALSE, nlam = 20)
reg.funs[["Ridge-SIS"]] = function(x,y) ridgeSIS(x, y, intercept = FALSE, nlam = 20)

foldername = paste0("results/rds/part2/lo/")
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

# Subset
n = 100; p = 10
fig_dir = sprintf("results/fig/part2/")
dir.create(file.path(fig_dir), showWarnings = TRUE, recursive=TRUE)
file.list = paste0("results/rds/part2/lo/", list.files("results/rds/part2/lo/"))
remove.list = c()
for(f in file.list){
  f.comp = unlist(strsplit(f, ".", fixed = TRUE))
  if("beta3" %in% f.comp | "00" %in% f.comp | "90" %in% f.comp){
    remove.list = c(remove.list, f)
  }
}
file.list = setdiff(file.list, remove.list)

method.nums = c(1:4, 8, 5:7, 11)
method.names = c("Best subset","Lasso","Forward stepwise",
                 "Relaxed lasso","LS-SIS",
                 "LS-GD","LS-CRI","LS-CRI.Z", "Ridge-CRI.Z")
pch = c(rep(19, 5), rep(15, 4))

plot.from.file.mod(file.list, what="error", rel.to=NULL, tuning="val",
                   method.nums=method.nums, method.names=method.names,
                   main=paste0("n=",n,", p=",p,", s=",s), pch=pch, make.pdf=TRUE,
                   fig.dir=fig_dir, w=6.5, h=7, lwd=0.75, subset=TRUE,
                   file.name=paste0("sim.n",n,".p",p,".val.err.rel.sub"))

plot.from.file.mod(file.list, what="F", tuning="val",
                   method.nums=method.nums, method.names=method.names,
                   main=paste0("n=",n,", p=",p,", s=",s), pch=pch, make.pdf=TRUE,
                   fig.dir=fig_dir, w=6.5, h=7, lwd=0.75, subset=TRUE,
                   file.name=paste0("sim.n",n,".p",p,".val.F.sub"))

# Full
n = 100; p = 10
fig_dir = sprintf("results/fig/part2/")
dir.create(file.path(fig_dir), showWarnings = TRUE, recursive=TRUE)
file.list = paste0("results/rds/part2/lo/", list.files("results/rds/part2/lo/"))

method.nums = c(1:4, 8, 12, 5:7, 9:11)
method.names = c("Best subset","Lasso","Forward stepwise",
                 "Relaxed lasso","LS-SIS","Ridge-SIS",
                 "LS-GD","LS-CRI","LS-CRI.Z",
                 "Ridge-GD","Ridge-CRI","Ridge-CRI.Z")
pch = c(rep(19, 6), rep(15, 6))

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