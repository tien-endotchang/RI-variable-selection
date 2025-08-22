## Medium-dimensional simulation, n=500 and p=100
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
message("Setup complete. Starting simulations...")
# --------------------------------------------------------------------------

# Set some overall simulation parameters
n = 500; p = 100 # Size of training set, and number of predictors
nval = n # Size of validation set
nrep = 100 # Number of repetitions for a given setting
seed = 42 # Random number generator seed
type.vec = c(1:3) # Simulation settings to consider
rho.vec = c(0.35,0.7,0.9) # Pairwise predictor correlations
snr.vec = exp(seq(log(0.05),log(6),length=4)) # Signal-to-noise ratios 
stem = paste0("sim.n",n,".p",p)
SIM.TYPE = "FL"

# Regression functions: lasso, forward stepwise, and best subset selection
select.funs = list()

select.funs[["CRI"]]   = function(x,y) order(cri(x, y), decreasing = T)
select.funs[["CRI.Z"]] = function(x,y) order(criz(x, y), decreasing = T)
select.funs[["SIS"]]   = function(x,y) order(abs(cor(x, y)), decreasing = T)

## NOTE: the loop below was not run in serial, it was in fact was split up
## and run on a Linux cluster

foldername = paste0("results/rds/part1/med/")
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
      
      sim.master.sel(n, p, nval, select.funs=select.funs, d=p, nrep=nrep, seed=seed, s=NULL,
                           verbose=TRUE, file=file, rho=rho, beta.type=beta.type, snr=snr, sim.type=SIM.TYPE)
      
      file.list = c(file.list, file)
      cat("\n")
    }
  }
}

##############################
# Run the code below to reproduce the figures without rerunning the sims

fig_dir = sprintf("results/fig/part1/")
dir.create(fig_dir, showWarnings = TRUE, recursive=TRUE)
w=12 
h=6  

problem = "med"
n = 500; p = 100 
foldername = paste0("results/rds/part1/",problem,"/")
dat = load_simulation_data(n, p, type.vec, rho.vec, snr.vec, SIM.TYPE, "minsize", foldername)
dat = dat[dat$method != "CAR", ]
gp_low_s = plot_simulation(dat, paste0("n=", n, ", p=", p), n, p, "minsize")
ggsave(sprintf("%s%s.pdf", fig_dir, paste0("sim.", problem, ".S")),
       height=h, width=w, device="pdf")

dat = load_simulation_data(n, p, type.vec, rho.vec, snr.vec, SIM.TYPE, "prop", foldername)
dat = dat[dat$method != "CAR", ]
gp_low_pr = plot_simulation(dat, paste0("n=", n, ", p=", p), n, p, "prop")
ggsave(sprintf("%s%s.pdf", fig_dir, paste0("sim.", problem, ".Pr")),
       height=h, width=w, device="pdf")