## Degrees of freedom and risk simulations
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
library(ggplot2)
# --------------------------------------------------------------------------
## Degrees of freedom simulation
# Note this is a "hand-made" simulation (rather than using built-in functions
# like above) because we want to keep x fixed for the df calculations
# Set some overall simulation parameters
n = 70; p = 30 # Size of training set, and number of predictors
nval = n # Size of validation set
nrep = 500
seed = 0 # Random number generator seed
s = 5 # Number of nonzero coefficients
beta.type = 2 # Coefficient type

set.seed(seed)
xy.obj = sim.xy.ext(n,p,nval,rho=0.35,s=s,beta.type=beta.type,snr=0.7,sim.type="HTT")
x = xy.obj$x
y = xy.obj$y
mu = as.numeric(x %*% xy.obj$beta)
sigma = xy.obj$sigma
nlam = 300
nrel = 9

ip.las = matrix(0,nrep,(p+1)*nrel)
ip.fs = ip.bs = ip.lscri = ip.lscriz = ip.lssis = matrix(0,nrep,p+1)

for (r in 1:nrep) {
  cat(r,"... \n")
  eps = rnorm(n)*sigma
  y = mu + eps
  beta.las = coef(lasso(x,y,intercept=FALSE,nlam=nlam,nrel=nrel))
  nzs.las = colSums(beta.las != 0)[0:(nlam-1)*nrel+1]
  j = nlam - rev(match(p:0, rev(nzs.las), NA)) # choose non-zeros with largest lambda
  ind = rep(j*nrel, each=nrel) + rep(1:nrel, p+1)
  yhat.las = (x %*% beta.las)[,ind]
  ip.las[r,] = colSums(yhat.las * eps)
  
  yhat.fs = predict(fs.mod(x,y,intercept=FALSE))[, 1:(p + 1)]
  yhat.bs = predict(bs(x,y,intercept=FALSE))[, 1:(p + 1)]
  yhat.lscri = predict(lscri(x,y,intercept=FALSE))[, 1:(p + 1)]
  yhat.lscriz = predict(lscriz(x,y,intercept=FALSE))[, 1:(p + 1)]
  yhat.lssis = predict(lsSIS(x,y,intercept=FALSE))[, 1:(p + 1)]
  
  ip.fs[r,] = colSums(yhat.fs * eps)
  ip.bs[r,] = colSums(yhat.bs * eps)
  ip.lscri[r,] = colSums(yhat.lscri * eps)
  ip.lscriz[r,] = colSums(yhat.lscriz * eps)
  ip.lssis[r,] = colSums(yhat.lssis * eps)
}

df.las = colMeans(ip.las, na.rm=TRUE) / sigma^2
df.las = matrix(df.las, p+1, nrel, byrow=TRUE)
df.fs = colMeans(ip.fs, na.rm=TRUE) / sigma^2
df.bs = colMeans(ip.bs, na.rm=TRUE) / sigma^2
df.lscri = colMeans(ip.lscri, na.rm=TRUE) / sigma^2
df.lscriz = colMeans(ip.lscriz, na.rm=TRUE) / sigma^2
df.lssis = colMeans(ip.lssis, na.rm=TRUE) / sigma^2

# save(list=ls(),file="sim.df.rda")

##############################
# Run the code below to reproduce the df figure without rerunning the sims
# library(bestsubset)
# load("rds/sim.df.rda")

# Plot the results
dat = data.frame(x=rep(0:p,4),
                 y=c(df.bs,df.fs,df.las[,1],df.lscri),
                 Method=factor(rep(c("Best subset","Forward stepwise","Lasso","LS-CRI"),
                                    rep(p+1,4))))

ggplot(dat, aes(x=x,y=y,color=Method)) +
  xlab("Number of nonzero coefficients") +
  ylab("Degrees of freedom") +
  geom_line(lwd=0.5, color="black", linetype=3, aes(x,x)) +
  geom_line(lwd=1) + geom_point(pch=19) +
    theme_bw() + theme(legend.just=c(1,0), legend.pos=c(0.95,0.05))
ggsave("results/fig/df1.pdf", height=6, width=6.4, device="pdf")

cbbPalette = c("Best subset" = "#F8766D75", 
               "Forward stepwise" = "#7CAE0075",
               "Lasso" = "#00BFC475", 
               "Relaxed lasso: 0.5" = "#C77CFF75",
               "Relaxed lasso: 0" = "#6600CC75",
               "LS-GD" = "#000000",
               "LS-CRI" = "#FF0000",
               "LS-CRI.Z" = "#0000FF",
               "LS-CAR" = "#3399FF",
               "LS-SIS" = "#FF990075",
               "Ridge-GD" = "#999999",
               "Ridge-CRI" = "#FF61CC",
               "Ridge-CRI.Z" = "#3366FF",
               "Ridge-CAR" = "#99CCFF",
               "Ridge-SIS" = "#FFCC0075"
)

methods.names = c("Best subset","Forward stepwise","Lasso",
                  "Relaxed lasso: 0.5","Relaxed lasso: 0","LS-SIS",
                  "LS-CRI", "LS-CRI.Z")
dat = data.frame(x=rep(0:p,8),
                 y=c(df.bs,df.fs,df.las[,1],df.las[,5],df.las[,9],df.lssis,
                     df.lscri,df.lscriz),
                 Method=factor(rep(methods.names,
                                    rep(p+1,8)), levels=methods.names))

shape_set = c(rep(19, 6), rep(15, 2))
ggplot(dat, aes(x=x,y=y,color=Method,shape=Method)) +
  xlab("Number of nonzero coefficients") +
  ylab("Degrees of freedom") +
  geom_line(lwd=0.5, color="black", linetype=3, aes(x,x)) +
  geom_line(lwd=1) + geom_point() +
  theme_bw() + theme(legend.just=c(1,0), legend.pos=c(0.95,0.05)) +
  scale_colour_manual(values = cbbPalette) +
  scale_shape_manual(values=shape_set)

ggsave("results/fig/df2.pdf", height=4, width=6, device="pdf")

