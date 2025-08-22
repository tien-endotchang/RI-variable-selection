# This function is an extended version of the original `plot.from.file` function
# in `plot.R` in Hastie et al. (2020).
# The following modifications and features have been added:
# 1. Added parameters to allow for fine-grained control over plot appearance, 
#    such as axis limits, line types, and point shapes.
# 2. A standardized color palette is now used to ensure that specific methods 
#    (e.g., "lasso", "LS-CRI-Z") are represented by the same color across all 
#    generated plots for consistency.
# 3. The function now utilizes the `ggh4x` package to create more advanced 
#    grid layouts (e.g., with independent axis scales), 
#    requiring this package as a new dependency.

plot.from.file.mod = function(file.list,
                          row=c("beta","rho","snr"), col=c("rho","beta","snr"),
                          method.nums=NULL, method.names=NULL,
                          what=c("error","risk","prop","F","nonzero"), rel.to=NULL,
                          tuning=c("validation","oracle"), type=c("ave","med"),
                          std=TRUE, lwd=1, pch=19, main=NULL, ylim=NULL,
                          legend.pos=c("bottom","right","top","left","none"),
                          make.pdf=FALSE, fig.dir=".", file.name="sim",
                          w=8, h=10, subset=FALSE) {
  
  # Check for ggplot2 package
  if (!require("ggplot2",quietly=TRUE)) {
    stop("Package ggplot2 not installed (required here)!")
  }
  
  if (!require("ggh4x",quietly=TRUE)) {
    stop("Package ggh4x not installed (required here)!")
  }

  row = match.arg(row)
  col = match.arg(col)
  if (row==col) stop("row and col must be different")
  
  what = match.arg(what)
  tuning = match.arg(tuning)
  type = match.arg(type)
  legend.pos = match.arg(legend.pos)
  
  # Set the method numbers and names
  sim.obj = readRDS(file.list[1])
  if (is.null(method.nums)) method.nums = 1:length(sim.obj$err.test)
  if (is.null(method.names)) method.names =
    names(sim.obj$err.test[method.nums])
  N = length(method.nums)
  
  # Set the base number and name
  if (is.null(rel.to)) {
    base.num = 0
    base.name = ifelse(what=="error","Bayes","null model")
  }
  else {
    base.num = which(method.nums==rel.to)
    base.name = tolower(method.names[base.num])
  }
  
  # Set the y-label
  ylab = switch(what,
                error=paste0("Relative test error (to ",base.name,")"),
                risk=paste0("Relative risk (to ",base.name,")"),
                prop="Proportion of variance explained",
                F="F classification of nonzeros",
                nonzero="Number of nonzeros")
  
  # Collect the y-variable from the file list
  yvec = ybar = beta.vec = rho.vec = snr.vec = c()
  for (i in 1:length(file.list)) {
    sim.obj = readRDS(file.list[i])
    beta.vec = c(beta.vec,rep(sim.obj$beta.type,N))
    rho.vec = c(rho.vec,rep(sim.obj$rho,N))
    snr.vec = c(snr.vec,rep(sim.obj$snr,N))
    
    z = sim.obj[[switch(what,
                        error="err.test",
                        risk="risk",
                        prop="prop",
                        F="F1",
                        nonzero="nzs")]]
    res = tune.and.aggregate(sim.obj, z)
    
    # For prop, F  and nonzero we ignore any request for a relative metric
    if (what=="prop" || what=="F" || what=="nonzero") {
      yvec = c(yvec,res[[paste0("z.",substr(tuning,1,3),".",type)]][method.nums])
      ybar = c(ybar,res[[paste0("z.",substr(tuning,1,3),".",
                                ifelse(type=="ave","std","mad"))]][method.nums])
    }
    
    # For err and risk we respect the request for a relative metric
    else {
      # First build the relative metric
      met = res[[paste0("z.",substr(tuning,1,3))]]#[method.nums]
      if (base.num == 0 && what=="error") denom = sim.obj$sigma^2
      else if (base.num == 0 && what=="risk") denom = sim.obj$risk.null
      else denom = met[[base.num]]
      z.rel = lapply(met, function(v) v / denom)
      # Now aggregate the relative metric
      res2 = tune.and.aggregate(sim.obj, z.rel, tune=FALSE)
      yvec = c(yvec,unlist(res2[[paste0("z.",type)]])[method.nums])
      ybar = c(ybar,unlist(res2[[paste0("z.",ifelse(type=="ave",
                                                    "std","mad"))]])[method.nums])
    }
  }
  # Set the x-variable and x-label
  xvec = snr.vec
  xlab = "Signal-to-noise ratio"
  
  # Set the y-limits
  if (is.null(ylim)) ylim = range(yvec-ybar, yvec+ybar)
  # Produce the plot
  beta.vec = factor(beta.vec + 3)
  rho.vec = factor(rho.vec)
  snr.vec = factor(snr.vec)
  # levels(beta.vec) = paste("Beta-type", levels(beta.vec)) ### HERE ###
  levels(beta.vec) = paste("Example", levels(beta.vec)) ### HERE ###
  levels(rho.vec) = paste("Correlation", levels(rho.vec))
  
  dat = data.frame(x=xvec, y=yvec, se=ybar,
                   beta=beta.vec, rho=rho.vec, snr=snr.vec,
                   Method=factor(rep(method.names, length=length(xvec))))
  
  cbbPalette = c("Best subset" = "#F8766D", 
                 "Forward stepwise" = "#7CAE00",
                 "Lasso" = "#00BFC4", 
                 "Relaxed lasso" = "#C77CFF",
                 "RELAXED lasso" = "#6600CC",
                 "LS-GD" = "#000000",
                 "LS-CRI" = "#FF0000",
                 "LS-CRI.Z" = "#0000FF",
                 "LS-CAR" = "#3399FF",
                 "LS-SIS" = "#FF9900",
                 "Ridge-GD" = "#999999",
                 "Ridge-CRI" = "#FF61CC",
                 "Ridge-CRI.Z" = "#3366FF",
                 "Ridge-CAR" = "#99CCFF",
                 "Ridge-SIS" = "#FFCC00"
  )
  
  gp = ggplot(dat, aes(x=x,y=y,color=Method,linetype=Method,shape=Method)) +
    xlab(xlab) + ylab(ylab) +
    geom_line(lwd=lwd) + geom_point() + 
    theme_bw() + theme(legend.position=legend.pos) +
    scale_colour_manual(values = cbbPalette)
  
  if(subset){
    shape_set = c(rep(pch, length(method.names)-1), 20)
    linetype_set = c(rep(1, length(method.names)-1), 2)
    gp = gp + scale_linetype_manual(values=linetype_set) +
      scale_shape_manual(values=shape_set)
  }else{
    shape_set = rep(pch, length(method.names))
    linetype_set = rep(1, length(method.names))
    gp = gp + scale_linetype_manual(values=linetype_set) +
      scale_shape_manual(values=shape_set)
  }
  
  if(what=="error") gp = gp + facet_grid2(formula(paste(row,"~",col)),
                                          scales = "free_y",
                                          independent = "y") 
  else gp = gp + facet_grid(formula(paste(row,"~",col)))
  
  if (!("snr" %in% c(row,col))) {
    # If SNR is being plotted on the x-axis in each plot, then define special
    # x-axis ticks and put the x-axis on a log scale
    snr.breaks = round(exp(seq(from=min(log(xvec)),
                               to=max(log(xvec)),length=4)),2)
    gp = gp + scale_x_continuous(trans="log", breaks=snr.breaks)
  }
  if (std) gp = gp + geom_errorbar(aes(ymin=y-se,ymax=y+se), width=0.02)
  # if (what=="error") gp = gp + geom_line(aes(x=x, y=1+x), lwd=0.5,
  #                                        linetype=3, color="black")
  if (what=="prop") gp = gp + geom_line(aes(x=x, y=x/(1+x)), lwd=0.5,
                                        linetype=3, color="black")
  if (what =="nonzero") gp = gp + geom_line(aes(x=x, y=sim.obj$s), lwd=0.5,
                                            linetype=3, color="black")
  if (!is.null(main)) gp = gp + ggtitle(main)
  # if (!is.null(ylim)) gp = gp + coord_cartesian(ylim=ylim)
  if (make.pdf) ggsave(sprintf("%s/%s.pdf",fig.dir,file.name),
                       height=h, width=w, device="pdf")
  else gp
}
