select_k = function(x, y) {
  pos = match(y, x)
  if (any(is.na(pos))) {
    stop("Some elements of y are not found in x.")
  }
  k = max(pos)
  return(k)
}

coverage_ratio <- function(x, y, d) {
  if (d < 0 || d > length(x)) {
    stop("d must be between 0 and length(x).")
  }
  if (length(y) == 0) {
    return(rep(0, d))
  }
  ratios <- vapply(1:d,
                   FUN.VALUE = numeric(1),
                   function(k) sum(x[1:k] %in% y) / length(y))
  return(ratios)
}

load_simulation_data = function(n, p, type.vec, rho.vec, snr.vec, sim.type, 
                                metric=c("minsize","prop"), foldername){
  stem = paste0("sim.n",n,".p",p)
  yvecs = beta.vecs = rho.vecs = snr.vecs = method.vecs = c()
  
  for (beta.type in type.vec) {
    for (rho in rho.vec) {
      name = paste0(stem, ".beta", beta.type, sprintf(".rho%0.2f", rho))
      for (snr in snr.vec) {
        file = paste0(foldername, name, ".snr", round(snr,2), ".rds")
        
        sim.obj = readRDS(file)
        
        if(metric == "minsize"){
          nrep = length(sim.obj$minsize$CRI)
          yvecs = c(yvecs, 
                    unlist(sim.obj$minsize, use.names = F))
          beta.vecs = c(beta.vecs,
                        rep(sim.obj$beta.type, nrep*length(sim.obj$minsize)))
          rho.vecs = c(rho.vecs,
                       rep(sim.obj$rho, nrep*length(sim.obj$minsize)))
          snr.vecs = c(snr.vecs,
                       rep(sim.obj$snr, nrep*length(sim.obj$minsize)))
          method.vecs = c(method.vecs,
                          rep(names(sim.obj$minsize), each=nrep))
        }else{
          yvecs = c(yvecs, 
                    as.vector(sapply(sim.obj$prop, colMeans)))
          beta.vecs = c(beta.vecs,
                        rep(sim.obj$beta.type, p*length(sim.obj$prop)))
          rho.vecs = c(rho.vecs,
                       rep(sim.obj$rho, p*length(sim.obj$prop)))
          snr.vecs = c(snr.vecs,
                       rep(sim.obj$snr, p*length(sim.obj$prop)))
          method.vecs = c(method.vecs,
                          rep(names(sim.obj$prop), each=p))
        }
      }
    }
  }
  beta.vecs = factor(beta.vecs)
  rho.vecs = factor(rho.vecs)
  snr.vecs = factor(round(snr.vecs, 2))
  levels(beta.vecs) = paste("Example", levels(beta.vecs))
  levels(rho.vecs) = paste("Correlation", levels(rho.vecs))
  levels(snr.vecs) = paste("SNR", levels(snr.vecs))
  xvec = rep(1:p, length(type.vec)*length(rho.vec)*length(snr.vec)*length(sim.obj$prop))
  res = data.frame(xvec=xvec, yvec=yvecs, beta=beta.vecs, rho=rho.vecs, snr=snr.vecs, 
                   method=factor(method.vecs, levels=c("GD","CRI","CRI.Z","CAR","SIS")))
  return( res )
}

plot_simulation = function(data, title, n, p, metric){
  if (!require("ggplot2",quietly=TRUE)) {
    stop("Package ggplot2 not installed (required here)!")
  }
  
  if (!require("ggh4x",quietly=TRUE)) {
    stop("Package ggh4x not installed (required here)!")
  }
  cbbPalette = c("GD" = "#00000090",
                 "CRI" = "#FF000090",
                 "CRI.Z" = "#0000FF90",
                 "CAR" = "#3399FF90",
                 "SIS" = "#FF990090"
  )
  
  legend.pos = "bottom"
  
  if(metric=="minsize"){
    gp = ggplot(data = data, aes(x = method, y = yvec, fill = method)) + # x = variable
      geom_boxplot() + ylab(switch(metric, minsize="Minimal Model Size (S)",
                                   prop="Selected Proprotion (Pr(k))")) +
      facet_nested(beta ~ rho + snr) + 
      scale_fill_manual(values = cbbPalette) +
      theme_bw() + theme(legend.position=legend.pos, 
                         axis.text.x=element_blank(), 
                         axis.ticks.x=element_blank(),
                         axis.title.x=element_blank()) +
      ggtitle(title)
    if(p > n){
      gp = gp + scale_y_continuous(trans='log10')
    }
  }else{
    if(p>10) xlim = c(0,50) else xlim = c(0,p)
    gp = ggplot(data = data, aes(x = xvec, y = yvec)) + # x = variable
      geom_line(aes(color = method, linetype = method), linewidth=0.75) +
      facet_nested(beta ~ rho + snr) + 
      scale_color_manual(values = cbbPalette) +
      theme_bw() + theme(legend.position=legend.pos) +
      ggtitle(title) + coord_cartesian(xlim=xlim) +
      xlab("Number of selected predictors k") + ylab(switch(metric, minsize="Minimal Model Size S",
                                                            prop="Selected Proprotion Pr(k)"))
    
    if(p==10){
      gp = gp + scale_x_continuous(breaks=c(1, 5, 10))
    }else{
      gp = gp + scale_x_continuous(breaks=c(1, 25, 50))
    }
  }
  return( gp )
}