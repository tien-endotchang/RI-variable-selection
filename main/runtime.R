## Compute runtime tables
# --------------------------------------------------------------------------
# This script should be run with the R working directory set to the root
# of the project folder (e.g., '~/RI-variable-selection/').
# setwd("G:/其他電腦/我的筆記型電腦/PhD/Journal Paper/CRI-feature-selection/RI-variable-selection")
# --------------------------------------------------------------------------
rm(list=ls())

base_dir = "G:/其他電腦/我的筆記型電腦/PhD/Journal Paper/CRI-feature-selection/RI-variable-selection_backup/results/rds"  
output_dir = "results/tab"
parts = c("part1", "part2")
levels = c("lo", "med", "hi50", "hi100")
methods = list(c("GD", "CRI", "CRI.Z", "CAR", "SIS"),
               c("Best subset", "Lasso", "Forward stepwise", "Relaxed lasso", 
                 "LS-GD", "LS-CRI", "LS-CRI.Z", "LS-CAR", "LS-SIS",
                 "Ridge-GD", "Ridge-CRI", "Ridge-CRI.Z", "Ridge-CAR", "Ridge-SIS"))
for(i in 1:length(parts)){
  res = data.frame(matrix(0, length(levels), length(methods[[i]])))
  colnames(res) = methods[[i]]
  row.names(res) = levels
  
  for(j in 1:length(levels)){
    
    folder_path = file.path(base_dir, parts[i], levels[j])
    rds_files = list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)
    
    count = 0
    for(k in 1:length(rds_files)){
      obj = readRDS(rds_files[k])
      count = count + colMeans(matrix(unlist(obj$runtime), ncol=length(obj$runtime))) / length(rds_files)
    }
    res[j, names(obj$runtime)] = count
  }
  
  write.csv(res, paste0(output_dir, "/", parts[i], ".csv"))
}
