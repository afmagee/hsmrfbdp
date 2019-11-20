source("1_simulation_study/0_setup/summarization/rank_based_convergence_diagnostics.R")

library(coda)

grid.sizes <- c(10,20,50,100,200)

hsrf.convergence <- vector("list",length(grid.sizes))
hsrf.ess <- vector("list",length(grid.sizes))

for (i in 1:length(grid.sizes)) {
  n <- grid.sizes[i]
  # Get HSRF analyses and combine
  hsrf1 <- read.table(paste0("2_empirical_analyses/output/HSMRFBDP_grid_size_",n,"_run_1.log"),header=TRUE)
  hsrf2 <- read.table(paste0("2_empirical_analyses/output/HSMRFBDP_grid_size_",n,"_run_2.log"),header=TRUE)
  
  # PSRF
  hsrf.convergence[[i]] <- diagnoseConvergence(list(hsrf1,hsrf2),FALSE)[-1,]
  hsrf.ess[[i]] <- rankESS(list(hsrf1,hsrf2))[-1]
}


gmrf.convergence <- vector("list",length(grid.sizes))
gmrf.ess <- vector("list",length(grid.sizes))

for (i in 1:length(grid.sizes)) {
  n <- grid.sizes[i]
  # Get gmrf analyses and combine
  gmrf1 <- read.table(paste0("2_empirical_analyses/output/GMRFBDP_grid_size_",n,"_run_1.log"),header=TRUE)
  gmrf2 <- read.table(paste0("2_empirical_analyses/output/GMRFBDP_grid_size_",n,"_run_2.log"),header=TRUE)
  
  # PSRF
  gmrf.convergence[[i]] <- diagnoseConvergence(list(gmrf1,gmrf2),FALSE)[-1,]
  gmrf.ess[[i]] <- rankESS(list(gmrf1,gmrf2))[-1]
}

# Look at worst rank-PSRF from each analysis
max(hsrf.convergence[[1]])
max(hsrf.convergence[[2]])
max(hsrf.convergence[[3]])
max(hsrf.convergence[[4]])
max(hsrf.convergence[[5]])

max(gmrf.convergence[[1]])
max(gmrf.convergence[[2]])
max(gmrf.convergence[[3]])
max(gmrf.convergence[[4]])
max(gmrf.convergence[[5]])

# Look at worst rank-ESS from each analysis
min(hsrf.ess[[1]])
min(hsrf.ess[[2]])
min(hsrf.ess[[3]])
min(hsrf.ess[[4]])
min(hsrf.ess[[5]])

min(gmrf.ess[[1]])
min(gmrf.ess[[2]])
min(gmrf.ess[[3]])
min(gmrf.ess[[4]])
min(gmrf.ess[[5]])
