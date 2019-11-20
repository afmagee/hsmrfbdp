source("1_simulation_study/0_setup/summarization/rank_based_convergence_diagnostics.R")

library(coda)

# Get HSRF analyses
hsrf1 <- read.table("2_empirical_analyses/output/HIV_env_HSMRFBDP_1.log",header=TRUE,row.names=1)
hsrf2 <- read.table("2_empirical_analyses/output/HIV_env_HSMRFBDP_2.log",header=TRUE,row.names=1)
hsrf3 <- read.table("2_empirical_analyses/output/HIV_env_HSMRFBDP_3.log",header=TRUE,row.names=1)
hsrf4 <- read.table("2_empirical_analyses/output/HIV_env_HSMRFBDP_4.log",header=TRUE,row.names=1)

# Get GMRF analyses
gmrf1 <- read.table("2_empirical_analyses/output/HIV_env_GMRFBDP_1.log",header=TRUE,row.names=1)
gmrf2 <- read.table("2_empirical_analyses/output/HIV_env_GMRFBDP_2.log",header=TRUE,row.names=1)
gmrf3 <- read.table("2_empirical_analyses/output/HIV_env_GMRFBDP_3.log",header=TRUE,row.names=1)
gmrf4 <- read.table("2_empirical_analyses/output/HIV_env_GMRFBDP_4.log",header=TRUE,row.names=1)

# PSRF
hsrf.convergence <- diagnoseConvergence(list(hsrf1,hsrf2,hsrf4),FALSE)
gmrf.convergence <- diagnoseConvergence(list(gmrf1,gmrf2),FALSE)

# exclude branch rate parameters
hsrf.names <- rownames(hsrf.convergence)
hsrf.is.branchrate <- grepl("branch_rates",rownames(hsrf.convergence))
hsrf.convergence <- hsrf.convergence[!hsrf.is.branchrate]
names(hsrf.convergence) <- hsrf.names[!hsrf.is.branchrate]
hsrf.convergence <- hsrf.convergence[!names(hsrf.convergence) == "Iteration"]
max(hsrf.convergence)

gmrf.names <- rownames(gmrf.convergence)
gmrf.is.branchrate <- grepl("branch_rates",rownames(gmrf.convergence))
gmrf.convergence <- gmrf.convergence[!gmrf.is.branchrate]
names(gmrf.convergence) <- gmrf.names[!gmrf.is.branchrate]
gmrf.convergence <- gmrf.convergence[!names(gmrf.convergence) == "Iteration"]
max(gmrf.convergence)

# ESS
hsrf.ess <- rankESS(list(hsrf1,hsrf2,hsrf4))
gmrf.ess <- rankESS(list(gmrf1,gmrf2))

# exclude branch rate parameters
hsrf.is.branchrate <- grepl("branch_rates",names(hsrf.ess))
hsrf.ess <- hsrf.ess[!hsrf.is.branchrate]
min(hsrf.ess)

gmrf.is.branchrate <- grepl("branch_rates",names(gmrf.ess))
gmrf.ess <- gmrf.ess[!gmrf.is.branchrate]
min(gmrf.ess)
