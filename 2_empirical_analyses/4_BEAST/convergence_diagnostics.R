source("1_simulation_study/0_setup/summarization/rank_based_convergence_diagnostics.R")

library(coda)

# Get beast analyses
beast1 <- read.table("2_empirical_analyses/output/env_subtA_run_1.log",header=TRUE,row.names=1)
beast2 <- read.table("2_empirical_analyses/output/env_subtA_run_2.log",header=TRUE,row.names=1)
beast3 <- read.table("2_empirical_analyses/output/env_subtA_run_3.log",header=TRUE,row.names=1)
beast4 <- read.table("2_empirical_analyses/output/env_subtA_run_4.log",header=TRUE,row.names=1)

# Discard burnin (20%)
beast1 <- beast1[-c(1:5001),]
beast2 <- beast2[-c(1:5001),]
beast3 <- beast3[-c(1:5001),]
beast4 <- beast4[-c(1:5001),]

# PSRF
beast.convergence <- diagnoseConvergence(list(beast1,beast2,beast3,beast4),FALSE)

max(beast.convergence)

# ESS
beast.ess <- rankESS(list(beast1,beast2,beast3,beast4))

min(beast.ess)