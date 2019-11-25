library(TESS)

## Broad information for target tree and number of datasets, we'll use this for all simulations
rho <- 1
tree.age <- 100
ntaxa <- 200
n.sims <- 100

## Get the helper functions we'll need
source("1_simulation_study/0_setup/simulation/helper_functions.R")



## We'll have constant extinction = 0.01 in all simulations
mu <- 0.01
m <- function(x) {
  return(mu)
}

## Simulations of a constant-rate tree
set.seed(47)

# Name simulation directory as below
sim.dir <- "1_simulation_study/1_constant_rate/"

# Make directories if needed
if ( !dir.exists(sim.dir)) {
  dir.create(sim.dir)
}
if ( !dir.exists(paste0(sim.dir,"data"))) {
  dir.create(paste0(sim.dir,"data"))
}
if ( !dir.exists(paste0(sim.dir,"output"))) {
  dir.create(paste0(sim.dir,"output"))
}
if ( !dir.exists(paste0(sim.dir,"summaries"))) {
  dir.create(paste0(sim.dir,"summaries"))
}

# Constant-rate is a 1-fold change under our setup, we can find our speciation rate with the above function
lambda <- findParams(t.shift=50,shift.radius=0,fold.change=1,extinction=mu,target.taxa=ntaxa,sampling.fraction=rho,tree.age=tree.age)

# Speciation-rate function for TESS
l <- function(x) {
  speciationRate(x,t.shift=50,shift.radius=0,rate.old=lambda[1],rate.new=lambda[2])
}

# Simulate and write trees
trees <- tess.sim.age(n.sims,age=tree.age,lambda=l,mu=m,samplingProbability=rho)
for (j in 1:n.sims) {
  write.tree(trees[[j]],paste0(sim.dir,"/data/",j,".tre"))
}
