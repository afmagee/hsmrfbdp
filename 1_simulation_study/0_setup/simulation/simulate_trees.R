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

# Get and write tree sizes
tree.sizes <- unlist(lapply(trees,function(phy){length(phy$tip.label)}))

cat(tree.sizes,sep="\n",file=paste0(sim.dir,"tree_sizes.txt"))
cat(quantile(tree.sizes,c(0,0.5,1)),"\n")



## Part 2, simulations of a 2-fold change

# Simulation directories are named with a numeric prefix (1:10) for convenience, 
#   then by the proportion of the tree the shift takes, 
#   then by the center of the shift as a proportion of the age of the tree
# So 3_0.5_0.5 is a shift that takes 1/2 the age of the tree, centered halfway through the age of the tree
# 1_0.0_0.5 thus fits into both the spectrum of shift duration and the spectrum of shift placement
sim.dirs <- c("1_simulation_study/2_two_fold/1_0.0_0.5/",
              "1_simulation_study/2_two_fold/2_0.25_0.5/",
              "1_simulation_study/2_two_fold/3_0.5_0.5/",
              "1_simulation_study/2_two_fold/4_0.75_0.5/",
              "1_simulation_study/2_two_fold/5_1.0_0.5/",
              "1_simulation_study/2_two_fold/6_0.0_0.4/",
              "1_simulation_study/2_two_fold/7_0.0_0.6/",
              "1_simulation_study/2_two_fold/8_0.0_0.7/",
              "1_simulation_study/2_two_fold/9_0.0_0.8/",
              "1_simulation_study/2_two_fold/10_0.0_0.9/")

# First we simulate 5 radius values that allow us to span from instant change to constant change
shift.radii <- c(0,12.5,25,37.5,50,
                 0,0,0,0,0)

shift.times <- c(0.5,0.5,0.5,0.5,0.5,
                 0.4,0.6,0.7,0.8,0.9) * 100

# Find parameters
piecewise.params <- vector("list",10)
piecewise.params[[1]]  <- findParams(t.shift=shift.times[1], shift.radius=shift.radii[1], fold.change=0.5,target.taxa=ntaxa,extinction=mu)
piecewise.params[[2]]  <- findParams(t.shift=shift.times[2], shift.radius=shift.radii[2], fold.change=0.5,target.taxa=ntaxa,extinction=mu)
piecewise.params[[3]]  <- findParams(t.shift=shift.times[3], shift.radius=shift.radii[3], fold.change=0.5,target.taxa=ntaxa,extinction=mu)
piecewise.params[[4]]  <- findParams(t.shift=shift.times[4], shift.radius=shift.radii[4], fold.change=0.5,target.taxa=ntaxa,extinction=mu)
piecewise.params[[5]]  <- findParams(t.shift=shift.times[5], shift.radius=shift.radii[5], fold.change=0.5,target.taxa=ntaxa,extinction=mu)
piecewise.params[[6]]  <- findParams(t.shift=shift.times[6], shift.radius=shift.radii[6], fold.change=0.5,target.taxa=ntaxa,extinction=mu)
piecewise.params[[7]]  <- findParams(t.shift=shift.times[7], shift.radius=shift.radii[7], fold.change=0.5,target.taxa=ntaxa,extinction=mu)
piecewise.params[[8]]  <- findParams(t.shift=shift.times[8], shift.radius=shift.radii[8], fold.change=0.5,target.taxa=ntaxa,extinction=mu)
piecewise.params[[9]]  <- findParams(t.shift=shift.times[9], shift.radius=shift.radii[9], fold.change=0.5,target.taxa=ntaxa,extinction=mu)
piecewise.params[[10]] <- findParams(t.shift=shift.times[10],shift.radius=shift.radii[10],fold.change=0.5,target.taxa=ntaxa,extinction=mu)


# Simulate and write trees, record tree sizes for posterity

set.seed(8472)

for (i in 1:10) {
  # Make directories if needed
  if ( !dir.exists(paste0(sim.dirs[i],"data"))) {
    dir.create(paste0(sim.dirs[i],"data"))
  }
  if ( !dir.exists(paste0(sim.dirs[i],"output"))) {
    dir.create(paste0(sim.dirs[i],"output"))
  }
  if ( !dir.exists(paste0(sim.dirs[i],"summaries"))) {
    dir.create(paste0(sim.dirs[i],"summaries"))
  }
  
  # Speciation-rate function for TESS
  l <- function(x) {
    speciationRate(x,t.shift=shift.times[i],shift.radius=shift.radii[i],rate.old=piecewise.params[[i]][1],rate.new=piecewise.params[[i]][2])
  }
  
  # Simulate and write trees
  trees <- tess.sim.age(n.sims,age=tree.age,lambda=l,mu=m,samplingProbability=rho)
  for (j in 1:n.sims) {
    write.tree(trees[[j]],paste0(sim.dirs[i],"/data/",j,".tre"))
  }
  
  # Get and write tree sizes
  tree.sizes <- unlist(lapply(trees,function(phy){length(phy$tip.label)}))
  
  cat(tree.sizes,sep="\n",file=paste0(sim.dirs[i],"tree_sizes.txt"))
  cat(quantile(tree.sizes,c(0,0.5,1)),"\n")
  
}

## Part 3, simulations of a 4-fold change

sim.dirs <- c("1_simulation_study/3_four_fold/1_0.0_0.5/",
              "1_simulation_study/3_four_fold/2_0.25_0.5/",
              "1_simulation_study/3_four_fold/3_0.5_0.5/",
              "1_simulation_study/3_four_fold/4_0.75_0.5/",
              "1_simulation_study/3_four_fold/5_1.0_0.5/",
              "1_simulation_study/3_four_fold/6_0.0_0.4/",
              "1_simulation_study/3_four_fold/7_0.0_0.6/",
              "1_simulation_study/3_four_fold/8_0.0_0.7/",
              "1_simulation_study/3_four_fold/9_0.0_0.8/",
              "1_simulation_study/3_four_fold/10_0.0_0.9/")

# First we simulate 5 radius values that allow us to span from instant change to constant change
shift.radii <- c(0,12.5,25,37.5,50,
                 0,0,0,0,0)

shift.times <- c(0.5,0.5,0.5,0.5,0.5,
                 0.4,0.6,0.7,0.8,0.9) * 100

# Find parameters
piecewise.params <- vector("list",10)
piecewise.params[[1]]  <- findParams(t.shift=shift.times[1], shift.radius=shift.radii[1], fold.change=0.25,target.taxa=ntaxa,extinction=mu)
piecewise.params[[2]]  <- findParams(t.shift=shift.times[2], shift.radius=shift.radii[2], fold.change=0.25,target.taxa=ntaxa,extinction=mu)
piecewise.params[[3]]  <- findParams(t.shift=shift.times[3], shift.radius=shift.radii[3], fold.change=0.25,target.taxa=ntaxa,extinction=mu)
piecewise.params[[4]]  <- findParams(t.shift=shift.times[4], shift.radius=shift.radii[4], fold.change=0.25,target.taxa=ntaxa,extinction=mu)
piecewise.params[[5]]  <- findParams(t.shift=shift.times[5], shift.radius=shift.radii[5], fold.change=0.25,target.taxa=ntaxa,extinction=mu)
piecewise.params[[6]]  <- findParams(t.shift=shift.times[6], shift.radius=shift.radii[6], fold.change=0.25,target.taxa=ntaxa,extinction=mu)
piecewise.params[[7]]  <- findParams(t.shift=shift.times[7], shift.radius=shift.radii[7], fold.change=0.25,target.taxa=ntaxa,extinction=mu)
piecewise.params[[8]]  <- findParams(t.shift=shift.times[8], shift.radius=shift.radii[8], fold.change=0.25,target.taxa=ntaxa,extinction=mu)
piecewise.params[[9]]  <- findParams(t.shift=shift.times[9], shift.radius=shift.radii[9], fold.change=0.25,target.taxa=ntaxa,extinction=mu)
piecewise.params[[10]] <- findParams(t.shift=shift.times[10],shift.radius=shift.radii[10],fold.change=0.25,target.taxa=ntaxa,extinction=mu)


# Simulate and write trees, record tree sizes for posterity

set.seed(42)

for (i in 1:10) {
  # Make directories if needed
  if ( !dir.exists(paste0(sim.dirs[i],"data"))) {
    dir.create(paste0(sim.dirs[i],"data"))
  }
  if ( !dir.exists(paste0(sim.dirs[i],"output"))) {
    dir.create(paste0(sim.dirs[i],"output"))
  }
  if ( !dir.exists(paste0(sim.dirs[i],"summaries"))) {
    dir.create(paste0(sim.dirs[i],"summaries"))
  }
  
  # Speciation-rate function for TESS
  l <- function(x) {
    speciationRate(x,t.shift=shift.times[i],shift.radius=shift.radii[i],rate.old=piecewise.params[[i]][1],rate.new=piecewise.params[[i]][2])
  }
  
  # Simulate and write trees
  trees <- tess.sim.age(n.sims,age=tree.age,lambda=l,mu=m,samplingProbability=rho)
  for (j in 1:n.sims) {
    write.tree(trees[[j]],paste0(sim.dirs[i],"/data/",j,".tre"))
  }
  
  # Get and write tree sizes
  tree.sizes <- unlist(lapply(trees,function(phy){length(phy$tip.label)}))
  
  cat(tree.sizes,sep="\n",file=paste0(sim.dirs[i],"tree_sizes.txt"))
  cat(quantile(tree.sizes,c(0,0.5,1)),"\n")
  
}


