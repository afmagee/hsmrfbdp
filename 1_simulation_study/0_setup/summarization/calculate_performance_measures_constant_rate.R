library(coda)
library(parallel)
library(TESS)

source("1_simulation_study/0_setup/simulation/helper_functions.R")

original.wd <- getwd()

## Information we used in simulating the trees
rho <- 1
tree.age <- 100
ntaxa <- 200
n.sims <- 100
n.windows <- 100

mu <- 0.01

interval.times <- tree.age * seq(1/n.windows,1-(1/n.windows),1/n.windows)
interval.midpoints <- tree.age * seq(0 + 1/(2*n.windows),1-(1/(2*n.windows)),1/n.windows)

# Point to high-level directory and note the fold-change, let looping handle the subsets of the simulation study
this.dir <- "1_simulation_study/1_constant_rate/"
this.fold.change <- 1.0

# Get true speciation rates for each time in the interval
shift.radius <- 0
shift.center <- 50
simulating.params <- findParams(t.shift=shift.center,shift.radius=shift.radius,fold.change=this.fold.change,extinction=mu,target.taxa=ntaxa,tree.age=tree.age,sampling.fraction=rho)

# Read in simulating values
true.spn <- simulating.params[1]

setwd(this.dir)

# Read in results
constant.psrf <- read.csv("summaries/CR_rank_psrf.csv",row.names=1)
constant.ess <- read.csv("summaries/CR_ess.csv",row.names=1)
constant.median <- read.csv("summaries/CR_median.csv")
constant.q.05 <- read.csv("summaries/CR_quantile_05.csv")
constant.q.95 <- read.csv("summaries/CR_quantile_95.csv")

# Process convergence diagnostics on model parameters
constant.psrf.failures <- which(apply(constant.psrf,1,function(psrf){any(psrf > 1.01)}))
length(constant.psrf.failures)
constant.psrf.successes <- which(!c(1:100) %in% constant.psrf.failures)
n.constant.successful <- length(constant.psrf.successes)

# Trim out poor performers
constant.ess <- constant.ess[constant.psrf.successes,]
constant.median <- constant.median[constant.psrf.successes,]
constant.q.05 <- constant.q.05[constant.psrf.successes,]
constant.q.95 <- constant.q.95[constant.psrf.successes,]

# Find appropriate columns in files
constant.speciation.columns <- (grepl("^speciation",names(constant.median)) & !grepl("_",names(constant.median),fixed=TRUE))

## Some analysis of performance, CR

# We use 90% CIs instead of 95% CIs, because 90% is more stable (and here the tails are fat enough we might worry), and 95% are often so big as to be almost meaningless
constant.spn.in.90 <- matrix(NA,n.constant.successful,n.windows)
constant.spn.width.90 <- matrix(NA,n.constant.successful,n.windows)

for (i in 1:n.constant.successful) {
  this.spn.05 <- as.numeric(constant.q.05[i,constant.speciation.columns])
  this.spn.95 <- as.numeric(constant.q.95[i,constant.speciation.columns])
  # Rev has speciation from present to past, we have past to present
  this.spn.05 <- rev(this.spn.05)
  this.spn.95 <- rev(this.spn.95)
  constant.spn.width.90[i,] <- this.spn.95 - this.spn.05
  # True param is in interval if 95% is > param and 5% is less
  # subtract true value and 95%>0, 5%<0, if both > or both <, this product is positive
  this.spn.05 <- this.spn.05 - true.spn
  this.spn.95 <- this.spn.95 - true.spn
  constant.spn.in.90[i,] <- as.numeric(this.spn.95*this.spn.05 < 0)
}

constant.coverage <- rowSums(constant.spn.in.90)/100

constant.width <- rowMeans(constant.spn.width.90)

# It is more interpretable to look at the widths relative to the true rate than the raw widths
constant.relative.width <- constant.width/mean(true.spn)

constant.speciation.mse <- rep(NA,n.constant.successful)
for (i in 1:n.constant.successful) {
  this.spn <- as.numeric(constant.median[i,constant.speciation.columns])
  constant.speciation.mse[i] <- (this.spn - true.spn)^2
}

constant.speciation.mad <- rep(NA,n.constant.successful)
for (i in 1:n.constant.successful) {
  this.spn <- as.numeric(constant.median[i,constant.speciation.columns])
  constant.speciation.mad[i] <- abs(this.spn - true.spn)
}

constant.minimum.ess <- apply(constant.ess,1,min)

# Write to files
cat(constant.coverage,sep="\n",file="summaries/CR_coverage.txt")  
cat(constant.width,sep="\n",file="summaries/CR_precision.txt")  
cat(constant.relative.width,sep="\n",file="summaries/CR_relative_precision.txt")  
cat(constant.speciation.mad,sep="\n",file="summaries/CR_MAD.txt")  
cat(constant.speciation.mse,sep="\n",file="summaries/CR_MSE.txt")  
cat(constant.minimum.ess,sep="\n",file="summaries/CR_min_ESS.txt")  


setwd(original.wd)
