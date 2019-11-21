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

# Point to high-level directory and note the fold-change
parent <- "1_simulation_study/1_constant_rate//"
this.fold.change <- 1.0

dirs <- list.dirs(parent,recursive=FALSE)


for (this.dir in dirs) {
  
  setwd(this.dir)
  
  # Read in results
  gaussian.psrf <- read.csv("summaries/GMRF_rank_psrf.csv",row.names=1)
  gaussian.ess <- read.csv("summaries/GMRF_ess.csv",row.names=1)
  gaussian.median <- read.csv("summaries/GMRF_median.csv")
  gaussian.q.05 <- read.csv("summaries/GMRF_quantile_05.csv")
  gaussian.q.95 <- read.csv("summaries/GMRF_quantile_95.csv")
  
  horseshoe.psrf <- read.csv("summaries/HSMRF_rank_psrf.csv",row.names=1)
  horseshoe.ess <- read.csv("summaries/HSMRF_ess.csv",row.names=1)
  horseshoe.median <- read.csv("summaries/HSMRF_median.csv")
  horseshoe.q.05 <- read.csv("summaries/HSMRF_quantile_05.csv")
  horseshoe.q.95 <- read.csv("summaries/HSMRF_quantile_95.csv")
  
  # Get true speciation rates for each time in the interval
  tmp <- as.numeric(strsplit(basename(this.dir),"_")[[1]][-1])
  shift.radius <- (tmp[1] * tree.age)/2
  shift.center <- tmp[2] * tree.age
  if ( this.fold.change == 1 ) {
    shift.center <- 50 # Doesn't techincally matter where we say the "shift" is, but this avoids TESS issues
  }
  simulating.params <- findParams(t.shift=shift.center,shift.radius=shift.radius,fold.change=this.fold.change,extinction=mu,target.taxa=ntaxa,tree.age=tree.age)
  
  # Read in simulating values
  true.spn <- sapply(interval.midpoints,speciationRate,t.shift=shift.center,shift.radius=shift.radius,rate.old=simulating.params[1],rate.new=simulating.params[2])
  
  # Process convergence diagnostics on model parameters
  gaussian.psrf.failures <- which(apply(gaussian.psrf,1,function(psrf){any(psrf > 1.01)}))
  length(gaussian.psrf.failures)
  gaussian.psrf.successes <- which(!c(1:100) %in% gaussian.psrf.failures)
  n.gaussian.successful <- length(gaussian.psrf.successes)

  horseshoe.psrf.failures <- which(apply(horseshoe.psrf,1,function(psrf){any(psrf > 1.01)}))
  length(horseshoe.psrf.failures)
  horseshoe.psrf.successes <- which(!c(1:100) %in% horseshoe.psrf.failures)
  n.horseshoe.successful <- length(horseshoe.psrf.successes)

  # Trim out poor performers
  gaussian.ess <- gaussian.ess[gaussian.psrf.successes,]
  gaussian.median <- gaussian.median[gaussian.psrf.successes,]
  gaussian.q.05 <- gaussian.q.05[gaussian.psrf.successes,]
  gaussian.q.95 <- gaussian.q.95[gaussian.psrf.successes,]
  
  horseshoe.ess <- horseshoe.ess[horseshoe.psrf.successes,]
  horseshoe.median <- horseshoe.median[horseshoe.psrf.successes,]
  horseshoe.q.05 <- horseshoe.q.05[horseshoe.psrf.successes,]
  horseshoe.q.95 <- horseshoe.q.95[horseshoe.psrf.successes,]

  
  # Find appropriate columns in files
  gaussian.speciation.columns <- (grepl("^speciation",names(gaussian.median)) & !grepl("_",names(gaussian.median),fixed=TRUE))
  gaussian.delta.columns <- grep("^delta",names(gaussian.median))

  horseshoe.speciation.columns <- (grepl("^speciation",names(horseshoe.median)) & !grepl("_",names(horseshoe.median),fixed=TRUE))
  horseshoe.delta.columns <- grep("^delta",names(horseshoe.median))

  ## Some analysis of performance, HSMRF

  # We use 90% CIs instead of 95% CIs, because 90% is more stable (and here the tails are fat enough we might worry), and 95% are often so big as to be almost meaningless
  horseshoe.spn.in.90 <- matrix(NA,n.horseshoe.successful,n.windows)
  horseshoe.spn.width.90 <- matrix(NA,n.horseshoe.successful,n.windows)

  for (i in 1:n.horseshoe.successful) {
    this.spn.05 <- as.numeric(horseshoe.q.05[i,horseshoe.speciation.columns])
    this.spn.95 <- as.numeric(horseshoe.q.95[i,horseshoe.speciation.columns])
    # Rev has speciation from present to past, we have past to present
    this.spn.05 <- rev(this.spn.05)
    this.spn.95 <- rev(this.spn.95)
    horseshoe.spn.width.90[i,] <- this.spn.95 - this.spn.05
    # True param is in interval if 95% is > param and 5% is less
    # subtract true value and 95%>0, 5%<0, if both > or both <, this product is positive
    this.spn.05 <- this.spn.05 - true.spn
    this.spn.95 <- this.spn.95 - true.spn
    horseshoe.spn.in.90[i,] <- as.numeric(this.spn.95*this.spn.05 < 0)
  }

  horseshoe.coverage <- rowSums(horseshoe.spn.in.90)/100
  
  horseshoe.width <- rowMeans(horseshoe.spn.width.90)
  
  # The average credible interval width across the whole trajectory is only part of the story
  # How much wider/narrower they are at the beginning/end is also important
  horseshoe.width.1.to.50 <- rowMeans(horseshoe.spn.width.90[,1:50])
  horseshoe.width.51.to.100 <- rowMeans(horseshoe.spn.width.90[,51:100])

  # It is more interpretable to look at the widths relative to the true rate than the raw widths
  horseshoe.relative.width <- horseshoe.width/mean(true.spn)

  horseshoe.relative.width.1.to.50 <- horseshoe.width.1.to.50/mean(true.spn)
  horseshoe.relative.width.51.to.100 <- horseshoe.width.51.to.100/mean(true.spn)

  horseshoe.speciation.mse <- rep(NA,n.horseshoe.successful)
  for (i in 1:n.horseshoe.successful) {
    this.spn <- as.numeric(horseshoe.median[i,horseshoe.speciation.columns])
    this.spn <- rev(this.spn)
    horseshoe.speciation.mse[i] <- 1/n.windows * sum(c(this.spn - true.spn)^2)
  }

  horseshoe.speciation.mad <- rep(NA,n.horseshoe.successful)
  for (i in 1:n.horseshoe.successful) {
    this.spn <- as.numeric(horseshoe.median[i,horseshoe.speciation.columns])
    this.spn <- rev(this.spn)
    horseshoe.speciation.mad[i] <- 1/n.windows * sum(abs(this.spn - true.spn))
  }
  
  horseshoe.speciation.tvn <- rep(NA,n.horseshoe.successful)
  for (i in 1:n.horseshoe.successful) {
    this.spn <- as.numeric(horseshoe.median[i,horseshoe.speciation.columns])
    this.spn <- rev(this.spn)
    horseshoe.speciation.tvn[i] <- sum(abs(this.spn[2:n.windows] - this.spn[1:(n.windows-1)]))
  }
  
  # This may be more interpretable than just the TVN
  horseshoe.relative.speciation.tvn <- horseshoe.speciation.tvn / abs(simulating.params[1] - simulating.params[2])
  
  horseshoe.minimum.ess <- apply(horseshoe.ess,1,min)
  
  horseshoe.fold.change <- horseshoe.median$speciation.100./horseshoe.median$speciation.1.
  
  # Write to files
  cat(horseshoe.coverage,sep="\n",file="summaries/HSMRF_coverage.txt")  
  cat(horseshoe.width,sep="\n",file="summaries/HSMRF_precision.txt")  
  cat(horseshoe.relative.width,sep="\n",file="summaries/HSMRF_relative_precision.txt")  
  cat(horseshoe.relative.width.1.to.50,sep="\n",file="summaries/HSMRF_relative_precision_oldest_half.txt")  
  cat(horseshoe.relative.width.51.to.100,sep="\n",file="summaries/HSMRF_relative_precision_youngest_half.txt")  
  cat(horseshoe.speciation.mad,sep="\n",file="summaries/HSMRF_MAD.txt")  
  cat(horseshoe.speciation.mse,sep="\n",file="summaries/HSMRF_MSE.txt")  
  cat(horseshoe.speciation.tvn,sep="\n",file="summaries/HSMRF_TVN.txt")  
  cat(horseshoe.relative.speciation.tvn,sep="\n",file="summaries/HSMRF_relative_TVN.txt")  
  cat(horseshoe.minimum.ess,sep="\n",file="summaries/HSMRF_min_ESS.txt")  
  cat(horseshoe.fold.change,sep="\n",file="summaries/HSMRF_fold_change.txt")  

  ## Some analysis of performance, GMRF
  
  # We use 90% CIs instead of 95% CIs, because 90% is more stable (and here the tails are fat enough we might worry), and 95% are often so big as to be almost meaningless
  gaussian.spn.in.90 <- matrix(NA,n.gaussian.successful,n.windows)
  gaussian.spn.width.90 <- matrix(NA,n.gaussian.successful,n.windows)
  
  for (i in 1:n.gaussian.successful) {
    this.spn.05 <- as.numeric(gaussian.q.05[i,gaussian.speciation.columns])
    this.spn.95 <- as.numeric(gaussian.q.95[i,gaussian.speciation.columns])
    # Rev has speciation from present to past, we have past to present
    this.spn.05 <- rev(this.spn.05)
    this.spn.95 <- rev(this.spn.95)
    gaussian.spn.width.90[i,] <- this.spn.95 - this.spn.05
    # True param is in interval if 95% is > param and 5% is less
    # subtract true value and 95%>0, 5%<0, if both > or both <, this product is positive
    this.spn.05 <- this.spn.05 - true.spn
    this.spn.95 <- this.spn.95 - true.spn
    gaussian.spn.in.90[i,] <- as.numeric(this.spn.95*this.spn.05 < 0)
  }
  
  gaussian.coverage <- rowSums(gaussian.spn.in.90)/100
  
  gaussian.width <- rowMeans(gaussian.spn.width.90)
  
  # The average credible interval width across the whole trajectory is only part of the story
  # How much wider/narrower they are at the beginning/end is also important
  gaussian.width.1.to.50 <- rowMeans(gaussian.spn.width.90[,1:50])
  gaussian.width.51.to.100 <- rowMeans(gaussian.spn.width.90[,51:100])
  gaussian.width.1.to.25 <- rowMeans(gaussian.spn.width.90[,1:25])
  gaussian.width.76.to.100 <- rowMeans(gaussian.spn.width.90[,76:100])
  
  # It is more interpretable to look at the widths relative to the true rate than the raw widths
  gaussian.relative.width <- gaussian.width/mean(true.spn)
  
  gaussian.relative.width.1.to.50 <- gaussian.width.1.to.50/mean(true.spn)
  gaussian.relative.width.51.to.100 <- gaussian.width.51.to.100/mean(true.spn)
  gaussian.relative.width.1.to.25 <- gaussian.width.1.to.25/mean(true.spn)
  gaussian.relative.width.76.to.100 <- gaussian.width.76.to.100/mean(true.spn)
  
  gaussian.speciation.mse <- rep(NA,n.gaussian.successful)
  for (i in 1:n.gaussian.successful) {
    this.spn <- as.numeric(gaussian.median[i,gaussian.speciation.columns])
    this.spn <- rev(this.spn)
    gaussian.speciation.mse[i] <- 1/n.windows * sum(c(this.spn - true.spn)^2)
  }
  
  gaussian.speciation.mad <- rep(NA,n.gaussian.successful)
  for (i in 1:n.gaussian.successful) {
    this.spn <- as.numeric(gaussian.median[i,gaussian.speciation.columns])
    this.spn <- rev(this.spn)
    gaussian.speciation.mad[i] <- 1/n.windows * sum(abs(this.spn - true.spn))
  }
  
  gaussian.speciation.tvn <- rep(NA,n.gaussian.successful)
  for (i in 1:n.gaussian.successful) {
    this.spn <- as.numeric(gaussian.median[i,gaussian.speciation.columns])
    this.spn <- rev(this.spn)
    gaussian.speciation.tvn[i] <- sum(abs(this.spn[2:n.windows] - this.spn[1:(n.windows-1)]))
  }
  
  # This may be more interpretable than just the TVN
  gaussian.relative.speciation.tvn <- gaussian.speciation.tvn / abs(simulating.params[1] - simulating.params[2])
  
  gaussian.minimum.ess <- apply(gaussian.ess,1,min)
  
  gaussian.fold.change <- gaussian.median$speciation.100./gaussian.median$speciation.1.
  
  # Write to files
  cat(gaussian.coverage,sep="\n",file="summaries/GMRF_coverage.txt")  
  cat(gaussian.width,sep="\n",file="summaries/GMRF_precision.txt")  
  cat(gaussian.relative.width,sep="\n",file="summaries/GMRF_relative_precision.txt")  
  cat(gaussian.relative.width.1.to.50,sep="\n",file="summaries/GMRF_relative_precision_oldest_half.txt")  
  cat(gaussian.relative.width.51.to.100,sep="\n",file="summaries/GMRF_relative_precision_youngest_half.txt")  
  cat(gaussian.relative.width.1.to.25,sep="\n",file="summaries/GMRF_relative_precision_oldest_quarter.txt")  
  cat(gaussian.relative.width.76.to.100,sep="\n",file="summaries/GMRF_relative_precision_youngest_quarter.txt")  
  cat(gaussian.speciation.mad,sep="\n",file="summaries/GMRF_MAD.txt")  
  cat(gaussian.speciation.mse,sep="\n",file="summaries/GMRF_MSE.txt")  
  cat(gaussian.speciation.tvn,sep="\n",file="summaries/GMRF_TVN.txt")  
  cat(gaussian.relative.speciation.tvn,sep="\n",file="summaries/GMRF_relative_TVN.txt")  
  cat(gaussian.minimum.ess,sep="\n",file="summaries/GMRF_min_ESS.txt")  
  cat(gaussian.fold.change,sep="\n",file="summaries/GMRF_fold_change.txt")  
  
  setwd(original.wd)
  
}
