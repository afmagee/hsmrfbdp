library(coda)
library(parallel)

source("1_simulation_study/0_setup/summarization/rank_based_convergence_diagnostics.R")

n.sims <- 100

# Point to high-level directory, let looping handle the subsets of the simulation study

this.dir <- "1_simulation_study/1_constant_rate/"

# Start reading file names
all.files <- list.files(paste0(this.dir,"/output/"),full.names=TRUE)

# Logs only
all.logs <- all.files[grepl(".log",all.files,fixed=TRUE)]

# To check convergence, we'll want the 4 runs for each
two.chains <- all.logs[grepl("_run",all.logs)]

# Create directory for summaries
if ( !dir.exists(paste0(this.dir,"/summaries")) ) {
  dir.create(paste0(this.dir,"/summaries"))
}

# Let's run convergence diagnostics and assemble parameter estimates
# We will take the mean, median, and 95% quantiles, and 90% quantiles for everything
# We will parallelize the calculations of diagnostics and summaries, then stitch the results into a table later
for (prior in c("CR")) {
  this.prior <- two.chains[grepl(prior,two.chains)]
  res <- mclapply(1:n.sims,function(i){
    # Read logs
    this.chain.logs <- this.prior[grepl(paste0("_",i,"_run_"),this.prior)]
    this.chain.logs <- lapply(this.chain.logs,function(f){
      try(read.table(f,sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE))
    })
    if ( any(unlist(lapply(this.chain.logs,class)) == "try-error") ) {
      est.q.025 <- NA
      est.q.975 <- NA
      est.q.05 <- NA
      est.q.95 <- NA
      est.mean <- NA
      est.median <- NA
      rank.psrf <- NA
      ess <- NA
    } else {
      # Concatenate
      this.chain.log <- do.call(rbind,this.chain.logs)
      # Grab estimates
      est.q.025 <- apply(this.chain.log,2,quantile,probs=0.025)
      est.q.975 <- apply(this.chain.log,2,quantile,probs=0.975)
      est.q.05 <- apply(this.chain.log,2,quantile,probs=0.05)
      est.q.95 <- apply(this.chain.log,2,quantile,probs=0.95)
      est.mean <- apply(this.chain.log,2,mean)
      est.median <- apply(this.chain.log,2,median)
      rank.psrf <- diagnoseConvergence(chains=this.chain.logs,return.both=FALSE)[,1]
      ess <- effectiveSize(this.chain.log)
    }
    return(list(index=i,
                est.q.025=est.q.025,
                est.q.975=est.q.975,
                est.q.05=est.q.05,
                est.q.95=est.q.95,
                est.mean=est.mean,
                est.median=est.median,
                rank.psrf=rank.psrf,
                ess=ess
    ))
  },mc.cores=8,mc.preschedule=FALSE)
  
  # Empty lists (will be made into tables later)
  est.q.025 <- vector("list",n.sims)
  est.q.975 <- vector("list",n.sims)
  est.q.05 <- vector("list",n.sims)
  est.q.95 <- vector("list",n.sims)
  est.mean <- vector("list",n.sims)
  est.median <- vector("list",n.sims)
  rank.psrf <- vector("list",n.sims)
  ess <- vector("list",n.sims)
  
  for (i in 1:n.sims) {
    index <- res[[i]]$index
    est.q.025[[index]] <-  res[[i]]$est.q.025
    est.q.975[[index]] <-  res[[i]]$est.q.975
    est.q.05[[index]] <-  res[[i]]$est.q.05
    est.q.95[[index]] <-  res[[i]]$est.q.95
    est.mean[[index]] <-  res[[i]]$est.mean
    est.median[[index]] <-  res[[i]]$est.median
    rank.psrf[[index]] <-  res[[i]]$rank.psrf
    ess[[index]] <-  res[[i]]$ess
  }
  
  est.q.025 <- do.call(rbind,est.q.025)
  est.q.975 <- do.call(rbind,est.q.975)
  est.q.05 <- do.call(rbind,est.q.05)
  est.q.95 <- do.call(rbind,est.q.95)
  est.mean <- do.call(rbind,est.mean)
  est.median <- do.call(rbind,est.median)
  rank.psrf <- do.call(rbind,rank.psrf)
  ess <- do.call(rbind,ess)
  
  write.csv(est.q.025,file=paste0(this.dir,"/summaries/",prior,"_quantile_025.csv"))
  write.csv(est.q.975,file=paste0(this.dir,"/summaries/",prior,"_quantile_975.csv"))
  write.csv(est.q.05,file=paste0(this.dir,"/summaries/",prior,"_quantile_05.csv"))
  write.csv(est.q.95,file=paste0(this.dir,"/summaries/",prior,"_quantile_95.csv"))
  write.csv(est.mean,file=paste0(this.dir,"/summaries/",prior,"_mean.csv"))
  write.csv(est.median,file=paste0(this.dir,"/summaries/",prior,"_median.csv"))
  write.csv(rank.psrf,file=paste0(this.dir,"/summaries/",prior,"_rank_psrf.csv"))
  write.csv(ess,file=paste0(this.dir,"/summaries/",prior,"_ess.csv"))
}
