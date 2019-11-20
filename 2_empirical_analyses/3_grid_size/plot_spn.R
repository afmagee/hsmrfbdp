library(latex2exp)

source("1_simulation_study/0_setup/simulation/zeta_helper_functions.R")
source("1_simulation_study/0_setup/simulation/set_zeta.R")

mylightblue <- "#4280f450"
mylightorange <- "#ffa03550"

myblue <- "#4280f490"
myorange <- "#ffa03590"
mygreen <- "#238B4590"

mysolidblue <- "#4280f4"
mysolidorange <- "#ffa035"
mysolidgreen <- "#238B45"

grid.end <- 29

makeTimes <- function(n) {
  start <- 0.5/n * grid.end
  end <- (n-0.5)/n * grid.end
  return(-seq(start,end,grid.end/n))
}

grid.sizes <- c(10,20,50,100,200)

hsrf <- vector("list",length(grid.sizes))

for (i in 1:length(grid.sizes)) {
  n <- grid.sizes[i]
  # Get HSRF analyses and combine
  hsrf1 <- read.table(paste0("2_empirical_analyses/output/HSMRFBDP_grid_size_",n,"_run_1.log"),header=TRUE)
  hsrf2 <- read.table(paste0("2_empirical_analyses/output/HSMRFBDP_grid_size_",n,"_run_2.log"),header=TRUE)
  
  spn.columns <- which(grepl("^speciation",names(hsrf1)))[1:n]
  
  hsrf[[i]] <- rbind(hsrf1[,spn.columns],hsrf2[,spn.columns])
  
  rm(hsrf1,hsrf2)
}


gmrf <- vector("list",length(grid.sizes))

for (i in 1:length(grid.sizes)) {
  n <- grid.sizes[i]
  # Get gmrf analyses and combine
  gmrf1 <- read.table(paste0("2_empirical_analyses/output/GMRFBDP_grid_size_",n,"_run_1.log"),header=TRUE)
  gmrf2 <- read.table(paste0("2_empirical_analyses/output/GMRFBDP_grid_size_",n,"_run_2.log"),header=TRUE)
  
  spn.columns <- which(grepl("^speciation",names(gmrf1)))[1:n]
  
  gmrf[[i]] <- rbind(gmrf1[,spn.columns],gmrf2[,spn.columns])
  
  rm(gmrf1,gmrf2)
}


# Plots!
pdf("2_empirical_analyses/3_grid_size/grid_size_comparison.pdf",width=4,height=1.5*length(grid.sizes))
par(mfrow=c(length(grid.sizes),2),mai=c(0.01,0.01,0.01,0.01),omi=c(0.65,0.6,0.01,0.01),lend=2)

for (i in 1:length(grid.sizes)) {
  spn.hsrf.50 <- apply(hsrf[[i]],2,median)
  
  spn.hsrf.10 <- apply(hsrf[[i]],2,quantile,probs=0.10)
  spn.hsrf.90 <- apply(hsrf[[i]],2,quantile,probs=0.90)
  
  spn.gmrf.50 <- apply(gmrf[[i]],2,median)
  
  spn.gmrf.10 <- apply(gmrf[[i]],2,quantile,probs=0.10)
  spn.gmrf.90 <- apply(gmrf[[i]],2,quantile,probs=0.90)
  
  times <- makeTimes(grid.sizes[i])
  
  plot(NULL,NULL,ylim=c(0.01,0.45),xlim=range(times),log="",ylab="",xlab="",xaxt="n",yaxt="n")
  polygon(x=c(times,rev(times)),y=c(spn.hsrf.10,rev(spn.hsrf.90)),col=myorange,border=NA)
  lines(times,spn.hsrf.50,lwd=3,col=mysolidorange)
  
  axis(side=2,cex.axis=1.4)
  
  mtext(TeX("$\\lambda(t)$"),2,line=2,cex=1.5)
  
  if (i == 1) {
    legend("topright",fill=myorange,legend="HSMRF",border=NA,bty="n")
  }
  
  if (i == length(grid.sizes)) {
    axis(side=1,at=seq(0,-28,-4),labels=seq(0,28,4),cex.axis=1.4)
    mtext("time",1,line=2.5,cex=1.4)
  }
  
  plot(NULL,NULL,ylim=c(0.01,0.45),xlim=range(times),log="",ylab="",xlab="",xaxt="n",yaxt="n")
  polygon(x=c(times,rev(times)),y=c(spn.gmrf.10,rev(spn.gmrf.90)),col=myblue,border=NA)
  lines(times,spn.gmrf.50,lwd=3,col=mysolidblue)

  if (i == 1) {
    legend("topright",fill=myblue,legend="GMRF",border=NA,bty="n")
  }
  
  if (i == length(grid.sizes)) {
    axis(side=1,at=seq(0,-28,-4),labels=seq(0,28,4),cex.axis=1.4)
    mtext("time",1,line=2.5,cex=1.4)
  }
  
}


dev.off()

