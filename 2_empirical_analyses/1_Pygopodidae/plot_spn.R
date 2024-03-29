library(latex2exp)
library(ape)

source("1_simulation_study/0_setup/simulation/zeta_helper_functions.R")
source("1_simulation_study/0_setup/simulation/set_zeta.R")

mylightblue <- "#4280f450"
mylightorange <- "#ffa03550"

myblue <- "#4280f490"
myorange <- "#ffa03590"

mysolidblue <- "#4280f4"
mysolidorange <- "#ffa035"

grid.end <- 29

times <- -seq(0.5 * grid.end/100,grid.end*(99.5/100),grid.end/100)
interval.times <- -grid.end * seq(1,99,1)/100
bin.breaks <- grid.end * seq(0,50,1)/50

#######
# Trees for heatmaps
#######
# Get HSRF analyses and combine (run 2 failed)
hsrf1 <- read.tree("2_empirical_analyses/output/Pygopodidae_HSMRFBDP_1.trees")
# hsrf2 <- read.tree("2_empirical_analyses/output/Pygopodidae_HSMRFBDP_2.trees")
hsrf3 <- read.tree("2_empirical_analyses/output/Pygopodidae_HSMRFBDP_3.trees")
hsrf4 <- read.tree("2_empirical_analyses/output/Pygopodidae_HSMRFBDP_4.trees")

hsrf <- c(hsrf1,hsrf3,hsrf4)

rm(hsrf1,hsrf3,hsrf4)

hsrf.bt <- lapply(hsrf,function(phy){sort(branching.times(phy))})
hsrf.bt <- do.call(rbind,hsrf.bt)

hsrf.bt.flat <- as.numeric(hsrf.bt)

hsrf.bt.counts <- sapply(2:length(bin.breaks),function(i) {
  sum(hsrf.bt.flat < bin.breaks[i] & hsrf.bt.flat > bin.breaks[i-1])
})

hsrf.bt.freqs <- hsrf.bt.counts/prod(dim(hsrf.bt))

hsrf.bt.cols <- rev(grey.colors(length(hsrf.bt.freqs)))
hsrf.bt.cols <- hsrf.bt.cols[cut(hsrf.bt.freqs,breaks=50,labels=FALSE,include.lowest=TRUE)]

# Get GMRF analyses and combine
gmrf1 <- read.tree("2_empirical_analyses/output/Pygopodidae_GMRFBDP_1.trees")
gmrf2 <- read.tree("2_empirical_analyses/output/Pygopodidae_GMRFBDP_2.trees")
gmrf3 <- read.tree("2_empirical_analyses/output/Pygopodidae_GMRFBDP_3.trees")
gmrf4 <- read.tree("2_empirical_analyses/output/Pygopodidae_GMRFBDP_4.trees")

gmrf <- c(gmrf1,gmrf2,gmrf3,gmrf4)

rm(gmrf1,gmrf2,gmrf3,gmrf4)

gmrf.bt <- lapply(gmrf,function(phy){sort(branching.times(phy))})
gmrf.bt <- do.call(rbind,gmrf.bt)

gmrf.bt.flat <- as.numeric(gmrf.bt)

gmrf.bt.counts <- sapply(2:length(bin.breaks),function(i) {
  sum(gmrf.bt.flat < bin.breaks[i] & gmrf.bt.flat > bin.breaks[i-1])
})

gmrf.bt.freqs <- gmrf.bt.counts/prod(dim(gmrf.bt))

gmrf.bt.cols <- rev(grey.colors(length(gmrf.bt.freqs)))
gmrf.bt.cols <- gmrf.bt.cols[cut(gmrf.bt.freqs,breaks=50,labels=FALSE,include.lowest=TRUE)]

#######
# Rates over time
#######

# Get HSRF analyses and combine (run 2 failed)
hsrf1 <- read.table("2_empirical_analyses/output/Pygopodidae_HSMRFBDP_1.log",header=TRUE)
# hsrf2 <- read.table("2_empirical_analyses/output/Pygopodidae_HSMRFBDP_2.log",header=TRUE)
hsrf3 <- read.table("2_empirical_analyses/output/Pygopodidae_HSMRFBDP_3.log",header=TRUE)
hsrf4 <- read.table("2_empirical_analyses/output/Pygopodidae_HSMRFBDP_4.log",header=TRUE)

spn.columns <- grepl("speciation_pygo",names(hsrf1))

hsrf <- rbind(hsrf1[,spn.columns],hsrf3[,spn.columns],hsrf4[,spn.columns])

rm(hsrf1,hsrf3,hsrf4)

# Get GSRF analyses and combine
gmrf1 <- read.table("2_empirical_analyses/output/Pygopodidae_GMRFBDP_1.log",header=TRUE)
gmrf2 <- read.table("2_empirical_analyses/output/Pygopodidae_GMRFBDP_2.log",header=TRUE)
gmrf3 <- read.table("2_empirical_analyses/output/Pygopodidae_GMRFBDP_3.log",header=TRUE)
gmrf4 <- read.table("2_empirical_analyses/output/Pygopodidae_GMRFBDP_4.log",header=TRUE)

spn.columns <- grepl("speciation_pygo",names(gmrf1))

gmrf <- rbind(gmrf1[,spn.columns],gmrf2[,spn.columns],gmrf3[,spn.columns],gmrf4[,spn.columns])

rm(gmrf1,gmrf2,gmrf3,gmrf4)

# Summaries for plotting
spn.hsrf.50 <- apply(hsrf,2,median)

spn.hsrf.05 <- apply(hsrf,2,quantile,probs=0.05)
spn.hsrf.95 <- apply(hsrf,2,quantile,probs=0.95)

spn.gmrf.50 <- apply(gmrf,2,median)

spn.gmrf.05 <- apply(gmrf,2,quantile,probs=0.05)
spn.gmrf.95 <- apply(gmrf,2,quantile,probs=0.95)

#########
# Plots!
#########
y1 <- -0.025
y2 <- 0.0
pdf("2_empirical_analyses/1_Pygopodidae/figures/speciation.pdf",width=6,height=2)
  par(mfrow=c(1,2),mai=c(0.65,0.65,0.01,0.01),omi=c(0.01,0.01,0.01,0.01),lend=2)
  plot(NULL,NULL,ylim=c(-0.025,0.65),xlim=range(times),log="",ylab="",xlab="",xaxt="n")
  axis(side=1,at=seq(0,-28,-7),labels=seq(0,28,7))
  polygon(x=c(times,rev(times)),y=c(spn.hsrf.05,rev(spn.hsrf.95)),col=myorange,border=NA)
  lines(times,spn.hsrf.50,lwd=3,col=mysolidorange)
  legend("topright",fill=myorange,legend="HSMRF",border=NA,bty="n")
  mtext(TeX("$\\lambda(t)$"),2,line=2,cex=1.25)
  mtext("time",1,line=2,cex=1.25)
  for (i in 1:length(hsrf.bt.freqs)) {
    polygon(c(-bin.breaks[i+1],-bin.breaks[i+1],-bin.breaks[i],-bin.breaks[i]),c(y1,y2,y2,y1),col=hsrf.bt.cols[i],border=NA)
  }
  
  plot(NULL,NULL,ylim=c(-0.025,0.65),xlim=range(times),log="",ylab="",xlab="",xaxt="n")
  axis(side=1,at=seq(0,-28,-7),labels=seq(0,28,7))
  polygon(x=c(times,rev(times)),y=c(spn.gmrf.05,rev(spn.gmrf.95)),col=myblue,border=NA)
  lines(times,spn.gmrf.50,lwd=3,col=mysolidblue)
  legend("topright",fill=myblue,legend="GMRF",border=NA,bty="n")
  mtext(TeX("$\\lambda(t)$"),2,line=2,cex=1.25)
  mtext("time",1,line=2,cex=1.25)
  for (i in 1:length(gmrf.bt.freqs)) {
    polygon(c(-bin.breaks[i+1],-bin.breaks[i+1],-bin.breaks[i],-bin.breaks[i]),c(y1,y2,y2,y1),col=gmrf.bt.cols[i],border=NA)
  }
  
dev.off()

#########
# Bayes factors
#########

# Shift is from approximately 12 million years ago to approximately 3 million years ago
# Our time is backwards, so it "starts" at 2 and "ends" at 12

end <- 43 # first interval on the other side of the break at ~12 million years ago
start <- 7 # first interval on this side of the break at ~2 million years ago

quantile(hsrf[,end]/hsrf[,start],prob=c(0.05,0.5,0.95))
quantile(gmrf[,end]/gmrf[,start],prob=c(0.05,0.5,0.95))

shift.size <- 1.0

hsrf.posterior.p.shift <- sum(hsrf[,end]/hsrf[,start] > shift.size)/dim(hsrf)[1]

hsrf.shift.bf <- (hsrf.posterior.p.shift/(1 - hsrf.posterior.p.shift))

# 2 ln BF for a shift between 2 and 12 ma with the HSMRF
2*log(hsrf.shift.bf)


gmrf.posterior.p.shift <- sum(gmrf[,end]/gmrf[,start] > shift.size)/dim(gmrf)[1]

gmrf.shift.bf <- (gmrf.posterior.p.shift/(1 - gmrf.posterior.p.shift))

# 2 ln BF for a shift between 2 and 12 ma with the GMRF
2*log(gmrf.shift.bf)


