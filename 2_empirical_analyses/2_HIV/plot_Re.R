library(latex2exp)
library(TreeSim)

mylightblue <- "#4280f450"
mylightorange <- "#ffa03550"

myblue <- "#4280f490"
myorange <- "#ffa03590"

mysolidblue <- "#4280f4"
mysolidorange <- "#ffa035"

grid.end <- 29.1

times <- -seq(0.5 * grid.end/100,grid.end*(99.5/100),grid.end/100)
interval.times <- -grid.end * seq(1,99,1)/100
bin.breaks <- grid.end * seq(0,50,1)/50

#######
# Trees for heatmaps
#######
# Get HSRF analyses and combine (run3 3 and 4 failed)
hsrf1 <- read.tree("2_empirical_analyses/output/HIV_env_HSMRFBDP_1.trees")
hsrf2 <- read.tree("2_empirical_analyses/output/HIV_env_HSMRFBDP_2.trees")
# hsrf3 <- read.tree("2_empirical_analyses/output/HIV_env_HSMRFBDP_3.trees")
# hsrf4 <- read.tree("2_empirical_analyses/output/HIV_env_HSMRFBDP_4.trees")

hsrf <- c(hsrf1,hsrf2)

rm(hsrf1,hsrf2)

hsrf.bt <- lapply(hsrf,function(phy){
  bt <- getx(phy,TRUE)
  bt <- bt[bt[,2] == 1,1]
  return(sort(bt))
})
hsrf.bt <- do.call(rbind,hsrf.bt)

hsrf.bt.flat <- as.numeric(hsrf.bt)

hsrf.bt.counts <- sapply(2:length(bin.breaks),function(i) {
  if (i == length(bin.breaks)) {
    sum(hsrf.bt.flat > bin.breaks[i-1])
  } else {
    sum(hsrf.bt.flat < bin.breaks[i] & hsrf.bt.flat > bin.breaks[i-1])
  }
})

hsrf.bt.freqs <- hsrf.bt.counts/prod(dim(hsrf.bt))

hsrf.bt.cols <- rev(grey.colors(length(hsrf.bt.freqs)))
hsrf.bt.cols <- hsrf.bt.cols[cut(hsrf.bt.freqs,breaks=50,labels=FALSE,include.lowest=TRUE)]

# Get gmrf analyses and combine (run 4 failed)
gmrf1 <- read.tree("2_empirical_analyses/output/HIV_env_GMRFBDP_1.trees")
gmrf2 <- read.tree("2_empirical_analyses/output/HIV_env_GMRFBDP_2.trees")
gmrf3 <- read.tree("2_empirical_analyses/output/HIV_env_GMRFBDP_3.trees")
# gmrf4 <- read.tree("2_empirical_analyses/output/HIV_env_GMRFBDP_4.trees")

gmrf <- c(gmrf1,gmrf2)

rm(gmrf1,gmrf2)

gmrf.bt <- lapply(gmrf,function(phy){
  bt <- getx(phy,TRUE)
  bt <- bt[bt[,2] == 1,1]
  return(sort(bt))
})
gmrf.bt <- do.call(rbind,gmrf.bt)

gmrf.bt.flat <- as.numeric(gmrf.bt)

gmrf.bt.counts <- sapply(2:length(bin.breaks),function(i) {
  if (i == length(bin.breaks)) {
    sum(gmrf.bt.flat > bin.breaks[i-1])
  } else {
    sum(gmrf.bt.flat < bin.breaks[i] & gmrf.bt.flat > bin.breaks[i-1])
  }
})

gmrf.bt.freqs <- gmrf.bt.counts/prod(dim(gmrf.bt))

gmrf.bt.cols <- rev(grey.colors(length(gmrf.bt.freqs)))
gmrf.bt.cols <- gmrf.bt.cols[cut(gmrf.bt.freqs,breaks=50,labels=FALSE,include.lowest=TRUE)]


# Get HSRF analyses and combine 
hsrf1 <- read.table("2_empirical_analyses/output/HIV_env_HSMRFBDP_1.log",header=TRUE)
hsrf2 <- read.table("2_empirical_analyses/output/HIV_env_HSMRFBDP_2.log",header=TRUE)
# hsrf3 <- read.table("2_empirical_analyses/output/HIV_env_HSMRFBDP_3.log",header=TRUE)
# hsrf4 <- read.table("2_empirical_analyses/output/HIV_env_HSMRFBDP_4.log",header=TRUE)

hsrf <- do.call(rbind,list(hsrf1,hsrf2))

lambda.columns <- grepl("^speciation",names(hsrf1)) & !grepl("_",names(hsrf1))

hsrf.lambda <- hsrf[,grepl("^birth_rate_env",names(hsrf))]

hsrf.mu <- c(hsrf$death_rate_env)

hsrf.phi <- c(hsrf$serial_sampling_rate)

hsrf.r <- 1.0

hsrf.delta <- hsrf.mu + hsrf.r * hsrf.phi

hsrf.Re <- hsrf.lambda/hsrf.delta

Re.hsrf.50 <- apply(hsrf.Re,2,median)

Re.hsrf.05 <- apply(hsrf.Re,2,quantile,probs=0.05)
Re.hsrf.95 <- apply(hsrf.Re,2,quantile,probs=0.95)

# Get gmrf analyses and combine 
gmrf1 <- read.table("2_empirical_analyses/output/HIV_env_GMRFBDP_1.log",header=TRUE)
gmrf2 <- read.table("2_empirical_analyses/output/HIV_env_GMRFBDP_2.log",header=TRUE)
gmrf3 <- read.table("2_empirical_analyses/output/HIV_env_GMRFBDP_3.log",header=TRUE)
# gmrf4 <- read.table("2_empirical_analyses/output/HIV_env_GMRFBDP_4.log",header=TRUE)

gmrf <- do.call(rbind,list(gmrf1,gmrf2,gmrf3))

lambda.columns <- grepl("^speciation",names(gmrf1)) & !grepl("_",names(gmrf1))

gmrf.lambda <- gmrf[,grepl("^birth_rate_env",names(gmrf))]

gmrf.mu <- c(gmrf$death_rate_env)

gmrf.phi <- c(gmrf$serial_sampling_rate)

gmrf.r <- 1.0

gmrf.delta <- gmrf.mu + gmrf.r * gmrf.phi

gmrf.Re <- gmrf.lambda/gmrf.delta

Re.gmrf.50 <- apply(gmrf.Re,2,median)

Re.gmrf.05 <- apply(gmrf.Re,2,quantile,probs=0.05)
Re.gmrf.95 <- apply(gmrf.Re,2,quantile,probs=0.95)

r <- range(c(unlist(Re.hsrf.10),unlist(Re.hsrf.90),unlist(Re.gmrf.10),unlist(Re.gmrf.90)),na.rm=TRUE)
y1 <- -0.5
y2 <- 0.0
pdf("2_empirical_analyses/2_HIV/figures/hiv_re.pdf",width=6,height=2)
  par(mfrow=c(1,2),mai=c(0.65,0.65,0.01,0.01),omi=c(0.01,0.01,0.01,0.01),lend=2)
  plot(NULL,NULL,ylim=c(-0.5,13),xlim=range(times),log="",ylab="",xlab="",xaxt="n",main="")
  axis(side=1,at=seq(0,-29,-4),labels=seq(2011,1982,-4))
  polygon(x=c(times,rev(times)),y=c(Re.hsrf.10,rev(Re.hsrf.90)),col=myorange,border=NA)
  lines((times),Re.hsrf.50,lwd=3,col=mysolidorange)
  abline(h=1,lty=2)
  mtext(TeX("$R_e(t)$"),2,line=2,cex=1.25)
  mtext("time",1,line=2,cex=1.25)
  legend("topright",fill=myorange,legend="HSMRF",border=NA,bty="n")
  for (i in 1:length(hsrf.bt.freqs)) {
    polygon(c(-bin.breaks[i+1],-bin.breaks[i+1],-bin.breaks[i],-bin.breaks[i]),c(y1,y2,y2,y1),col=hsrf.bt.cols[i],border=NA)
  }
  
  plot(NULL,NULL,ylim=c(-0.5,r[2]),xlim=range(times),log="",ylab="",xlab="",xaxt="n",main="")
  axis(side=1,at=seq(0,-29,-4),labels=seq(2011,1982,-4))
  polygon(x=c(times,rev(times)),y=c(Re.gmrf.10,rev(Re.gmrf.90)),col=myblue,border=NA)
  lines((times),Re.gmrf.50,lwd=3,col=mysolidblue)
  abline(h=1,lty=2)
  mtext(TeX("$R_e(t)$"),2,line=2,cex=1.25)
  mtext("time",1,line=2,cex=1.25)
  legend("topright",fill=myblue,legend="GMRF",border=NA,bty="n")
  for (i in 1:length(gmrf.bt.freqs)) {
    polygon(c(-bin.breaks[i+1],-bin.breaks[i+1],-bin.breaks[i],-bin.breaks[i]),c(y1,y2,y2,y1),col=gmrf.bt.cols[i],border=NA)
  }
  
dev.off()

###########
# Bayes Factors
###########

# Shift at the end of the 1990s
start <- 34 # First interval this side of 2001
end <- 42 # First interval on the other side of 1999

quantile(hsrf.lambda[,end]/hsrf.lambda[,start],c(0.05,0.95))

shift.size <- 1.0

hsrf.posterior.p.shift <- sum(hsrf.lambda[,end]/hsrf.lambda[,start] > shift.size)/dim(hsrf.lambda)[1]
hsrf.shift.bf <- (hsrf.posterior.p.shift/(1 - hsrf.posterior.p.shift))
2*log(hsrf.shift.bf)

gmrf.posterior.p.shift <- sum(gmrf.lambda[,end]/gmrf.lambda[,start] > shift.size)/dim(gmrf.lambda)[1]
gmrf.shift.bf <- (gmrf.posterior.p.shift/(1 - gmrf.posterior.p.shift))
2*log(gmrf.shift.bf)


# Shift at the beginning of the 1990s
start <- 57 # First interval this side of mid-1994
end <- 66 # First interval on the other side of 1999

shift.size <- 1.0

hsrf.posterior.p.shift <- sum(hsrf.lambda[,end]/hsrf.lambda[,start] < shift.size)/dim(hsrf.lambda)[1]
hsrf.shift.bf <- (hsrf.posterior.p.shift/(1 - hsrf.posterior.p.shift))
2*log(hsrf.shift.bf)

gmrf.posterior.p.shift <- sum(gmrf.lambda[,end]/gmrf.lambda[,start] < shift.size)/dim(gmrf.lambda)[1]
gmrf.shift.bf <- (gmrf.posterior.p.shift/(1 - gmrf.posterior.p.shift))
2*log(gmrf.shift.bf)


# summaries and CI on elevated Re
elevated.rate <- c(42:56) # Range of birth rate intervals where the rate is elevated, between the shift times specified

min(Re.hsrf.05[elevated.rate])
max(Re.hsrf.95[elevated.rate])
mean(Re.hsrf.50[elevated.rate])

min(Re.gmrf.05[elevated.rate])
max(Re.gmrf.95[elevated.rate])
mean(Re.gmrf.50[elevated.rate])


