library(latex2exp)
library(TreePar)

grid.end <- 29.1

times <- -seq(0,29.1,length.out=10)
interval.times <- -c(0,8:21,29.1)
bin.breaks <- c(0,8:21,29.1)

mygreen <- "#238B4590"
mysolidgreen <- "#238B45"

#######
# Trees for heatmaps
#######
beast <- read.nexus("2_empirical_analyses/output/HIV_BEAST_combined_trees.trees")

beast.bt <- lapply(beast,function(phy){
  bt <- getx(phy,TRUE)
  bt <- bt[bt[,2] == 1,1]
  return(sort(bt))
})
beast.bt <- do.call(rbind,beast.bt)

beast.bt.flat <- as.numeric(beast.bt)

beast.bt.counts <- sapply(2:length(bin.breaks),function(i) {
  if (i == length(bin.breaks)) {
    sum(beast.bt.flat > bin.breaks[i-1])
  } else {
    sum(beast.bt.flat < bin.breaks[i] & beast.bt.flat > bin.breaks[i-1])
  }
})

beast.bt.freqs <- beast.bt.counts/prod(dim(beast.bt))

beast.bt.cols <- rev(grey.colors(length(beast.bt.freqs)))
beast.bt.cols <- beast.bt.cols[cut(beast.bt.freqs,breaks=15,labels=FALSE,include.lowest=TRUE)]

# Get beast analyses and combine 
beast1 <- read.table("2_empirical_analyses/output/env_subtA_run_1.log",header=TRUE)
beast2 <- read.table("2_empirical_analyses/output/env_subtA_run_2.log",header=TRUE)
beast3 <- read.table("2_empirical_analyses/output/env_subtA_run_3.log",header=TRUE)
beast4 <- read.table("2_empirical_analyses/output/env_subtA_run_4.log",header=TRUE)

beast <- do.call(rbind,list(beast1,beast2,beast3,beast4))

re.cols <- which(grepl("reproductiveNumber",names(beast)))

beast.Re <- beast[,rev(re.cols)]

Re.beast.50 <- apply(beast.Re,2,median)

Re.beast.10 <- apply(beast.Re,2,quantile,probs=0.10)
Re.beast.90 <- apply(beast.Re,2,quantile,probs=0.90)

Re.beast.05 <- apply(beast.Re,2,quantile,probs=0.05)
Re.beast.95 <- apply(beast.Re,2,quantile,probs=0.95)

r <- range(c(unlist(Re.beast.10),unlist(Re.beast.90)*1.1),na.rm=TRUE)

plot(NULL,NULL,ylim=c(-0.5,r[2]),xlim=range(times),log="",ylab=TeX("$R_e(t)$"),xlab="time",xaxt="n",main="")
axis(side=1,at=seq(0,-26,-2),labels=seq(2011,1985,-2))
for (i in 1:(length(interval.times)-1)) {
  lines(x=c(interval.times[i],interval.times[i+1]),y=c(Re.beast.50[i],Re.beast.50[i]),col=mysolidgreen,lwd=3)
  lines(x=c(interval.times[i+1],interval.times[i+1]),y=c(Re.beast.50[i],Re.beast.50[i+1]),col=mysolidgreen,lwd=3)
}
for (i in 1:(length(interval.times)-1)) {
  polygon(x=c(interval.times[i],interval.times[i],interval.times[i+1],interval.times[i+1]),y=c(Re.beast.05[i],Re.beast.95[i],Re.beast.95[i],Re.beast.05[i]),col=mygreen,border=NA)
}
abline(h=1,lty=2)

#######
# Plots!
#######

r <- range(c(unlist(Re.beast.10),unlist(Re.beast.90)*1.1),na.rm=TRUE)
y1 <- -0.5
y2 <- 0.0
pdf("2_empirical_analyses/4_BEAST/figures/beast_hiv_re.pdf",width=3.5,height=2)
par(mfrow=c(1,1),mai=c(0.65,0.65,0.01,0.01),omi=c(0.01,0.01,0.01,0.01),lend=2)
plot(NULL,NULL,ylim=c(-0.5,r[2]),xlim=range(times),log="",ylab="",xlab="",xaxt="n",main="")
axis(side=1,at=seq(0,-26,-2),labels=seq(2011,1985,-2))
for (i in 1:(length(interval.times)-1)) {
  polygon(x=c(interval.times[i],interval.times[i],interval.times[i+1],interval.times[i+1]),y=c(Re.beast.05[i],Re.beast.95[i],Re.beast.95[i],Re.beast.05[i]),col=mygreen,border=NA)
}
for (i in 1:(length(interval.times)-1)) {
  lines(x=c(interval.times[i],interval.times[i+1]),y=c(Re.beast.50[i],Re.beast.50[i]),col=mysolidgreen,lwd=3)
  lines(x=c(interval.times[i+1],interval.times[i+1]),y=c(Re.beast.50[i],Re.beast.50[i+1]),col=mysolidgreen,lwd=3)
}
abline(h=1,lty=2)
mtext(TeX("$R_e(t)$"),2,line=2,cex=1.25)
mtext("time",1,line=2,cex=1.25)
legend("topright",fill=mygreen,legend="BEAST (IID)",border=NA,bty="n",cex=0.8)

for (i in 1:length(beast.bt.freqs)) {
  polygon(c(-bin.breaks[i+1],-bin.breaks[i+1],-bin.breaks[i],-bin.breaks[i]),c(y1,y2,y2,y1),col=beast.bt.cols[i],border=NA)
}

dev.off()

#######
# Bayes Factors
#######

beast.fold.change.end.90s <- beast.Re[,6]/beast.Re[,5]
quantile(beast.fold.change.end.90s,c(0.05,0.95))
beast.pr.post.shift.end.90s <- sum(beast.fold.change.end.90s > 1)/length(beast.fold.change.end.90s)
beast.BF.shift.end.90s <- beast.pr.post.shift.end.90s/(1-beast.pr.post.shift.end.90s)
2 * log(beast.BF.shift.end.90s)

beast.fold.change.start.90s <- beast.Re[,10]/beast.Re[,12]
quantile(beast.fold.change.start.90s,c(0.05,0.95))
beast.pr.post.shift.start.90s <- sum(beast.fold.change.start.90s > 1)/length(beast.fold.change.start.90s)
beast.BF.shift.start.90s <- beast.pr.post.shift.start.90s/(1-beast.pr.post.shift.start.90s)
2 * log(beast.BF.shift.start.90s)


beast.fold.change.early.mid.90s <- beast.Re[,10]/beast.Re[,9]
quantile(beast.fold.change.early.mid.90s,c(0.05,0.95))
beast.pr.post.shift.early.mid.90s <- sum(beast.fold.change.early.mid.90s > 1)/length(beast.fold.change.early.mid.90s)
beast.BF.shift.early.mid.90s <- beast.pr.post.shift.early.mid.90s/(1-beast.pr.post.shift.early.mid.90s)
2 * log(beast.BF.shift.early.mid.90s)

beast.fold.change.late.mid.90s <- beast.Re[,6]/beast.Re[,7]
quantile(beast.fold.change.late.mid.90s,c(0.05,0.95))
beast.pr.post.shift.late.mid.90s <- sum(beast.fold.change.late.mid.90s > 1)/length(beast.fold.change.late.mid.90s)
beast.BF.shift.late.mid.90s <- beast.pr.post.shift.late.mid.90s/(1-beast.pr.post.shift.late.mid.90s)
2 * log(beast.BF.shift.late.mid.90s)
