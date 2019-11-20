library(phangorn)
library(treedater)

aln <- as.character(read.phyDat("2_empirical_analyses/2_HIV/data/env_subtA.fas","fasta"))
n <- dim(aln)[1]
m <- dim(aln)[2]

# Make text file with ages
taxa <- row.names(aln)
ages <- as.numeric(do.call(rbind,strsplit(taxa,".",fixed=TRUE))[,1])
ages.backwards.time <- max(ages) - ages

rev.ages <- cbind(taxa,ages.backwards.time)
colnames(rev.ages) <- c("taxon","age")

write.table(rev.ages,"2_empirical_analyses/2_HIV/data/env_subtA_taxa.txt",quote=FALSE,row.names=FALSE)

# Guesstimate of age of tree
# We already have a RAxML tree, we just need to ultrametricize it
raxml <- read.tree("2_empirical_analyses/2_HIV/data/raxml/RAxML_bestTree.UkraineENV")
raxml <- unroot(raxml)

sts <- raxml$tip.label
sts <- as.numeric(do.call(rbind,strsplit(sts,".",fixed=TRUE))[,1])
names(sts) <- raxml$tip.label

ultrametric <- dater(raxml,sts=sts,s=m,searchRoot=10,strictClock=TRUE)

write.tree(ultrametric,"2_empirical_analyses/2_HIV/data/treedater.tre")

# Guesstimate birth rate for prior
estimateNetDivFromSerialSamples <- function(n.serial.samples,age,MRCA=FALSE) {
  n_births <- n.serial.samples - 1 - MRCA
  r <- (log(n_births + 1 + MRCA) - MRCA*log(2))/age
  return(r)
}

estimateNetDivFromSerialSamples(n,max(getx(ultrametric,TRUE)),TRUE)

