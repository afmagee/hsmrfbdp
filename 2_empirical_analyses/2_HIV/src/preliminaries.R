library(phangorn)
setwd("")

# Get BEAST tree
chloro <- read.nexus("data/Chloranthaceae.tre")

# 4th tree is BEAST
chloro <- chloro[[4]]

### Get a rough estimate of a clock rate
its <- as.character(read.phyDat("data/subsets/ITS.nex","nexus"))
rbcl <- as.character(read.phyDat("data/subsets/rbcL.nex","nexus"))
rps16 <- as.character(read.phyDat("data/subsets/rps16.nex","nexus"))

weights <- c(dim(its)[2],dim(rbcl)[2],dim(rps16)[2])
weights <- weights/sum(weights)

# No gap-only sequences
its <- its[lengths(apply(its,1,table)) > 1,]
rbcl <- rbcl[lengths(apply(rbcl,1,table)) > 1,]
rps16 <- rps16[lengths(apply(rps16,1,table)) > 1,]

# What is seen at each site
its.obs <- apply(its,2,table)
rbcl.obs <- apply(rbcl,2,table)
rps16.obs <- apply(rps16,2,table)

# Remove all gaps before considering sites invariant
its.obs <- lapply(its.obs,function(obs){obs[names(obs) != "?"]})
rbcl.obs <- lapply(rbcl.obs,function(obs){obs[names(obs) != "?"]})
rps16.obs <- lapply(rps16.obs,function(obs){obs[names(obs) != "?"]})

# Number of observed invariant sites
its.inv <- sum(lengths(its.obs) == 1)
rbcl.inv <- sum(lengths(rbcl.obs) == 1)
rps16.inv <- sum(lengths(rps16.obs) == 1)

its.tl <- -log(its.inv/dim(its)[2])
rbcl.tl <- -log(rbcl.inv/dim(rbcl)[2])
rps16.tl <- -log(rps16.inv/dim(rps16)[2])

partition.substitution.tl <- c(its.tl,rbcl.tl,rps16.tl)

substitution.tl <- sum(partition.substitution.tl * weights)

clock.guestimate <- substitution.tl/sum(chloro$edge.length)
