library(phangorn)

# Re-assembler Carpho alignment into PartitionFinder subsets, drop gaps

pygo <- read.phyDat("2_empirical_analyses/1_Pygopodidae/data/PartitionFinder/Pygopodidae.phy")

pygo <- as.character(pygo)

subset1 <- c(1:3756)[-seq(2,1087,3)]
subset2 <- seq(2,1087,3)

subset1 <- pygo[,subset1]
subset2 <- pygo[,subset2]

table1 <- apply(subset1,2,table)
table2 <- apply(subset2,2,table)

allgap1 <- unlist(lapply(table1,function(x){
  (length(x) == 1) && (names(x)[1] == "-")
}))
allgap2 <- unlist(lapply(table2,function(x){
  (length(x) == 1) && (names(x)[1] == "-")
}))

subset1 <- as.phyDat(subset1[,!allgap1])
subset2 <- as.phyDat(subset2[,!allgap2])

write.phyDat(subset1,"2_empirical_analyses/1_Australia/data/Pygopodidae_subset_1.nex","nexus")
write.phyDat(subset2,"2_empirical_analyses/1_Australia/data/Pygopodidae_subset_2.nex","nexus")
