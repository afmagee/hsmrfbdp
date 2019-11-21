require(parallel)

rb.path <- "rb"
n.cores <- 4

# List of all scripts for all 2 chains * 2 priors * 100 simulations analyses
revscripts <- list.files("../../0_setup/analysis_scripts",full.names=TRUE)

x <- mclapply(revscripts,function(this.script){
  system2(rb.path,args=this.script,stdout=paste0("output/",basename(this.script),".stdout"),stderr=paste0("output/",basename(this.script),".stderr"))
},mc.cores=n.cores)

