cr_template <- '
# This will analyze the simulated datasets with a CRBDP

seed(THISSEED)

#######################
# Reading in the Data #
#######################

### Read in the "observed" tree
psi <- readTrees("data/" + REPLICATE + ".tre")[1]

# Get some useful variables from the data. We need these later on.
taxa <- psi.taxa()
num_species <- psi.ntips()
root_height <- psi.rootAge()

# set my move index
mvi = 0

# set number of diversification regimes
NUM_INTERVALS = 100

NUM_BREAKS = NUM_INTERVALS - 1

# Chop up duration from root to present into equally sized intervals
interval_times <- abs(root_height * seq(1, NUM_BREAKS, 1)/NUM_INTERVALS)

rho <- 1.0

####################
# Create the rates #
####################

speciation_logmean <- ln(ln(num_species/rho/2)/root_height)
speciation_logsd <- 2*0.587405

speciation ~ dnLognormal(speciation_logmean,speciation_logsd)
extinction ~ dnExponential(10)

moves[++mvi] = mvScaleBactrian(speciation,weight=2)
moves[++mvi] = mvScaleBactrian(extinction,weight=2)

timetree ~ dnEpisodicBirthDeath(rootAge          = psi.rootAge(), 
lambdaRates                                      = speciation, 
muRates                                          = extinction, 
rho                                              = rho, 
samplingStrategy                                 = "uniform", 
condition                                        = "time", 
taxa                                             = taxa)

### clamp the model with the "observed" tree
timetree.clamp(psi)

#############
# The Model #
#############


### workspace model wrapper ###
mymodel = model(rho)

### set up the monitors that will output parameter values to file and screen 
monitors[1] = mnModel(filename="output/CRBDP_" + REPLICATE + "_run_" + THISRUN + ".log",printgen=200, separator = TAB)


################
# The Analysis #
################

### workspace mcmc ###
mymcmc = mcmc(mymodel, monitors, moves, nruns=1)

### Use burnin to tune the MH proposals ###
mymcmc.burnin(generations=25000,tuningInterval=500)

### run the MCMC ###
mymcmc.run(generations=200000)


q()
'

# This script requires src and src/analysis, src/simulation, src/hyak to exist
# This script assumes the data are already simulated
# For ease of running on hyak, we print one copy of each analysis script to a directory for each simulation set we run

# Make sure we know all the seeds
seed <- 8472

for (i in 1:100) {
  this.rb <- gsub("REPLICATE",i,cr_template)
  for (j in 1:2) {
    this.chain.rb <- gsub("THISRUN",j,this.rb)
    this.chain.rb <- gsub("THISSEED",seed,this.chain.rb)
    for (k in 1:4) {
      cat(this.chain.rb,file=paste0("1_simulation_study/0_setup/analysis_scripts/CRBDP_batch_",i,"_chain_",j,".Rev"))
    }
    seed <- seed + 1
  }
}
