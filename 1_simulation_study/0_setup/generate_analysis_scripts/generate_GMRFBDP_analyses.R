gmrf_template <- '
# This will analyze the simulated datasets with an GMRFBDP

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

speciation_at_present ~ dnLognormal(speciation_logmean,speciation_logsd)
extinction ~ dnExponential(10)

moves[++mvi] = mvScaleBactrian(speciation_at_present,weight=2)
moves[++mvi] = mvScaleBactrian(extinction,weight=2)

# Prior on global scale of process (diffusion rate)
speciation_global_scale_hyperprior <- 0.0094
speciation_sd ~ dnHalfCauchy(0,1)
speciation_sd.setValue(runif(1,0.005,0.1)[1])

for (i in 1:(NUM_INTERVALS-1)) {
  index = i+1
  
  # specify normal priors (= Brownian motion) on the differences between log of the rates
  delta_log_speciation[i] ~ dnNormal( mean=0, sd=speciation_sd*speciation_global_scale_hyperprior )
  delta_log_speciation[i].setValue(runif(1,-0.1,0.1)[1])
  
}

# Assemble first-order differences and speciation at present into the random field
speciation := fnassembleContinuousMRF(speciation_at_present,delta_log_speciation,initialValueIsLogScale=FALSE,order=1)

# Move all field parameters in one go
moves[++mvi] = mvEllipticalSliceSamplingSimple(delta_log_speciation,weight=5,tune=FALSE)

# Gibbs sampler on global scale
moves[++mvi] = mvGMRFHyperpriorGibbs(speciation_sd, delta_log_speciation, speciation_global_scale_hyperprior, weight=5)


timetree ~ dnEpisodicBirthDeath(rootAge          = psi.rootAge(), 
lambdaRates                                      = speciation, 
lambdaTimes                                      = interval_times, 
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
monitors[1] = mnModel(filename="output/GMRFBDP_" + REPLICATE + "_run_" + THISRUN + ".log",printgen=200, separator = TAB)


################
# The Analysis #
################

### workspace mcmc ###
mymcmc = mcmc(mymodel, monitors, moves, nruns=1)

### pre-burnin to tune the proposals ###
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
  this.rb <- gsub("REPLICATE",i,gmrf_template)
  for (j in 1:2) {
    this.chain.rb <- gsub("THISRUN",j,this.rb)
    this.chain.rb <- gsub("THISSEED",seed,this.chain.rb)
    for (k in 1:4) {
      cat(this.chain.rb,file=paste0("1_simulation_study/0_setup/analysis_scripts//GMRFBDP_batch_",i,"_chain_",j,".Rev"))
    }
    seed <- seed + 1
  }
}
