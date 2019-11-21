# Simulations

## Software
Simulating trees requires the R package TESS (version 3.0.0 or higher).
Using the R script to manage all analyses of simulated data requires the R package parallel (version 3.4.4 or higher, this comes with the base R installation).
Post-processing the simulation results requires the R package coda (version 0.19-3 or higher).

## Running a simulation for a single grid cell
The following section describes the process of going from a desired birth-rate function through time to obtaining a set of summary statistics describing the behavior of the MRF models.
The steps must be performed in the order they are written here.
Plotting and summarizing the summary statistics is left to the user.

#### Simulating trees
All simulations are from a piecewise-linear model, where the rate is lambda_1 from the root until some time t_1, lambda_2 from t_2 until the present, and a linear interpolation between t_1 and t_2.
This includes as a special case piecewise-constant models (t1=t2), constant-rate models (lambda_1 = lambda_2), and fully linear models (t1 = tree age, t2 = present = 0).

The following is largely a blow-by-blow walkthrough of the script /1_simulation_study/simulation/simulate_trees.R.

1. Open R and make sure the working directory is /hsmrfbdp.

1. Open the script /1_simulation_study/simulation/simulate_trees.R.
Choose a tree age (`tree.age`), sampling fraction at present (`rho`), number of simulations (`n.sims`), a target number of taxa per tree (`ntaxa`), and the extinction rate (`mu`).
The currently specified values are those we used in the study.
Changing the number of trees is possible, but will require changing other files (such as the script to generate the analysis scripts) where the value 100 is hard-coded.

1. Choose and set a seed.

1. Set the simulation (sub)directory `sim.dir` where all data and analyses will go. This is currently set up to simulate to /hsmrf/1_simulation_study/1_constant_rate.

1. Create the appropriate subdirectories within `sim.dir` (/data, /output, /summaries)

1. Use the function `findParams()` to find simulating parameters lambda_1 and lambda_2 that will produce an expected tree size (`ntaxa`), given a desired shift location (`t.shift`), a desired shift radius (`shift.radius`, the radius is one half of the duration), a desired fold change (`fold.change`), and the extinction rate (`extinction`), sampling fraction (`sampling.fraction`), and tree age (`tree.age`).
The function returns a named vector, "rate.old" is the (constant) rate in the oldest part of the tree, "rate.new" is the (constant) rate in the newest portion of the tree.

1. Use these lambda values and the function `speciationRate()` to make a continuous function TESS can use to simulate trees.

1. Use the code block labeled `# Simulate and write trees` to simulate the trees and write them to the /data subdirectory.

The script is currently set up to simulate constant-rate trees.
To simulate time-varying trees, ensure that `shift.magnitude` in `findParams()` is not equal to 1, and use the following guidelines.
Note that the parameter `shift.magnitude` is equal to lambda_1/lambda_2, meaining it is the fold change from past to present.
- To simulate an _increase_ in speciation, set `shift.magnitude` > 1.
- To simulate an _decrease_ in speciation, set 0 < `shift.magnitude` < 1.
- To simulate piecewise-constant models, keep `shift.radius` set to 0.
- To simulate piecewise-linear models, set `shift.radius` larger than 0.
- To simulate one-piece linear models, set `shift.radius` to be `tree.age`/2.

#### Running MRF analyses
Before running analyses using the MRF models, you must populate the directory /0_setup/analysis_scripts with the analysis scripts.
To do this, open R, ensure the working directory is /hsmrfbdp, and run 0_setup/generate_analysis_scripts/generate_GMRFBDP_analyses.R and 0_setup/generate_analysis_scripts/generate_HSMRFBDP_analyses.R.
To optionally be able to compare to the results of a constant-rate model, generate the constant-rate analysis scripts using 0_setup/generate_analysis_scripts/generate_CRBDP_analyses.R.
Any desired changes to the RevBayes analyses, such as changing the grid size, zeta, or the length of the MCMC can be accomplished simply by editing the template scripts in these files before generating the analysis scripts.

Once /analysis_scripts is populated, use the script 0_setup/automation/run_simulation_cell.R to automate running every analysis script in 0_setup/analysis_scripts.
This script requires specifying a path to the RevBayes executable (default assumes it is in the $PATH) and a number of cores to use with parallel::mclapply (default 4).
From the directory for this grid cell, call `Rscript run_simulation_cell.R` and all analyses will be run.
Individual analyses in a grid cell may be run by setting the working directory to that cell and calling RevBayes on a script in /analysis_scripts.
For example, if RevBayes is in you $PATH, you could run the second replicate of the HSMRF-based analysis of the first simulated dataset with the command `rb ../../0_setup/analysis_scripts/HSMRFBDP_batch_1_chain_1.Rev`.

#### Summarizing parameter estimates and diagnosing convergence
Open R and ensure the working directory is /hsmrfbdp.
Open the script /0_setup/summarization/summarize_convergence_and_estimates.R and ensure that `this.dir` is pointing to the correct simulation grid cell (currently its set to /1_constant_rate) and the `n.sims` is the number of simulations performed (currently set to 100).
You can run all lines from inside R or if you prefer call it with `Rscript`.
Running this script will summarize convergence diagnostics (rank-ESS and rank-PSRF) for both the GMRF-based and HSMRF-based models for all model parameters.
It additionally calculates the posterior mean and median for all parameters, as well as the posterior 2.5%, 5%, 95%, and 97.5% quantiles.
For each summary and each model (GMRF or HSMRF), a table is written in /summaries where row i is the analysis of simulation i and each column is the summary for that parameter.
Columns and rows are labeled.

As not all grid cells will use the CRBDP analysis scripts, summarizing convergence and parameter estimates for these analyses is done with the script /0_setup/summarization/summarize_convergence_and_estimates_constant_rate.R, which works as above.

#### Computing summary statistics
To compute all the summaries mentioned in the paper, use the script /0_setup/summarization/calculate_performance_measures.R.
As with summarizing parameter estimates and computing convergence diagnostics, the script must be pointed at the correct directory `this.dir`.
To compute measures of accuracy, the script will need to be able to compute the true values, and thus the simulating parameters must be specified as in simulating the trees.
The number of grid cells, `n.windows`, is also needed as an argument.
Unless any of these have been changed, the current values will suffice to run the script.
Once the values are correct, run the script, either in R or with Rscript.

For each summary statistic, two text files will be added to the /summaries directory (one for the HSMRF-based model and one for the GMRF-based model).
These files will not contain all 100 values, they will only contain summaries for analyses which pass convergence checks (PSRF <= 1.01).
In practice, these files will contain 97-100 values as most analyses pass convergence checks.

As not all grid cells will use the CRBDP analysis scripts, computing summary statistics for these analyses is done with the script /0_setup/summarization/calculate_performance_measures_constant_rate.R, which works as above.

## Directory structure
    /0_setup
        /automation/run_simulation_cell.R This script automates running RevBayes on a simulation grid cell.
        /analysis_scripts This directory will contain Rev scripts to perform 2 independent replicate analyses of all 100 simulated trees.
        /generate_analysis_scripts This directory contains R scripts that populate /analysis_scripts.
            /generate_GMRFBDP_analyses.R This script contains a template Rev script for running analyses using the GMRF-based model and a loop to generate 2 replicate analyses for each of 100 simulated datasets.
            /generate_HSMRFBDP_analyses.R This script contains a template Rev script for running analyses using the HSMRF-based model and a loop to generate 2 replicate analyses for each of 100 simulated datasets.
            /generate_CRBDP_analyses.R This script contains a template Rev script for running constant-rate analyses and a loop to generate 2 replicate analyses for each of 100 simulated datasets.
        /simulation This directory contains scripts for simulating trees and MRF trajectories, and calculating the global scale parameter.
            /simulate_trees.R This script will simulate 100 trees for a user-specified birth-rate-through-time function.
            /set_zeta.R This script contains R functions for calculating zeta, the global scale parameter, for both GMRF-based and HSMRF-based models.
            /helper_functions.R This script contains functions for finding simulation parameters given the location, duration, and size of a shift and the target number of tips in the tree.
            /zeta_helper_functions.R This script contains helper functions for calculating zeta and functions for simulating from MRF models.
        /summarization This directory contains R scripts for summarizing the analyses in a grid cell and diagnosing MCMC convergence.
            /summarize_convergence_and_estimates.R This script computes convergence diagnostics and summarizes parameter estimates for all HSMRF and GMRF analyses in a grid cell.
            /calculate_performance_measures.R This script computes summary statistics of performance for all HSMRF and GMRF analyses that passes convergence checks.
            /summarize_convergence_and_estimates_constant_rate.R This script computes convergence diagnostics and summarizes parameter estimates for all analyses performed using a constant-rate BDP model.
            /calculate_performance_measures_constant_rate.R This script computes summary statistics of performance for all constant-rate BDP analyses that passes convergence checks.
            /rank_based_convergence_diagnostics.R This file contains R functions to compute rank-based PSRF and ESS as MCMC convergence diagnostics.
