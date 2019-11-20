# Empirical analyses
This directory This contains code for running full inference of phylogeny and diversification through time on the macroevolutionary and phylodynamic datasets.
There is also code for running fixed-tree analyses over a variety of grid sizes, and an XML file for running BEAST on the phylodynamic dataset.
Note that many of these scripts depend on being able to source files in /hsmrfbdp/1_simulation_study, so changing the top-level directory structure is not advised.

## Software
The following R packages are required:

  - latex2exp (version 0.4.0 or higher)
  - phangorn (version 2.5.5 or higher)
  - treedater (version 0.3.0 or higher)
  - TreeSim (version 2.3 or higher)

To reproduce the [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html) analysis, RAxML version 8.2.12 or higher is required (see the RAxML_info file for more details).
To run partitioning on Pygopodidae dataset, PartitionFinder version 1.1.1 (or higher, but not version 2.X) is required.
Running BEAST requires [BEAST2](https://www.beast2.org/) version 2.6.0 or higher and the BEAST package BDSKY version 1.4.5 or higher.

## Running analyses
All Rev analyses will be placed into the folder /hsmrfbdp/2_empirical_analyses/output.
The first analysis run will create this folder.
A detailed description of the files and filestructure is available at the end of this document.

#### Analyses of Pygopodidae and HIV
Set your working directory to /hsmrfbdp/2_empirical_analyses.
Use RevBayes to run all 8 master RevScripts in the relevant /src subdirectory, these are scripts named `run_analysis_<n>_<pygo|env>.Rev`.
For example, if the RevBayes executable is in your $PATH, running the 8th analysis for Pygopodidae would look like `rb 1_Pygopodidae/src/run_analysis_1_pygo.Rev`.
The first 4 analyses are the independent replicates of the GMRF-based model analyses and the last four are the independent replicates of the HSMRF-based analyses.
The replicates are split across separate run scripts to enable running in parallel as they each take several days.

#### Sensitivity to grid size
Set your working directory to /hsmrfbdp/2_empirical_analyses.
Use RevBayes to run all RevScripts (`.Rev` files) in /src.
For example, if the RevBayes executable is in your $PATH, running the GMRF-based model analysis using a grid size of 10 would look like `rb 3_grid_size/src/GMRFBDP_10.Rev`
These scripts run relatively quickly and so the multiple replicate chains per analysis are handled at the Rev level.

#### BEAST analysis
Use BEAST2 to run the XML file /4_BEAST/src/env_SubtA_BDSKY.xml.
We recommend using the BEAST GUI.
This will produce a single run, to obtain additional runs, edit the xml to change the name suffix for the parameter and tree logfiles from `_run_1.<log|trees>` and call BEAST again.
The downstream scripts assume there are 4 chains.
Each analysis takes several days.
The location of these output files will depend on how BEAST was run, the downstream scripts require you to place all BEAST output in /hsmrfbdp/2_empirical_analyses/output.

## Diagnosing MCMC convergence
This works the same for all of the analyses.
Each directory has a script named `convergence_diagnostics.R`.
Open R, make sure the working directory in R is /hsmrfbdp, and run all lines to get summaries of the maximum rank-PSRF (across all model parameters) and the minimum rank-ESS (across all model parameters).
More detailed summaries can be seen be more thoroughly examining the output of the calls to `diagnoseConvergence()` and `rankESS()`.

## Making plots
In each directory there is an R script named `plot_<something>.R`.
To make the plots, open this script in R and run all code through the call to `dev.off()`.
Code after this point is for computing summaries of the trajectories such as Bayes Factors.
Plots will be made as PDFs in the subdirectory for the given analysis.
Alternately, the plot function can be called with `Rscript`, and the summaries will be harmlessly printed to stdout.
In either case, the working directory must be /hsmrfbdp.

Before plotting Re for the BEAST analyses, a combined treefile must be created.
Use the logcombiner program to combine the trees from all separate BEAST runs into the file /hsmrfbdp/2_empirical_analyses/output/HIV_BEAST_combined_trees.trees.
When doing this, we recommend downsampling to every 400000th tree.
The trees are only used for the branching time heatmaps, which does not require a huge number of tree samples.


## Directory structure
    /1_Pygopodidae This contains everything for running the macroevolutionary example.
        /data This contains the Rev input files.
            /PartitionFinder This contains the inputs to and outputs from PartitionFinder.
            /NodeCalibration.Rev This is a file of fossil calibrations.
            /Pygopodidae_subset_<1,2>.nex These are the two data subsets identified by PartitionFinder.
            /pygo_starting_tree.tre This is the Pygopodidae subtree from MCC tree of Brennan and Oliver (2017), to ensure tree is within node age constraints.
        /src This contains the key Rev scripts.
            /run_analysis<1-8>.Rev These are the master scripts to run 4 HSMRF and 4 GMRF analyses.
            /Pygopodidae.Rev This is the top-level script that reads in all model components. It does not specify a seed, a replicate ID, or the tree prior.
            /HSMRFBDP.Rev This script sets up the HSMRF tree prior.
            /HSMRFBDP.Rev This script sets up the GMRF tree prior.
            /clock_model.Rev This script sets up the clock model
            /substitution_model.Rev This script sets up the substitution model.
      convergence_diagnostics.R This script is used to calculate rank-based PSRF and ESS for all model parameters excluding branch rate parameters (due to the changing tree topology convergence here cannot be meaningfully addressed).
      plot_spn.R this script will create speciation-through-time plots using the Rev output and computes Bayes Factors for the significance of shifts.
    /2_HIV This contains everything for running the phylodynamic example.
        /data This contains the Rev input files.
            /raxml This contains the inputs to and outputs from running RAxML.
            /NodeCalibration.Rev This is a file of fossil calibrations.
            /env_subtA.fas This is the multiple sequence alignment.
            /env_subtA_taxa.txt This is a Rev input file specifying the ages (in time since 2011) of all the taxa in the tree.
            /treedater.tre This is the RAxML tree made time-calibrated with treedater and used as a starting tree.
        /src This contains the key Rev scripts.
            /run_analysis<1-8>.Rev These are the master scripts to run 4 HSMRF and 4 GMRF analyses.
            /Pygopodidae.Rev This is the top-level script that reads in all model components. It does not specify a seed, a replicate ID, or the tree prior.
            /HSMRFBDP.Rev This script sets up the HSMRF tree prior.
            /HSMRFBDP.Rev This script sets up the GMRF tree prior.
            /clock_model.Rev This script sets up the clock model
            /substitution_model.Rev This script sets up the substitution model.
      convergence_diagnostics.R This script is used to calculate rank-based PSRF and ESS for all model parameters excluding branch rate parameters (due to the changing tree topology convergence here cannot be meaningfully addressed).
      plot_Re.R This script will create Re-through-time plots using the Rev output and computes Bayes Factors for the significance of shifts.
      preliminary_analyses.R This script runs treedater on the RAxML tree and uses that to determine the prior median on the birth rate at the present.

    /3_grid_size This contains everything for assessing the sensitivity of the model to the grid size.
        /data/pygopodidae_hsmrf_run_4_tree_500000.tre This is the last tree sampled in the posterior from run 4 of the Pygopodidae analyses using the HSMRF-based model.
        /src Contains all scripts used to analyze the Pygopodidae tree.
            GMRFBDP_<n>.Rev This script runs a GMRF-based model on the fixed tree using a grid of n cells.
            HSMRFBDP_<n>.Rev This script runs a HSMRF-based model on the fixed tree using a grid of n cells.
        calc_zetas.R This script is used to calculate the global scale parameters for all the different grid sizes used.
        convergence_diagnostics.R This script is used to calculate rank-based PSRF and ESS for all model parameters excluding branch rate parameters (due to the changing tree topology convergence here cannot be meaningfully addressed).

    /4_BEAST This contains BEAST code for analyzing the phylodynamic dataset.
        /src/env_SubtA_BDSKY.xml This is the XML used to run an analysis. To run multiple analyses, change the name suffix of the logfile and treefile.
        convergence_diagnostics.R This script is used to calculate rank-based PSRF and ESS for all model parameters excluding branch rate parameters (due to the changing tree topology convergence here cannot be meaningfully addressed).
        plot_Re_beast.R this script will create Re-through-time plots using the BEAST output and computes Bayes Factors for the significance of shifts.
