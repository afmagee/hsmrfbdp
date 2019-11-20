# Empirical analyses
This directory This contains code for running full inference of phylogeny and diversification through time on the macroevolutionary and phylodynamic datasets.
There is also code for running fixed-tree analyses over a variety of grid sizes, and an XML file for running BEAST on the phylodynamic dataset.

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


## Software
To make any plots, the R package latex2exp (version 0.4.0 or higher) is required.
To reproduce preliminary analyses of the HIV dataset, the R packages phangorn (version 2.5.5 or higher) and treedater (version 0.3.0 or higher) are necessary.
To reproduce the [RAxML ](https://cme.h-its.org/exelixis/web/software/raxml/index.html) analysis, RAxML version 8.2.12 or higher is required (see the RAxML_info file for more details).
To run partitioning on Pygopodidae dataset, PartitionFinder version 1.1.1 (or higher, but not version 2.X) is required.
Running BEAST requires [BEAST2](https://www.beast2.org/) version 2.6.0 or higher and the BEAST package BDSKY version 1.4.5 or higher.
