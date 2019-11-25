# hsmrfbdp
This repository contains code for analyses from "Locally adaptive Bayesian birth-death model successfully detects slow and rapid rate shifts" (Magee et al. 2019).
In the directory 1_simulation_study is code for simulating a number of trees under a time-dependent birth-death process, inferring diversification rates under the HSMRF-based and GMRF-based models, and processing those results.
In the directory 2_empirical_analyses is code for the empirical analyses presented in the manuscript.
Both directories contain more extensive readme files.

## Software requirements
Running any analyses in this repository requires an installation of [RevBayes](https://revbayes.github.io/), version 1.0.11 or higher.
Pre- and post-processing steps require R 3.4.4 or higher.
All R scripts are designed to be called from the top level directory.
Further requirements are discussed in the relevant readme files.
