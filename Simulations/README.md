# Simulations

This folder contains scripts to perform simulations under the null and under
the alternative. Parameters for simulation were taken from the MCMC run on the
microglia cells of the Zeisel dataset (Simulations.R). Simulations_FPR.R and
Simulations_TPR.R are cluster scripts to start the MCMCs for different sample 
sizes, iterations and conditions. FPR.R and TPR.R perform differential testing 
and compute the false positive rate (FPR) and true positive rate (TPR) as well 
as the expected false discovery rate (EFDR).