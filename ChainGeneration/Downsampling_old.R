#!/usr/bin/env Rscript

# This file reads in all MCMC chain objects from subsetted CA1 pyramidal neurons
# and extracts mu and delta parameters. MCMC chains were generated 
# using the non-regression approach.

library(BASiCS)

setwd("/nfs/research2/marioni/Nils/BASiCS/Tdist/Results/Testing/")

list_files = list.files("Downsampling_Zeisel_old/", full.names = TRUE)

first <- readRDS(list_files[1])

output <- matrix(data=NA, ncol=2, nrow=ncol(first@parameters$mu))
output[,1] <- colMedians(first@parameters$mu)
output[,2] <- colMedians(first@parameters$delta)

for(n in list_files[-1]){
	cur_MCMC = readRDS(n)
	cur_mat = matrix(data=NA, ncol=2, nrow=ncol(first@parameters$mu))
	cur_mat[,1] <- colMedians(cur_MCMC@parameters$mu)
	cur_mat[,2] <- colMedians(cur_MCMC@parameters$delta)
	output <- cbind(output, cur_mat)
}

names.down <- list.files("Downsampling_Zeisel_old/", full.names = FALSE)
colnames(output) <- paste(rep(names.down, each=2), c("mu", "delta"), sep="_")

saveRDS(output, "Downsampling_Zeisel_old.rds")

