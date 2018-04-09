#!/usr/bin/env Rscript

################################################
#### Script to simulat data for TPR testing ####
################################################

args = commandArgs(trailingOnly=TRUE)

condition = as.character(args[1])
nocells = as.character(args[2])
iter = as.numeric(args[3])

library(BASiCS)

setwd("/nfs/research1/marioni/Nils/BASiCS/")

load("Data/SimData/DataForSimulations.RData")

if(condition == "B"){
	add.values <- rep(5, 1000)

	# Change the sign for 50% of values
	a <- as.numeric(Sys.time())
	set.seed(a)
	position <- sample(1:1000, 500)
	add.values[position] <- -add.values[position]

	# Turn log2 fold change in fold change
	add.values <- 2^add.values

	# Sample genes to be replaced
	position <- sample(1:length(parameters$delta), 1000)
	delta <- parameters$delta
	delta[position] <- delta[position]*add.values
	saveRDS(delta, paste("Revisions/Simulations/parameters/", "Deltas_B_", as.character(nocells), "_", as.character(iter), ".rds", sep = ""))
} else {
	delta = parameters$delta
}

a <- as.numeric(Sys.time())
set.seed(a)
sam <- sample(1:200, nocells)

phi = c(parameters$phi, parameters$phi, parameters$phi)[sam]
phi = phi*(length(phi)/sum(phi))

# Simulate Data
a <- as.numeric(Sys.time())
set.seed(a)
Test_Data <- BASiCS::BASiCS_Sim(Mu = parameters$mu,
                                Mu_spikes = mu_spikes,
                                Delta = delta,
                                Phi = phi,
                                S =  c(parameters$s, parameters$s, parameters$s)[sam],
                                Theta = parameters$theta)

Counts <- assays(Test_Data)$Counts
rownames(Counts) <- c(rownames(Counts)[1:length(parameters$mu)], names(mu_spikes))

ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num/2500000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[names(mu_spikes),1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

Data <- newBASiCS_Data(Counts = Counts, Tech = c(rep(FALSE, length(parameters$mu)), rep(TRUE, length(mu_spikes))), SpikeInfo = SpikeInput.1)

a <- as.numeric(Sys.time())
set.seed(a)
MCMC <- BASiCS_MCMC(Data = Data, 40000, 20, 20000, PrintProgress=FALSE, Regression=TRUE, PriorDelta="log-normal")

saveRDS(MCMC, paste("Revisions/Simulations/MCMCs/TPR/", "MCMC_", condition, "_", nocells, "_", iter, ".rds", sep=""))

