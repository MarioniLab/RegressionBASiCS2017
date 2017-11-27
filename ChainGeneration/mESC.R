#!/usr/bin/env Rscript

#######################################################
#### Script to run the model on mESC cells ###########
#######################################################

setwd("/nfs/research2/marioni/Nils/BASiCS/")

library(BASiCS)

args = commandArgs(trailingOnly=TRUE)

k = as.numeric(args[1])

Var = as.numeric(args[2])

eta = as.numeric(args[3])

# Read in ERCCs 
ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

#### mESC data ####
input.mESC <- read.table("Data/Test_Data/mESC.txt", sep = "\t")
# Select lif condition
input.mESC <- input.mESC[,grepl("_lif_", colnames(input.mESC))]

ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

ERCC.num.final <- ERCC.num/100
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input.mESC)[grepl("ERCC", rownames(input.mESC))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)
	
Data.mESC <- newBASiCS_Data(Counts = input.mESC, Tech = grepl("ERCC", rownames(input.mESC)), SpikeInfo = SpikeInput.1)

MCMC.mESC <- BASiCS_MCMC(Data = Data.mESC, N=40000, Thin = 20, Burn = 20000, prior = "log-normal", Regression = TRUE, k=k, Var=Var, eta=eta)

#saveRDS(MCMC.mESC, paste("Tdist/Results/Testing/Gridsearch/MCMC_mESC_", k, "_", Var, "_", eta, ".rds", sep=""))

saveRDS(MCMC.mESC, paste("Tdist/Results/Testing/Datasets/MCMC_mESC_", k, "_", Var, "_", eta, "_reg.rds", sep=""))

MCMC.mESC.old <- BASiCS_MCMC(Data = Data.mESC, N=40000, Thin = 20, Burn = 20000, prior = "log-normal")

saveRDS(MCMC.mESC.old, paste("Tdist/Results/Testing/Datasets/MCMC_mESC_old.rds", sep=""))


# Do Posterior predictive Checks
#source("/nfs/research2/marioni/Nils/BASiCS/Tdist/Scripts/PostPredTest.R")

#postpred <- PostPredTest(Data.mESC, MCMC.mESC, variance=1.2)

#saveRDS(postpred, paste("Tdist/Results/Testing/PostPred_DF/PPC_mESC_", eta, ".rds", sep=""))
