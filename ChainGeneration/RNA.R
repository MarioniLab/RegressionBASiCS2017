#!/usr/bin/env Rscript

#######################################################
#### Script to run the model on split RNA ###########
#######################################################

setwd("/nfs/research2/marioni/Nils/BASiCS/")

args = commandArgs(trailingOnly=TRUE)

k = as.numeric(args[1])

Var = as.numeric(args[2])

eta = as.numeric(args[3])

library(BASiCS)

# Read in ERCCs 
ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

#### Pool and split data ####
input.ps <- read.table("Data/Test_Data/PoolSplit.txt", sep = "\t")
# Select split condition 
input.ps <- input.ps[,grepl("RNA", colnames(input.ps))]
chips <- sapply(colnames(input.ps), function(n){unlist(strsplit(n, "_"))[1]})

ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

ERCC.num.final <- ERCC.num/2500000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input.ps)[grepl("ERCC", rownames(input.ps))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

Data.ps <- newBASiCS_Data(Counts = input.ps, Tech = grepl("ERCC", rownames(input.ps)), SpikeInfo = SpikeInput.1, BatchInfo = chips)

MCMC.ps <- BASiCS_MCMC(Data = Data.ps, N=40000, Thin = 20, Burn = 20000, prior = "log-normal", Regression = TRUE, k=k, Var=Var, eta=eta)

#saveRDS(MCMC.ps, paste("Tdist/Results/Testing/Gridsearch/MCMC_RNA_", k, "_", Var, "_", eta, ".rds", sep=""))

saveRDS(MCMC.ps, paste("Tdist/Results/Testing/Datasets/MCMC_RNA_", k, "_", Var, "_", eta, "_reg.rds", sep=""))

MCMC.ps.old <- BASiCS_MCMC(Data = Data.ps, N=40000, Thin = 20, Burn = 20000, prior = "log-normal")

saveRDS(MCMC.ps.old, paste("Tdist/Results/Testing/Datasets/MCMC_RNA_old.rds", sep=""))


# Do Posterior predictive Checks
#source("/nfs/research2/marioni/Nils/BASiCS/Tdist/Scripts/PostPredTest.R")

#postpred <- PostPredTest(Data.ps, MCMC.ps, variance=1.2)

#saveRDS(postpred, paste("Tdist/Results/Testing/PostPred_DF/PPC_ps_", eta, ".rds", sep=""))


