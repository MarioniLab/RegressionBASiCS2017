#!/usr/bin/env Rscript

#######################################################
#### Script to run the model on CD4 T cells ###########
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

#### CD4 T cell data ####
input.CD4 <- read.table("Data/Test_Data/CD4_NaiveActiveYoungB6.txt", sep = "\t")
# Select naive CD4 T cells from young B6 
input.CD4 <- input.CD4[,grepl("SS51_naive", colnames(input.CD4)) | grepl("SS52_naive", colnames(input.CD4))]
chips <- sapply(colnames(input.CD4), function(n){unlist(strsplit(n, "_"))[1]})

ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

ERCC.num.final <- ERCC.num/50000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input.CD4)[grepl("ERCC", rownames(input.CD4))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

Data.CD4 <- newBASiCS_Data(Counts = input.CD4, Tech = grepl("ERCC", rownames(input.CD4)), SpikeInfo = SpikeInput.1, BatchInfo = chips)

MCMC.CD4 <- BASiCS_MCMC(Data = Data.CD4, N=40000, Thin = 20, Burn = 20000, prior = "log-normal", Regression = TRUE, k=k, Var=Var, eta=eta)

#saveRDS(MCMC.CD4, paste("Tdist/Results/Testing/Gridsearch/MCMC_CD4_", k, "_", Var, "_", eta, ".rds", sep=""))

saveRDS(MCMC.CD4, paste("Tdist/Results/Testing/Datasets/MCMC_CD4_", k, "_", Var, "_", eta, "_reg.rds", sep=""))

MCMC.CD4.old <- BASiCS_MCMC(Data = Data.CD4, N=40000, Thin = 20, Burn = 20000, prior = "log-normal")

saveRDS(MCMC.CD4.old, paste("Tdist/Results/Testing/Datasets/MCMC_CD4_old.rds", sep=""))



# Do Posterior predictive Checks
#source("/nfs/research2/marioni/Nils/BASiCS/Tdist/Scripts/PostPredTest.R")

#postpred <- PostPredTest(Data.CD4, MCMC.CD4, variance=1.2)

#saveRDS(postpred, paste("Tdist/Results/Testing/PostPred_DF/PPC_CD4_", eta, ".rds", sep=""))

