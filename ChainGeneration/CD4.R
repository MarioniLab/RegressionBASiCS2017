#!/usr/bin/env Rscript

#######################################################
#### Script to run the model on CD4 T cells ###########
#######################################################

# The regression and non-regression model is run on naive CD4 T cells from
# B6 animals (Martinez et al). The script takes the number of GRBFs, their scale parameter
# and the degrees of freedom as input.

setwd("/nfs/research2/marioni/Nils/BASiCS/")

args = commandArgs(trailingOnly=TRUE)

k = as.numeric(args[1])

Var = as.numeric(args[2])

eta = as.numeric(args[3])

library(BASiCS)

# Read in ERCCs and calculate the number of spike-in molecules
ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", 
                        fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

#### CD4 T cell data ####
input.CD4 <- read.table("Data/Test_Data/CD4_NaiveActiveYoungB6.txt", sep = "\t")
# Select naive CD4 T cells from young B6 
input.CD4 <- input.CD4[,grepl("Unstimulated", colnames(input.CD4))]
chips <- sapply(colnames(input.CD4), function(n){unlist(strsplit(n, "_"))[1]})

ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", 
                        fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

ERCC.num.final <- ERCC.num/50000

# Here, the volume of the well needs to be incorporated
# In the original publication, we missed to scale the final ERCC counts by a factor of 
# ERCC.num.final <- ERCC.num.final * 0.009

rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input.CD4)[grepl("ERCC", rownames(input.CD4))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

# Generate the Data object
Data.CD4 <- newBASiCS_Data(Counts = input.CD4, 
                           Tech = grepl("ERCC", rownames(input.CD4)), 
                           SpikeInfo = SpikeInput.1, BatchInfo = chips)

# Run the regression model
MCMC.CD4 <- BASiCS_MCMC(Data = Data.CD4, 
                        N=40000, Thin = 20, Burn = 20000,
                        Regression = TRUE, k=k, Var=Var, eta=eta)

saveRDS(MCMC.CD4, paste("Tdist/Results/Testing/Datasets/MCMC_CD4_", 
                        k, "_", Var, "_", eta, "_reg.rds", sep=""))

# Run the non-regression model
MCMC.CD4.old <- BASiCS_MCMC(Data = Data.CD4, N=40000, Thin = 20, Burn = 20000, 
                            prior = "log-normal")

saveRDS(MCMC.CD4.old, paste("Tdist/Results/Testing/Datasets/MCMC_CD4_old.rds", sep=""))

