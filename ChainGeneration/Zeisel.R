#!/usr/bin/env Rscript

##########################################################
#### Script to run the model on microglia cells ##########
##########################################################

# The regression and non-regression model is run on microglia cells 
# from Zeisel et al. 
# The script takes the number of GRBFs, their scale 
# parameter and the degrees of freedom as input.

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

#### Zeisel data ####
input.Zeisel <- read.table("Data/Test_Data/microglia_Zeisel.txt", sep = "\t")

ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

ERCC.num.final <- ERCC.num/2500000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input.Zeisel)[grepl("ERCC", rownames(input.Zeisel))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

# Generate Data object
Data.Zeisel <- newBASiCS_Data(Counts = input.Zeisel, 
                              Tech = grepl("ERCC", rownames(input.Zeisel)), 
                              SpikeInfo = SpikeInput.1)

# Run regression model
MCMC.Zeisel <- BASiCS_MCMC(Data = Data.Zeisel, N=40000, Thin = 20, Burn = 20000, 
                           Regression = TRUE, k=k, Var=Var, eta=eta)

saveRDS(MCMC.Zeisel, paste("Tdist/Results/Testing/Datasets/MCMC_Zeisel_", 
                           k, "_", Var, "_", eta, "_reg.rds", sep=""))

# Run non-regression model
MCMC.Zeisel.old <- BASiCS_MCMC(Data = Data.Zeisel, N=40000, Thin = 20, 
                               Burn = 20000, prior = "log-normal")

saveRDS(MCMC.Zeisel.old, paste("Tdist/Results/Testing/Datasets/MCMC_Zeisel_old.rds", 
                               sep=""))
