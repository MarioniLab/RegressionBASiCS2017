#!/usr/bin/env Rscript

###############################################################
#### Script to run the model on Dictyostelium cells ###########
###############################################################

# The regression and non-regression model is run on Dictyostelium cells before 
# differentiation started (Antolović). The script takes the number of GRBFs, their scale 
# parameter and the degrees of freedom as input.

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

#### Dictyostelium data ####
input.dict <- read.table("Data/Test_Data/Dictyostelium.txt", sep = "\t")
# Select first time point
input.dict <- input.dict[,grepl("X0h", colnames(input.dict))]
chips <- sapply(colnames(input.dict), function(n){unlist(strsplit(n, "_"))[1]})

ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))

ERCC.num.final <- ERCC.num/1000

# Here, the volume of the well needs to be incorporated
# In the original publication, we missed to scale the final ERCC counts by a factor of 
# ERCC.num.final <- ERCC.num.final * 0.009

rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input.dict)[grepl("ERCC", rownames(input.dict))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

# Generate Data object			   
Data.dict <- newBASiCS_Data(Counts = input.dict, 
                            Tech = grepl("ERCC", rownames(input.dict)), 
                            SpikeInfo = SpikeInput.1, BatchInfo = chips)

# Run the regression model 
MCMC.dict <- BASiCS_MCMC(Data = Data.dict, N=40000, Thin = 20, Burn = 20000, 
                         Regression = TRUE, k=k, Var=Var, eta=eta)

saveRDS(MCMC.dict, paste("Tdist/Results/Testing/Datasets/MCMC_Dict_", 
                         k, "_", Var, "_", eta, "_reg.rds", sep=""))

# Run the non-regression model 
MCMC.dict.old <- BASiCS_MCMC(Data = Data.dict, N=40000, Thin = 20, Burn = 20000, 
                             prior = "log-normal")

saveRDS(MCMC.dict.old, paste("Tdist/Results/Testing/Datasets/MCMC_Dict_old.rds", 
                             sep=""))


