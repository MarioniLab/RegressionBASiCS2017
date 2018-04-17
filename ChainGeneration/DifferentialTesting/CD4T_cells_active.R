########################################################################
#### Script to run the model on CD4 T cells of different conditions ####
########################################################################

# This script runs the regression model on CD4 T cells after 3h of in vitro activation.
# Data was taken from Martinez et al.

library(BASiCS)
setwd("/nfs/research2/marioni/Nils/BASiCS/")

#### Read in data

input <- read.table("Data/Test_Data/CD4_NaiveActiveYoungB6.txt", sep = "\t")

#### Read in Spike-ins
ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num/50000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input)[grepl("ERCC", rownames(input))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

#### Create Data objects for each condition

# Young active B6
input <- input[,grepl("Active", colnames(input))]
chips <- sapply(colnames(input), function(n){unlist(strsplit(n, "_"))[1]})

Data.YoungActiveB6 <- newBASiCS_Data(Counts = input,
                                     Tech = grepl("ERCC", rownames(input)), 
                                     SpikeInfo = SpikeInput.1, BatchInfo=chips)

#### Run MCMC on these conditions

# Young active B6

MCMC.YoungActiveB6 <- BASiCS_MCMC(Data.YoungActiveB6, 40000, 20, 20000, 
                                  Regression = TRUE, k = 12, Var = 1.2, 
                                  PrintProgress=FALSE)

saveRDS(MCMC.YoungActiveB6, "Tdist/Results/Differential_testing/MCMC_Active_B6.rds")
