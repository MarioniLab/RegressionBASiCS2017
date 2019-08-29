# This script runs the regression model on CD4 T cells 4 days after malaria infection.
# Data was taken from Loennberg et al.

library(BASiCS)
setwd("/nfs/research2/marioni/Nils/BASiCS/")

#### Read in data

input <- read.table("Data/Test_Data/CD4_diff.txt", sep = "\t")

#### Read in Spike-ins
ERCC.conc <- read.table("Data/Test_Data/ERCC_malaria.txt", header=TRUE, sep = "\t")

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,2]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num

# Here, the volume of the well needs to be incorporated
# In the original publication, we missed to scale the final ERCC counts by a factor of 
# ERCC.num.final <- ERCC.num.final * 0.009

rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,1]

SpikeInput <- ERCC.num.final[rownames(input)[grepl("ERCC", rownames(input))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
			   stringsAsFactors = FALSE)


input <- input[,grepl("4_", colnames(input))]
chips <- sapply(colnames(input), function(n){unlist(strsplit(n, "\\."))[1]})

# Generate data object
Data.4day <- newBASiCS_Data(Counts = input, Tech = grepl("ERCC", rownames(input)), 
                            SpikeInfo = SpikeInput.1, BatchInfo=chips)

#### Run MCMC on these condition

MCMC.4day <- BASiCS_MCMC(Data.4day, 40000, 20, 20000, Regression = TRUE, 
                         k = 12, Var = 1.2, PrintProgress=FALSE)

saveRDS(MCMC.4day, "Tdist/Results/Differential_testing/MCMC_CD4diff_4day.rds")
