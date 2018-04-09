################################################################
#### Script to simualte data based on the Zeisel brain data ####
################################################################

library(BASiCS)
# Read in MCMC run on the microglia cells of the Zeisel dataset using the regression
model
MCMC <- readRDS("Google Drive File Stream/My Drive/BASiCS_add-on/Analysis/Results/Testing/Tdist/All_datasets.rds")
Zeisel_MCMC <- MCMC$MCMC_Zeisel_12_1.2_5_reg.rds
rm(MCMC)

# Read in the microglia Zeisel data
Zeisel_data <- read.table("Google Drive File Stream/My Drive/BASiCS_add-on/Data/Test_Data/microglia_Zeisel.txt")

# Save parameters
parameters <- lapply(Zeisel_MCMC@parameters, colMedians)

mu_spikes <-rowMeans(Zeisel_data[grepl("ERCC", rownames(Zeisel_data)),])
rm(Zeisel_data)
rm(Zeisel_MCMC)

save.image("Google Drive File Stream/My Drive/BASiCS_add-on/Analysis/Revisions/Simulations/DataForSimulations.RData")
