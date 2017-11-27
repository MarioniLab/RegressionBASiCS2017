##########################################################
#### Script to preprocess datasets used in this study ####
##########################################################

setwd("/Users/eling01/Google Drive/BASiCS_add-on/")

#### CD4 T cell dataset ####
input <- read.table("Data/Raw_data/CD4/QC_B6_CAST_all.txt", header = TRUE, sep = "\t")

# Select young B6 cells
input.CD4 <- input[,which(grepl("SS51", colnames(input)) | grepl("SS52", colnames(input)))]

# Filtering of bological and technical genes
g.bio <- input.CD4[grepl("ENS", rownames(input.CD4)),]
ERCC <- input.CD4[grepl("ERCC", rownames(input.CD4)),]
ERCC <- ERCC[rowSums(ERCC) > 0,]
g.bio <- g.bio[rowMeans(g.bio) > 1,]

write.table(rbind(g.bio, ERCC), "Data/Test_Data/CD4_NaiveActiveYoungB6.txt", sep = "\t")

#### Zeisel data ####
Zeisel.bio <- read.table("Data/Raw_data/Zeisel/Zeisel_data.txt", sep = "\t")
rownames(Zeisel.bio) <- Zeisel.bio[,1]
Zeisel.bio <- Zeisel.bio[,-c(1,2)]
Zeisel.meta <- read.table("Data/Raw_data/Zeisel/Zeisel_meta.txt", sep = "\t", stringsAsFactors = FALSE)
rownames(Zeisel.meta) <- Zeisel.meta$V1
Zeisel.meta <- Zeisel.meta[,-1]
Zeisel.ERCC <- read.table("Data/Raw_data/Zeisel/Zeisel_ERCC.txt", sep = "\t")
rownames(Zeisel.ERCC) <- Zeisel.ERCC[,1]
Zeisel.ERCC <- Zeisel.ERCC[,-c(1,2)]

# Filter biological and technical genes
Zeisel.bio <- Zeisel.bio[rowMeans(Zeisel.bio) > 0.1,]
colnames(Zeisel.bio) <- Zeisel.meta[8,]
colnames(Zeisel.ERCC) <- Zeisel.meta[8,]

# Look at sizes of cell groups 
table(as.character(Zeisel.meta[9,]))

# Filter out one cell type - microglia
microglia <- Zeisel.bio[,Zeisel.meta[2,] == 5]
ERCC <- Zeisel.ERCC[,colnames(microglia)]
ERCC <- ERCC[rowSums(ERCC) > 0,]

write.table(rbind(microglia, ERCC), "Data/Test_Data/microglia_Zeisel.txt", sep = "\t")

# Another cell type for downsampling - CA1
CA1 <- Zeisel.bio[,Zeisel.meta[2,] == 3]
ERCC <- Zeisel.ERCC[,colnames(CA1)]

write.table(rbind(CA1, ERCC), "Data/Test_Data/CA1_Zeisel.txt", sep = "\t")

#### Pool and split data ####
input <- read.table("Data/Raw_data/PoolSplit/GSE54695_data_transcript_counts.txt", sep = "\t", header = TRUE)
rownames(input) <- input$GENENAME
input <- input[,-1]

# Select 2i conditions
input.2i <- input[,grepl("2i", colnames(input))]

# Batch info
Batch <- c(rep(1, 40), rep(2, 40), rep(1, 40), rep(2, 40))

# Generate raw UMI counts
UMICount <- function(MoleculeCount, UMILength)
{
  # MoleculeCount is the normalized count
  M = 4^UMILength
  UMICount = M*(1-exp(-MoleculeCount/M))
  return(UMICount)
}
CountsUMI = round(UMICount(input.2i, 4))

# Quality control - remove low Nanog cells
Pou5f1.per.cell <- as.numeric(CountsUMI["Pou5f1",])

CountsUMI_1 <- CountsUMI[,Pou5f1.per.cell >= 10]
CountsUMI.bio <- CountsUMI_1[!grepl("ERCC", rownames(CountsUMI_1)), ]
ERCC <- CountsUMI_1[grepl("ERCC", rownames(CountsUMI_1)), ]

CountsUMI.bio <- CountsUMI.bio[rowMeans(CountsUMI.bio) > 0.5,]
ERCC <- ERCC[rowMeans(ERCC) > 0,]
Batch <- Batch[Pou5f1.per.cell >= 10]

Data <- rbind(CountsUMI.bio, ERCC)
colnames(Data) <- paste(sapply(colnames(Data), function(n){unlist(strsplit(n, "_"))[1]}), Batch, "_", sapply(colnames(Data), function(n){unlist(strsplit(n, "_"))[3]}), sep = "")

write.table(Data, "Data/Test_Data/PoolSplit.txt", sep = '\t')

#### Dictyostelium data ####
input.Dict <- read.csv("Data/Raw_data/Chubb/mmc2.csv", header = TRUE, skip = 2)
rownames(input.Dict) <- input.Dict[,1]
input.Dict <- input.Dict[,-1]

# QC
input.Dict <- input.Dict[,!grepl("\\...", colnames(input.Dict))]
input.Dict.bio <- input.Dict[!(grepl("ERCC", rownames(input.Dict)) | grepl("__", rownames(input.Dict))),]
ERCC <- input.Dict[grepl("ERCC", rownames(input.Dict)),]

# Look at batch effects
pca <- prcomp(log10(t(input.Dict.bio) + 1))
library(ggplot2)

ggplot(data.frame(PC1=pca$x[,1], PC2=pca$x[,2], 
                  day=sapply(colnames(input.Dict.bio), function(n){substr(n, 1,3)}),
                  batch=sapply(colnames(input.Dict.bio), function(n){substr(n, 4,4)}))) +
  geom_point(aes(x=PC1, y=PC2, colour=day, shape=batch))

# Batch 3 appears to cluster separatlety - remove this one
input.Dict.bio <- input.Dict.bio[,!(grepl("X0h3", colnames(input.Dict.bio)) | 
                                      grepl("X3h3", colnames(input.Dict.bio)) | 
                                      grepl("X6h3", colnames(input.Dict.bio)))]

input.Dict.bio <- input.Dict.bio[rowMeans(input.Dict.bio) > 1,]
ERCC <- ERCC[rowMeans(ERCC) > 0, colnames(input.Dict.bio)]

write.table(rbind(input.Dict.bio, ERCC), "Data/Test_Data/Dictyostelium.txt", sep = "\t")

#### CD4 T cell differentiation ####
input <- read.table("Data/Raw_data/CD4_Diff/Tcell_all.txt", sep = "\t", header = TRUE)
colnames(input) <- sapply(colnames(input), function(n){unlist(strsplit(n, "\\."))[1]})
nam <- data.frame(ERS = as.character(c(paste("ERS", seq(1064048,1064799,1), sep = ""), paste("ERS", seq(1616102,1616947,1), sep = ""))),
                  ERR = as.character(c(paste("ERR", seq(1276858,1277609,1), sep = ""), paste("ERR", seq(1892405,1893250,1), sep = ""))), stringsAsFactors = FALSE)
meta <- read.table("Data/Raw_data/CD4_Diff/Malaria_meta.txt", sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
colnames(input) <- nam$ERS[match(colnames(input), nam$ERR)]

# Select CD4 T cells
meta <- meta[meta$Characteristics.cell.type. == "CD4+ T cell",]

input <- input[,unique(meta$Comment.ENA_SAMPLE.)]

meta <- meta[match(unique(meta$Comment.ENA_SAMPLE.), meta$Comment.ENA_SAMPLE.),]

input <- input[-c((nrow(input)-4):nrow(input)),]

# renames colnames
colnames(input) <- paste(meta$Characteristics.time., ifelse(meta$Characteristics.infect. != 'none', "infect", "none"), sapply(meta$Source.Name, function(n){unlist(strsplit(n, "-"))[3]}), sep = "_")

# Filter cells
day_0 <- input[,grepl("0_", colnames(input))]
day_2 <- input[,grepl("2_", colnames(input))]
day_3 <- input[,grepl("3_", colnames(input))]
day_4 <- input[,grepl("4_", colnames(input))]
day_7 <- input[,grepl("7_none", colnames(input))]
day_7_infect <- input[,grepl("7_infect", colnames(input))]

plot(log10(colSums(day_0[grepl("ENS", rownames(day_0)),])), 
     log10(colSums(day_0[grepl("ERCC", rownames(input)),])), pch=16)

plot(log10(colSums(day_2[grepl("ENS", rownames(day_2)),])), 
     log10(colSums(day_2[grepl("ERCC", rownames(input)),])), pch=16)

day_2 <- day_2[,log10(colSums(day_2[grepl("ENS", rownames(day_2)),])) > 5.5 &
                 log10(colSums(day_2[grepl("ERCC", rownames(day_2)),])) > 5.8]

plot(log10(colSums(day_3[grepl("ENS", rownames(day_3)),])), 
     log10(colSums(day_3[grepl("ERCC", rownames(input)),])), pch=16)

plot(log10(colSums(day_4[grepl("ENS", rownames(day_4)),])), 
     log10(colSums(day_4[grepl("ERCC", rownames(input)),])), pch=16)

day_4 <- day_4[,log10(colSums(day_4[grepl("ENS", rownames(day_4)),])) > 6]

plot(log10(colSums(day_7[grepl("ENS", rownames(day_7)),])), 
     log10(colSums(day_7[grepl("ERCC", rownames(input)),])), pch=16)

plot(log10(colSums(day_7_infect[grepl("ENS", rownames(day_7_infect)),])), 
     log10(colSums(day_7_infect[grepl("ERCC", rownames(input)),])), pch=16)

day_7_infect <- day_7_infect[,log10(colSums(day_7_infect[grepl("ENS", rownames(day_7_infect)),])) > 6]

input <- cbind(day_0, day_2, day_3, day_4, day_7, day_7_infect)

# Filter genes 
genes <- rownames(input)[rowMeans(day_2) > 1 & rowMeans(day_3) > 1 & rowMeans(day_4) > 1 & rowMeans(day_7_infect) > 1]
input <- input[genes,]

# Test if MCMC works
ERCC.conc <- read.table("Data/Raw_data/CD4_Diff/ERCC_malaria.txt", header=TRUE, sep = "\t")

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,2]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,1]

SpikeInput <- ERCC.num.final[rownames(input)[grepl("ERCC", rownames(input))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

Data <- newBASiCS_Data(Counts = input[,grepl("2_", colnames(input))],
                            Tech = grepl("ERCC", rownames(input)), SpikeInfo = SpikeInput.1)

#### Run MCMC on these conditions

MCMC <- BASiCS_MCMC(Data, 40000, 20, 20000, Regression = TRUE, k = 12, Var = 1.2)


write.table(input, "Data/Test_Data/CD4_diff.txt", sep = "\t")

#### Simulated data ####

# Simulate data for testing of one population

library(BASiCSplus)

# Select naive T cells from young and old B6 
input <- read.table("Data/Raw_data/CD4/QC_B6_CAST_all.txt", header = TRUE, sep = "\t")
input <-  input[rowMeans(input) > 1,]

input <- input[,which(grepl("SS51_naive", colnames(input)) | grepl("SS52_naive", colnames(input)) | 
                        grepl("SS63_naive", colnames(input)) | grepl("SS64_naive", colnames(input)))]

length(which(rowSums(input) == 0))

# Run BASiCS with and without regression
ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num/50000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input)[grepl("ERCC", rownames(input))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

Data <- newBASiCS_Data(Counts = input, Tech = grepl("ERCC", rownames(input)), SpikeInfo = SpikeInput.1)

# Without regression
MCMC <- BASiCS_MCMC(Data = Data, N=40000, Thin = 20, Burn = 20000, prior = "log-normal")

rm(list=setdiff(ls(), c("Data", "MCMC")))
save.image("Analysis/Results/Sim_data/Sample_MCMCs.RData")
load("/Users/eling01/Google Drive/BASiCS_add-on/Analysis/Results/Sim_data/Sample_MCMCs.RData")

# Simulate data from these runs
Sim.Data.old <- BASiCSplus::BASiCS_Sim(mu = colMedians(MCMC@parameters$mu), mu_spikes = metadata(Data)$SpikeInput,
                                       delta = colMedians(MCMC@parameters$delta), phi = colMedians(MCMC@parameters$phi), s = colMedians(MCMC@parameters$s), theta = median(MCMC@parameters$theta))

rownames(Sim.Data.old)[grepl("ERCC", rownames(input))] <- rownames(input)[grepl("ERCC", rownames(input))]

write.table(assay(Sim.Data.old), "Data/Test_Data/Sim_Data_old.txt", sep = "\t")
saveRDS(MCMC, "Analysis/Results/Sim_data/Sim_MCMC.rds")
saveRDS(Data, "Analysis/Results/Sim_data/Sim_Data.rds")

# Simulated data for 2 conditions

library(BASiCSplus)

# Select naive and activated T cells from young B6 

input <- read.table("Data/Raw_data/CD4/QC_B6_CAST_all.txt", header = TRUE, sep = "\t")

naive <- input[,which(grepl("SS51_naive", colnames(input)) | grepl("SS52_naive", colnames(input)) |
                        grepl("SS63_naive", colnames(input)) | grepl("SS64_naive", colnames(input)))]
set.seed(123)
naive <- naive[,sample(1:ncol(naive), 100)]

active <- input[,which(grepl("SS51_active", colnames(input)) | grepl("SS52_active", colnames(input)) |
                        grepl("SS63_active", colnames(input)) | grepl("SS64_active", colnames(input)))]
set.seed(123)
active <- active[,sample(1:ncol(naive), 100)]

input <- cbind(naive, active)
input <- input[rowMeans(input) > 1,]

# Read in the spike-in information

ERCC.conc <- read.table("Data/ERCC_conc.txt", header=TRUE, sep = "\t", fill = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.final <- ERCC.num/50000
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(input)[grepl("ERCC", rownames(input))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

# Create the SummarizedExperiment object

Data.naive <- newBASiCS_Data(Counts = input[,grepl("naive", colnames(input))], 
                             Tech = grepl("ERCC", rownames(input)), SpikeInfo = SpikeInput.1) 

Data.active <- newBASiCS_Data(Counts = input[,grepl("active", colnames(input))], 
                             Tech = grepl("ERCC", rownames(input)), SpikeInfo = SpikeInput.1) 

# Run the MCMC on these two conditions

MCMC.naive <- BASiCS_MCMC(Data.naive, 40000, 20, 20000, Regression = TRUE, k=12, Var=1.2, eta=5)

MCMC.active <- BASiCS_MCMC(Data.active, 40000, 20, 20000, Regression = TRUE, k=12, Var=1.2, eta=5)

# Correct for offset 

OffSetCorrection <- function(MCMC1, MCMC2){
  median(rowSums(MCMC1@parameters$mu)/rowSums(MCMC2@parameters$mu)) 
}

OffSet <- OffSetCorrection(MCMC.naive, MCMC.active)

MCMC.naive@parameters$mu <- MCMC.naive@parameters$mu/OffSet
MCMC.naive@parameters$phi <- MCMC.naive@parameters$phi*OffSet

rm(list=setdiff(ls(), c("Data.naive", "Data.active", "MCMC.naive", "MCMC.active", "SpikeInput.1")))

save.image("Analysis/Results/Sim_Data/SimForDiffTesting.RData")
load("Analysis/Results/Sim_Data/SimForDiffTesting.RData")
# Simulate data from these runs
Sim.Data.A <- BASiCSplus::BASiCS_Sim(mu = colMedians(MCMC.naive@parameters$mu), mu_spikes = metadata(Data.naive)$SpikeInput,
                                       delta = colMedians(MCMC.naive@parameters$delta), 
                                     phi = colMedians(MCMC.active@parameters$phi), 
                                     s = colMedians(MCMC.naive@parameters$s), 
                                     theta = median(MCMC.naive@parameters$theta))

rownames(Sim.Data.A)[grepl("ERCC", rownames(input))] <- rownames(input)[grepl("ERCC", rownames(input))]

Sim.Data.B <- BASiCSplus::BASiCS_Sim(mu = colMedians(MCMC.active@parameters$mu), mu_spikes = metadata(Data.active)$SpikeInput,
                                     delta = colMedians(MCMC.active@parameters$delta), 
                                     phi = colMedians(MCMC.active@parameters$phi), 
                                     s = colMedians(MCMC.active@parameters$s), 
                                     theta = median(MCMC.active@parameters$theta))

rownames(Sim.Data.B)[grepl("ERCC", rownames(input))] <- rownames(input)[grepl("ERCC", rownames(input))]

write.table(assay(Sim.Data.A), "Data/Test_Data/Sim_Data_A.txt", sep = "\t")
write.table(assay(Sim.Data.B), "Data/Test_Data/Sim_Data_B.txt", sep = "\t")


