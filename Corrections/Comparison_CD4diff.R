library(BASiCS)

mypath <- "~/Dropbox (Cambridge University)//BASiCSplus/Comparisons/CellSystems/"

#### Read in data
input <- read.table(file.path(mypath, "Data/CD4_diff.txt"), sep = "\t")

# Annotate genes
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", 
                      mart = ensembl)
genenames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = ensembl)
rownames(genenames) <- genenames$ensembl_gene_id

# Filter out pseudogenes
#cur_genenames <- genenames[rownames(input),]
#input <- input[!grepl("^Gm", cur_genenames[,2]),]
#cur_genenames <- genenames[rownames(input),]
#input <- input[!grepl("-ps", cur_genenames[,2]),]

ERCC.conc <- read.table(file.path(mypath, "Data/ERCC_malaria.txt"), header=TRUE, sep = "\t")

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,2]*(10^(-18)))*(6.0221417*(10^23))

# Here, we add the well volume
ERCC.num.final <- ERCC.num * 0.009
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,1]

SpikeInput <- ERCC.num.final[rownames(input)[grepl("ERCC", rownames(input))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput, stringsAsFactors = FALSE)

#### Create Data object for each condition
# Day 2
input.day2 <- input[,grepl("2_", colnames(input))]
chips.day2 <- sapply(colnames(input.day2), function(n){unlist(strsplit(n, "\\."))[1]})

Data.2day <- newBASiCS_Data(Counts = input.day2, 
                            Tech = grepl("ERCC", rownames(input.day2)), 
                            SpikeInfo = SpikeInput.1, 
                            BatchInfo=chips.day2)
saveRDS(Data.2day, file.path(mypath, "Data/Data_day2.rds"))

# Day 4
input.day4 <- input[,grepl("4_", colnames(input))]
chips.day4 <- sapply(colnames(input.day4), function(n){unlist(strsplit(n, "\\."))[1]})

Data.4day <- newBASiCS_Data(Counts = input.day4, 
                            Tech = grepl("ERCC", rownames(input.day4)), 
                            SpikeInfo = SpikeInput.1, 
                            BatchInfo=chips.day4)
saveRDS(Data.4day, file.path(mypath, "Data/Data_day4.rds"))

# Day 7
input.day7 <- input[,grepl("7_infect", colnames(input))]
chips.day7 <- sapply(colnames(input.day7), function(n){unlist(strsplit(n, "\\."))[1]})

Data.7day <- newBASiCS_Data(Counts = input.day7, 
                            Tech = grepl("ERCC", rownames(input.day7)), 
                            SpikeInfo = SpikeInput.1, 
                            BatchInfo=chips.day7)
saveRDS(Data.7day, file.path(mypath, "Data/Data_day7.rds"))

chains.path <- file.path(mypath, "/MCMCs/")

# Running the MCMC
MCMC.day2 <- BASiCS_MCMC(Data.2day, N = 100000, Thin = 50, Burn = 50000,
                         Regression = TRUE, WithSpikes = TRUE,
                         StoreChains = TRUE, RunName = "day2_long",
                         StoreDir = chains.path)

MCMC.day4 <- BASiCS_MCMC(Data.4day, N = 100000, Thin = 50, Burn = 50000,
                         Regression = TRUE, WithSpikes = TRUE,
                         StoreChains = TRUE, RunName = "day4_long",
                         StoreDir = chains.path)

MCMC.day7 <- BASiCS_MCMC(Data.7day, N = 100000, Thin = 50, Burn = 50000,
                         Regression = TRUE, WithSpikes = TRUE,
                         StoreChains = TRUE, RunName = "day7_long",
                         StoreDir = chains.path)
