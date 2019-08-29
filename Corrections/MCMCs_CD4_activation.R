library(BASiCS)

mypath <- "~/Dropbox (Cambridge University)//BASiCSplus/Comparisons/CellSystems/"

CD4.cells <- read.table(file.path(mypath, "/Data/CD4_NaiveActiveYoungB6.txt"), 
                        sep = "\t", header = TRUE)

#### Read in Spike-ins
website <- "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/"
file <- "cms_095046.txt"
ERCC.conc <- read.table(url(paste(website, file, sep = "")), 
                        sep = "\t", header = TRUE)

ERCC.num <- matrix(data=NA, nrow=nrow(ERCC.conc), ncol=1)
ERCC.num[,1] <- (ERCC.conc[,4]*(10^(-18)))*(6.0221417*(10^23))
ERCC.num.count <- ERCC.num/50000
ERCC.num.final <- ERCC.num.count*0.009
rownames(ERCC.num) <- rownames(ERCC.num.final) <- ERCC.conc[,2]

SpikeInput <- ERCC.num.final[rownames(CD4.cells)[grepl("ERCC", rownames(CD4.cells))],1]
SpikeInput.1 <- data.frame("Name" = names(SpikeInput),
                           "Molecules" = SpikeInput,
                           stringsAsFactors = FALSE)

#### Create Data objects for each condition

# Young active B6
input.active <- CD4.cells[,grepl("Active", colnames(CD4.cells))]
chips <- as.numeric(sapply(colnames(input.active), 
                           function(n){unlist(strsplit(n, "\\."))[3]}))

Data.YoungActiveB6 <- newBASiCS_Data(Counts = input.active,
                                     Tech = grepl("ERCC", rownames(input.active)), 
                                     SpikeInfo = SpikeInput.1, BatchInfo=chips)
saveRDS(Data.YoungActiveB6, file.path(mypath, "/Data/CD4_active.rds"))

# Young naive B6
input.naive <- CD4.cells[,grepl("Unstimulated", colnames(CD4.cells))]
chips <- as.numeric(sapply(colnames(input.naive), 
                           function(n){unlist(strsplit(n, "\\."))[3]}))

Data.YoungNaiveB6 <- newBASiCS_Data(Counts = input.naive,
                                     Tech = grepl("ERCC", rownames(input.naive)), 
                                     SpikeInfo = SpikeInput.1, BatchInfo=chips)
saveRDS(Data.YoungNaiveB6, file.path(mypath, "/Data/CD4_naive.rds"))


#### Run MCMC on these conditions

chains.path <- file.path(mypath, "/MCMCs/")

# Young active B6
#MCMC.YoungActiveB6 <- BASiCS_MCMC(Data.YoungActiveB6, 160000, 80, 80000, 
#                                  Regression = TRUE, WithSpikes = TRUE,
#                                  StoreChains = TRUE, RunName = "B6active_cata_long",
#                                  StoreDir = chains.path)

#MCMC.YoungActiveB6.noReg <- BASiCS_MCMC(Data.YoungActiveB6, 40000, 20, 20000, 
#                                       Regression = FALSE, WithSpikes = TRUE, 
#                                       StoreChains = TRUE, RunName = "B6active_noReg_cata", 
#                                       StoreDir = chains.path)

# Young naive B6
input.naive <- CD4.cells[,grepl("Unstimulated", colnames(CD4.cells))]
chips <- sapply(colnames(input.naive), function(n){unlist(strsplit(n, "\\."))[3]})

Data.YoungNaiveB6 <- newBASiCS_Data(Counts = input.naive,
                                     Tech = grepl("ERCC", rownames(input.naive)), 
                                     SpikeInfo = SpikeInput.1, BatchInfo=chips)
saveRDS(Data.YoungNaiveB6, file.path(mypath, "Data/CD4_naive.rds"))

#### Run MCMC on these conditions

# Young naive B6
MCMC.YoungNaiveB6 <- BASiCS_MCMC(Data.YoungNaiveB6, 160000, 80, 80000,
                                 Regression = TRUE, WithSpikes = TRUE,
                                 StoreChains = TRUE, RunName = "B6naive_cata_long",
                                 StoreDir = chains.path)

#MCMC.YoungNaiveB6.noReg <- BASiCS_MCMC(Data.YoungNaiveB6, 40000, 20, 20000, 
#                                 Regression = FALSE, WithSpikes = TRUE,
#                                 StoreChains = TRUE, RunName = "B6naive_noReg_cata",
#                                 StoreDir = chains.path)


#### Mixture experiment
set.seed(123456)
sampled.active <- input.active[,sample(1:ncol(input.active), 5)]
input.mixed <- cbind(input.naive, sampled.active)

Data.mixture <- newBASiCS_Data(Counts = input.mixed,
                                    Tech = grepl("ERCC", rownames(input.mixed)), 
                                    SpikeInfo = SpikeInput.1)
saveRDS(Data.mixture, file.path(mypath, "/Data/CD4_mixed.rds"))

MCMC.mixture <- BASiCS_MCMC(Data.mixture, 160000, 80, 80000,
                            Regression = TRUE, WithSpikes = TRUE,
                            StoreChains = TRUE, RunName = "B6mixed_cata_long",
                            StoreDir = chains.path)

