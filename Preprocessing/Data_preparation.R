##########################################################
#### Script to preprocess datasets used in this study ####
##########################################################

setwd("/Users/nils/Google Drive File Stream/My Drive/BASiCS_add-on/")
library(data.table)

#### CD4 T cell dataset ####

# Downloaded from:
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4888/

# Read in the data from Martinez et al. - E-MTAB-4888.processed.2.zip
input <- read.table("Data/Raw_data/CD4/raw_data.txt", header = TRUE, sep = "\t")

# Metadata file - E-MTAB-4888.additional.1.zip
meta <- read.table("Data/Raw_data/CD4/metadata_file.txt", header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)

# Select young B6 cells
input.CD4 <- input[,meta$X[meta$Individuals == "B6 young 1" |
                             meta$Individuals == "B6 young 2"]]
meta <- meta[meta$Individuals == "B6 young 1" |
               meta$Individuals == "B6 young 2",] 

# Rename columns
colnames(input.CD4) <- paste(meta$Individuals, meta$Stimulus, rownames(meta), sep = " ")

# Filtering of bological and technical genes
g.bio <- input.CD4[grepl("ENS", rownames(input.CD4)),]
ERCC <- input.CD4[grepl("ERCC", rownames(input.CD4)),]
ERCC <- ERCC[rowSums(ERCC) > 0,]
g.bio <- g.bio[rowMeans(g.bio) > 1,]

write.table(rbind(g.bio, ERCC), "Data/Test_Data/CD4_NaiveActiveYoungB6.txt", sep = "\t")

#### Mouse brain cells data ####
# Data from Zeisel et al.
# Downloaded from: http://linnarssonlab.org/cortex/

# Raw counts
Zeisel <- read.table("Data/Raw_data/Zeisel/expression_mRNA_17-Aug-2014.txt", 
                    header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                    fill = TRUE)

# Meta information
Zeisel.meta <- Zeisel[1:11,]
Zeisel.meta <- Zeisel.meta[,-c(1,2)]

# Collect biological counts
Zeisel.bio <- Zeisel[12:nrow(Zeisel),]
rownames(Zeisel.bio) <- as.character(Zeisel.bio[,1])
Zeisel.bio <- Zeisel.bio[,-c(1,2)]
colnames(Zeisel.bio) <- as.character(Zeisel.meta[8,])
Zeisel.bio <- data.matrix(Zeisel.bio[rowMeans(data.matrix(Zeisel.bio)) > 0.1,])

# Spike-in counts
Zeisel.ERCC <- read.table("Data/Raw_data/Zeisel/expression_spikes_17-Aug-2014.txt", 
                          header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                          fill = TRUE)
Zeisel.ERCC <- Zeisel.ERCC[12:nrow(Zeisel.ERCC),]
rownames(Zeisel.ERCC) <- as.character(Zeisel.ERCC[,1])
Zeisel.ERCC <- Zeisel.ERCC[,-c(1,2)]
colnames(Zeisel.ERCC) <- as.character(Zeisel.meta[8,])

# Look at sizes of cell groups 
table(as.character(Zeisel.meta[9,]))

# Filter out one cell type for model testing - microglia
microglia <- Zeisel.bio[,Zeisel.meta[2,] == 5]
ERCC <- Zeisel.ERCC[,colnames(microglia)]
ERCC <- ERCC[rowSums(ERCC) > 0,]

write.table(rbind(microglia, ERCC), "Data/Test_Data/microglia_Zeisel.txt", sep = "\t")

# Another cell type for downsampling - CA1
CA1 <- Zeisel.bio[,Zeisel.meta[2,] == 3]
ERCC <- Zeisel.ERCC[,colnames(CA1)]

write.table(rbind(CA1, ERCC), "Data/Test_Data/CA1_Zeisel.txt", sep = "\t")

#### Pool and split data ####

# Downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54695

# Data from Gruen et al.
input <- read.table("Data/Raw_data/PoolSplit/GSE54695_data_transcript_counts.txt", sep = "\t", header = TRUE)
rownames(input) <- input$GENENAME
input <- input[,-1]

# Batch info
Batch <- c(rep(1, 40), rep(2, 40), rep(1, 40), rep(2, 40), rep(1, 40), rep(2, 40), rep(1, 40), rep(2, 40))

# Generate raw UMI counts
UMICount <- function(MoleculeCount, UMILength)
{
  # MoleculeCount is the normalized count
  M = 4^UMILength
  UMICount = M*(1-exp(-MoleculeCount/M))
  return(UMICount)
}
CountsUMI = round(UMICount(input, 4))

# Quality control - remove low Nanog cells
Pou5f1.per.cell <- as.numeric(CountsUMI["Pou5f1",])

CountsUMI_1 <- CountsUMI[,Pou5f1.per.cell >= 10]
CountsUMI.bio <- CountsUMI_1[!grepl("ERCC", rownames(CountsUMI_1)), ]
ERCC <- CountsUMI_1[grepl("ERCC", rownames(CountsUMI_1)), ]

CountsUMI.bio <- CountsUMI.bio[rowMeans(CountsUMI.bio) > 0.5,]
ERCC <- ERCC[rowMeans(ERCC) > 0,]
Batch <- Batch[Pou5f1.per.cell >= 10]

Data <- rbind(CountsUMI.bio, ERCC)
colnames(Data) <- paste(sapply(colnames(Data), function(n){unlist(strsplit(n, "_"))[1]}),
                        sapply(colnames(Data), function(n){unlist(strsplit(n, "_"))[2]}),
                        Batch, "_", 
                        sapply(colnames(Data), function(n){unlist(strsplit(n, "_"))[3]}), sep = "")

write.table(Data, "Data/Test_Data/PoolSplit.txt", sep = '\t')

#### Dictyostelium data ####

# Downloaded from: http://www.cell.com/action/showMethods?pii=S0960-9822%2817%2930564-X

# Data from Antolović et al.
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

# Batch 3 appears to cluster separatlety 
# BASiCS will capture this batch effect with the parameter theta

# Filtering on biological and technical genes
input.Dict.bio <- input.Dict.bio[rowMeans(input.Dict.bio) > 1,]
ERCC <- ERCC[rowMeans(ERCC) > 0, colnames(input.Dict.bio)]

write.table(rbind(input.Dict.bio, ERCC), "Data/Test_Data/Dictyostelium.txt", sep = "\t")

#### CD4 T cell differentiation ####

# Data downlodaed from: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4388/
# And mapped against mouse genome using gsnap

# Data from Loennberg et al.
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

# Renames colnames
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

write.table(input, "Data/Test_Data/CD4_diff.txt", sep = "\t")
