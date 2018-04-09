library(BASiCS)

files.A <- list.files("/nfs/research1/marioni/Nils/BASiCS/Revisions/Simulations/MCMCs/FPR/", full.names = TRUE)[1:25]
files.B <- list.files("/nfs/research1/marioni/Nils/BASiCS/Revisions/Simulations/MCMCs/FPR/", full.names = TRUE)[26:50]

file.names <- list.files("/nfs/research1/marioni/Nils/BASiCS/Revisions/Simulations/MCMCs/FPR/", full.names = FALSE)[1:25]

FPR.mu <- vector(length = 25)
names(FPR.mu) <- file.names
FPR.delta <- vector(length = 25)
names(FPR.delta) <- file.names
FPR.epsilon <- vector(length = 25)
names(FPR.epsilon) <- file.names

EFDR_mu_NULL <- vector(length = 25)
names(EFDR_mu_NULL) <- file.names
EFDR_delta_NULL <- vector(length = 25)
names(EFDR_delta_NULL) <- file.names
EFDR_eps_NULL <- vector(length = 25)
names(EFDR_eps_NULL) <- file.names

for(i in 1:length(files.A)){
	cur_A <- readRDS(files.A[i])
	cur_B <- readRDS(files.B[i])

	cur_Test <- BASiCS_TestDE(Chain1 = cur_A, Chain2 = cur_B,
			         Plot = FALSE, PlotOffset = FALSE,
				 GroupLabel1 = "A", GroupLabel2 = "B")

	FPR.mu[i] <- sum(cur_Test$TableMean$ResultDiffMean != "NoDiff")/
			nrow(cur_Test$TableMean)
	FPR.delta[i] <- sum(cur_Test$TableDisp$ResultDiffDisp != "NoDiff" & cur_Test$TableDisp$ResultDiffDisp != "ExcludedFromTesting")/
                        sum(cur_Test$TableDisp$ResultDiffDisp != "ExcludedFromTesting")
	FPR.epsilon[i] <- sum(cur_Test$TableResDisp$ResultDiffResDisp != "NoDiff" & cur_Test$TableResDisp$ResultDiffResDisp != "ExcludedFromTesting")/
                        sum(cur_Test$TableResDisp$ResultDiffResDisp != "ExcludedFromTesting")

	EFDR_mu_NULL[i] <- sum(1 - cur_Test$TableMean$ProbDiffMean[cur_Test$TableMean$ResultDiffMean != "NoDiff"])/
						sum(cur_Test$TableMean$ResultDiffMean != "NoDiff")
	EFDR_delta_NULL[i] <- sum(1 - cur_Test$TableDisp$ProbDiffDisp[cur_Test$TableDisp$ResultDiffDisp != "ExcludedFromTesting" &
                                          cur_Test$TableDisp$ResultDiffDisp != "NoDiff"])/
  			   sum(cur_Test$TableDisp$ResultDiffDisp != "ExcludedFromTesting" &
                                          cur_Test$TableDisp$ResultDiffDisp != "NoDiff")
	EFDR_eps_NULL[i] <- sum(1 - cur_Test$TableResDisp$ProbDiffResDisp[cur_Test$TableResDisp$ResultDiffResDisp != "ExcludedFromTesting" &
                                          cur_Test$TableResDisp$ResultDiffResDisp != "NoDiff"])/
  		    sum(cur_Test$TableResDisp$ResultDiffResDisp != "ExcludedFromTesting" &
                                          cur_Test$TableResDisp$ResultDiffResDisp != "NoDiff")
}

saveRDS(FPR.mu, "FPR_mu.rds")
saveRDS(FPR.delta, "FPR_delta.rds")
saveRDS(FPR.epsilon, "FPR_eps.rds")

saveRDS(EFDR_mu_NULL, "EFDR_mu_NULL.rds")
saveRDS(EFDR_delta_NULL, "EFDR_delta_NULL.rds")
saveRDS(EFDR_eps_NULL, "EFDR_eps_NULL.rds")

