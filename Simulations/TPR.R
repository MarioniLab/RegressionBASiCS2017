library(BASiCS)

files.A <- list.files("/nfs/research1/marioni/Nils/BASiCS/Revisions/Simulations/MCMCs/TPR/", full.names = TRUE)[1:25]
files.B <- list.files("/nfs/research1/marioni/Nils/BASiCS/Revisions/Simulations/MCMCs/TPR/", full.names = TRUE)[26:50]
deltas <- list.files("/nfs/research1/marioni/Nils/BASiCS/Revisions/Simulations/parameters", full.names = TRUE)

file.names <- list.files("/nfs/research1/marioni/Nils/BASiCS/Revisions/Simulations/MCMCs/TPR/", full.names = FALSE)[1:25]

TPR.delta <- vector(length = 25)
names(TPR.delta) <- file.names

EFDR.delta <- vector(length = 25)
names(EFDR.delta) <- file.names


TPR.eps <- vector(length = 25)
names(TPR.eps) <- file.names

EFDR.eps <- vector(length = 25)
names(EFDR.eps) <- file.names

load("/nfs/research1/marioni/Nils/BASiCS/Data/SimData/DataForSimulations.RData")
delta.ref <- parameters$delta

for(i in 1:length(files.A)){
	cur_A <- readRDS(files.A[i])
	cur_B <- readRDS(files.B[i])
	
	cur_delta <- readRDS(deltas[i])

	cur_test <- BASiCS_TestDE(Chain1 = cur_A, Chain2 = cur_B,
			         Plot = FALSE, PlotOffset = FALSE,
				 GroupLabel1 = "A", GroupLabel2 = "B")
	
	changes <- colnames(cur_A@parameters$delta)[cur_delta/delta.ref != 1]

	TPR.delta[i] <- sum(cur_test$TableDisp$ResultDiffDisp != "NoDiff" &
                   		cur_test$TableDisp$ResultDiffDisp != "ExcludedFromTesting" &
                   		cur_test$TableDisp$GeneName %in% changes)/
             		sum(cur_test$TableDisp$GeneName %in% changes &
                   		cur_test$TableDisp$ResultDiffDisp != "ExcludedFromTesting")
	
	EFDR.delta[i] <- sum(1 - cur_test$TableDisp$ProbDiffDisp[cur_test$TableDisp$ResultDiffDisp != "ExcludedFromTesting" &
                                          cur_test$TableDisp$ResultDiffDisp != "NoDiff"])/
  		         sum(cur_test$TableDisp$ResultDiffDisp != "ExcludedFromTesting" &
                                          cur_test$TableDisp$ResultDiffDisp != "NoDiff")

	TPR.eps[i] <- sum(cur_test$TableResDisp$ResultDiffResDisp != "NoDiff" & 
                   	cur_test$TableResDisp$ResultDiffResDisp != "ExcludedFromTesting" & 
                   	cur_test$TableResDisp$GeneName %in% changes)/
             	   sum(cur_test$TableResDisp$GeneName %in% changes &
                   	cur_test$TableResDisp$ResultDiffResDisp != "ExcludedFromTesting")

	EFDR.eps[i] <- sum(1 - cur_test$TableResDisp$ProbDiffResDisp[cur_test$TableResDisp$ResultDiffResDisp != "ExcludedFromTesting" &
                    	cur_test$TableResDisp$ResultDiffResDisp != "NoDiff"])/
  		    sum(cur_test$TableResDisp$ResultDiffResDisp != "ExcludedFromTesting" &
                        cur_test$TableResDisp$ResultDiffResDisp != "NoDiff")	

}


saveRDS(TPR.delta, "TPR_delta_alternative.rds")

saveRDS(TPR.eps, "TPR_epsilon_alternative.rds")

saveRDS(EFDR.delta, "EFDR_delta_alternative.rds")

saveRDS(EFDR.eps, "EFDR_eps_alternative.rds")
