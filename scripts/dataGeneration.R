### run the script
# module load R/4.4.0-gfbf-2023a
# Rscript scripts/dataGeneration.R 
setwd('/home/s4162315/cjica/ClusterwiseJICA_v2/')
library("mclust")
library("multichull")
source("./functions/data_simulation_functions.R")
source("./functions/ClusterwiseJICA_varyQ.R")
source("./functions/ClusterwiseJICA.R")
source("./functions/CJICA_asist_functions.R")
source("./functions/ILS_CJICA.R")

#### Generate simulation data
M <- 2
Vm <- 500
Nk <- 100
Qvect_list <- list(c(2, 2, 2), c(2,2,6), c(2,5,8), c(2,2,2,2,2), c(2,4,4,6,8), c(2,3,6,7,8))
E <- c(0.2, 0.6)
VAF <- c(0.7, 1)
rep <- 1:20

# Create the simulation grid
base_grid <- expand.grid(Qid=seq_along(Qvect_list), E = E, VAF=VAF, rep=rep)
base_grid$Qvect <- Qvect_list[base_grid$Qid]
# Remove Qid if not needed
grid <- base_grid[, !(names(base_grid) %in% "Qid")]

# Create directory to store simulation results if it doesn't exist
if (!dir.exists("simulation_data")) {
    dir.create("simulation_data")
}

print(paste("Starting", nrow(grid), "simulations..."))

for(sim in 1:nrow(grid)){
    cat("Running simulation", sim, "of", nrow(grid), "\n")
    
    # Extract current simulation parameters
    current_K <- length(grid$Qvect[[sim]])
    current_Qvect <- grid$Qvect[[sim]]
    current_E <- grid[sim, ]$E
    current_VAF <- grid[sim, ]$VAF
    current_rep <- grid[sim, ]$rep
    
    # Generate simulation data
    simdata <- Simulate_CJICA(
        Nk = Nk, 
        Vm = Vm,
        K = current_K,
        Qvect = current_Qvect,
        E = current_E,
        M = M,
        type = 1,
        VAF = current_VAF
    )
    
    # Create filename for individual simulation
    Qvect_str <- paste(current_Qvect, collapse = "_")
    filename <- paste0("simulation_data/sim_", sim, 
                      "_K", current_K,
                      "_Q", Qvect_str,
                      "_E", current_E, 
                      "_VAF", current_VAF, 
                      "_rep", current_rep, ".RData")
    
    # Save only the simulation data (not wrapped in additional structure)
    save(simdata, file = filename)
    
}

# Save the parameter grid for reference
save(grid, file = "simulation_data/simulation_parameters.RData")

print("All simulations completed!")
print(paste("Total simulations:", nrow(grid)))
print("Each simulation saved as separate .RData file in 'simulation_data' directory")
print("Parameter grid saved as 'simulation_parameters.RData'")