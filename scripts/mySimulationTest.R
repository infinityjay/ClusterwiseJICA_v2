## # load R module
## module load R/4.4.0-gfbf-2023a
## # run script
## Rscript my_script.R 1

setwd('/home/s4162315/cjica/ClusterwiseJICA_v2/')

library("mclust")
library("multichull")
source("./functions/data_simulation_functions.R")
source("./functions/ClusterwiseJICA_varyQ.R")
source("./functions/ClusterwiseJICA.R")
source("./functions/CJICA_asist_functions.R")
source("./functions/ILS_CJICA.R")

#### generate simulation data

Nk <- c(50, 100)
# Qvect <- list(c(3,3), c(8,8), c(2,8), c(2,2,2,2), c(2,8,3,8))
Qvect_list <- list(c(3,3), c(8,8))
E <- c(0.2, 0.6)
cor <- c(0, .60, .94)
VAF <- c(0.8, 1)
rep <- 1:5

base_grid <- expand.grid(Nk=Nk, Qid=seq_along(Qvect_list), E=E, cor=cor, VAF=VAF, rep=rep)
base_grid$Qvect <- Qvect_list[base_grid$Qid]

# Remove Qid if not needed
grid <- base_grid[, !(names(base_grid) %in% "Qid")]

# Print total grid information
cat(sprintf("Total grid size: %d rows\n", nrow(grid)))

# Get command line argument
args <- commandArgs(TRUE)
if (length(args) == 0) {
  stop("Please provide a split number as argument")
}
sp <- as.numeric(args[1])

# Split grid into chunks of 10 rows each
splits <- split(1:nrow(grid), ceiling(seq_along(1:nrow(grid))/10))
cat(sprintf("Total number of splits: %d\n", length(splits)))
cat(sprintf("Processing split: %d\n", sp))

# Check if split number is valid
if (sp > length(splits)) {
  cat(sprintf("Split %d does not exist. Maximum split number is %d\n", sp, length(splits)))
  quit(status = 0)  # Exit gracefully, not an error
}

rows <- splits[[sp]]
cat(sprintf("Processing rows: %s\n", paste(rows, collapse = ", ")))

#### run simulation and c-jica function

Vm <- 2500
complex <- c(1,3,4,6,7)
grid <- grid[rows,]

# Create output directory
output_dir <- './result/sim1/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

for(sim in 1:nrow(grid)){
  cat(sprintf("\n--- Running simulation %d of %d (row index %d) at %s ---\n", sim, nrow(grid), rows[sim]), Sys.time())
  cat(sprintf("Params: Nk=%d, Qvect=%s, E=%.2f, cor=%.2f, VAF=%.2f\n", 
              grid[sim, ]$Nk, 
              paste(grid$Qvect[[sim]], collapse = ","), 
              grid[sim, ]$E, 
              grid[sim, ]$cor, 
              grid[sim, ]$VAF))
  seed <- 123
  set.seed(seed)
  simdata <- Simulate_CJICA(
    Nk = grid[sim, ]$Nk, 
    Vm = Vm,
    K = length(grid$Qvect[[sim]]),
    Qvect = grid$Qvect[[sim]],
    E = grid[sim, ]$E,
    M = 2,
    type = 1,
    cor = grid[sim, ]$cor,
    VAF = grid[sim, ]$VAF
  )

  for(comp in complex) {
    # run cjica and calculate executing time
    cat(sprintf("  > Running ClusterwiseJICA with complex=%d...\n", comp))
    ptm <- proc.time()
    cjica <- ClusterwiseJICA_varyQ(
      X = simdata$Xe,
      k = length(grid$Qvect[[sim]]), 
      useInputQ = F, 
      starts = 100, 
      scale = F, 
      VAF = grid[sim, ]$VAF, 
      complex = comp, 
      useChull = T, 
      Vm = Vm
    )
    time <- proc.time() - ptm
    cat(sprintf("    Completed in %.2f seconds\n", time[3]))
    output <- list()
    output$time <- time
  }
  
  # Save results
  filename <- paste('CJICA_sim1_', rows[sim], '.Rdata', sep = '')
  filepath <- file.path(output_dir, filename)
  save(output, file = filepath)
  cat(sprintf("    Result saved to %s\n", filepath))
}
cat("=== CJICA Simulation Script Completed ===\n")