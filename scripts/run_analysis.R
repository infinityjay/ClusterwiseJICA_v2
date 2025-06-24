# run_analysis.R
# Script to process a single simulation file

# Load required libraries
library("mclust")
library("multichull")
source("./functions/data_simulation_functions.R")
source("./functions/ClusterwiseJICA_varyQ.R")
source("./functions/ClusterwiseJICA.R")
source("./functions/CJICA_asist_functions.R")
source("./functions/ILS_CJICA.R")

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript process_single_simulation.R <filename>")
}

filename <- args[1]
cat("Processing file:", filename, "\n")

# Parse parameters from filename
# Expected format: sim_X_KY_QZ_EW_VAFV_repR.RData
parse_filename <- function(fname) {
    # Remove .RData extension and split by underscores
    parts <- strsplit(gsub("\\.RData$", "", fname), "_")[[1]]
    
    # Extract parameters
    sim_id <- as.numeric(gsub("sim", "", parts[2]))
    K <- as.numeric(gsub("K", "", parts[3]))
    
    # Extract Q vector (everything between Q and E)
    q_start <- which(grepl("^Q", parts))
    e_pos <- which(grepl("^E", parts))
    q_parts <- parts[(q_start+1):(e_pos-1)]
    Qvect <- as.numeric(q_parts)
    
    E <- as.numeric(gsub("E", "", parts[e_pos]))
    VAF <- as.numeric(gsub("VAF", "", parts[e_pos + 1]))
    rep <- as.numeric(gsub("rep", "", parts[e_pos + 2]))
    
    return(list(
        sim_id = sim_id,
        K = K,
        Qvect = Qvect,
        E = E,
        VAF = VAF,
        rep = rep
    ))
}

# Parse parameters from filename
params <- parse_filename(filename)
cat("Parsed parameters:\n")
cat("  sim_id:", params$sim_id, "\n")
cat("  K:", params$K, "\n")
cat("  Qvect:", paste(params$Qvect, collapse = ", "), "\n")
cat("  E:", params$E, "\n")
cat("  VAF:", params$VAF, "\n")
cat("  rep:", params$rep, "\n")

# Set simulation parameters
Vm <- 500
Nk <- 100
M <- 2

# Load simulation data from file
cat("Loading simulation data from file...\n")
data_filepath <- file.path("simulation_data", filename)

if (!file.exists(data_filepath)) {
    stop(paste("Data file not found:", data_filepath))
}

# create results folder
results_dir <- "results"
if (!dir.exists(results_dir)) {
    dir.create(results_dir)
}

simdata <- load(data_filepath)
cat("Data loaded successfully from:", data_filepath, "\n")

# Run ClusterwiseJICA_varyQ analysis
cat("Running ClusterwiseJICA_varyQ analysis...\n")

complex <- c(0,1,3,4,6,7)

for(comp in complex) {
    ## use Chull, not use input Q
    start_time <- Sys.time()
    cjica_chull <- ClusterwiseJICA_varyQ(
        X = simdata$Xe,
        k = params$K, 
        useInputQ = F,                                                   
        starts = 100, 
        scale = F, 
        VAF = params$VAF, 
        complex = comp, 
        useChull = T, 
        Vm = Vm
    )
    end_time <- Sys.time()
    cat("Analysis (use Chull) completed in", difftime(end_time, start_time, units = "mins"), "minutes\n")
    # calculate ARI
    loss100 <- sapply(seq_along(cjica_chull), function(anom) tail(cjica_chull[[anom]]$aiciter, n = 1))
    optimal <- cjica[[which.min(loss100)]]
    cjica_chull_ARI <- adjustedRandIndex(simdata$P, optimal$p)
    # calculate turker
    tucker_cor_lap <- unlist(TuckCheck(simdata$S))
    clusper <- FindOptimalClusPermut(optimal$p, simdata$P)
    
    if(params$K == 3){
        tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov))
        
    }else{
        tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[4]]], Strue = simdata$S[[4]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[5]]], Strue = simdata$S[[5]])$BestRecov))
    }
    # save result

    ## use VAF=0.7

    

}
## original algorithm, use inpute Q vec
nc1 <- c(2,2,2) 
nc2 <- c(8,8,8)
if(params$K == 5) {
    nc1 <- c(2,2,2,2,2)
    nc2 <- c(8,8,8,8,8)
}
## min, all nc = 2
start_time <- Sys.time()
cjica_origin_min <- ClusterwiseJICA_varyQ(
    X = simdata$Xe,
    k = params$K,
    nc = nc1
    useInputQ = T,                                                   
    starts = 100, 
    scale = F, 
    VAF = params$VAF, 
    complex = 0, 
    useChull = F, 
    Vm = Vm
)
end_time <- Sys.time()
cat("Analysis (use Chull) completed in", difftime(end_time, start_time, units = "mins"), "minutes\n")
# calculate ARI
loss100 <- sapply(seq_along(cjica_origin_min), function(anom) tail(cjica_origin_min[[anom]]$lossiter, n = 1))
optimal <- cjica_origin_min[[which.min(loss100)]]
cjica_origin_min_ARI <- adjustedRandIndex(simdata$P, optimal$p)
# calculate turker
tucker_cor_lap <- unlist(TuckCheck(simdata$S))
clusper <- FindOptimalClusPermut(optimal$p, simdata$P)

if(params$K == 3){
    tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
    FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
    FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov))
    
}else{
    tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
    FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
    FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov,
    FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[4]]], Strue = simdata$S[[4]])$BestRecov,
    FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[5]]], Strue = simdata$S[[5]])$BestRecov))
}

# Create result filename
result_filename <- paste0(results_dir, "/result_", gsub("\\.RData$", "", filename), ".RData")

# Save analysis results along with parameters (simdata already loaded from file)
analysis_result <- list(
    result = result,
    parameters = params,
    filename = filename,
    runtime = difftime(end_time, start_time, units = "mins")
)

save(analysis_result, file = result_filename)
cat("Results saved to:", result_filename, "\n")

cat("Job completed successfully!\n")