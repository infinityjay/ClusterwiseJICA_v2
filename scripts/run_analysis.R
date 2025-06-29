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

# Helper function to add timestamp to log messages
log_with_time <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat("[", timestamp, "] ", message, "\n", sep = "")
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript process_single_simulation.R <filename>")
}

filename <- args[1]
log_with_time(paste("Processing file:", filename))

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
log_with_time("Parsed parameters:")
log_with_time(paste("  sim_id:", params$sim_id))
log_with_time(paste("  K:", params$K))
log_with_time(paste("  Qvect:", paste(params$Qvect, collapse = ", ")))
log_with_time(paste("  E:", params$E))
log_with_time(paste("  VAF:", params$VAF))
log_with_time(paste("  rep:", params$rep))

# set seed for reproduce
set.seed(params$sim_id)

# Set simulation parameters
Vm <- 500
Nk <- 100
M <- 2

# Load simulation data from file
log_with_time("Loading simulation data from file...")
data_filepath <- file.path("simulation_data", filename)

if (!file.exists(data_filepath)) {
    stop(paste("Data file not found:", data_filepath))
}

# create results folder
results_dir <- "results"
if (!dir.exists(results_dir)) {
    dir.create(results_dir)
}

load(data_filepath)
log_with_time(paste("Data loaded successfully from:", data_filepath))

# Run ClusterwiseJICA_varyQ analysis
log_with_time("Running ClusterwiseJICA_varyQ analysis...")

complex <- c(0,1,3,4,6,7)
for(comp in complex) {
    log_with_time(paste("Starting analysis for complexity:", comp))
    
    ### use Chull, not use input Q ###
    log_with_time("Running analysis with Chull method...")
    start_time <- Sys.time()
    cjica_chull <- ClusterwiseJICA_varyQ(
        X = simdata$Xe,
        k = params$K, 
        useInputQ = FALSE,                                                   
        starts = 100, 
        scale = FALSE, 
        VAF = params$VAF, 
        complex = comp, 
        useChull = TRUE, 
        Vm = Vm
    )
    end_time <- Sys.time()
    runtime_chull <- difftime(end_time, start_time, units = "mins")
    log_with_time(paste("Analysis (use Chull) completed in", round(runtime_chull, 2), "minutes"))
    
    # calculate ARI
    loss100 <- sapply(seq_along(cjica_chull), function(anom) tail(cjica_chull[[anom]]$aiciter, n = 1))
    optimal <- cjica_chull[[which.min(loss100)]]
    cjica_chull_ARI <- adjustedRandIndex(simdata$P, optimal$p)
    log_with_time(paste("Chull ARI calculated:", round(cjica_chull_ARI, 4)))
    
    # calculate tucker
    tucker_cor_lap <- unlist(TuckCheck(simdata$S))
    clusper <- FindOptimalClusPermut(optimal$p, simdata$P)
    log_with_time(paste("FindOptimalClusPermut calculated: ", clusper$BestPerm))
    
    if(params$K == 3){
        log_with_time(" calculate k = 3 tucker_S")
        tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov))
        
    }else{
        log_with_time(" calculate k = 5 tucker_S")
        tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[4]]], Strue = simdata$S[[4]])$BestRecov,
        FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[5]]], Strue = simdata$S[[5]])$BestRecov))
    }
    log_with_time(paste("Chull Tucker S calculated:", round(tucker_S, 4)))
    
    # save result
    # Create result filename
    result_filename <- paste0(results_dir, "/result_chull_", comp, "_", gsub("\\.RData$", "", filename), ".RData")
    # Save analysis results along with parameters
    analysis_result <- list(
        params = params,
        seed = params$sim_id,
        complex = comp,
        ari = cjica_chull_ARI,
        tucker_S = tucker_S,
        Araw = optimal$ica$Mr,
        runtime = runtime_chull,
        optimal = optimal,
        aicloss = loss100,
        tuckercheck = tucker_cor_lap
    )
    save(analysis_result, file = result_filename)
    log_with_time(paste("Chull results saved to:", result_filename))

    ### use VAF, not use chull ###
    log_with_time("Running analysis with VAF method...")
    start_time <- Sys.time()
    cjica_vaf <- ClusterwiseJICA_varyQ(
        X = simdata$Xe,
        k = params$K, 
        useInputQ = FALSE,                                                   
        starts = 100, 
        scale = FALSE, 
        VAF = params$VAF, 
        complex = comp, 
        useChull = FALSE, 
        Vm = Vm
    )
    end_time <- Sys.time()
    runtime_vaf <- difftime(end_time, start_time, units = "mins")
    log_with_time(paste("Analysis (use VAF) completed in", round(runtime_vaf, 2), "minutes"))
    
    # calculate ARI
    loss100 <- sapply(seq_along(cjica_vaf), function(anom) tail(cjica_vaf[[anom]]$aiciter, n = 1))
    optimal <- cjica_vaf[[which.min(loss100)]] 
    cjica_vaf_ARI <- adjustedRandIndex(simdata$P, optimal$p)
    log_with_time(paste("VAF ARI calculated:", round(cjica_vaf_ARI, 4)))
    
    # calculate tucker
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
    log_with_time(paste("VAF Tucker S calculated:", round(tucker_S, 4)))
    
    # save result
    # Create result filename
    result_filename <- paste0(results_dir, "/result_vaf_", comp, "_", gsub("\\.RData$", "", filename), ".RData")

    # Save analysis results along with parameters
    analysis_result <- list(
        params = params,
        seed = params$sim_id,
        complex = comp,
        ari = cjica_vaf_ARI,
        tucker_S = tucker_S,
        Araw = optimal$ica$Mr,
        runtime = runtime_vaf,
        optimal = optimal,
        aicloss = loss100,
        tuckercheck = tucker_cor_lap
    )
    save(analysis_result, file = result_filename)
    log_with_time(paste("VAF results saved to:", result_filename))
    
    log_with_time(paste("Completed analysis for complexity:", comp))
}

log_with_time("Starting original algorithm analysis...")

### original algorithm, use input Q vec ###
nc1 <- c(2,2,2) 
nc2 <- c(8,8,8)
if(params$K == 5) {
    nc1 <- c(2,2,2,2,2)
    nc2 <- c(8,8,8,8,8)
}

### min, all nc = 2 ###
log_with_time("Running original algorithm with minimum nc (all = 2)...")
start_time <- Sys.time()
cjica_origin_min <- ClusterwiseJICA_varyQ(
    X = simdata$Xe,
    k = params$K,
    nc = nc1,
    useInputQ = TRUE,                                                   
    starts = 100, 
    scale = FALSE, 
    VAF = params$VAF, 
    complex = 0, 
    useChull = FALSE, 
    Vm = Vm
)
end_time <- Sys.time()
runtime_min <- difftime(end_time, start_time, units = "mins")
log_with_time(paste("Analysis (original min nc) completed in", round(runtime_min, 2), "minutes"))

# calculate ARI
loss100 <- sapply(seq_along(cjica_origin_min), function(anom) tail(cjica_origin_min[[anom]]$lossiter, n = 1))
optimal <- cjica_origin_min[[which.min(loss100)]]
cjica_origin_min_ARI <- adjustedRandIndex(simdata$P, optimal$p)
log_with_time(paste("Original min ARI calculated:", round(cjica_origin_min_ARI, 4)))

# calculate tucker
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
log_with_time(paste("Original min Tucker S calculated:", round(tucker_S, 4)))

# save result
# Create result filename
result_filename <- paste0(results_dir, "/result_origin_min_", gsub("\\.RData$", "", filename), ".RData")  # Fixed: unique filename
# Save analysis results along with parameters
analysis_result <- list(
    params = params,
    seed = params$sim_id,
    method = "original_min",
    ari = cjica_origin_min_ARI,
    tucker_S = tucker_S,
    Araw = optimal$ica$Mr,
    runtime = runtime_min,
    optimal = optimal,
    ssloss = loss100,
    tuckercheck = tucker_cor_lap
)
save(analysis_result, file = result_filename)
log_with_time(paste("Original min results saved to:", result_filename))

### max, all nc = 8 ###
log_with_time("Running original algorithm with maximum nc (all = 8)...")
start_time <- Sys.time()
cjica_origin_max <- ClusterwiseJICA_varyQ(
    X = simdata$Xe,
    k = params$K,
    nc = nc2,
    useInputQ = TRUE,                                                   
    starts = 100, 
    scale = FALSE, 
    VAF = params$VAF, 
    complex = 0, 
    useChull = FALSE, 
    Vm = Vm
)
end_time <- Sys.time()
runtime_max <- difftime(end_time, start_time, units = "mins")
log_with_time(paste("Analysis (original max nc) completed in", round(runtime_max, 2), "minutes"))

# calculate ARI
loss100 <- sapply(seq_along(cjica_origin_max), function(anom) tail(cjica_origin_max[[anom]]$lossiter, n = 1))
optimal <- cjica_origin_max[[which.min(loss100)]]
cjica_origin_max_ARI <- adjustedRandIndex(simdata$P, optimal$p)
log_with_time(paste("Original max ARI calculated:", round(cjica_origin_max_ARI, 4)))

# calculate tucker
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
log_with_time(paste("Original max Tucker S calculated:", round(tucker_S, 4)))

# save result
# Create result filename
result_filename <- paste0(results_dir, "/result_origin_max_", gsub("\\.RData$", "", filename), ".RData")  # Fixed: unique filename
# Save analysis results along with parameters
analysis_result <- list(
    params = params,
    seed = params$sim_id,
    method = "original_max",
    ari = cjica_origin_max_ARI,
    tucker_S = tucker_S,
    Araw = optimal$ica$Mr,
    runtime = runtime_max,
    optimal = optimal,
    ssloss = loss100,
    tuckercheck = tucker_cor_lap
)
save(analysis_result, file = result_filename)
log_with_time(paste("Original max results saved to:", result_filename))

log_with_time("Job completed successfully!")