# collect_results.R
# Script to collect all analysis results into a single comprehensive dataset

library(data.table)

# Helper function to add timestamp to log messages
log_with_time <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat("[", timestamp, "] ", message, "\n", sep = "")
}

# Function to parse filename and extract parameters
parse_result_filename <- function(fname) {
    # Remove path if present
    fname <- basename(fname)
    
    # Extract method from filename
    if (grepl("^result_origin_min_", fname)) {
        method <- "origin_min"
        remaining <- gsub("^result_origin_min_", "", fname)
    } else if (grepl("^result_origin_max_", fname)) {
        method <- "origin_max"
        remaining <- gsub("^result_origin_max_", "", fname)
    } else if (grepl("^result_chull_", fname)) {
        method <- "chull"
        # Extract complexity number
        complex_match <- regmatches(fname, regexpr("^result_chull_[0-9]+_", fname))
        complex_num <- as.numeric(gsub("^result_chull_|_$", "", complex_match))
        remaining <- gsub("^result_chull_[0-9]+_", "", fname)
    } else if (grepl("^result_vaf_", fname)) {
        method <- "vaf"
        # Extract complexity number
        complex_match <- regmatches(fname, regexpr("^result_vaf_[0-9]+_", fname))
        complex_num <- as.numeric(gsub("^result_vaf_|_$", "", complex_match))
        remaining <- gsub("^result_vaf_[0-9]+_", "", fname)
    } else {
        stop(paste("Unknown filename pattern:", fname))
    }
    
    # Remove .RData extension
    remaining <- gsub("\\.RData$", "", remaining)
    
    # Parse the simulation parameters from remaining part
    # Expected format: sim_X_KY_QZ_EW_VAFV_repR
    parts <- strsplit(remaining, "_")[[1]]
    
    # Extract parameters
    sim_id <- as.numeric(gsub("sim", "", parts[2]))
    K <- as.numeric(gsub("K", "", parts[3]))
    
    # Extract Q vector (everything between Q and E)
    q_start <- which(grepl("^Q", parts))
    e_pos <- which(grepl("^E", parts))
    
    # Get the first Q value from the Q part itself
    q_first <- as.numeric(gsub("Q", "", parts[q_start]))
    # Get the remaining Q values from parts between Q and E
    if ((e_pos - q_start) > 1) {
        q_remaining <- as.numeric(parts[(q_start+1):(e_pos-1)])
        Qvect <- c(q_first, q_remaining)
    } else {
        Qvect <- q_first
    }
    
    E <- as.numeric(gsub("E", "", parts[e_pos]))
    VAF <- as.numeric(gsub("VAF", "", parts[e_pos + 1]))
    rep <- as.numeric(gsub("rep", "", parts[e_pos + 2]))
    
    result <- list(
        method = method,
        sim_id = sim_id,
        K = K,
        Qvect = Qvect,
        E = E,
        VAF = VAF,
        rep = rep
    )
    
    # Add complexity if it exists
    if (exists("complex_num")) {
        result$complex <- complex_num
    } else {
        result$complex <- NA
    }
    
    return(result)
}

# Function to safely extract nested list elements
safe_extract <- function(obj, path, default = NA) {
    tryCatch({
        # Navigate through nested list structure
        result <- obj
        for (element in path) {
            result <- result[[element]]
        }
        return(result)
    }, error = function(e) {
        return(default)
    })
}

# Function to handle tuckercheck data properly
handle_tuckercheck <- function(tuckercheck_data) {
    if (is.null(tuckercheck_data) || length(tuckercheck_data) == 0) {
        return(list(
            tuckercheck_1 = NA,
            tuckercheck_2 = NA,
            tuckercheck_mean = NA,
            tuckercheck_min = NA,
            tuckercheck_max = NA,
            tuckercheck_length = 0
        ))
    }
    
    # Convert to numeric if it's not already
    if (!is.numeric(tuckercheck_data)) {
        tuckercheck_data <- as.numeric(tuckercheck_data)
    }
    
    # Remove any NA values for calculations
    valid_data <- tuckercheck_data[!is.na(tuckercheck_data)]
    
    result <- list(
        tuckercheck_1 = if(length(tuckercheck_data) >= 1) tuckercheck_data[1] else NA,
        tuckercheck_2 = if(length(tuckercheck_data) >= 2) tuckercheck_data[2] else NA,
        tuckercheck_mean = if(length(valid_data) > 0) mean(valid_data) else NA,
        tuckercheck_min = if(length(valid_data) > 0) min(valid_data) else NA,
        tuckercheck_max = if(length(valid_data) > 0) max(valid_data) else NA,
        tuckercheck_length = length(tuckercheck_data)
    )
    
    return(result)
}

# Function to calculate summary statistics for loss vectors
calc_loss_stats <- function(loss_vector) {
    if (is.null(loss_vector) || length(loss_vector) == 0) {
        return(list(
            loss_min = NA,
            loss_max = NA,
            loss_mean = NA,
            loss_sd = NA,
            loss_length = 0
        ))
    }
    
    return(list(
        loss_min = min(loss_vector, na.rm = TRUE),
        loss_max = max(loss_vector, na.rm = TRUE),
        loss_mean = mean(loss_vector, na.rm = TRUE),
        loss_sd = sd(loss_vector, na.rm = TRUE),
        loss_length = length(loss_vector)
    ))
}

# Main collection function
collect_all_results <- function(results_dir = "results", output_dir = "collect_result") {
    
    log_with_time("Starting result collection process...")
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        log_with_time(paste("Created output directory:", output_dir))
    }
    
    # Find all result files
    result_files <- list.files(results_dir, pattern = "^result_.*\\.RData$", full.names = TRUE)
    
    if (length(result_files) == 0) {
        stop("No result files found in the results directory!")
    }
    
    log_with_time(paste("Found", length(result_files), "result files"))
    
    # Initialize list to store all results
    all_results <- list()
    failed_files <- character()
    
    # Process each result file
    for (i in seq_along(result_files)) {
        file_path <- result_files[i]
        file_name <- basename(file_path)
        
        if (i %% 50 == 0 || i == length(result_files)) {
            log_with_time(paste("Processing file", i, "of", length(result_files), ":", file_name))
        }
        
        tryCatch({
            # Load the result file
            load(file_path)  # This should load 'analysis_result'
            
            # Parse filename to get method and parameters
            parsed_info <- parse_result_filename(file_name)
            
            # Handle tuckercheck data properly
            tuckercheck_raw <- safe_extract(analysis_result, "tuckercheck")
            tuckercheck_processed <- handle_tuckercheck(tuckercheck_raw)
            
            # Create a comprehensive result record
            result_record <- list(
                # File information
                filename = file_name,
                filepath = file_path,
                
                # Method information
                method = parsed_info$method,
                complex = parsed_info$complex,
                
                # Simulation parameters
                sim_id = parsed_info$sim_id,
                K = parsed_info$K,
                E = parsed_info$E,
                VAF = parsed_info$VAF,
                rep = parsed_info$rep,
                
                # Q vector information
                Q1 = if(length(parsed_info$Qvect) >= 1) parsed_info$Qvect[1] else NA,
                Q2 = if(length(parsed_info$Qvect) >= 2) parsed_info$Qvect[2] else NA,
                Q3 = if(length(parsed_info$Qvect) >= 3) parsed_info$Qvect[3] else NA,
                Q4 = if(length(parsed_info$Qvect) >= 4) parsed_info$Qvect[4] else NA,
                Q5 = if(length(parsed_info$Qvect) >= 5) parsed_info$Qvect[5] else NA,
                Qvect_string = paste(parsed_info$Qvect, collapse = "_"),
                
                # Performance metrics
                ari = safe_extract(analysis_result, "ari"),
                tucker_S = safe_extract(analysis_result, "tucker_S"),
                
                # Runtime information
                runtime = safe_extract(analysis_result, "runtime"),
                
                # Seed information
                seed = safe_extract(analysis_result, "seed")
            )
            
            # Add processed tuckercheck data
            result_record <- c(result_record, tuckercheck_processed)
            
            # Handle loss information (could be ssloss or aicloss)
            ssloss <- safe_extract(analysis_result, "ssloss")
            aicloss <- safe_extract(analysis_result, "aicloss")
            
            if (!is.na(ssloss[1]) && !is.null(ssloss)) {
                loss_stats <- calc_loss_stats(ssloss)
                result_record$loss_type <- "ssloss"
            } else if (!is.na(aicloss[1]) && !is.null(aicloss)) {
                loss_stats <- calc_loss_stats(aicloss)
                result_record$loss_type <- "aicloss"
            } else {
                loss_stats <- calc_loss_stats(NULL)
                result_record$loss_type <- NA
            }
            
            # Add loss statistics to result record
            result_record <- c(result_record, loss_stats)
            
            # Add to results list
            all_results[[i]] <- result_record
            
        }, error = function(e) {
            log_with_time(paste("ERROR processing file", file_name, ":", e$message))
            failed_files <- c(failed_files, file_name)
        })
    }
    
    # Remove NULL entries (failed files)
    all_results <- all_results[!sapply(all_results, is.null)]
    
    log_with_time(paste("Successfully processed", length(all_results), "files"))
    if (length(failed_files) > 0) {
        log_with_time(paste("Failed to process", length(failed_files), "files"))
        log_with_time("Failed files:")
        for (failed_file in failed_files) {
            log_with_time(paste("  -", failed_file))
        }
    }
    
    # Convert to data frame
    log_with_time("Converting to data frame...")
    results_df <- rbindlist(all_results, fill = TRUE)
    
    # Convert to regular data frame
    results_df <- as.data.frame(results_df)
    
    # Add collection timestamp
    results_df$collection_timestamp <- Sys.time()
    
    # Sort by sim_id, method, complex for better organization
    results_df <- results_df[order(results_df$sim_id, results_df$method, results_df$complex, na.last = TRUE), ]
    
    # Reset row names
    rownames(results_df) <- NULL
    
    log_with_time(paste("Final dataset dimensions:", nrow(results_df), "rows,", ncol(results_df), "columns"))
    
    # Save the comprehensive results
    output_file <- file.path(output_dir, "comprehensive_results.RData")
    save(results_df, file = output_file)
    log_with_time(paste("Comprehensive results saved to:", output_file))
    
    # Also save as CSV for easier viewing
    csv_file <- file.path(output_dir, "comprehensive_results.csv")
    write.csv(results_df, csv_file, row.names = FALSE)
    log_with_time(paste("Results also saved as CSV:", csv_file))
    
    # Create summary statistics
    log_with_time("Creating summary statistics...")
    
    # Summary by method
    method_summary <- aggregate(cbind(ari, tucker_S, runtime) ~ method, 
                               data = results_df, 
                               FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                  sd = sd(x, na.rm = TRUE), 
                                                  count = sum(!is.na(x))))
    
    # Summary by simulation parameters
    param_summary <- aggregate(cbind(ari, tucker_S, runtime) ~ K + E + VAF, 
                              data = results_df, 
                              FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                 sd = sd(x, na.rm = TRUE), 
                                                 count = sum(!is.na(x))))
    
    # Save summaries
    summary_file <- file.path(output_dir, "result_summaries.RData")
    save(method_summary, param_summary, file = summary_file)
    log_with_time(paste("Summary statistics saved to:", summary_file))
    
    # Print some basic statistics
    log_with_time("=== COLLECTION SUMMARY ===")
    log_with_time(paste("Total results collected:", nrow(results_df)))
    log_with_time(paste("Unique simulation IDs:", length(unique(results_df$sim_id))))
    log_with_time(paste("Methods found:", paste(unique(results_df$method), collapse = ", ")))
    log_with_time(paste("K values:", paste(sort(unique(results_df$K)), collapse = ", ")))
    log_with_time(paste("E values:", paste(sort(unique(results_df$E)), collapse = ", ")))
    log_with_time(paste("VAF values:", paste(sort(unique(results_df$VAF)), collapse = ", ")))
    
    # Check for missing combinations
    expected_methods <- c("origin_min", "origin_max", "chull", "vaf")
    expected_complex <- c(0, 1, 3, 4, 6, 7)  # For chull and vaf methods
    
    missing_results <- 0
    for (method in expected_methods) {
        method_count <- sum(results_df$method == method, na.rm = TRUE)
        if (method %in% c("chull", "vaf")) {
            expected_count <- length(unique(results_df$sim_id)) * length(expected_complex)
        } else {
            expected_count <- length(unique(results_df$sim_id))
        }
        log_with_time(paste(method, "results:", method_count, "/ expected ~", expected_count))
        missing_results <- missing_results + max(0, expected_count - method_count)
    }
    
    if (missing_results > 0) {
        log_with_time(paste("WARNING: Approximately", missing_results, "results may be missing"))
    }
    
    log_with_time("Result collection completed successfully!")
    
    return(results_df)
}

# Run the collection
if (!interactive()) {
    # If script is run from command line
    results <- collect_all_results()
} else {
    # If run interactively, you can call:
    # results <- collect_all_results()
    log_with_time("Script loaded. Run collect_all_results() to start collection.")
}