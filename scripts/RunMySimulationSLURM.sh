#!/bin/bash
#SBATCH --job-name=test_R_parallel
#SBATCH -N 1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition="cpu-short"
#SBATCH --time=3:59:00
#SBATCH --mem=1G

# Print job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Load R module
echo "Loading R module..."
module load R/4.4.0-gfbf-2023a

# Verify R is loaded
echo "R version:"
R --version | head -1

# Create scripts directory if it doesn't exist
mkdir -p scripts

# Check if the R script exists
if [ ! -f "./scripts/mySimulationTest.R" ]; then
    echo "ERROR: mySimulationTest.R not found in ./scripts/"
    exit 1
fi

# Create output directories
mkdir -p job_logs
mkdir -p result/sim1

# Calculate total grid size and number of splits
echo "Calculating grid dimensions..."

# Grid parameters from R script:
# Nk <- c(50, 100)  # 2 values
# Qvect_list <- list(c(3,3), c(8,8))  # 2 values  
# E <- c(0.2, 0.6)  # 2 values
# cor <- c(0, .60, .94)  # 3 values
# VAF <- c(0.8, 1)  # 2 values
# rep <- 1:5  # 5 values

TOTAL_GRID_SIZE=$((2 * 2 * 2 * 3 * 2 * 5))  # = 240
ROWS_PER_SPLIT=10
NUM_SPLITS=$(((TOTAL_GRID_SIZE + ROWS_PER_SPLIT - 1) / ROWS_PER_SPLIT))  # Ceiling division

echo "Total grid combinations: $TOTAL_GRID_SIZE"
echo "Rows per split: $ROWS_PER_SPLIT"
echo "Number of splits needed: $NUM_SPLITS"

# Track background process PIDs
pids=()

echo "Starting parallel R jobs..."

# Run R scripts in parallel using background processes
for i in $(seq 1 $NUM_SPLITS); do
    echo "Starting job $i at $(date)"
    
    # Run with output redirection and error handling
    (
        echo "Job $i started on PID $$"
        Rscript ./scripts/mySimulationTest.R $i 2>&1
        exit_code=$?
        echo "Job $i finished with exit code: $exit_code"
        exit $exit_code
    ) > job_logs/job_${i}.log 2>&1 &
    
    # Store the PID
    pids+=($!)
    echo "Job $i launched with PID: ${pids[-1]}"
    
    # Small delay to prevent overwhelming the system
    sleep 1
done

echo "All $NUM_SPLITS jobs launched. Waiting for completion..."
echo "Active PIDs: ${pids[@]}"

# Wait for all background processes and check their exit status
failed_jobs=0
for i in "${!pids[@]}"; do
    job_num=$((i + 1))
    pid=${pids[$i]}
    
    echo "Waiting for job $job_num (PID: $pid)..."
    
    if wait $pid; then
        echo "Job $job_num completed successfully"
    else
        exit_code=$?
        echo "Job $job_num failed with exit code: $exit_code"
        ((failed_jobs++))
    fi
done

# Print summary
echo "="*50
echo "Job completion summary:"
echo "Total jobs: $NUM_SPLITS"
echo "Failed jobs: $failed_jobs"
echo "Successful jobs: $((NUM_SPLITS - failed_jobs))"

# Print individual job logs for failed jobs if any
if [ $failed_jobs -gt 0 ]; then
    echo "="*50
    echo "Failed job logs:"
    for i in $(seq 1 $NUM_SPLITS); do
        if [ -f "job_logs/job_${i}.log" ]; then
            last_line=$(tail -1 "job_logs/job_${i}.log")
            if [[ $last_line == *"exit code: 0"* ]]; then
                continue
            else
                echo "--- Job $i log (last 10 lines) ---"
                tail -10 "job_logs/job_${i}.log"
                echo ""
            fi
        fi
    done
fi

echo "All R scripts processing completed at: $(date)"

# Exit with error if any jobs failed
if [ $failed_jobs -gt 0 ]; then
    echo "Exiting with error due to $failed_jobs failed jobs"
    exit 1
else
    echo "All jobs completed successfully!"
    exit 0
fi