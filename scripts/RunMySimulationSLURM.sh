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

# Create output directory for individual job logs
mkdir -p job_logs

# Track background process PIDs
pids=()

echo "Starting parallel R jobs..."

# Run R scripts in parallel using background processes
for i in {1..10}; do
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

echo "All jobs launched. Waiting for completion..."
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
echo "Total jobs: 10"
echo "Failed jobs: $failed_jobs"
echo "Successful jobs: $((10 - failed_jobs))"

# Print individual job logs for failed jobs if any
if [ $failed_jobs -gt 0 ]; then
    echo "="*50
    echo "Failed job logs:"
    for i in {1..10}; do
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