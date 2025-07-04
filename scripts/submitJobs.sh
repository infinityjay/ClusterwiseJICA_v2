#!/bin/bash

# submitJob.sh - Submit individual jobs for each simulation file

module load R/4.4.0-gfbf-2023a

# Set the path to simulation data folder
SIMULATION_DATA_PATH="/home/s4162315/cjica/ClusterwiseJICA_v2/simulation_data"
# # use test data
# SIMULATION_DATA_PATH="/home/s4162315/cjica/ClusterwiseJICA_v2/data_test"

# remove the previous logs
rm -rf ./logs

# Create directories for logs if they don't exist
mkdir -p logs
mkdir -p results

# Get all .RData files in the simulation_data folder
DATA_FILES=($(ls ${SIMULATION_DATA_PATH}/sim_*.RData 2>/dev/null))

if [ ${#DATA_FILES[@]} -eq 0 ]; then
    echo "Error: No simulation files found in ${SIMULATION_DATA_PATH}"
    exit 1
fi

echo "Found ${#DATA_FILES[@]} simulation files"
echo "Submitting jobs..."

# Submit a job for each data file
for FILE_PATH in "${DATA_FILES[@]}"; do
    # Extract just the filename from the full path
    FILENAME=$(basename "$FILE_PATH")
    JOB_NAME="cjica_${FILENAME%.RData}"
    
    echo "Submitting job for: $FILENAME"
    
    # Submit job with proper Slurm parameters
    sbatch --job-name="$JOB_NAME" \
           --nodes=1 \
           --ntasks=1 \
           --cpus-per-task=1 \
           --mem=2G \
           --time=23:59:00 \
           --partition="cpu-long" \
           --output="logs/%x_%j.out" \
           --error="logs/%x_%j.err" \
           --mail-type=FAIL \
           --wrap="module load R/4.4.0-gfbf-2023a; cd /home/s4162315/cjica/ClusterwiseJICA_v2/; Rscript ./scripts/run_analysis.R $FILENAME"
    
    # Small delay to avoid overwhelming the scheduler
    sleep 0.1
done

echo "All ${#DATA_FILES[@]} jobs submitted!"
echo "Job status: squeue -u s4162315"
echo "Check running job numbers: squeue --user=s4162315 | grep ' R ' | wc -l"
echo "Cancel all jobs: scancel --user=s4162315"