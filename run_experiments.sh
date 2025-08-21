#!/bin/bash

# Configuration
EXECUTABLE="./cmake-build-release/BISAM"  # Replace with your actual executable path
MAX_RETRIES=10  # Maximum number of retries per dataset
TOTAL_RUNS=20  # Total runs per dataset
TIMING_CSV="./simulations/step_size_fixed_timing.csv"  # Path to timing CSV
RUN_NAME_PREFIX="jenkins"

DATASETS=(
    "rootfind_stepsize_050"
    "rootfind_stepsize_075"
    "rootfind_stepsize_100"
    "rootfind_stepsize_150"
    "rootfind_stepsize_300"
    "rootfind_stepsize_500"
)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Log file
LOG_FILE="dataset_runs_$(date +%Y%m%d_%H%M%S).log"

# Function to get run name from dataset name and prefix
get_run_name() {
    local dataset_name="$1"
    local step_size="000"  # Default fallback

    # Extract step size from dataset name (e.g., "rootfind_stepsize_050" -> "050")
    local pos=$(echo "$dataset_name" | awk -F'_' '{print $NF}')
    if [ ! -z "$pos" ]; then
        step_size="$pos"
    fi

    echo "${RUN_NAME_PREFIX}_${step_size}"
}

# Function to count completed runs for a dataset
count_completed_runs() {
    local dataset_name="$1"
    local run_name=$(get_run_name "$dataset_name")

    if [ -f "$TIMING_CSV" ]; then
        # Count lines in CSV for this specific run_name AND dataset combination
        # CSV format assumed: something,run_name,dataset_name,...
        local count=$(grep "^[^,]*,${run_name},${dataset_name}," "$TIMING_CSV" 2>/dev/null | wc -l)
        echo $count
    else
        echo 0
    fi
}

# Function to get the next run number to start from
get_next_run() {
    local dataset_name="$1"
    local completed=$(count_completed_runs "$dataset_name")
    local next_run=$((completed + 1))

    if [ $next_run -gt $TOTAL_RUNS ]; then
        echo 0  # All runs completed
    else
        echo $next_run
    fi
}

echo -e "${BLUE}=== BISAM Dataset Runner with Resume ===${NC}"
echo "Log file: $LOG_FILE"
echo "Executable: $EXECUTABLE"
echo "Timing CSV: $TIMING_CSV"
echo "Run name prefix: $RUN_NAME_PREFIX"
echo "Total runs per dataset: $TOTAL_RUNS"
echo "Max retries per dataset: $MAX_RETRIES"
echo "Datasets to run: ${#DATASETS[@]}"
echo ""

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo -e "${RED}Error: Executable '$EXECUTABLE' not found!${NC}"
    echo "Please update the EXECUTABLE variable in this script to point to your compiled program."
    exit 1
fi

# Function to run a single dataset with retries
run_dataset() {
    local dataset_name="$1"
    local start_run="$2"
    local attempt=1

    while [ $attempt -le $MAX_RETRIES ]; do
        echo -e "${BLUE}--- Running dataset: $dataset_name from run $start_run (attempt $attempt/$MAX_RETRIES) ---${NC}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting $dataset_name from run $start_run (attempt $attempt)" >> "$LOG_FILE"

        # Capture output to a temporary file while also showing it real-time
        local temp_output=$(mktemp)
        local start_time=$(date +%s)

        # Run the program with timeout, showing output in real-time AND capturing to temp file
        echo -e "${YELLOW}Starting execution (output will appear below):${NC}"
        if timeout 3600 $EXECUTABLE "$dataset_name" "$start_run" "$RUN_NAME_PREFIX" 2>&1 | tee "$temp_output"; then
            local exit_code=${PIPESTATUS[0]}  # Get exit code of the actual program, not tee
            local end_time=$(date +%s)
            local run_time=$((end_time - start_time))

            # Copy captured output to log file as well
            cat "$temp_output" >> "$LOG_FILE"

            # Check if run actually completed successfully
            local success_markers=(
                "Completed multithreading performance test"
                "Final save of all results"
                "Successfully completed"
            )

            local completion_found=false
            for marker in "${success_markers[@]}"; do
                if grep -q "$marker" "$temp_output"; then
                    completion_found=true
                    break
                fi
            done

            # Count timing measurements in this run
            local timing_count=$(grep -o '[0-9]\+\.[0-9]\+ms' "$temp_output" | wc -l)
            local expected_runs=$((TOTAL_RUNS - start_run + 1))

            if [ "$completion_found" = true ] && [ "$timing_count" -ge "$expected_runs" ]; then
                echo -e "${GREEN}✓ Dataset $dataset_name completed successfully (${timing_count} runs, ${run_time}s)${NC}"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] $dataset_name completed successfully (${timing_count} runs, ${run_time}s)" >> "$LOG_FILE"
                echo "" >> "$LOG_FILE"
                rm "$temp_output"
                return 0
            else
                echo -e "${RED}✗ Dataset $dataset_name incomplete (exit: $exit_code, runs: ${timing_count}/${expected_runs}, time: ${run_time}s)${NC}"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] $dataset_name incomplete (exit: $exit_code, runs: ${timing_count}/${expected_runs}, time: ${run_time}s)" >> "$LOG_FILE"

                # Show last few lines of output for debugging
                echo -e "${YELLOW}Last 5 lines of captured output:${NC}"
                tail -5 "$temp_output" | while read line; do
                    echo "  $line"
                done
            fi
        else
            local exit_code=${PIPESTATUS[0]}  # Get exit code of the actual program, not tee
            local end_time=$(date +%s)
            local run_time=$((end_time - start_time))

            # Copy captured output to log file
            cat "$temp_output" >> "$LOG_FILE"

            # Handle specific exit codes
            if [ $exit_code -eq 124 ]; then
                echo -e "${RED}✗ Dataset $dataset_name timed out after 1 hour${NC}"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] $dataset_name timed out after 1 hour" >> "$LOG_FILE"
            elif [ $exit_code -eq 137 ]; then
                echo -e "${RED}✗ Dataset $dataset_name killed by system (likely out of memory - SIGKILL)${NC}"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] $dataset_name killed by system (OOM - SIGKILL)" >> "$LOG_FILE"
            elif [ $exit_code -eq 138 ]; then
                echo -e "${RED}✗ Dataset $dataset_name terminated (likely out of memory - SIGTERM)${NC}"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] $dataset_name terminated (OOM - SIGTERM)" >> "$LOG_FILE"
            elif [ $exit_code -eq 139 ]; then
                echo -e "${RED}✗ Dataset $dataset_name segmentation fault${NC}"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] $dataset_name segmentation fault" >> "$LOG_FILE"
            else
                echo -e "${RED}✗ Dataset $dataset_name failed with exit code $exit_code (${run_time}s)${NC}"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] $dataset_name failed with exit code $exit_code (${run_time}s)" >> "$LOG_FILE"
            fi

            # Show last few lines for debugging
            echo -e "${YELLOW}Last 5 lines of captured output:${NC}"
            tail -5 "$temp_output" | while read line; do
                echo "  $line"
            done
        fi

        rm "$temp_output"

        # Check how many runs were actually completed before the failure
        local new_completed=$(count_completed_runs "$dataset_name")
        local new_next_run=$(get_next_run "$dataset_name")

        if [ $new_next_run -eq 0 ]; then
            echo -e "${GREEN}Dataset $dataset_name actually completed all runs despite error${NC}"
            return 0
        elif [ $new_completed -gt $((start_run - 1)) ]; then
            echo -e "${YELLOW}Partial progress made: completed $new_completed/$TOTAL_RUNS runs${NC}"
            echo -e "${YELLOW}Will resume from run $new_next_run${NC}"
            start_run=$new_next_run
        fi

        if [ $attempt -lt $MAX_RETRIES ]; then
            echo -e "${YELLOW}Retrying in 10 seconds...${NC}"
            sleep 10
            ((attempt++))
        else
            echo -e "${RED}Max retries reached for dataset $dataset_name${NC}"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Max retries reached for $dataset_name" >> "$LOG_FILE"
            echo "" >> "$LOG_FILE"
            return 1
        fi
    done
}

# Main execution
start_time=$(date +%s)
successful_datasets=()
failed_datasets=()
skipped_datasets=()

echo -e "${BLUE}Starting dataset runs...${NC}"
echo ""

# Check existing progress first
echo -e "${BLUE}=== CHECKING EXISTING PROGRESS ===${NC}"
for dataset in "${DATASETS[@]}"; do
    completed=$(count_completed_runs "$dataset")
    next_run=$(get_next_run "$dataset")
    run_name=$(get_run_name "$dataset")

    if [ $next_run -eq 0 ]; then
        echo -e "${GREEN}✓ Dataset $dataset (${run_name}): All $TOTAL_RUNS runs completed${NC}"
    else
        echo -e "${YELLOW}○ Dataset $dataset (${run_name}): $completed/$TOTAL_RUNS runs completed, will start from run $next_run${NC}"
    fi
done
echo ""

for dataset in "${DATASETS[@]}"; do
    completed=$(count_completed_runs "$dataset")
    next_run=$(get_next_run "$dataset")
    run_name=$(get_run_name "$dataset")

    if [ $next_run -eq 0 ]; then
        echo -e "${GREEN}Skipping $dataset (${run_name}) - all runs already completed${NC}"
        skipped_datasets+=("$dataset")
        continue
    fi

    echo -e "${BLUE}Processing $dataset (${run_name}): $completed/$TOTAL_RUNS completed, starting from run $next_run${NC}"

    if run_dataset "$dataset" "$next_run"; then
        # Verify final completion after successful run
        final_completed=$(count_completed_runs "$dataset")
        if [ $final_completed -eq $TOTAL_RUNS ]; then
            successful_datasets+=("$dataset")
            echo -e "${GREEN}✓ Dataset $dataset (${run_name}) fully completed ($final_completed/$TOTAL_RUNS runs)${NC}"
        else
            failed_datasets+=("$dataset")
            echo -e "${YELLOW}⚠ Dataset $dataset (${run_name}) partially completed ($final_completed/$TOTAL_RUNS runs)${NC}"
        fi
    else
        failed_datasets+=("$dataset")
        final_completed=$(count_completed_runs "$dataset")
        echo -e "${RED}✗ Dataset $dataset (${run_name}) failed ($final_completed/$TOTAL_RUNS runs completed)${NC}"
    fi

    # Add a pause between datasets to let system recover
    if [ ${#DATASETS[@]} -gt 1 ]; then
        echo -e "${YELLOW}Pausing 10 seconds before next dataset...${NC}"
        sleep 10
        echo ""
    fi
done

# Summary
end_time=$(date +%s)
total_time=$((end_time - start_time))

echo -e "${BLUE}=== EXECUTION SUMMARY ===${NC}"
echo "Total execution time: ${total_time} seconds"
echo "Run name prefix: $RUN_NAME_PREFIX"
echo "Successful datasets: ${#successful_datasets[@]}"
echo "Failed datasets: ${#failed_datasets[@]}"
echo "Skipped datasets (already complete): ${#skipped_datasets[@]}"
echo ""

if [ ${#skipped_datasets[@]} -gt 0 ]; then
    echo -e "${GREEN}Skipped (already complete):${NC}"
    for dataset in "${skipped_datasets[@]}"; do
        run_name=$(get_run_name "$dataset")
        echo "  ✓ $dataset ($run_name)"
    done
    echo ""
fi

if [ ${#successful_datasets[@]} -gt 0 ]; then
    echo -e "${GREEN}Successful:${NC}"
    for dataset in "${successful_datasets[@]}"; do
        run_name=$(get_run_name "$dataset")
        echo "  ✓ $dataset ($run_name)"
    done
    echo ""
fi

if [ ${#failed_datasets[@]} -gt 0 ]; then
    echo -e "${RED}Failed/Incomplete:${NC}"
    for dataset in "${failed_datasets[@]}"; do
        completed=$(count_completed_runs "$dataset")
        run_name=$(get_run_name "$dataset")
        echo "  ✗ $dataset ($run_name) ($completed/$TOTAL_RUNS runs completed)"
    done
    echo ""
    echo -e "${YELLOW}Check $LOG_FILE for detailed error information.${NC}"
    echo -e "${YELLOW}You can re-run this script to resume from where it left off.${NC}"
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Execution completed. Success: ${#successful_datasets[@]}, Failed: ${#failed_datasets[@]}, Skipped: ${#skipped_datasets[@]}" >> "$LOG_FILE"

# Exit with error code if any datasets failed
if [ ${#failed_datasets[@]} -gt 0 ]; then
    exit 1
else
    exit 0
fi