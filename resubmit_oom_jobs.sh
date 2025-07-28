#!/bin/bash

# Script to monitor and resubmit jobs that fail with OOM errors
# Usage: ./resubmit_oom_jobs.sh <job_id>

JOB_ID=$1
SCRIPT_PATH="/home/iyellan/tfbs_age_estimation/run_tfbs_liftover_parse.sh"
MAX_RETRIES=3

if [ -z "$JOB_ID" ]; then
    echo "Usage: $0 <job_id>"
    echo "Example: $0 12345"
    exit 1
fi

check_and_resubmit() {
    local job_id=$1
    local retry_count=${2:-0}
    
    echo "Monitoring job $job_id (retry $retry_count)..."
    
    # Wait for job to complete
    while squeue -j $job_id &>/dev/null; do
        sleep 60
    done
    
    # Check exit codes for all array tasks
    failed_tasks=()
    while IFS= read -r line; do
        task_id=$(echo "$line" | cut -d'|' -f1)
        exit_code=$(echo "$line" | cut -d'|' -f2)
        
        # Check for OOM error (exit code 125) or memory-related failures
        if [[ "$exit_code" == *"125"* ]] || [[ "$exit_code" == *"137"* ]] || [[ "$exit_code" == *"9"* ]] || [[ "$exit_code" == *"143"* ]]; then
            failed_tasks+=("$task_id")
            echo "Task $task_id failed with memory-related error (exit code: $exit_code)"
        elif [[ "$exit_code" != "0:0" ]] && [[ "$exit_code" != "" ]]; then
            echo "Task $task_id failed with non-memory error (exit code: $exit_code) - NOT resubmitting"
        fi
    done < <(sacct -j $job_id --format=JobID,ExitCode --noheader --parsable2 | grep -E "_[0-9]+\|")
    
    if [ ${#failed_tasks[@]} -eq 0 ]; then
        echo "All tasks completed successfully!"
        return 0
    fi
    
    if [ $retry_count -ge $MAX_RETRIES ]; then
        echo "Maximum retries ($MAX_RETRIES) reached. Failed tasks: ${failed_tasks[*]}"
        return 1
    fi
    
    # Calculate new memory allocation (increase by 50% each retry)
    base_mem=30
    new_mem=$(echo "$base_mem * (1.5^($retry_count + 1))" | bc -l | cut -d. -f1)
    
    echo "Resubmitting ${#failed_tasks[@]} failed tasks with ${new_mem}G memory..."
    
    # Create array string from failed tasks
    array_string=$(IFS=','; echo "${failed_tasks[*]}")
    
    # Resubmit only the failed tasks
    new_job_id=$(sbatch --array="$array_string" --mem="${new_mem}G" --parsable "$SCRIPT_PATH")
    
    if [ $? -eq 0 ]; then
        echo "Resubmitted as job $new_job_id"
        # Recursively monitor the new job
        check_and_resubmit "$new_job_id" $((retry_count + 1))
    else
        echo "Failed to resubmit job"
        return 1
    fi
}

check_and_resubmit "$JOB_ID"
