#!/bin/bash

# Enhanced monitoring script optimized for screen sessions
# Usage: ./resubmit_oom_jobs_screen.sh <job_id> [log_file]

JOB_ID=$1
LOG_FILE=${2:-"/home/iyellan/scratch/tfbs_ages/job_monitor_${JOB_ID}.log"}
SCRIPT_PATH="/home/iyellan/tfbs_age_estimation/run_tfbs_liftover_parse.sh"
MAX_RETRIES=4

# Create log directory if it doesn't exist
mkdir -p "$(dirname "$LOG_FILE")"

# Function to log with timestamp
log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to show progress
show_progress() {
    local job_id=$1
    local total_tasks=136  # Based on your array 0-135
    
    while squeue -j $job_id &>/dev/null; do
        # Get job status
        pending=$(squeue -j $job_id -t PD -h | wc -l)
        running=$(squeue -j $job_id -t R -h | wc -l)
        completed=$(sacct -j $job_id --state=CD,CA,F,TO,NF --noheader | wc -l)
        
        log_msg "Status - Pending: $pending, Running: $running, Completed: $completed/$total_tasks"
        sleep 300  # Check every 5 minutes
    done
}

if [ -z "$JOB_ID" ]; then
    echo "Usage: $0 <job_id> [log_file]"
    echo "Example: $0 12345"
    echo "Log will be written to: $LOG_FILE"
    exit 1
fi

log_msg "Starting enhanced monitoring for job $JOB_ID"
log_msg "Log file: $LOG_FILE"
log_msg "Screen session: ${STY:-'Not in screen'}"

check_and_resubmit() {
    local job_id=$1
    local retry_count=${2:-0}
    
    log_msg "Monitoring job $job_id (retry $retry_count)..."
    
    # Show progress while job is running
    show_progress $job_id
    
    log_msg "Job $job_id finished. Analyzing results..."
    
    # Check exit codes for all array tasks
    failed_tasks=()
    oom_tasks=()
    
    while IFS= read -r line; do
        if [[ -z "$line" ]]; then continue; fi
        
        task_id=$(echo "$line" | cut -d'|' -f1)
        exit_code=$(echo "$line" | cut -d'|' -f2)
        state=$(echo "$line" | cut -d'|' -f3)
        
        # Extract just the task number
        task_num=$(echo "$task_id" | sed 's/.*_//')
        
        # Check for OOM error (exit code 125) or memory-related failures
        if [[ "$exit_code" == *"125"* ]] || [[ "$exit_code" == *"137"* ]]; then
            failed_tasks+=("$task_num")
            oom_tasks+=("$task_num")
            log_msg "Task $task_num failed with OOM (exit code: $exit_code) - WILL RESUBMIT"
        elif [[ "$exit_code" == *"9"* ]] || [[ "$exit_code" == *"143"* ]]; then
            # Signal-based kills (often memory related)
            failed_tasks+=("$task_num")
            oom_tasks+=("$task_num")
            log_msg "Task $task_num killed by signal (exit code: $exit_code) - WILL RESUBMIT"
        elif [[ "$state" == "FAILED" ]] && [[ "$exit_code" != "0:0" ]]; then
            # Log other failures but don't resubmit them
            if [[ "$exit_code" == *"1"* ]]; then
                log_msg "Task $task_num failed with general error (exit code: $exit_code) - NOT RESUBMITTING (likely not memory-related)"
            else
                log_msg "Task $task_num failed (exit code: $exit_code, state: $state) - NOT RESUBMITTING (non-memory failure)"
            fi
        fi
    done < <(sacct -j $job_id --format=JobID,ExitCode,State --noheader --parsable2 | grep -E "_[0-9]+\|")
    
    log_msg "Analysis complete. Failed tasks: ${#failed_tasks[@]}, OOM tasks: ${#oom_tasks[@]}"
    
    if [ ${#failed_tasks[@]} -eq 0 ]; then
        log_msg "🎉 All tasks completed successfully!"
        return 0
    fi
    
    if [ $retry_count -ge $MAX_RETRIES ]; then
        log_msg "❌ Maximum retries ($MAX_RETRIES) reached. Failed tasks: ${failed_tasks[*]}"
        log_msg "Consider manually investigating these tasks or increasing memory further."
        return 1
    fi
    
    # Calculate new memory allocation (increase by 50% each retry)
    base_mem=15
    new_mem=$(echo "$base_mem * (1.5^($retry_count + 1))" | bc -l | cut -d. -f1)
    
    log_msg "🔄 Resubmitting ${#failed_tasks[@]} failed tasks with ${new_mem}G memory..."
    
    # Create array string from failed tasks
    array_string=$(IFS=','; echo "${failed_tasks[*]}")
    
    # Resubmit only the failed tasks
    new_job_id=$(sbatch --array="$array_string" --mem="${new_mem}G" --parsable "$SCRIPT_PATH")
    
    if [ $? -eq 0 ]; then
        log_msg "✅ Resubmitted as job $new_job_id with tasks: $array_string"
        log_msg "Memory allocation: ${new_mem}G"
        
        # Recursively monitor the new job
        check_and_resubmit "$new_job_id" $((retry_count + 1))
    else
        log_msg "❌ Failed to resubmit job"
        return 1
    fi
}

# Trap signals to log when script is interrupted
trap 'log_msg "Monitoring interrupted by user"; exit 130' INT TERM

# Start monitoring
check_and_resubmit "$JOB_ID"
final_status=$?

if [ $final_status -eq 0 ]; then
    log_msg "🎉 Monitoring completed successfully!"
else
    log_msg "❌ Monitoring completed with errors (exit code: $final_status)"
fi

log_msg "Final log location: $LOG_FILE"
echo ""
echo "Monitoring completed. Check log file for details:"
echo "  tail -f $LOG_FILE"
echo ""
echo "If running in screen, press Ctrl+A then D to detach, or any key to exit."
read -n 1 -s
