#!/bin/bash

# Script to submit job and automatically start monitoring in screen session
# Usage: ./submit_and_monitor.sh [screen_session_name]

SCRIPT_PATH="/home/iyellan/tfbs_age_estimation/run_tfbs_liftover_parse.sh"
MONITOR_SCRIPT="/home/iyellan/tfbs_age_estimation/resubmit_oom_jobs.sh"
SESSION_NAME=${1:-"tfbs_monitor_$(date +%Y%m%d_%H%M%S)"}

echo "Submitting job..."
job_id=$(sbatch --parsable "$SCRIPT_PATH")

if [ $? -eq 0 ]; then
    echo "Job submitted with ID: $job_id"
    echo "Starting monitoring in screen session: $SESSION_NAME"
    
    # Check if screen is available
    if ! command -v screen &> /dev/null; then
        echo "Error: screen is not available. Installing or loading screen module..."
        # Try to load screen module (common on HPC systems)
        module load screen 2>/dev/null || {
            echo "Please install screen or run manually:"
            echo "$MONITOR_SCRIPT $job_id"
            exit 1
        }
    fi
    
    # Start screen session and run monitoring script
    screen -dmS "$SESSION_NAME" bash -c "
        echo 'Starting monitoring for job $job_id at $(date)'
        echo 'Session: $SESSION_NAME'
        echo 'Use: screen -r $SESSION_NAME to reattach'
        echo '----------------------------------------'
        $MONITOR_SCRIPT $job_id
        echo 'Monitoring completed at $(date)'
        echo 'Press any key to exit...'
        read -n 1
    "
    
    echo "Screen session '$SESSION_NAME' started."
    echo "To check on progress:"
    echo "  screen -r $SESSION_NAME"
    echo "To list all sessions:"
    echo "  screen -ls"
    echo "To detach from session: Ctrl+A, then D"
    
else
    echo "Failed to submit job"
    exit 1
fi
