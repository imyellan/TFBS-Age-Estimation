#!/bin/bash

# SLURM Exit Code Reference for TFBS Age Estimation Jobs
# =====================================================

echo "SLURM Exit Code Reference for Memory-Related Failures"
echo "======================================================"
echo ""
echo "MEMORY-RELATED (Will be automatically resubmitted):"
echo "  0:125  - Out of Memory (OOM) killer"
echo "  0:137  - Process killed (SIGKILL, often due to memory)"
echo "  0:9    - Process killed (SIGKILL)"
echo "  0:143  - Process terminated (SIGTERM, sometimes memory-related)"
echo ""
echo "NON-MEMORY RELATED (Will NOT be resubmitted):"
echo "  0:0    - Success"
echo "  0:1    - General error (script/command failure)"
echo "  0:2    - Misuse of shell builtins"
echo "  0:126  - Command not executable"
echo "  0:127  - Command not found"
echo "  0:128  - Invalid argument to exit"
echo "  0:130  - Script terminated by Ctrl+C"
echo ""
echo "SLURM SPECIFIC:"
echo "  CANCELLED  - Job was cancelled by user or admin"
echo "  TIMEOUT    - Job exceeded time limit"
echo "  NODE_FAIL  - Node failure"
echo "  PREEMPTED  - Job was preempted by higher priority job"
echo ""
echo "To check exit codes for a job:"
echo "  sacct -j <job_id> --format=JobID,ExitCode,State"
echo ""
echo "To see only failed tasks:"
echo "  sacct -j <job_id> --format=JobID,ExitCode,State | grep -v '0:0'"
echo ""
echo "Example usage:"
echo "  ./check_exit_codes.sh 12345"

if [ ! -z "$1" ]; then
    echo ""
    echo "Exit codes for job $1:"
    echo "====================="
    sacct -j $1 --format=JobID,ExitCode,State,MaxRSS --noheader | \
    while IFS='|' read -r jobid exitcode state maxrss; do
        if [[ "$jobid" == *"_"* ]]; then  # Array task
            task_num=$(echo "$jobid" | sed 's/.*_//')
            printf "Task %3s: Exit=%s, State=%s" "$task_num" "$exitcode" "$state"
            
            # Add interpretation
            if [[ "$exitcode" == "0:0" ]]; then
                echo " ✅ Success"
            elif [[ "$exitcode" == *"125"* ]] || [[ "$exitcode" == *"137"* ]] || [[ "$exitcode" == *"9"* ]]; then
                echo " 🔄 Memory-related (would resubmit)"
            elif [[ "$exitcode" == *"1"* ]]; then
                echo " ❌ General error (would NOT resubmit)"
            else
                echo " ❓ Other failure (would NOT resubmit)"
            fi
        fi
    done
fi
