#!/bin/bash

# Script to estimate memory usage from previous job runs
# Usage: ./estimate_memory.sh <job_id>

JOB_ID=$1

if [ -z "$JOB_ID" ]; then
    echo "Usage: $0 <job_id>"
    echo "This will show memory usage statistics for the specified job"
    exit 1
fi

echo "Memory usage statistics for job $JOB_ID:"
echo "========================================="

# Get memory usage for all tasks in the array job
sacct -j $JOB_ID --format=JobID,MaxRSS,ReqMem,State,ExitCode -P | \
while IFS='|' read -r jobid maxrss reqmem state exitcode; do
    if [[ "$jobid" == *"_"* ]]; then  # Array task
        task_num=$(echo "$jobid" | sed 's/.*_//')
        
        # Convert MaxRSS to GB if it's not empty
        if [ -n "$maxrss" ] && [ "$maxrss" != "" ]; then
            if [[ "$maxrss" == *"K" ]]; then
                mem_gb=$(echo "$maxrss" | sed 's/K//' | awk '{print $1/1024/1024}')
            elif [[ "$maxrss" == *"M" ]]; then
                mem_gb=$(echo "$maxrss" | sed 's/M//' | awk '{print $1/1024}')
            elif [[ "$maxrss" == *"G" ]]; then
                mem_gb=$(echo "$maxrss" | sed 's/G//')
            else
                mem_gb=$(echo "$maxrss" | awk '{print $1/1024/1024/1024}')
            fi
            
            printf "Task %3s: Used %6.2f GB, Requested %s, Status: %s" "$task_num" "$mem_gb" "$reqmem" "$state"
            
            if [[ "$exitcode" == *"125"* ]] || [[ "$exitcode" == *"137"* ]]; then
                echo " (OOM KILL)"
            else
                echo ""
            fi
        fi
    fi
done

echo ""
echo "Recommendations:"
echo "==============="

# Calculate 95th percentile memory usage
max_mem=$(sacct -j $JOB_ID --format=MaxRSS --noheader -P | \
    grep -v "^$" | \
    sed 's/K$//' | sed 's/M$/*1024/' | sed 's/G$/*1024*1024/' | \
    bc -l 2>/dev/null | \
    sort -n | \
    tail -n 1)

if [ -n "$max_mem" ] && [ "$max_mem" != "0" ]; then
    recommended_gb=$(echo "$max_mem / 1024 / 1024 * 1.2" | bc -l | cut -d. -f1)
    echo "- Recommended memory allocation: ${recommended_gb}G (20% buffer)"
    echo "- For safety, consider: $((recommended_gb + 10))G"
else
    echo "- No memory usage data available (job may have failed immediately)"
    echo "- Try increasing from current allocation by 50-100%"
fi
