#!/bin/bash
# waits for the stage-A queue to finish, then runs the main queue (2 processes)
cd "$(dirname "$0")"
while pgrep -f "queue.sh jobs_r2_stageA.txt" > /dev/null; do sleep 15; done
P=2 ./queue.sh jobs_r2_main.txt
