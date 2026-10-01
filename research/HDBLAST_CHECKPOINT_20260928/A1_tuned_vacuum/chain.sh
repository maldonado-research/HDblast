#!/bin/bash
cd "$(dirname "$0")"
while pgrep -f "queue.sh jobs_pre.txt" > /dev/null; do sleep 10; done
./queue.sh jobs_main.txt > queue_main.out 2>&1
