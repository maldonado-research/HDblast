#!/bin/bash
cd "$(dirname "$0")"
while pgrep -f "chain.sh" > /dev/null || pgrep -f "queue.sh jobs_main.txt" > /dev/null; do sleep 15; done
./queue.sh jobs_ctl.txt > queue_ctl.out 2>&1
