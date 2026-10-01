#!/bin/bash
cd "$(dirname "$0")"
while pgrep -f "chain2.sh" > /dev/null || pgrep -f "queue.sh jobs_ctl.txt" > /dev/null; do sleep 15; done
./queue.sh jobs_d3a.txt > queue_d3a.out 2>&1
./queue.sh jobs_d3b.txt > queue_d3b.out 2>&1
