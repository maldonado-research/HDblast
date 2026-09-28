#!/bin/sh
# Full reproduction (about 20-25 min on one core; the four scans can run in parallel on 2 cores).
set -e
cd "$(dirname "$0")"
mkdir -p runs logs
python3 series_expansion.py 8            > logs/series_run.log      # exact series, ~50 s
python3 gegenbauer_structure.py                                      # exact identities, ~20 s
python3 run_bvp_scan.py reg reg 50 50 0.125 0.25 0.0001,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064,0.0003,0.001,0.003,0.01,0.03,0.1 > logs/scan_reg.log
python3 run_bvp_scan.py cm04 -0.4 50 50 0.125 0.25 0.0001,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064 > logs/scan_cm04.log
python3 run_bvp_scan.py cp13 1.3 50 50 0.125 0.25 0.0001,0.0002,0.0004,0.0008,0.0016,0.0032,0.0064 > logs/scan_cp13.log
python3 run_bvp_scan.py c0 0 50 50 0.125 0.25 0.0001,0.001,0.01,0.1 > logs/scan_c0.log
python3 run_bvp_scan.py reg_conv reg 65 64 0.1 0.2 0.0001,0.0008,0.0064,0.1 > logs/scan_reg_conv.log
python3 analyze_series_vs_bvp.py         > logs/analysis.log
python3 resonance_probe.py               > logs/resonance.log
python3 spectrum_leading_order.py        > logs/spectrum.log
python3 make_manifest.py
