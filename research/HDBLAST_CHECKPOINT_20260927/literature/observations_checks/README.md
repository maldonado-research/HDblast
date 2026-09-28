# observations_checks: scripts and outputs for literature sweep 2 (observations)

Full description, results table, limitations and reproduction commands: `../README_OBSERVATIONS.md`. The deliverable is `../observations.md`.

Quick reproduction, from `research/HDBLAST_CHECKPOINT_20260927/literature/`:

```
python3 observations_checks/knee_and_scale_checks.py
python3 observations_checks/project_scan_1yr.py
python3 observations_checks/build_observations_md.py
python3 observations_checks/build_observations_md.py --selftest
```

All source entries are search-snippet-only (no paper pages could be opened). Every computed number is in `knee_and_scale_checks.json`.
