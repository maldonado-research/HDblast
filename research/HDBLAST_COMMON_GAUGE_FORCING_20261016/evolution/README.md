# Full formal second response

`FULL_SECOND_RESPONSE_CONSTRUCTION.md` defines the complete second metric/scalar
coefficient by a convergent mixed-boundary integral series on the retained causal
collar. It keeps all287 angular section modes and every physical/gauge boundary
row. The registered contraction factor is0.002109. Its quantitative bounds cover
values and first scaled derivatives; smoothness and subsidiary completion are
proved separately without claiming numerical higher-derivative bounds.

The construction is conditional on the retained background and actual first
solution. It does not give an exact nonlinear family, an evaluated waveform,
a finite-amplitude remainder, or empirical cosmological evidence.

Run the exact arithmetic join checker from a portable research-input tree:

```sh
python check_second_response.py --research-root retained-inputs --output fresh.json
python -O check_second_response.py --research-root retained-inputs --output fresh-optimized.json
```

`INPUT_PINS.json` specifies the relative paths and hashes below that root. All
scientific inputs are read only. The checker refuses to overwrite its output.
The accepted golden reports are `SECOND_RESPONSE_FINAL_NORMAL.json` and
`SECOND_RESPONSE_FINAL_OPTIMIZED.json`, which are byte identical. The earlier
`SECOND_RESPONSE_NORMAL.json` and `SECOND_RESPONSE_OPTIMIZED.json` are preserved
workspace drafts and are excluded from the publication inventory.

Six direct rejection controls are exposed by `--inject`: `omit_TT`,
`coordinate_Robin_in_frame`, `drop_scalar_DC`, `wrong_metric_force_sign`,
`wrong_lift_sign`, and `omit_trace_weight`. Each must exit nonzero for its named
reason and leave no report. `DIRECT_CONTROL_REPLAY.json` records their twelve
normal/optimized checks. These finite checks supplement the analytic proof;
they do not replace it or replay every inherited certificate.

The authoritative selected publication files and SHA256 values are in
`EVOLUTION_PRODUCER_FREEZE.json`. Internal independent scientific acceptance is
recorded separately by the fresh referee.
