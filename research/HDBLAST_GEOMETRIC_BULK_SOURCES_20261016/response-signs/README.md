# Evolved formal scalar-mean signs

The accepted C1 response bound quantitatively preserves two complete initial
scalar-mean sign cells for `10^-9<=t<=5*10^-9`. The mean exceeds `1.85e30`
on the positive cell and is below `-1.85e30` on the negative cell. These are
conditional statements about the selected formal coefficient in its registered
gauge and initial prescription, not a computed waveform or invariant force.

See `EVOLVED_SCALAR_MEAN_SIGNS.md` for the proof and scope. Replay from this
directory, using unused output paths:

```sh
python check_evolved_scalar_signs.py --output EVOLVED_SCALAR_SIGNS_NORMAL.json
python -O check_evolved_scalar_signs.py --output EVOLVED_SCALAR_SIGNS_OPTIMIZED.json
python run_sign_controls.py --output EVOLVED_SIGN_CONTROLS.json
```

All runtime inputs are local relative files authenticated by `INPUT_PINS.json`.
No third-party dependency, external service, or parent directory is required.
