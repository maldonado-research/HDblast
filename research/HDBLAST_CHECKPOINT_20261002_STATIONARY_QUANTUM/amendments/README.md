# Mechanical validator correction after registered integration

The original 48 registered solves and three wrong-model controls completed with
producer exit0 and no solver failures. The frozen validator exited1 and reports
FAIL: its diagnostic metadata keyword `condition=` collided with the acceptance
helper's `condition` argument, producing 48 TypeErrors. Original code, results,
checks and logs remain unchanged and are included.

`code/validate_stationary_v2.py` changes exactly that metadata keyword to
`condition_number=`. All numerical equations, inputs, cases, gates, thresholds
and producer data remain byte-identical. The amendment pins failed and corrected
sources. It is published before the corrected recheck, after the numerical
values already existed. It is a disclosed implementation repair, not a new
prospective scientific protocol or an originally passing validator.
