# Independent audit of reported postcontact arithmetic

The report-data algebra passes 202 exact rational and directed-rounding checks. The reconstructed screen interval is `[-0.4043111633297405, -0.014601815851683174]` (decimal display only), which excludes zero. Producer normal and optimized reports are canonically identical. Running this separate audit normally and with `python -O` also produces canonically identical output.

The audit independently implements rational interval operations at the reported denominator `2^112`. It checks every global source-cover rectangle, all 33 characteristic segment geometries and lengths, reflection at `x=9/8`, both endpoint characteristic contributions, their accumulated integral, and the `-G(9/8)` Neumann contribution. It reconstructs the complete feedback cap algebra, then recomputes the fixed quadratic screen expression from the reported apex jet intervals and combines it with the reported full-feedback derivative using the second-parameter-derivative normalization.

The feedback contraction parameter is about `5.4265e-5`; the complete derivative tail is about `5.6063e-6`. All declared interval widths and charged feedback error are reproduced exactly.

Five altered report-data controls are rejected with exact failed witnesses: changed Neumann endpoint data, a negative-side unreflected segment, a dropped fixed quadratic, disabled feedback, and an extra factor of one half in the linear observable. These controls test report consistency only; they do not exercise altered scientific source files.

**Scope:** this is an independent audit of report aggregation and observable arithmetic. It imports neither the postcontact scientific producer nor numerical model-source functions. It accepts each reported source-cell interval and apex first-field interval as input data, so it does not independently certify those values, the model-source closure, the analytic tail proof, or a physical claim. Source numerical recomputation is documented in `../PUBLIC_NUMERIC_REVIEW.md`.

Artifacts: `audit_report.py`, `ALGEBRA_NORMAL.json`, and `ALGEBRA_OPTIMIZED.json`. Each JSON records input report hashes before and after the run, the audit source hash, all checks, and exact negative-control witnesses.

Portable replay uses `python audit_report.py --producer-dir PRODUCER_REPORTS --output FRESH_OUTPUT.json`; repeat with `python -O` and a different fresh output. `PRODUCER_REPORTS` contains `POSTCONTACT_NORMAL.json` and `POSTCONTACT_OPTIMIZED.json`. Python 3.9 or newer and its standard library are sufficient. Five distinct altered report-data controls run in process in each mode, separately from scientific implementation and authentication mutation subprocesses.
