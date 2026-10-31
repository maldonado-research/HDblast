# Separate exact Cauchy-packet initial experiment

`W1_exact_Cauchy_packet` passes the frozen initial and local gates at 128-bit
outward precision. The complete [0,7/4] velocity value/radial-derivative
certificate uses 49 exact leaves, zero splits and zero unresolved cells;
the frozen 4096-leaf/depth 16 budgets were respected. The 192-bit fallback
was unnecessary. This is a distinct velocity candidate; all older
candidate functions and published artifacts were preserved.

| Complete velocity quantity in full23 H | Certified upper |
| --- | ---: |
| Absolute nominal candidate v | 1.2828974160481017e24 |
| Absolute nominal candidate v_r | 7.671822936222669e25 |
| Paired actual-minus-nominal v error | 3.73804140889402e11 |
| Paired actual-minus-nominal v_r error | 2.2352381383293816e13 |

The two new error certificates are each below half their respective old
C2 velocity error certificates. Exact q/E functional differences vanish;
both background coefficient and actual ODE remainder charges remain.
The complete absolute candidate bounds are separate quantities and need
not improve when velocity changes.

The fresh local v_r certificate at (0,65/64) has target H-norm lower
D=2.250277124262796e-12, new error upper U=6.555352734793914e-25,
and old smoothed-P1 error lower L=6.194797618620037e-6.
Thus U/D is at most approximately 2.9132e-13 and U<L/2. The old error
lower uses the exact half-sum of adjacent P1 secants, freshly compared
with the same actual analytic target in all 23 Gram components. This is
a proved local true-error improvement. The whole-domain comparisons
concern upper certificates. Relative quotients at the flat support edge
are undefined.

`FREEZE.json` and `FREEZE_RECEIPT.json` predate the first benchmark.
`producer.py` authenticates inherited sources and their manifests before
project imports and again after evaluation. `INITIAL_NORMAL.json` stores
all 23 component leaf caps, separated error channels, exact coverage,
denominator separation, local contribution rows and endpoint traces.
`REGULARITY_AND_LOCAL_PROOF.md` gives the completed-component affine
velocity exterior, inherited parity, exact flat collar, C1 velocity and
C2 displacement proof, derivative inventory and local convolution proof.

`CANONICAL_IDENTITY.json` records byte-identical normal/optimized initial,
control and combined reports. `CONTROLS_NORMAL.json` contains 14 exact
mathematical or proof-ledger rejection witnesses. All 30 subprocesses in
`CONTROL_SUBPROCESSES.json` reject the corresponding corruption under
normal and optimized Python, including source authentication failure
before imports, without producing a success artifact. These controls
do not replace the independent evaluator in the separate directory
`/workspace/hdblast-exact-velocity-independent-31c`.

`COMBINED_INITIAL_NORMAL.json` combines only the final separately proved
absolute displacement refinement with this fresh velocity certificate.
The published displacement function and its paired errors are inherited
exactly. The source reports remain separate and unmodified. This assembly
has its own source hashes, exact field provenance and optimized identity.
It is an initial ledger; response calculations belong to the separate
response experiment.

No Git operation, publication, natural residual, invariant observable,
nonlinear shell-contact solution, physical evidence or novelty claim is
part of this experiment.
