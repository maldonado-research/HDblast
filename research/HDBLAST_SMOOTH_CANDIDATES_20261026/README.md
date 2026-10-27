# Smooth approximate full23 candidates

The [previous audit](../HDBLAST_SHELL_DERIVATIVE_C2_20261026/) found a
second-derivative jump in the existing approximate response. This checkpoint
constructs two distinct smooth nominal replacements and accounts for their
changes to the inherited C¹ error estimate. They are approximate candidates;
neither is a complete second-order solution or a numerical C² error certificate.

## 1. Executable P1-based reconstruction

[The construction](construction/NOTE.md) convolves the frozen piecewise-linear
displacement and velocity profiles with a compact even kernel, yielding C⁴
profiles, and defines the radial profile as the actual derivative of the new
displacement. The original mixed reflections and zero shell collar are retained.
The script evaluates the new free term plus first feedback, checks boundary
conditions at a shell point, and an independent exact calculation finds zero
second-derivative jump at the earlier kink point. The candidate still omits the
affine bulk and shell sources; its full boundary defect is not established to vanish.

[Its conditional error transfer](transfer/P1_TRANSFER.md) gives the conservative
whole-domain C¹ bound `X(U−V_s)<9.711e31`, relying on the previous certificates.
This is substantially looser than the earlier `1e29` target. It is an upper
bound, not a measured error; it does not imply a useful sign or physical signal.

## 2. A tighter integral-defined reconstruction

[A separate construction](transfer/NOMINAL_CONVOLUTION.md) smooths the frozen
expression-defined nominal displacement directly with a narrower kernel, while
smoothing the piecewise-linear velocity. Its conditional C¹ transfer gives
`X(U−V_s)<2.184e28`, below the earlier `1e29` target. This figure applies to
this second candidate only; it is not the error bound for the executable
P1-based candidate above.

Validated quadrature checks enclose one component of the smoothed nominal
initial displacement and its derivative, and [the free-wave part at the old
kink point](construction/FREEWAVE_NOTE.md). They do not evaluate the complete
integral-defined candidate, its feedback, the actual field or the Einstein
residual. A complete validated evaluator remains future work.

## One localized bulk derivative

[The scalar potential note](bulk/POTENTIAL_NOTE.md) uses an inherited
actual-first-field cell certificate to bound the derivative of one isolated
term in the scalar bulk source. Other terms may cancel it; the result does not
bound or determine the sign of the complete bulk source derivative.

## Reproduction

All scripts use Python 3's standard library, authenticated prior inputs and
explicit failure checks. From the repository root:

```sh
C=research/HDBLAST_SMOOTH_CANDIDATES_20261026
python "$C/construction/smooth_candidate.py" --repo-root .
python "$C/checks/check_transfer.py" --repo-root .
python "$C/checks/check_smoothed_kink.py" --candidate "$C/construction/smooth_candidate.py" --repo-root .
python "$C/checks/check_nominal_transfer.py" --repo-root .
python "$C/checks/check_nominal_point.py" --repo-root .
python "$C/checks/check_smoothed_freewave.py" --repo-root .
python "$C/checks/check_potential_derivative.py" --repo-root .
```

Run with `-O` as well; normal and optimized output must agree. The point check
uses the pinned earlier interval producer, which authenticates its upstream
roots. The checks do not repeat the full upstream field certificates. Internal
AI-assisted reviews are not external specialist review or proof-assistant
verification. [Focused literature](RELATED.md) provides context and does not
establish originality or supply the missing coupled-system estimate.

The next steps are a validated evaluator for the tighter candidate, bounds for
all differentiated bulk terms and numerical second-derivative errors, and a
complete coupled shell/bulk residual audit. The exact nonlinear spacetime and a
higher-dimensional origin of the Big Bang remain unestablished.
