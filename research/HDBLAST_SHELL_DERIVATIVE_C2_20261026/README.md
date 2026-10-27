# Shell source derivatives and the route to C² response control

This checkpoint derives a numerical bound for the derivative of the **actual
registered affine shell source** and identifies an obstruction in the saved
approximate response. Both are conditional on the cited frozen inputs and
normalizations. They are not a complete second-order solution.

[The full shell-source calculation](shell/FULL23_NOTE.md) differentiates the
existing eight-profile factorization of the 23-mode shell source. With x=t/ell,
it gives the uniform formal-norm bound

    ||G_x||_H < 58,807,538       on 0<=x<=3 at the shell.

This is a conservative absolute upper bound in the registered coefficient
norm. It is not a bound on a physical force, the entire boundary residual, or
the full response error. The saved approximation omits the affine G, so the
same figure bounds the derivative of that particular source-only omission;
feedback and bulk terms require separate estimates. For each of the two TT
rows, [a sharper calculation](shell/TT_NOTE.md) gives `|G_T,x| < 3.01e-19`.
These TT rows have zero linear Robin feedback and the saved TT candidate is
zero, so this also bounds their candidate boundary-defect derivatives, subject
to the inherited first-field and shell conditions.

[The candidate audit](candidate/NOTE.md) establishes a nonzero jump in one
second coordinate derivative of the existing C¹ candidate at an interior mesh
characteristic. The audit includes both its free term and first feedback; two
independent exact implementations agree. The jump belongs to the **approximate
candidate**, not the exact field, and does not invalidate its existing C¹
certificate. A classical pointwise second-derivative residual for this candidate
is undefined at that characteristic. A smooth replacement or a precisely
specified piecewise/weak treatment is needed for the stronger residual claim.

[A conditional error gate](gate/NOTE.md) shows how time-differentiating the
full 23-mode system can reuse the inherited mixed-boundary contraction to
control the second derivative of an error between suitably smooth fields. It
specifies the required bulk/source derivatives, coefficient derivatives,
initial jets and corner assumptions. The present P1-based candidate does not
satisfy the smoothness assumptions. No numerical full-solution C² or Einstein
residual bound is computed here.

## Reproduction

Python 3 and its standard library are sufficient. From the repository root:

```sh
python research/HDBLAST_SHELL_DERIVATIVE_C2_20261026/checks/check_full23_shell_derivative.py --repo-root .
python research/HDBLAST_SHELL_DERIVATIVE_C2_20261026/checks/check_tt_shell_derivative.py --repo-root .
python research/HDBLAST_SHELL_DERIVATIVE_C2_20261026/checks/check_candidate_kink.py --repo-root .
python research/HDBLAST_SHELL_DERIVATIVE_C2_20261026/checks/audit_b1.py --repo-root .
python research/HDBLAST_SHELL_DERIVATIVE_C2_20261026/checks/check_gate.py
```

All programs use explicit failures and can also run with `python -O`; normal
and optimized output must agree. They authenticate the inherited source or
profile files where used and perform exact rational arithmetic. The full-shell
check rebuilds and compares all 64 entries of the source Gram matrix. The two
kink checks independently compute the same rational jump; one imports a pinned
frozen implementation, the other reconstructs P1 segment integrals separately.
The gate check is synthetic. These checks do not repeat the upstream numerical
field certification. Internal AI-assisted review is not external specialist
review or a proof-assistant verification.

[Related primary literature](RELATED_IBVP.md) provides context for compatible
initial and boundary differentiation; it does not establish novelty or prove
this model's estimates. Next, construct a smooth, derivative-consistent
candidate, bound the differentiated bulk forcing and initial-data errors,
then combine them with the shell bound in a full residual ledger. A complete
nonlinear spacetime and a higher-dimensional origin of the Big Bang remain
unestablished.
