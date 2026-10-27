# Conditional derivative bound for the two TT shell-source coefficients

This is a bound for the two affine shell-source coefficients in the existing
23-mode angular decomposition. It does not evaluate the second-order solution,
all boundary rows, or any pointwise Einstein residual. The bounds are
conditional on the registered smooth first field, its first physical shell
rows, and the previously authenticated continuous amplitude enclosures.

## What is bounded

The frozen `HDBLAST_FULL_ANGULAR_RESPONSE_20261026/structure/ANGULAR_STRUCTURE_PROOF.md`
defines the two TT coefficients in a fixed tensor-harmonic basis with
`<T+,T+>=<T-,T->=1/6`. In the source and derivative convention used there,

    G_T+ = G_T- = -(8/3) ell rho_b (a_b J + sigma'_b v/3) d,   (1)

where `J` is the shell amplitude of the first normal metric component
`H=J Y`, `v` is the first scalar amplitude, and `d` is the first angular
gauge amplitude; `a_b`, `sigma'_b` and `rho_b` are fixed shell background
constants. The factor `(8/3)` includes the exact `Y Hess Y` TT projections
and the physical coordinate-forcing adapter. Each displayed bound applies to
an individual registered TT coefficient; it is not the norm of the sum or a
physical force.

The derivative caps `J_t`, `v_t` and `d_t` in the inherited shell report use
ordinary physical coordinate time t, not shell proper time or the scaled x.
For `t=ell x`, the chain rule gives

    G_T,x = -(8/3) ell^2 rho_b [ (a_b J_t + sigma'_b v_t/3) d
                                 + (a_b J + sigma'_b v/3) d_t ].    (2)

The TT linear Robin coefficient is exactly zero in the frozen operator;
thus no homogeneous Robin-feedback term has been omitted from (1). This is
not true for the general set of 18 Robin rows.

## Exact arithmetic and inherited hypotheses

The source `HDBLAST_COMMON_GAUGE_FORCING_20261016/shell-producer/COMMON_HARMONIC_SHELL_NORMAL.json`
is pinned in `check_tt_shell_derivative.py` by SHA-256
`559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a`.
It certifies absolute upper bounds for `J,J_t,v,v_t,d,d_t` and shell constants
on `0<=x<=3,r=0`, conditional on its authenticated actual master/gauge
C3 and background inputs and physical first-shell equations. Let upper bounds
be written with bars and define

    B  = abar Jbar + sigmabar vbar/3,
    Bt = abar Jtbar + sigmabar vtbar/3.

Then the full product rule and triangle inequality give

    |G_T|   <= (8/3) ell rho_max B dbar < 5e-21,
    |G_T,x| <= (8/3) ell^2 rho_max (Bt dbar+B dtbar) < 3.01e-19. (3)

These are uniform absolute bounds; they do not imply the actual coefficient
or its derivative is nonzero. The checker uses exact fractions and includes
the three independent negative controls for the `d_t` product, extra factor
of ell, and sigma'/3. The large exact fractions in its JSON output—not the
rounded decimal—are the primary computed bounds. The report's C3 certification
is inherited, not repeated by this arithmetic check.

## Consequence for the existing nominal candidate

`HDBLAST_FULL_ANGULAR_RESPONSE_20261026/compute/full23_response.py` describes
the saved candidate as having both TT slots identically zero, by zero TT
initial tables and autonomous TT blocks, while retaining actual TT forcing
in its error budget. The exact TT Robin coefficient is zero and the target
boundary law is `tau_r=G_T`. Hence, in the frozen coordinate and sign
convention, the candidate's TT boundary residual is `-G_T`, and its x
variation satisfies `|(R_b)_x|<3.01e-19` for each TT row. This is the
derivative of that boundary defect, not the entire response error derivative.

For the isolated characteristic source lift
`L G=-integral_0^(x-r)_+ G(s) ds`, its nonzero-side second partials are
`(LG)_xx=(LG)_rr=-G_x(x-r)`, `(LG)_xr=+G_x(x-r)`. Under the inherited
smooth/flat zero-corner compatibility that makes the lift C2 at x=r, (3)
bounds all three of those dimensionless-coordinate second partials by
`3.01e-19`. This bound covers only this isolated affine shell contribution;
the bulk quadratic source, second-order response feedback and remaining
boundary rows remain to be bounded. It is not a bound on the total TT
response's second derivatives.

For context, (2) differentiates the on-shell source *tangentially*; it does
not differentiate the first shell equations normally or assign a normal
derivative to the source. Although physical |partial_t G_T| < 3.01e-13
follows from |partial_x G_T|/ell, no shell or gravitational-wave energy or
physical force interpretation is claimed.

## Reproduce

From the repository root, run

    python /path/to/check_tt_shell_derivative.py --repo-root .
    python -O /path/to/check_tt_shell_derivative.py --repo-root .

The outputs should agree byte for byte. The program reads and authenticates
the inherited report, checks its declared scope and domain, uses only Python's
standard library and exact fractions, and tests three deliberately wrong
formulas. It does not recompute the upstream actual-field C3 certificates,
verify existence of a nonlinear family or directly evaluate physical boundary
residuals. Existing full-angular and shell-source capsules remain unchanged.
