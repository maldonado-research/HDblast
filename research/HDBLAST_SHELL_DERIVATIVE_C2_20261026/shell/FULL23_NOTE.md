# Conditional x-derivative bound for the complete 23-mode affine shell source

This note bounds the actual registered 23-mode affine source vector `G_x`,
not an independently evaluated boundary residual or a full second-order
solution. The estimate is deliberately conservative and conditional on the
inherited first-field and background enclosures. The source and basis are the
frozen second-parameter-derivative conventions of
`HDBLAST_FULL_ANGULAR_RESPONSE_20261026/shell/FULL23_SHELL_AND_SUBSIDIARY.md`.

## Exact eight-profile structure

The frozen audit gives `G=B c` with constant rational 23-by-8 matrix B, a
positive comparison norm `||G||_H`, and `S=B^T H B` recorded exactly in
`FULL23_SHELL_NORMAL.json`. The matrix B and S are constant in shell time x;
their entries include negative signs. The eight source profiles, with h=J,
L=Z and V=V_gauge in the retained shell producer's notation, are

    c0 = -ell*rho*g0,       c1 = -ell*rho*g1, c2 = -ell*rho*g2,
    c3 = 2*ell*rho*(-alpha+2*beta*X), c4 = ell*rho*pi,
    c5 = -2*ell*rho*beta*V, c6 = 2*ell*rho*(alpha-2*beta*L),
    c7 = -4*ell*rho*beta*d,

where

    beta = a*h + sigma1*v/3,
    alpha = a*h^2/4 - sigma1*h*v/6 - sigma2*v^2/6,
    pi = p*h^2/4 + sigma2*h*v/2 + sigma3*v^2/2,
    tau = 2*X + 6*L - 16*d,
    g0 = 4*sigma2*v^2/3 + 4*sigma1*h*v - 10*a*h^2
         + 8*a*h*tau - 4*a*(4*X^2+12*L^2-64*L*d),
    g1 = 8*a*V^2,           g2 = -16*a*d^2.

The shell constants `rho,a,p,sigma1,sigma2,sigma3` have no t derivative in
this registered fixed-background problem. The report bounds the listed field
amplitudes and their **ordinary physical-coordinate t derivatives**; it does
not identify these derivatives with proper shell time. Because `t=ell*x`,
every c_x has an additional factor ell relative to its t derivative. In
particular, the c0 derivative must retain the derivatives of both h*tau and
L*d, and c3, c5, c6, c7 must retain both product-rule terms.

## Uniform exact bound

Let bars denote the positive upper bounds recorded in the inherited shell
report for the corresponding magnitude, and write derivatives with `t`
subscripts. Define

    beta_bar = a_bar*h_bar + sigma1_bar*v_bar/3,
    beta_t_bar = a_bar*ht_bar + sigma1_bar*vt_bar/3,
    alpha_t_bar = a_bar*h_bar*ht_bar/2
        + sigma1_bar*(ht_bar*v_bar+h_bar*vt_bar)/6
        + sigma2_bar*v_bar*vt_bar/3,
    pi_t_bar = p_bar*h_bar*ht_bar/2
        + sigma2_bar*(ht_bar*v_bar+h_bar*vt_bar)/2
        + sigma3_bar*v_bar*vt_bar,
    g0_t_bar = (8/3)*sigma2_bar*v_bar*vt_bar
        + 4*sigma1_bar*(ht_bar*v_bar+h_bar*vt_bar)
        + 20*a_bar*h_bar*ht_bar
        + 8*a_bar*(ht_bar*tau_bar+h_bar*taut_bar)
        + 4*a_bar*(8*X_bar*Xt_bar+24*L_bar*Lt_bar
                    +64*(Lt_bar*d_bar+L_bar*dt_bar)),
    g1_t_bar = 16*a_bar*V_bar*Vt_bar,
    g2_t_bar = 32*a_bar*d_bar*dt_bar.

The vector of absolute profile-derivative bounds is

    C = ell^2*rho_max * (
         g0_t_bar, g1_t_bar, g2_t_bar,
         2*(alpha_t_bar+2*beta_t_bar*X_bar+2*beta_bar*Xt_bar),
         pi_t_bar, 2*(beta_t_bar*V_bar+beta_bar*Vt_bar),
         2*(alpha_t_bar+2*beta_t_bar*L_bar+2*beta_bar*Lt_bar),
         4*(beta_t_bar*d_bar+beta_bar*dt_bar) ).

For every `0<=x<=3` at the shell, the triangle inequality gives

    |(G_x)_i| <= sum_j |B_ij| C_j,
    ||G_x||_H^2 = c_x^T S c_x <= sum_ij |S_ij| C_i C_j.

The exact rational calculations in `check_full23_shell_derivative.py` yield

    ||G_x||_H < 58,807,538,
    max_i |(G_x)_i| < 182,093,220,
    |(G_T+)_x|=|(G_T-)_x| < 3.01e-19.                    (3)

All five homogeneous Dirichlet rows of G vanish exactly. The apparent
contrast between the very small TT bound and the large full-norm bound reflects
the individual global upper bounds and the particular source factorization; neither inequality asserts a nonzero derivative. The norms and
23 amplitudes are the fixed formal decomposition's conventions, not a physical
radiation or energy scale. The exact squared norm bound and maximal modal-row
bound are printed as fractions by the checker; it takes an integer upper bound
on their positive square root without floating point.

The input reports are authenticated by SHA-256 before use:

* `COMMON_HARMONIC_SHELL_NORMAL.json`: `559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a`.
* `FULL23_SHELL_NORMAL.json`: `0b4a7fd639deec8bfbe65f7b2427e6b084f18195ef954d0add39471528b6f6fe`.

These hashes authenticate the recorded inputs; they do not independently
re-prove either source identities or first-field PDE certificates. Inherited
background, exact-shell, master/gauge C3 and first-shell equations remain
conditional assumptions.

## What the bound enables and does not establish

For the saved candidate, the affine source was omitted: it is compared with
zero and its omission was already included in the existing C1 error ledger.
Thus (3) is also a bound on `partial_x(G-0)` for that particular **source-only**
comparison. It is not the derivative of the entire finite-candidate Robin
residual: homogeneous feedback `M PN V` and other terms must be included.
For TT only, `M=0` and the saved candidate has TT coefficients identically
zero, giving the stronger per-row boundary-defect bound explained in NOTE.md.

For the isolated characteristic lift `L G`, on the nonzero side of `x=r`,
its `xx, xr, rr` derivatives are `-G_x,+G_x,-G_x`. Subject to the inherited
smooth zero-corner compatibility, (3) bounds this **one isolated term** in
the H norm and also applies to each modal row. It does not bound differentiated
bulk forcing, the coupled solution/feedback, the full Robin residual, or the
pointwise Einstein equations. Obtaining those requires further analysis and
quantitative error bounds. No new field evolution, physical detection, exact
nonlinear spacetime, or higher-dimensional Big Bang claim follows.

## Reproduction

    python /path/to/check_full23_shell_derivative.py --repo-root .
    python -O /path/to/check_full23_shell_derivative.py --repo-root .

Outputs must match byte for byte. The standard-library checker authenticates
the two reports, checks their scope and 23-by-8 shapes, explicitly uses the
signed B and Gram matrices, reconstructs all 64 entries of B^T H B from
the frozen section-space weights, independently differentiates all eight profile
expressions with exact rational dual numbers at a nondegenerate fixture, and
checks intentionally dropped `d_t`, `V_t`, and `L_t` terms and an incorrect
TT Gram weight. These algebraic
controls verify arithmetic and normalization for the displayed inputs; they
are not an independent full-field evolution or outside expert review.
