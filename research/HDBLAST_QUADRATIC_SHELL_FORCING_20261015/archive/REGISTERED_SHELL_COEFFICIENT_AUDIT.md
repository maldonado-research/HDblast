# Registered shell tension and quadratic coefficient audit

This note audits the actual near-phi=+1 regular branch at delta=1/1000 and
c=96714321425692/161839258930731. It uses the authenticated prior inertial
waveform and the registered full-contact junction conventions. It introduces no
alternate tension or background and executes no historical program.

## Authenticated public inputs

The inertial manifest is pinned externally to
`9b42295d819f0ed5caf94d0a46af38274a27fb11648dd6bff91f5d992385f140`.
All 75 source records, its manifest and all 76 ZIP members were checked against
their recorded sizes and SHA256 hashes. The ZIP itself matches its external
receipt. The full-contact manifest is pinned by that checkpoint's external
frozen receipt to
`0ee47e71ac4371a0f118e2c9a6bff835181a8300e2805423a444ffe700eac843`;
the same checks cover its 80 source records and 81 ZIP members.

Nine selected public files retain their original filenames in the flat
`inputs/prior/` and `inputs/reference/` directories. `inputs/INPUT_PINS.json`
records their bytes, hashes, original repository paths, manifest identities and
immutable published commit URLs. `SOURCE_AUTHENTICATION.json` records the full
source and ZIP authentication. These checks authenticate bytes; the retained
analytic background theorem remains a mathematical input. No private research
material or raw paper was copied.

## Exact tension and correlated shell bounds

Set eta=phi-1. The actual registered tension is

\[
\sigma(\eta)=\frac23+2\eta^2+\frac23\eta^3
 +\delta(1+c+c\eta).
\]

Consequently

\[
\sigma'=4\eta+2\eta^2+\delta c=-2p_b,
\qquad \sigma''=4+4\eta=4\phi_b,\qquad \sigma'''=4.
\]

Primes differentiate with respect to phi, equivalently eta. On the retained
eta interval, sigma'' is positive. Thus sigma' is strictly increasing and its
exact endpoint evaluations give a sharper correlated enclosure than separately
combining the eta and eta-squared terms. The full rational endpoints are in
`REGISTERED_SHELL_COEFFICIENT_AUDIT_NORMAL.json`; simple strict outward bounds
are

\[
\frac{248679003}{10^{12}}<\sigma'_b<\frac{274244855}{10^{12}},
\qquad
\frac{3999651068}{10^9}<\sigma''_b<\frac{3999676637}{10^9}.
\]

The original input box is preserved unchanged. Its sigma' interval equals the
negative twice of its p interval, exactly. The sharper interval also gives a
sharper p interval by p=-sigma'/2, without any new trajectory computation.

## Material source versus a coordinate source

The convention is the derivative coefficient

\[
g=g_0+\varepsilon h_1+\tfrac12\varepsilon^2h_2,
\qquad \phi=\phi_0+\varepsilon\chi_1+\tfrac12\varepsilon^2\chi_2.
\]

With one common tension and identical copies, the single-copy material scalar
junction is

\[
p_2+\frac{\sigma''}2\psi_2
=-\frac{\sigma'''}2\psi_1^2=-2\psi_1^2.
\]

The unsplit two-copy equation instead has source -4 psi_1 squared. The
physical first material trace is psi_1=S(t)Y in the retained straight-shell
observable construction. This does not identify it with the coordinate scalar
trace in every wave gauge. The material term -2 psi_1 squared is one exact
inhomogeneous Taylor coefficient, not the full coordinate scalar forcing, a
second-order solution, an energy, or an observable second-order waveform.

For the retained longitudinal first representative, the normalized unit normal
gives

\[
D^2K=(4uu_y+20au^2)Y^2h_0,\qquad
D^2p=(4uv_y+12pu^2)Y^2.
\]

Using the first shell rows

\[
u_y=-2au-pv/3,\qquad v_y=-2pu-\sigma''v/2,
\]

the complete Eulerian right sides are

\[
S_K=\frac{\sigma''}6v^2-12au^2,\qquad
S_p=-4pu^2+2\sigma''uv-2v^2.
\]

Replacing S_p by only -2v squared would discard normal and gauge-dependent
terms. The retained degree-four harmonic Q satisfies
Y squared=Q+F_2/2+1/6, so the projection coefficient of Y squared onto this
specific unnormalized Q is one. Coefficient signs have no independent physical
energy or cosmological interpretation.

## Fresh independent checks

`check_registered_shell_coefficients.py` reads the nine pinned copies and uses
exact rational arithmetic and symbolic differentiation. It checks the cubic
tension, monotone endpoint bounds, p correlation, normal normalization, Taylor
factorials, and both Eulerian sources. The eight embedded wrong-formula controls
cover the alternate quadratic tension, third-derivative signs/factors, material
sign/factor, omitted normal terms and Eulerian coefficient errors.

Normal and optimized Python outputs are retained separately; their mathematical
content must agree. Four separate injected mathematical errors are rejected in
each mode. The modes, process outcomes and log hashes are recorded in
`REGISTERED_SHELL_COEFFICIENT_REPLAY.json`. These finite checks supplement the
displayed derivation; they do not formally verify the inherited analytic PDE
hypotheses or establish a nonlinear family, a remainder bound or external
novelty.
