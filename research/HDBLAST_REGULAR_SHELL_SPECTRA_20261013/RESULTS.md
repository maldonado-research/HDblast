# Claim ledger — 13 October 2026

The exact registered model is

\[
 W(\eta)=\frac13+\eta^2+\frac13\eta^3,\qquad
 U=\frac12W_\eta^2-\frac23W^2,\qquad
 \sigma=2W+\delta(1+c+c\eta),\quad
 c=\frac{96714321425692}{161839258930731}.
\]

The unchanged static basis certifies a regular nonconstant branch for
`0 < delta <= 0.002`, with `|eta| <= 0.01`, positive warp expansion and negative
scalar derivative away from the cone. The divided scalar variable does not
cover `delta = 0`. Binary floating-point coefficient variants are distinguished
from this exact model.

The field shift `eta = phi - 1` takes the earlier potential
`W(phi)=1-phi+phi^3/3` to the displayed expression. The regular branch near
`phi=+1` and the original unstable hilltop branch are distinct solutions of
the same shell theory. No evolution connecting them is proved. The earlier
scalar-independent-detuning comparator is a separate boundary theory.

Let `alpha = 1/(rho^2 p^2)`, `q = 2/(3 rho^2)`,
`w = 1/(rho^4 p^2)`, `Lambda = mu^2 + 4`, and
`beta = rho_b^4 p_b^2 B > 0`. The scalar equation is
`-(alpha P')' + q P = Lambda w P`.

| Result | Assumptions and evidence | Boundary of the claim |
| --- | --- | --- |
| Even scalar realization | Closed energy form in `L2(w dy) + C`, with trace mass `1/beta`; endpoint equation `alpha P'=Lambda P/beta` | The finite-energy regular-cone phase space is specified, not derived for every global extension |
| Odd scalar realization | Full two-copy matching gives `alpha P'=-c_o P`, `c_o=1/(3a_b rho_b^2)>0`; bulk Hilbert space with positive Robin form | Different reflection parity is included; arbitrary junction actions are not |
| Entire scalar spectrum positive after subtracting 4 | Positive Jacobi solution, strict endpoint mismatch, closed Picone identity and weighted Hardy inequality | More than an eigenvalue scan; the functional gap is fixed-background, not a computed uniform detuning bound |
| Scalar essential threshold | Liouville potential tends to `25/4`, hence shifted threshold `9/4` | No blanket assertion about time evolution or all spectral types above threshold |
| Generic local scalar kinetic signs | Even `K=18 Lambda (integral w P^2+P_b^2/beta)`; odd `K=18 Lambda integral w P^2`; complete local boundary action independently checked | Global gauge, cone flux and physical-state completeness remain conditional |
| Exceptional scalar subspaces | Static-family determinants and unreduced translation-type constraints | Exceptions are harmonic subspaces, not whole mass values |
| Tensor gap and zero profile | Closed partner operators and exact `U<0` on the branch tube; even zero mode, all nonzero spectrum at least `9/4`, odd spectrum at least `9/4` | Standard techniques applied to this family; not unrestricted gravitational stability |

For the even scalar gap put `a_g=alpha P0^2`, `m_g=w P0^2`,
`b=P0(T)^2/beta`, and
`k=[alpha(T)P0'(T)/P0(T)-4/beta]P0(T)^2 > 0`.
With `alpha_*=inf(a_g/y^4)>0`, `m_*=sup(m_g/y^2)<infinity`, and
`M=integral m_g`, the proved full-spectrum bound is

\[
 \epsilon_\delta=
 \left[\frac{4m_*}{9\alpha_*}+\frac{M+b}{k}\right]^{-1}>0.
\]

For the odd problem take `b=0` and
`k=[alpha(T)P0'(T)/P0(T)+c_o]P0(T)^2>0`.
The proof establishes these background functionals' finiteness and positivity;
it does not assign them a numerically certified value.

The exact tube bound used for the tensor theorem is

\[
 U(\eta)\le-\frac{997874976899}{13500000000000}<0,
 \qquad |\eta|\le\frac1{100}.
\]

The standard partner potential is
`V_-=9/4+rho^2(9p^2/16-U/8)`. Both parities therefore have the claimed
nonzero tensor threshold without relying on a sampled eigenvalue search.

Earlier soft-boundary and brane-stability literature is addressed explicitly
in the literature report. The compact Minkowski single-brane instability
theorem's endpoint hypotheses differ from the regular de Sitter cone used
here; that distinction is an assumption comparison, not a general refutation
of the earlier theorem.

No result in this capsule establishes a higher-dimensional Big Bang origin,
observational preference, a new cosmological history, APS acceptance or a
Nobel-level discovery. The scientific advance is a more complete conditional
spectral and local-action calculation with reproducible evidence.
