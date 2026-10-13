# Final manuscript mathematical-review addendum

Date: 13 October 2026. Separate internal mathematical review; not external peer review. This addendum leaves the earlier `MATHEMATICAL_REVIEW.md` and its receipts unchanged.

## Reviewed version and decision

Reviewed file: `manuscript-authoring/REGULAR_CONE_EVOLUTION.md`.

SHA-256: `645b5889718e63cde8688b8314e52e435c17c9b64a806370a41edc7ba13cad35`.

Exact bytes: 22,813. Snapshot: `final_reviewed_REGULAR_CONE_EVOLUTION.md`.

**No genuine mathematical error or required correction was found in this version.** The earlier parity and positive-`s` support clarifications are present. Appendix A's boundary laws, physical Cauchy form, coefficients, positivity argument and kernel test are correct within the stated local smooth compact-data scope. The new abstract sentence is justified by that argument. The earlier assessment that the manuscript supports a rigorous conditional local theorem remains unchanged.

“Sufficiently small” in the new abstract sentence is understood as support in a sufficiently small ball, exactly as specified in Appendix A, rather than a nonlinear small-amplitude hypothesis. The result concerns the image of this master reconstruction and does not claim nondegeneracy or completeness of an unrestricted gravitational phase space.

## 1. Shell-law check

For `P=r⁴F`,

\[
P_y/r^4=nF+4F/\rho.
\]

The inherited even condition is `alpha P_y=(Box_gamma+4)P/beta`, with `beta=rho⁴p_bg² B` and `alpha=1/(rho²p_bg²)`. Since `Box_h=Box_gamma/rho²`, it becomes precisely

\[
nF-B^{-1}\Box_hF+(4/\rho-4/(\rho^2B))F=0.
\]

The inherited odd condition is `alpha P_y=−c_oP`, with `c_o=1/(3a rho²)`, and therefore becomes

\[
nF+(4/\rho+p_{\rm bg}^2/(3a))F=0.
\]

The outward-normal signs and all factors in (20) are correct. These identities do not constitute a boundary well-posedness theorem, as the manuscript explicitly states.

## 2. Actual Cauchy-current measure and the powers of r

At `T=t=0`, `R=r`, `dy=Omega dr`, and `partial_t=r partial_T`. Transforming the inherited physical time current with future Cauchy orientation gives the positive pairing convention with measure

\[
\rho^2\,dy\,\partial_t=\Omega^3r^3\,dr\,\partial_T.
\]

This is the covariant Einstein–scalar current, not the auxiliary canonical master's current. On the slice,

\[
\Psi={r^2\over\Omega^2}q,\quad
\Psi_T={r^2\over\Omega^2}p,
\]

\[
\chi=-{3\over2\Omega^2f}(D_r+4)q,\quad
\chi_T=-{3\over2\Omega^2f}(D_r+5)p.
\]

The prefactor's `T` derivative vanishes at `T=0`, since it is a function of `s=R²−T²`. The extra unit in `D_r+5` comes from `partial_T E F=E partial_TF+partial_TF`, not from a field equation. Thus the two coefficients in (21) follow independently:

\[
\Omega^3r^3\left({3\over2\Omega^2f}\right)^2
={9r^3\over4\Omega f^2}=m_r,
\]

\[
6\Omega^3r^3\left({r^2\over\Omega^2}\right)^2
={6r^7\over\Omega}.
\]

In particular `r³`, `r⁷`, and `D_r+5` are all correct. Antisymmetrizing the two data pairs gives exactly `K(q1,p2)−K(q2,p1)`. Reversing the global Cauchy orientation would reverse the complete symplectic form, without changing its nondegeneracy; no relative sign discrepancy was found here.

## 3. Integration by parts and the origin

Expand the scalar contribution on the diagonal:

\[
(D_rq+4q)(D_rq+5q)=(D_rq)^2+9qD_rq+20q^2.
\]

Because `D_r log m_r=3−2sL`, the exact density identity is

\[
m_r(D_rq+4q)(D_rq+5q)
=m_r[(D_rq)^2+(2+9sL)q^2]
+\partial_r\left({9\over2}r m_r q^2\right).
\]

The coefficient is therefore `20−(9/2)(4−2sL)=2+9sL`, confirming (22). At the outer boundary compact support removes the total derivative. At the origin, smooth Cartesian data are bounded uniformly in the angular variables and `r m_r=O(r⁴)`, so the boundary term also vanishes. The polar-coordinate origin causes no missing corner contribution.

`Omega>0` and `f!=0` imply `m_r>0` for `r>0`. Since `L` is bounded near zero, a ball can be chosen with `2+9sL≥1`. Then the diagonal integrand contains a strictly positive `m_r q²` term, as well as nonnegative derivative and `6r⁷q²/Omega` terms. A nonzero smooth function cannot be supported only at the single origin. Consequently `K(q,q)>0` for every nonzero smooth real compact `q` in that ball.

## 4. Kernel test and the abstract claim

Let `V` be all smooth real master-data pairs compactly supported in the chosen ball. If `(q,p)` lies in the kernel of `sigma` against every pair in `V`, then `(0,q)` is an admissible test and gives `K(q,q)=0`; hence `q=0`. Testing next against `(p,0)` gives `−K(p,p)=0`; hence `p=0`.

This argument does not assume that `K` is symmetric. It establishes weak algebraic nondegeneracy on the stated smooth compact-data space and therefore on its reconstructed solution image. There is no unproved completion, topology equivalence, Hamiltonian positivity or full quotient claim. The shell contribution vanishes because the data and their reconstructed derivatives vanish near the shell. For common reflection parity the doubled form is twice the one-copy form, which preserves this conclusion.

The positivity of the auxiliary diagonal `K(q,q)` used in the kernel test is not a claim that the antisymmetric symplectic form is positive. It is also not a proof that the physical Hamiltonian is bounded below. The manuscript keeps these distinctions explicit.

Accordingly, the new abstract sentence—nondegeneracy on a sufficiently small smooth compact master-data image—is supported. It should retain that restriction in any later summary. No assertion about all Einstein–scalar perturbations, all boundary gauges, or the global gauge quotient follows.

## 5. Executed independent checks and unchanged scope

The new independent checker `check_final_appendix_independently.py` has SHA-256 `57df92c0a317a7dfb96b764dfca7be45d1b9588f3996caab22215c2244cec000`. It checks the shell conversions, physical time measure, both radial weights, the arbitrary-function Euler commutator, `D_r+5`, the logarithmic weight derivative, and the exact integration-by-parts density identity.

Normal and optimized Python each pass nine exact identities and four wrong-formula controls. Their JSON outputs are byte-identical. An injected failed identity exits one in both modes. `FINAL_APPENDIX_CHECKS.json` records the actual commands, source hashes and output hashes. Smooth endpoint decay, positivity and the kernel argument are analytic reasoning, not claims proved by these symbolic tests.

The prior full-residual, local hyperbolic evolution, support, flux and relational-scalar review remains in force. Nothing in Appendix A upgrades the finite-time result to evolution after shell contact, establishes a completed global physical phase space, or certifies novelty. Within those unchanged assumptions and limits, this version is mathematically suitable to proceed to the author's final compilation and review.
