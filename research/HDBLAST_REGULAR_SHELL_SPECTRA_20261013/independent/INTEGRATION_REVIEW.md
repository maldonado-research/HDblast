# Independent integration review

13 October 2026 UTC. This review concerns the new proofs supplied in this
package and their dated static/bridge basis. It is an internal independent
check, not external peer review. Symbolic receipts validate equations and
arithmetic; they do not replace the functional-analytic or physical hypotheses.

The scalar Jacobi bridge, full two-copy matching, local even and odd kinetic
calculations, closed operator domains, fixed-background scalar gaps and tensor
partner threshold are mutually consistent under the stated regular-cone and
finite-energy interpretation. No blocking mathematical defect was found in the
reviewed corrected domain or tensor proofs. The unresolved scope is the global
physical state space, admissible global gauge transformations and possible
symplectic flux for other cone extensions.

## What is actually combined

| Sector | Shell realization | Spectral result at fixed positive detuning | Local physical norm |
| --- | --- | --- | --- |
| Generic even scalar | alpha P'=Lambda P/beta; positive endpoint mass | Entire shifted spectrum bounded below by epsilon_even>0 | K=18Lambda[ integral wSL P²+P(T)²/beta ] |
| Generic odd scalar | alpha P'=-P/(3a rho²); ordinary Robin | Entire shifted spectrum bounded below by epsilon_odd>0 | Full EH+GHY+moving-domain+shell calculation gives K=18Lambda integral wSL P² |
| Constant scalar harmonic | Separate static Jacobi constraints | No modulus in the complete regular static family: even det=H_b²M0>0, odd det=ah-pv>0 | Generic N=-2Psi formula is not used |
| Translation-type scalar harmonic | Separate unreduced bulk/gauge classification and scalar matching | No relative embedding under admissible smooth gauge quotient: two-copy scalar determinant=2p²B!=0 | Classified as gauge after matching; generic reduction is not imposed |
| Even transverse-traceless tensor | Neumann h'(T)=0 | One zero radial profile; nonzero spectrum starts at dimensionless 9/4 | Standard positive reduced tensor norm |
| Odd transverse-traceless tensor | Dirichlet h(T)=0 | Spectrum starts at dimensionless 9/4 | Standard positive reduced tensor norm |

Here Lambda=mu²+4 for scalars and alpha=1/(rho²p²). Each generic scalar
essential spectrum starts at mu²=9/4. The scalar fixed-background gap is an
explicit positive functional of that background, established by weighted Hardy
coercivity. It is not inferred merely from strict positivity on eigenfunctions.
The tensor nonzero threshold uses closed partner operators and polar
decomposition, so it concerns full nonzero spectra rather than sampled roots.

The time-independent and translation-type scalar harmonics are exceptional by
their tensor structure. Other harmonics at mu²=0 and mu²=-4 remain in the
generic radial problem. Excluding a whole mass value by its exceptional
representative would be invalid.

## Review corrections and genuine new calculation

The scalar domain proof's endpoint resolvent comparison was clarified to use
a Dirichlet operator with the same actual Liouville potential. The comparison
with the constant 25/4 potential is a separate compact perturbation. With
that clarification, density, trace identification, cone cutoff core, closed
form, endpoint action, limit-point classification, Picone closure and the
Hardy gap proof were accepted under their stated hypotheses.

The initial odd kinetic audit had only the bulk term and an isolated shell
volume term. That limitation has now been removed locally by an independently
checked graph calculation. The GHY graph correction has a lapse cross term;
the moving-domain EH tangential divergence cancels it. Per copy the remaining
graph term is +3a rho²zeta²(DQ)², and the moving intrinsic-curvature term is
+3rho²zeta psi(DQ)². With odd continuity psi=-a zeta they cancel. The common
shell leaves -3a rho²zeta²(DQ)² and hence K_boundary=18c_oddP(T)².
This restores K_odd=18Lambda integral wSLP² at actual odd spectral roots.
Compactly supported tangential test variations avoid external time corners.
The remaining global gauge and cone-flux qualifications still apply.

The dated static proof and bridge are available in [basis/static](../basis/static)
and [basis/bridge](../basis/bridge). The new scalar operator proof is
[SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md](../domain/SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md).
The exact-decimal coupling is retained. The binary64 coupling and stored
coefficients are explicitly recorded as distinct rational mathematical
models; none of these algebra checks validates rounding during historical ODE
trajectories.

## Portable check and limits

From any working directory, run `python independent/check_scalar.py` from
this package's path, or supply `--basis PATH` and `--output PATH`. The basis
defaults to the sibling basis directory, its required inputs have dated
SHA256 pins, and output parent directories are created. Normal and optimized
Python executions have explicit failure exits and equal JSON results. No
absolute authoring-workspace path or clock timestamp is required.

The independent checks reconstruct the Jacobi transformation, scalar boundary
orientation, odd matching, exceptional constraint determinants, local Picone
identity, half-line potential and full local odd graph kinetic coefficient.
Deliberate wrong formulas include reversed shell signs and omitted physical
terms. A separate inequality witness uses g(y)=y^(-7/5)-1 on (0,1]; its exact
weighted ratio 5/12 rejects a corrupted Hardy constant 2/5. This checks an
actual function rather than a tautological comparison of declared constants.

Historical numerical results are not used as rigorous proof in this package.
An isolated archive replay reproduced the old symbolic 55/55 equations and
four wrong controls, but its sampled numerical scans do not prove completeness.
The new algebra controls and reviewed closed-domain arguments are the relevant
evidence for the new result.

Every conclusion concerns the stated regular nonconstant branch at positive
detuning. The divided scalar variable does not include delta=0. A physical
mass is H_b² times its dimensionless squared mass. Since H_b tends to zero,
the essential physical thresholds approach zero: no uniform strictly positive
physical mass floor is claimed. No uniform dimensionless scalar lower bound
has been proved either. Vector perturbations, nonlinear attraction, vacuum
decay, cosmological evolution and all choices of nonnormalizable cone data
remain outside the claim.
