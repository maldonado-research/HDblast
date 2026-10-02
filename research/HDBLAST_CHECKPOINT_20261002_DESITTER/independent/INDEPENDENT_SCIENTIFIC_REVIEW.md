# Independent review of the fixed-reference massive de Sitter sources

Reviewed against public baseline `0205cc651bfb614c32e39dfe833d93d957264229`.
The frozen JUNCTIONS and SMOOTH_FRW checkpoints were read, not changed.
This is a separate internal derivation and implementation, not external peer
review or a physical HDBLAST parameter study.

**Outcome:** the specified S4 action, mass current and radius-variation stress
are consistent with the inherited finite-action convention. Independent
proper-time integration agrees with the producer's spectral-zeta computation
at all four mathematical points. A separate exact-mode calculation obtains
pressure, energy and variance directly from the frozen subtraction orders at
one curved point. No new correction to the conditional frozen static trace
formula is needed for this action. These results do not identify a physical
shell mass law, choose gravitational matching, or solve coupled backreaction.

## 1. Action and variation checks

Write z=H², h=x−r, c=h²/2−2zh+29z²/15. The full heat trace includes l=0,
with eigenvalues z*l(l+3)+x and degeneracies (l+1)(l+2)(2l+3)/6.
Its density begins (4πs)^−2[1+2zs+29z²s²/15+…]. Multiplication by e^−xs
and subtraction of the specified reference polynomial cancel through s².
The action integrand is finite at s=0. Positive x and r ensure the large-s
integrals converge. Neither conclusion is uniform at the excluded x=0
endpoint, where the homogeneous mode obstructs the massive Euclidean state.

For Gamma_E=Vol(S4)*W, varying the radius changes both Vol and z. It gives
rho=W−z*W_z/2, while the mass variation gives Q=2W_x. De Sitter invariance
then gives p=−rho. These definitions fix the signs of the scalar source
J_phi=x_phi*Q/2 and the curvature contribution. A static Ward identity alone
cannot determine those signs.

The reference derivative is exact:

    ∂r{e^−rs[1+(2z−h)s+c*s²]} = −c*s³ e^−rs,
    W_r = −c/(32π²r).

Together with z*W_z+x*W_x+r*W_r=2W this yields

    −4rho = −xQ−2r*W_r = −xQ+c/(16π²).

Thus the conditional A_r expression in JUNCTIONS is correct for this explicit
action. Holding r fixed in physical variations is essential. Treating r=x
during a mass derivative adds 2W_r to Q and changes the theory. The trace
identity follows from action homogeneity; a tiny residual generated from the
same truncated spectral formulas is not independent evidence of small
truncation error. The report must not describe the producer's approximately
10^−83 trace residual as an 83-digit accuracy certification.

## 2. Matching and an independent resolvent derivation

The first three heat coefficients give the finite local pieces

    V=[x²ln(x/r)−3x²/2+2rx−r²/2]/(64π²),
    F=[x ln(x/r)−x+r]/(192π²),
    W_curvature²=29z² ln(x/r)/(480π²).

Consequently V(r)=V_x(r)=V_xx(r)=F(r)=F_x(r)=0. The covariant heat-kernel
basis supplies the inherited R²,C²,E4 logarithms and their zero-at-r finite
convention. S4 alone cannot independently determine all covariant EFT
coefficients: C² vanishes, and constant curvature-squared/topological terms
have vanishing de Sitter metric stress. Matching is a declared convention,
not an observation or a unique determination of physical couplings. The full
curved W and Q need not vanish at x=r.

An independent partial-fraction finite-part calculation uses n=l+3/2,
u=x/z and nu²=9/4−u. Subtract n−(u−2)/n from the resolvent summand and sum the
convergent remainder with digamma functions. Keeping the zeta pole before
taking the finite part gives

    FP Z(1;u)=[(u−2)Psi−u+4/3]/6,
    Psi=psi(3/2+nu)+psi(3/2−nu),
    Q=[(x−2z)(Psi+ln(z/r))−x+r+4z/3]/(16π²).

This independently reproduces the theory agent's expression. The +4z/3
finite term matters. Expressions in other renormalization conventions must
be converted by a common local action before comparison.

## 3. Separate proper-time computation and its error limits

`independent_proper_time.py` imports no producer implementation. It uses the
explicit harmonic heat sum, half-integer Euler–Maclaurin heat coefficients,
and distinct convergent integrands for W, rho and Q. It evaluates the four
registered dimensionless points x/r in {1,2}, z/r in {1/2,1}, r=1 twice.

The review's `PREREGISTRATION.md` was saved locally before any review source
values were computed. Its numerical execution preceded the producer's public
registration commit. This is **local prospective registration only**, not a
claim that the review run was publicly preregistered. The producer's separate
public registration and execution ordering remain its own provenance.

All 12 primary/refinement comparisons pass the predeclared
1e−9+1e−7*abs(value) gate. The largest change is below 8.69e−24. The separate
cross-method audit compares all 12 refined values to the frozen producer and
passes, with maximum absolute difference below 5.99e−31. These are observed
differences, not total-error bounds. Current-sign, metric-sign and heat-a2
coefficient changes are rejected.

The harmonic truncation uses an analytic Gaussian polynomial envelope.
The large-proper-time tail uses analytic exponential-integral/incomplete-
gamma bounds, retaining the zero mode. The small-proper-time asymptotic
remainder has no evaluated rigorous enclosing constant in this implementation.
Precision/join/order refinement tests it empirically. Quadrature and mpmath
arithmetic are not interval-certified. No total numerical error enclosure is
claimed even though the available source accuracy evidence exceeds the gates.

The initial review aggregator recorded trace residuals but omitted them from
its aggregate PASS calculation. The original evidence is preserved.
`audit_independent_sources.py` explicitly enforces all eight trace gates after
the run; they pass. This is a post-run correction restoring enforcement of an
originally registered gate; the correction itself was not preregistered. The same audit verifies the exact producer gate
name set, finite numeric values, booleans, the negative-control set and
recomputable thresholds. It independently rejects five corrupted records
(renamed gate, infinite threshold, nonboolean pass, negative tail and duplicate
point), without changing the frozen primary validator. The saved pairing
residual is checked against the stricter absolute floor because its relative
scale is not separately saved.

## 4. Direct pressure-sensitive bridge to the prescribed subtraction

`exact_mode_bridge.py` is an analytic diagnostic added after the review
registration. At x=r=2H², the minimally coupled canonical mode equation is
u_k''+k²u_k=0 and the Euclidean/Bunch–Davies modes are exact plane waves.
The script derives the reference WKB frequency from its Riccati equation,
forms the stress bilinears, and separately integrates the energy, pressure and
variance differences through stress order four/current order two. The change
of variable v=p/sqrt(p²+r) converts each integral to a rational integral on
[0,1]. The exact answers are

    rho = 11H⁴/(960π²),
    p   = −11H⁴/(960π²),
    Q   = H²/(12π²).

Pressure is calculated directly, without using p=−rho or the trace equation
to obtain it. These answers agree with the source action at the corresponding
registered point x=r=1,z=1/2 after restoring dimensions. This provides an
independent curved-state subtraction check beyond comparing two evaluations
of one determinant. It is one exact point, not a general numerical test of
every Bunch–Davies momentum integral.

## 5. Literature and scope cautions

The locally supplied text of García-Consuegra–Rajantie, arXiv:2511.23076v2,
was inspected around Eqs.33–38. Eq.33 has a positive coefficient multiplying
log[(n+a)²−nu_d²]. Differentiating with respect to nu_d² gives a negative
resolvent, whereas the displayed intermediate sum in Eq.36 has a positive
coefficient. This supports the literature agent's carefully bounded report
of an apparent intermediate sign inconsistency. It does not establish that
the final renormalized potential or the paper's physics is wrong. Our source
formula does not depend on copying that intermediate sign.

The Euclidean determinant defines the regular massive invariant state and
its constant stress/current. Its local subtraction convention agrees with
the frozen covariant action. The exact-mode check independently supports the
source identification at a curved point. None of this establishes a causal
in-in response kernel for time-dependent geometry or state perturbations.
One must not transplant this Euclidean equilibrium action into arbitrary
time-dependent backreaction without that additional construction.

No physical mass law or coupling was selected, no 5D or shell evolution was
performed, and no reheating, thermalization, or hot-universe conclusion
follows from these static vacuum-polarization results.

## Reproduction

From this directory, with Python, mpmath and SymPy available:

    python independent_proper_time.py --output /absolute/path/proper_time.json
    python exact_mode_bridge.py --output /absolute/path/exact_mode.json
    python audit_independent_sources.py --producer /path/to/primary/results.json \
      --proper /absolute/path/proper_time.json --output /absolute/path/audit.json

The proper-time output option accepts either an absolute path or a filename
relative to the script directory. The exact-mode symbolic script defaults to
`EXACT_MODE_BRIDGE.json` beside itself and supports an explicit output path.
All scripts raise
explicit exceptions and retain checks under `python -O`.
