# Prospective four-point Euclidean de Sitter source benchmark

Prepared before evaluating any source values, 2 October 2026 UTC. This is a
dimensionless mathematical benchmark, not a choice of HDBLAST masses, couplings,
or gravitational matching. Source commit:
`0205cc651bfb614c32e39dfe833d93d957264229`.

## Fixed experiment

Use the minimally coupled massive scalar in the Euclidean/Bunch–Davies state.
Let y=H², h=x−r, r>0 held fixed under every metric and mass variation. The
benchmark action density is exactly

    W = −1/2 ∫₀∞ ds/s { e^(−xs) K_H(s)/Vol
              − e^(−rs)/(4πs)² [1+(2y−h)s+c s²] },
    c=h²/2−2hy+29y²/15,
    K_H(s)=Σₗ₌₀∞ (l+1)(l+2)(2l+3)e^(−y l(l+3)s)/6,
    Vol=8π²/(3y²).

This is the finite common-action convention in the inherited smooth-FRW
matching, extended to the Euclidean state by this explicitly convergent action.
Its low-proper-time subtraction through order s² and positive x,r make the
integral convergent. Local reference subtraction is inside W exactly once.

The grid is the four combinations x/r∈{1,2}, y/r∈{1/2,1}, with r=1 as the unit
choice. No other points, parameter fitting, response tuning or bulk solution
are authorized by this registration. The domain is chosen prospectively so
q=x/y−9/4 obeys |q|/(3/2)²≤7/9, allowing an absolutely convergent zeta expansion
with an explicit truncation bound, while containing x=r and x≠r and both signs
of q. All source values include the l=0 massive mode.

Compute W/r², rho/r²=[W−(y/2)W_y]/r², p=−rho and Q/r=2W_x/r. Do not use the
trace identity to compute rho. J_phi=x_phi Q/2 remains symbolic because no
physical mass law has been supplied.

## Analytic continuation and numerical method

Set a=3/2 and define

    B0=ζH(−3,a)−ζH(−1,a)/4+q/8+q²/4,
    B1=2ζH'(−3,a)−ζH'(−1,a)/2
       −q[ζH(−1,a)+ψ(a)/4]
       +q²[1/4−ψ(a)/2−ζH(3,a)/8]
       +Σₖ₌₃∞ (−q)ᵏ[ζH(2k−3,a)−ζH(2k−1,a)/4]/k.

The spectral zeta derivative is ζ_D'(0)=(B1−B0 ln y)/3, and the finite part of
the convergent proper-time difference gives

    A=x²/2−2xy+29y²/15,
    W=−y²(B1−B0 ln y)/(16π²)
      +(rx−r²/4−2ry−A ln r)/(32π²).

Metric and mass derivatives are obtained by analytically differentiating this
same convergent series. The special k=1,2 pole contributions displayed above
are retained. No spectral cutoff or regulator limit is fit from the answers.

At truncation k≤N, t=|q|/a²<1, define C=a³[1+a/(2N−2)]. Then rigorous absolute
series truncation bounds are

    |δB1| ≤ C t^(N+1)/[(N+1)(1−t)],
    |δB1_q| ≤ C t^(N+1)/[|q|(1−t)].

They follow from ζH(2k−3,a)≤a^(3−2k)[1+a/(2k−4)] and dropping the positive
subtracted ζH term. Propagate them through W, W_x, W_y and rho with triangle
inequalities. At q=0 the tail and its first derivative vanish. These bounds
control analytic series truncation; mpmath arithmetic and transcendental
evaluation are not interval-certified.

Run 50 decimal digits with N=128, then 80 digits with N=256. For independent
direct variation of the action, evaluate centered W differences at h=10⁻⁴ and
h/2 in each dimensionless variable and use fourth-order Richardson
extrapolation. The halved-step difference is retained as a diagnostic, without
calling it a rigorous derivative enclosure. For the common-action cross
derivative use centered differences of independently evaluated rho and Q.
An independent coincident Green-function check uses the digamma expression

    u=x/y, v=r/y, ν=sqrt(9/4−u), Ψ=ψ(3/2+ν)+ψ(3/2−ν),
    Q=y[(u−2)(Ψ−ln v)−u+v+4/3]/(16π²).

The independent reviewer separately registers a direct proper-time quadrature
with heat-trace UV expansion and explicit harmonic/IR tail bounds. Its UV
asymptotic and quadrature refinement are empirical convergence evidence, not a
claimed total-error enclosure. It imports no primary producer implementation.

## Frozen gates and failure policy

All thresholds apply to dimensionless quantities. Thresholds are fixed before
source values, using the unit reference scale and the expected one-loop factor,
not rescaled from observed discrepancies.

1. Analytic series tail bounds for W, rho and Q must each be <10⁻¹⁴ at both
   resolutions. Numerical 50/N128→80/N256 change must be ≤10⁻¹¹+10⁻⁹|refined|.
2. Independently differenced mass/metric derivatives must agree with analytic
   Q and rho within 10⁻¹⁰+10⁻⁸|analytic|. The source-pairing check
   rho_x=Q/2−y Q_y/4 uses the same absolute/relative rule.
3. The independent digamma Q must agree within 10⁻¹¹+10⁻⁹|Q|.
4. Independently computed rho,Q must satisfy −4rho=−xQ+c/(16π²) within
   10⁻¹¹+10⁻⁹ max(|4rho|,|xQ|,|c/(16π²)|). This checks a consequence of W_r;
   rho is never manufactured from it.
5. Direct proper-time cross-method comparison, separately implemented, targets
   ≤10⁻⁹+10⁻⁷|value| for W,rho,Q, contingent on its own refinement gates.
6. Negative controls must reject rho_wrong=W (omit metric variation),
   Q_wrong=−Q, and omission of the current paired to a local F=x curvature
   term. A local V=x²/2,F=x perturbation must obey the action-pairing identity.

Failing gates are retained as FAIL and stop acceptance. No tolerance inflation,
unregistered parameter changes or deletion of failures is permitted. A new
method after failure requires a separately dated prospective repair record.
The executable freezes hashes of this registration, producer, validator and
the inherited source documents before the first run; every run refuses a
hash mismatch and writes to a new output directory. Numerical execution
provenance is distinct from the source commit pin.

This benchmark does not solve the static junctions or coupled five-dimensional
problem, select physical matching, establish dynamic attraction, calculate
particle production or heating, or extrapolate beyond the four declared points.
