# Active-source Ward ledger from independently defined operators

Prepared October 2, 2026 (America/Los_Angeles). The completed source-free
checkpoint and all previous metric FAIL results remain unchanged. This is a
prospective mathematical diagnostic of the fixed interval `[-4.5,-3.5]`, not a
new theory or a verified trajectory. No actual source values, retained mode or
stress arrays, physical time integral, or physical momentum sum were evaluated
in preparing this derivation and its pure proof receipts.

## Direct action and unrestricted amplitudes

At fixed physical mass and subtraction reference `x=r=2`, `H=1`, `xi=0`, write
`a=a0(1+epsilon h)`, `a0=L=-1/eta`, and `L'=L²`. Incoming canonical modes are
`v0=exp(-ik eta)/sqrt(2k)`. With actual unprojected variations
`delta v=v0 u`, `delta v'=v0(w-ik u)`, the exact linear equation is

```
u'=w,  w'=i Omega w-epsilon g,  Omega=2k,
g=4L²h-2Lh'-h''.
```

`h` is the metric perturbation divided by epsilon; `g` has the same convention.
The physical density and pressure operators are varied separately:

```
D=v'-L v, delta D=delta v'-L delta v-epsilon h' v0,
rho=[|D|²+(k²+2a²)|v|²]/(2a⁴),
p=[|D|²-(k²/3+2a²)|v|²]/(2a⁴).
```

For comparisons with retained data, pi, epsilon, time nodes and momentum
weights preserve the producer's exact represented binary80 values. The pi in
measure and finite-band primitives is that same fixed represented constant,
not an independently rounded replacement. MP80/100 accumulation of native-LD
source/phase values does not imply 80/100-digit accuracy of those inputs.

The subtraction is the full fixed-reference W0/W2/W4 stress inventory.
Comoving k and its cutoff remain fixed; physical `a^-4` volume factors vary.
The pressure is never obtained by solving a Ward equation. Put
`X=Re(u)`, `Z=Re(w)`, `T=Im(w)`. The direct mode contributions before the
measure and momentum weights are

```
R_m=[(2k²+3L²)X-kT-LZ]/(2k epsilon),
P_m=[(2k²/3-L²)X-kT-LZ]/(2k epsilon).
```

The bare baseline and explicit metric contacts, including volume variation,
are

```
R0,b=(2k²+3L²)/(4k),  P0,b=(2k²/3-L²)/(4k),
C_R,b=Lh'/(2k)-2kh-2L²h/k,
C_P,b=Lh'/(2k)-2kh/3.
```

These expressions contain the variation of the metric D operator, the fixed
physical mass operator, and the physical volume. For subtraction numerators
`S_R,S_P` and their metric directions `delta S_R,delta S_P`, the renormalized
contacts and baselines are

```
C_R=C_R,b-delta S_R+4h S_R,
C_P=C_P,b-delta S_P+4h S_P,
R0=R0,b-S_R, P0=P0,b-S_P.
```

Every retained subtraction grade independently satisfies its fixed-mass Ward
identity. The new contact verifier rechecks these arbitrary-jet identities
from the preserved formal inventory before using its finite-K primitives.

## Exact forcing and source work

For arbitrary real forcing `f` in `w'=i Omega w-epsilon f`,

```
c_R=X-T/Omega,  c_R'=0,
G=(3L²X-LZ)/(2k epsilon), R_m=k c_R/epsilon+G,
G'=L(R_m-3P_m)+L f/(2k).
```

No Wronskian projection sets `c_R` to zero. Direct metric contacts obey

```
C_R'=L(C_R-3C_P)-3h'(R0+P0)-L g/(2k).
```

Therefore, with the independent direct stress integrand

```
F=L(R_m+C_R-3(P_m+C_P))-3h'(R0+P0),
(G+C_R)'=F+L(f-g)/(2k).
```

The last term is a signed source-definition defect. It is essential if an
independently fitted forcing is combined with metric jets from a different
interpolant. The approved current design uses the actual analytic bump B and
its consistent analytic derivatives, without a source interpolant. The native
long-double evaluations and fixed GL controls still have arithmetic/quadrature
uncertainty; increasing only the accumulator precision cannot remove it.

The full source-energy convention is distinct from physical heating. The
retained canonical free-energy diagnostic is `2kX-T=2k c_R`, which is constant
for arbitrary real forcing at first order. Including the Hamiltonian source
potential gives `(2kX-T)/(2epsilon)+g/(4k)`, whose derivative is `g'/(4k)`.
That explicit external work/contact is not a drift in the physical fixed-mass
Ward equation and does not establish particle yield or heating.

## Phase-aware forced propagation and a genuine integral check

From retained `u_a,w_a` at `a=-4.5`, propagate over the fixed interval using

```
E(t)=exp(i Omega(t-a)),
M0(t)=integral_a^t g(s) ds,
ME(t)=integral_a^t exp(i Omega(t-s)) g(s) ds,
w_pred=E w_a-epsilon ME,
u_pred=u_a+(E-1)w_a/(i Omega)-epsilon(ME-M0)/(i Omega).
```

Use `expm1` for the short-step phase difference. The continuous integral must
be formed from an independently expanded F using these reconstructed modes,
separate contact integration, and the fixed finite-K baseline work. A phase
primitive evaluated from predicted modes is a useful check, but it is not an
independent numerical integral if it merely recomputes the same endpoint
stress twice. The approved primary route integrates the mode primitive and
source work separately using GL24/32; the independent route integrates direct
F with GL16/24. Both retain their lower/higher-order results and fail the fixed
control threshold rather than choosing a new rule after seeing the data.

An optional exact phase construction is proved on fabricated polynomial
sources only: a polynomial h gives `g=polynomial+A/eta+B/eta²`. For
`lambda=-i Omega`, the necessary primitives are

```
J0=exp(lambda eta)/lambda,
Jn=eta^n exp(lambda eta)/lambda-n J_(n-1)/lambda,
J_-1=Ei(lambda eta),
J_-2=lambda Ei(lambda eta)-exp(lambda eta)/eta.
```

At lambda=0 use ordinary polynomial primitives, `log|eta|` and `-1/eta`.
The fixed study interval has eta<0 and does not cross the Ei branch or eta=0.
These polynomial/Ei identities and numerical artificial-source controls do not
select an interpolation method for the actual experiment.

## Exact finite-K contacts versus saved momentum quadrature

Let `mu=k²/(2pi²)`. The original stored trajectory and stresses use the saved
momentum nodes and weights. Their baseline contact is also the weighted sum:
`3h' a0⁴(rho0,K+p0,K)`. It is not the later raw-audit closed baseline primitive.
Define

```
M_d=sum_j weight_j mu_j/(2k_j),
M_A=integral_0^K mu/(2k) dk=K²/(8pi²),
J=integral_a^b L g.
```

The superscript A denotes the independently reduced analytic finite-band
contact. It uses `v=K/sqrt(K²+2L²)`,
`A=asinh(K/(sqrt(2)L))-v`, `A'=-Lv³`, and `q_sub=g A/(8pi²)`.
Writing `q_sub'` and `q_sub''` as actual derivatives, the contact corrections
matched to the unrestricted raw mode bilinears are

```
C_R^A=(3L²q_sub-Lq_sub')/2+L²Q0,K g/2+density_contact+rho_local,
C_P^A=(q_sub''-3Lq_sub'-3L²q_sub)/6-L²Q0,K g/6
      +pressure_contact+p_local-M_A g/3.
```

Here `rho_local=-4h R0,K+Er` and `p_local=-4h P0,K+Ep` are the established
full metric contact primitives, and the density/pressure contacts are the
separately derived scalar mass-subtraction reductions. The pressure term
`-M_A g/3` is necessary: the exact raw mode pressure contains `+M_d g/3`
when expressed through derivatives of its variance. The new symbolic proof
checks

```
(C_R^A)'=F_C^A-M_A Lg,
F_C^A=L(C_R^A-3C_P^A)-3h'(R0,K^A+P0,K^A).
```

Mixing the discrete modes with analytic finite-band contacts does not give an
exact discrete Ward identity. The correct independent ledger and primitive
comparison are

```
I_A=Delta G_d-M_d J+integral F_C^A,
I_A=Delta(G_d+C_R^A)+(M_A-M_d)J.
```

The approved raw primary target uses `-M_A g/3` in `C_P^A`. An alternative
explicit mixed target would use `C_P^mix=C_P^A+(M_A-M_d)g/3`; its direct
integrand is `F_mix=F_A+(M_d-M_A)Lg`, and its primitive is simply
`Delta(G_d+C_R^A)`. That shifted pressure must be derived from the direct
mode-variance relation and declared as a changed target; it cannot silently
replace `C_P^A` while retaining the original momentum-correction labels.
The contact proof includes these mixed-target identities as omission controls.

The primitive correction cannot be silently discarded, even if a sampled
momentum rule is accurate. Neither GL order agreement nor a precision check
proves momentum continuum convergence.

## Signed error decomposition and the smallest decisive test

Keep the original stored long-double Simpson global prefixes and subtract
`S_b-S_a` in long double before exact-ratio conversion. Subtract exact converted
stored density endpoints in MP. For all twelve cases (two sources, two saved
resolutions, K=64/128/256), on the sole fixed interval `[-4.5,-3.5]`, define

```
DeltaR=R_stored,b-R_stored,a,
D_S=DeltaR-S_ab, D_cont=DeltaR-I_A, E_Q=I_A-S_ab,
D_S=D_cont+E_Q.
```

`u_1,w_1` are the initial raw modes; `u_3,w_3` are the final measured modes.
The registered `u_2,w_2` midpoint is an additional fixed diagnostic, not a new
anchor chosen after the result. The endpoint flow defect projection is

```
delta u=u_saved,b-u_pred,b, delta w=w_saved,b-w_pred,b,
E_flow,k=[(2k²+3L_b²)Re(delta u)/k-Im(delta w)
          -L_b Re(delta w)/k]/(2epsilon).
```

Weight and sum the signed projection. Its triangle bound uses
`(2k²+3L_b²)|delta u|/k+(1+L_b/k)|delta w|`, with the same positive weights and
`1/(2epsilon)`. It only bounds the projection of these measured interval flow
defects, not inherited initial-state error or total physical trajectory error.

Reconstruct the stored contacts `C_R^d` separately from the direct action/WKB
inventory or retained component members. Put

```
delta C=C_R^d-C_R^A,
E_momentum=Delta(delta C)+(M_d-M_A)J,
E_operator=Delta[R_stored-(R_mode,saved^d+C_R^d)],
D_cont=E_flow+E_momentum+E_operator.
```

That three-term identity describes exact real forced flow and an exact contact
integral. Actual fixed GL/native-LD calculations require a fourth term rather
than assuming the exact Ward cancellation has occurred in floating point. Let
`J_num` be the computed source-work integral, `Q_C` the computed contact
integral, and `DeltaG_pred` the computed mode primitive. Define

```
E_reconstruction=DeltaR_mode,pred-DeltaG_pred
                 +DeltaC_R^A+M_A*J_num-Q_C,
D_cont=E_flow+E_momentum+E_operator+E_reconstruction.
```

When the computed mode density and G use the same direct arithmetic definition,
`DeltaR_mode,pred-DeltaG_pred` is the weighted change in
`C_inv=sum weight*mu*k*c_R/epsilon`. It is a canonical-invariant rounding/drift
term, not a Wronskian projection. The remaining term is the measured contact
primitive/quadrature residual. Both contribute with their signs. For exact
real propagation and an exact contact integral `E_reconstruction=0`; the actual
experiment reports it and applies the root's fixed empirical 2e-7 consistency
gate. The four-term signed arithmetic closure includes this term and retains
its separate 1e-12 threshold. Neither small `E_reconstruction` nor MP80/100
accumulation supplies an interval certificate for the native-LD modes/source.


Differences between represented old source jets and the consistent new source
jets belong in an explicitly labeled source/operator endpoint residual rather
than being attributed to momentum alone. If `delta C` is obtained by simply
subtracting measured endpoint R from predicted mode R, report it as a combined
contact/operator residual: that arithmetic definition does not independently
validate the subtraction operator. Retaining its sign keeps the bookkeeping
closure useful while preserving this limitation.

Thus there are separate questions: does the original Simpson history differ
from the direct continuous reconstructed F integral; do the stored endpoint
modes agree with reconstructed active forcing; and how much does analytic
finite-K contact reduction differ from the stored momentum quadrature/operator
rounding? No single residual proves all three. GL controls and native long-
double arithmetic are empirical checks, not rigorous integral enclosures.

Preserve every previous FAIL and every original metric gate. Where used in
this new diagnostic, retain the inherited 2e-7 profile/continuous/flow
consistency scale, 1e-12 arithmetic/serialization closure scale, and original
2e-6 attribution trigger. Raw full-profile analytic-contact differences must
not be called a verified momentum-equivalence test when independently
reconstructed discrete contacts are checked only at the three fixed points.
The final prospective protocol supplies its exact registered quantity universe
and acceptance criteria before physical evaluation. A positive result can only identify a ledger contribution on this
registered interval. It cannot explain all 59 historical failures, grant a
metric calibration PASS, establish coupled evolution, replace Einstein's
field equations, or support a higher-dimensional blast as the Big Bang origin.
