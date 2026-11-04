# Causal identifiability from the normalized shell scalar history

4 November 2026. Conditional result on the retained compact h-only preparation
fiber. This is an inverse stability estimate on its admissible image, not an
existence theorem for arbitrary prescribed histories or a physical state rule.

Fix the same background, first packet, parameter normalization, other second
canonical sectors, harmonic conditioning and physical shell/anchor as in the
preparation comparison. Let the initial conditioned scalar mean difference be
w(r)=-p(r)zeta(r). This is distinct from the evolved five-field scalar channel.
The area profile relation gives

```
D=3h/2-3a w/p.
```

Substitution in the previously verified initial constraint IVP gives

```
z_r=3b z-h,
h_r=3ell^2 z-bh+c w, c=2ell rho p/3,
z(1)=h(1)=0.                                            (1)
```

Thus the initial scalar profile w determines h causally through a two-component
ODE. This statement uses the fixed other canonical sectors; it is not true for
unrestricted second Cauchy data. The complete packet's flatness/moment conditions
remain restrictions on the image, not conditions automatically satisfied by
arbitrary w on the causal segment I=[1,T], T=5/4.

Define `Y(w)=max(sup_I|w|,sup_I|w_r|)`. On the initial interval, the retained
continuous background bounds imply

```
bM=917/62500000, cM=2ell*131*(1/7000)/3=131/10500000000,
Hw=(T-1)*(4/3)*cM=cM/3,
Hrw=(3ell^2+bM)Hw+cM, Vw=3ell^2 Hw+cM.                    (2)
```

The same matrix infinity-norm/exponential argument used for the former initial
IVP applies to (1): its row sums are bounded by20001/20000 and its exponential
on length1/4 by4/3. Consequently `|h|,|z|<=Hw Y(w)`,
`|h_r|<=Hrw Y(w)`. The mixed metric velocity is exactly
`V_L=h_r+bh=3ell^2 z+c w`, hence `|V_L|<=Vw Y(w)`.
The conditioned initial five means are

```
U0=(h,0,0,h,w), U0_x=(0,V_L,0,0,0), U0_r=(h_r,0,0,h_r,w_r).
```

With H=diag(1,2,1,3,1), Gram triangle bounds give

```
D0<=(1+2Hw)Y(w), Dr<=(1+2Hrw)Y(w), Dv<=2Vw Y(w),
N<=max((1+2Hw+2T Vw)/3,1+2Hrw+2Vw)Y(w)
 <Nguard Y(w), Nguard=20000001/20000000                    (3)
```

for Y(w)>0; the non-strict bound also includes Y=0. The precise coefficient in
(3) is `7875000393000961016393/7875000000000000000000`, strictly below Nguard.

Let the evolved homogeneous five-mean difference be W and its reflected free
wave B. The retained full coupled resolvent bound gives
`X(W-B)<=qN/(1-q)`, `q=1201/3200000`, in
`X=max(sup|U|H/3,sup|U_x|H,sup|U_r|H)`. The scalar free wave has exactly
`B_w(x,0)=w(x)`: its initial scalar velocity is zero and its reflection is even.
Both solutions vanish in the exact wedge x+r<=1. Therefore the trace error
`e(x)=W_w(x,0)-w(x)` satisfies

```
|e_x|<=qN/(1-q), |e(x)|<=(x-1)qN/(1-q),
Y(e)<=kappa Y(w),
kappa=q Nguard/(1-q)=24020001201/63975980000000
       <47/125000<1.                                    (4)
```

Every internal metric/scalar and Robin feedback remains in the resolvent.
The derivative of the boundary scalar trace is controlled because X includes
the actual x derivative; no coordinate force is inferred.

The history `m(x)=W_w(x,0)` has a physical difference interpretation. Specifically,
let `M_e(tau)=(integral_{Sigma_e(tau)} phi_e dV_e)/V_e(tau)`, meaning the scalar
integrated over the anchored physical shell section and divided by its proper
volume. For the two exact families sharing the first perturbation, m is the
difference of their second parameter derivatives of this **normalized scalar
mean** at the fixed proper clock `tau=rho_b ell x`. Background phi_b is constant
on the shell. Hence its term phi_b times the second volume variation cancels
between numerator and denominator; the first scalar/measure/clock quadratic
terms are common and cancel in the family difference. A second clock shift
multiplies the zero derivative of phi_b. The fixed shell harmonic representative
includes the required embedding pullback and has zero pure second normal gauge
component at the boundary. Thus the normalized physical mean difference is
exactly m. An unnormalized integral of phi would contain an additional
`phi_b Delta V2` and is not the observable used here.

Using (4) and the triangle inequality gives, including Y=0,

```
(1-kappa)Y(w)<=Y(m)<=(1+kappa)Y(w),
Y(w)<=Y(m)/(1-kappa).                                    (5)
```

The C1 norm of m uses x, while that of w uses r, with both normalizations fixed
by the same ell. They are not SI norms. Applying the estimate to a difference
of two admissible profiles proves causal injectivity: equal normalized physical
scalar histories on I force equal w on I, then equal h,z,zeta there by (1).
This does not determine exterior preparations beyond the causal segment, other
canonical sectors or higher nonlinear coefficients. A single scalar value or
one integrated volume response is not the complete history assumed here.

No onto/completion statement is made for arbitrary C1 histories; the proof uses
the retained regular solutions of the admissible smooth preparation image. A
retarded state still needs independent past data and a compatible development
before its history can be inserted into this estimate. The constants and
analytic/exact-family premises remain conditional, with no effective common
nonlinear parameter radius or cosmological-origin conclusion.


## Retained source links

- [Actual first packet and initial master](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/HIGHER_JETS.md).
- [Physical first Cauchy data and nonzero second constraint witness](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md).
- [Exact harmonic lapse/shift velocities](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md).
- [Canonical action/momentum normalization](../../HDBLAST_SHELL_NOETHER_20261013/manuscript/SHELL_NOETHER.tex).
- [Five mean physical boundary rows](../../HDBLAST_PREPARED_SHELL_VOLUME_20261101/gauge/MEAN_INITIAL_GAUGE.md).
- [Conditional continuous operator premises](../../HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json).
- [Retained coupled preparation comparison and proper scalar meaning](../../HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/DERIVATION.md).
- [Earlier direct H1 energy theorem](../../HDBLAST_H1_PREPARATION_RESPONSE_20261104/energy/PROOF.md).
