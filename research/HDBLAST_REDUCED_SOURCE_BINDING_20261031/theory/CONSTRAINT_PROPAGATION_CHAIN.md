# Exact reduced-source joins and their conditional composition

31 October 2026. The analytic Bianchi/subsidiary construction already exists.
This checkpoint checks its current registered source roles and covector
operator explicitly, and supplies an additive actual third-jet interface.
It does not construct a new evolution or finite-amplitude solution.

## Source conventions

Use the affine path `g(e)=g0+e q`, `phi(e)=phi0+e eta` and second parameter
derivatives, so a Hessian Q occurs as `e^2 Q/2`. The natural residuals are
`E=Ric-dphi*dphi-(2/3)Vg` and `S=Box(phi)-V'`.

Let P be the linear connection variation with its first index lowered by g0,
and let `C1_b=g0^{ac}P_bac`. For the relative harmonic gauge covector lowered
by the full metric, its affine Hessian obeys

```
QC_actual_b=-2q^{ac}P_bac,
QC_actual_b-QC_background_b=2q_bd C1^d,
PC1_ab=2P^d_ab C1_d,
gauge_gradient_ab=(nabla_a QC_actual_b+nabla_b QC_actual_a)/2,
QEred_general=QE-gauge_gradient+PC1,
QEred=QE-gauge_gradient.
```

`QEred` is the registered harmonic specialization. At generic first jets
`QEred_general` carries the PC1 correction. Natural QE, these reduced
Hessians and principal forcing `-2ell^2 rho^2 QEred` are distinct inputs.
Raising contracted indices uses the Lorentzian background signature.

Independent source/frame constructions check the full angular lowering,
PC1 and gauge-gradient identities, all 20 symmetric-tensor modal role joins,
and the scalar correction

```
Qphi_general-Qphi=2C1^b(q_b{}^a(dphi0)_a-(deta)_b).
```

The first field is the registered seven-profile metric quadrupole and scalar
profile. The source-role checks use its free coherent ordinary profile jets;
they are not a completeness claim for unrestricted metric perturbations.
The role identities use the stated warped/de Sitter kinematics. Reducing the
curvature and Bianchi identities below also uses the actual background
Einstein–scalar equations.

## Independent covector wave

For `B0(T)_b=nabla^a T_ab-(1/2)nabla_b tr(T)`, direct contractions give

```
B0(sym nabla Z)=(1/2)L0 Z,   L0 Z=Box0 Z+Ric0.Z,
Ric0.Z=(2/3)V Z+p^2 Z_y e^y.
```

The symmetric gradient is half-normalized. Two independent arithmetic/frame
implementations verify the operator on an arbitrary registered eight-mode
lower covector and 48 free ordinary jets, with full angular reconstruction.
The rough wave is contracted directly, not defined by rearranging B0.
Its connection acts on every covariant slot, including the first derivative
slot. The radial p^2 curvature term is retained; the background is not vacuum.

The eight modal coefficients consist of three radial, three temporal and two
spherical-gradient rows. The spacetime is five-dimensional; four redundant
ambient sphere slots represent three tangent dimensions. The background
constraint and differential relations, including the de Sitter quotient,
are explicit in the coefficient reduction.

## General reduced Bianchi identity

The preceding natural-source checkpoint verified the residual-corrected
identity before imposing the first equations. The new direct composition is

```
B0(QEred_general)+p Qphi_general e^y
 +2 deltaB_q(E1)+2S1 d eta+(1/2)L0 QC_actual-B0(PC1)=0,
deltaB_q(E1)_b=-q^{ac}nabla_a E1_cb
               +(1/2)q^{ac}nabla_b E1_ac-C1^d E1_db.
```

All eight rows vanish as collected exact coefficients, with the first
residuals independently constructed and left free. The wave helper was
independently derived and compared coefficient by coefficient before being
composed with the same QC_actual map. All needed angular polynomials are
reconstructed, rather than only compared after projection.

The homogeneous specialization needs the local first equations, including
their required prolongations: E1 and nabla E1 vanish, and S1=0. Harmonic
specialization separately needs C1=0 and its first derivative to remove
`B0(PC1)`. This gives

```
B0(QEred)=-p Qphi e^y-(1/2)L0 QC_actual.
```

Replacing QC_actual by QC_background inside L0 would additionally require
the appropriate second prolongation of C1. Keeping QC_actual avoids that
substitution. Equality of the QC values at a point does not suffice.
Separate coherent ordinary-jet examples demonstrate this distinction, and
also show why E1=0 at a point does not imply its derivative vanishes. They
are formal premise counterexamples, not new physical solutions or a new
obstruction for the retained on-shell state.

The preceding first-master/gauge checks supply the implementation joins
used for the local first equations and harmonicity, conditional on the
actual smooth solutions of their retained PDEs. Those are neighbourhood
identities; unrelated values chosen inside interval enclosures cannot
replace them.

## Connection to the second constraint

For second parameter coefficients `(k,zeta)`, the constraint covector is
`C2=DC[k]+QC_actual`. Combining the linearized Bianchi identity with the
second reduced Einstein and natural scalar equations and the displayed source identity
gives `L0 C2=0`, under the first/background premises above. The already
retained initial and physical-shell compatibility arguments supply the
homogeneous subsidiary data and boundary law. Uniqueness and regularity of
that problem remain identified analytic inputs. The finite coefficient
checks do not establish those PDE theorems anew, nor do they turn a numerical
response candidate into an exact constrained solution.

The B0 contraction annihilates transverse-traceless tensors. Its eight-row
test alone therefore cannot validate every natural tensor-source component.
The all-20 role joins retain the TT entries where effective; the older full
geometric source checks keep their separate responsibility.

## Actual third jets

The additive adapter preserves the original actual lower two-jet interface
and returns ordinary y/t derivatives through total order three: 70 metric
profile slots and 10 scalar profile slots. Of the 32 third-order entries,
24 match the earlier time-containing intervals exactly and eight pure yyy
entries are added. They enclose the same actual first fields on
`0<=x<=3, r>=0, x+r<=8`, with physical-side one-sided shell traces.

For any smooth composite profile F and `Lrad=1/(ell rho)`, the new radial
conversion includes

```
F_yyy=-Lrad^3 F_rrr-3a Lrad^2 F_rr+(a_y-a^2)Lrad F_r.
```

Raw actual master/gauge fourth-derivative certificates and coefficient jets
through third order are inherited authenticated inputs. The new operations
compose their derivative enclosures; they do not differentiate selected
interval values or unknown Taylor remainders. Independent nonlinear
coordinate substitution and separate product/inverse arithmetic check the
conversion and reconstruction. The third entries are coarse uniform
magnitude enclosures. Cartesian selections from those boxes need not be
coherent jets, and no smaller derivative approximation error is certified.

The numerical adapter and generic coefficient identity are separate results.
Neither supplies a new actual PDE/regularity theorem, an improved response
candidate, a nonlinear remainder, a physical force or observational evidence
for a higher-dimensional Big Bang origin.
