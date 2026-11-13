# All-mode momentum coercivity for the selected weighted inverse

13 November 2026 UTC. Reviewed mathematical continuation. This note concerns the
actual selected background momentum operator on its entire invariant
one-form space. It does not replace that operator by the three angular
channels used for the actual second source. The completed quantitative
claim below is an all-mode weighted H1 coercivity and inverse bound.
An evaluated inverse constant from H12 into H14 is not asserted here.

## 1. Operator, density, norm, and domain audit

All retained source paths in this note lie in
`research/ in the pinned base tree`, at the supplied base
`d26f530122f8bfa186291ce601e4876e3646d760`. In particular use
`HDBLAST_WEIGHTED_REGULARITY_20261106/{README,MOMENTUM_COERCIVITY,
MOMENTUM_REGULARITY,MOMENTUM_EXTENSION}.md` and the actual localized
initial-data theorem pinned there. No frozen file is modified.

The exact numerical member and physical background are

```
ell=10^-6, kappa_w=ell/100=10^-8,
A={ell/2<z<5ell/2} x S3,
r=(z-ell/2)(5ell/2-z)/ell^2,
q=kappa_w^2 r^2, w=exp(-1/r),
gamma0=dy^2+rho(y)^2 kappa_S3, dy=-rho dz, rho_y=a rho,
129<=rho<=131, 111/1000<=a<=112/1000.
```

The symbol `kappa_w` denotes the fixed weight scale, whereas `kappa_S3`
denotes the unit round metric. Write `rho0=sqrt(det gamma0)` for the
fixed density. The ordinary contravariant correction tensor is tau,
and its canonical density is pi=rho0 tau. With the fixed physical
pairings the linear operators are

```
(P_m tau)_i=-2 gamma0_ik nabla_j tau^jk,
P_m^*Z=(Lie_(Zsharp) gamma0)^ij,
L_M Z=w^-2 P_m(w^2 q^2 P_m^*Z).
```

Both raised tensor indices in P_m^* use gamma0. Lowering them commutes
with nabla and preserves their full Frobenius norm. Consequently this
is the symmetric-gradient energy, with its factor two retained; it is
not a trace-free conformal Killing operator. Multiplication by the
fixed density rho0 does not change the normalized operator or pairing.

The actual weighted norms are

```
||Z||_Hk(q,w)^2=sum_(j=0)^k integral_A w^2 q^(2j)
                                |nabla^j Z|_gamma0^2 dmu_gamma0,
||Z||_w=||Z||_H0(q,w), E(Z)=||q P_m^* Z||_w.
```

Every tensor and derivative index uses gamma0. The dot means closure
of compactly supported smooth sections in this precise norm. The
range norm is the same weighted norm on one-forms. Thus q is the
derivative scale: it multiplies each derivative, not its square root.
The principal coefficient q^2 of this second-order normal operator
has the correct scaling for

```
L_M: dot H14_(q,w)(T*A)^Gamma -> dot H12_(q,w)(T*A)^Gamma,
Gamma=(Z2)^4.
```

For intermediate homogeneous inequalities displayed with a strict sign,
the strict statement applies when the right-hand norm or pointwise
energy is nonzero; at zero both sides vanish and the sign is read as
non-strict. The final operator bounds below use non-strict signs.

The retained weighted isomorphism theorem applies to these spaces on
the complement of the adjoint Killing kernel; Gamma has no invariant
part in that kernel. This is a qualitative existence/uniqueness premise,
not a numerical high-order inverse constant. In particular there is no
norm or differential-order mismatch that would invalidate the requested
mapping. Its high-order constant must still be quantified.

Probability angular measure is convenient in intermediate formulas:
write `<f>=(2pi^2)^-1 integral_S3 f dOmega`. Changing to full area
measure multiplies every norm and energy below by sqrt(2pi^2), so all
homogeneous operator constants are unchanged. Neither an endpoint
Dirichlet trace nor an ordinary unweighted potential extension is
imposed. Physical corrections, rather than potentials, have the retained
zero-extension property.

## 2. Complete angular decomposition and excluded kernel

Decompose a one-form into its radial scalar part and its angular exact
and coexact parts on the round S3. There are no harmonic one-forms on
S3. Scalar harmonics obey

```
-Delta_S3 E=lambda E, lambda=l(l+2), l=0,1,2,... .
```

For lambda>0 and `<E^2>=1`, the radial/exact block is

```
Z=u(y)E dy+(rho/sqrt(lambda))v(y)dE,
X=(u,v), <|Z|^2>=|X|^2.
```

The scalar monopole has only u. The l=1 scalar eigenspace consists of
the four ambient coordinate restrictions; every coordinate reflection
changes the sign of its corresponding component. Hence an invariant
scalar has no l=1 part. Every nonconstant scalar block that needs to
be estimated therefore has lambda>=8. No upper bound is placed on l.

For a coexact one-form alpha, normalize `<|alpha|_kappa^2>=1` and use

```
delta_S3 alpha=0, Delta_H alpha=mu alpha,
mu=(l+1)^2, l>=1,
Z=rho v(y)alpha, <|Z|^2>=v^2.
```

This includes both coexact helicity families. At mu=4 these are exactly
the six round Killing one-forms. The ambient rotations K_ab transform
under a coordinate reflection by the character epsilon_a epsilon_b.
Each such character is nontrivial on Gamma, so their entire invariant
subspace vanishes, including radially varying coefficients. Every
remaining coexact block has mu>=9. The present proof is deliberately
valid for all mu>=9, whether or not a particular eigenspace actually
has an invariant vector.

The complete round vector spectrum used here is independently stated in
Lindblom, Taylor and Zhang, [*Scalar, Vector and Tensor Harmonics on the
Three-Sphere*, section 3, equations (12)--(23) and Table 1](https://arxiv.org/html/1709.08020v1#S3).
Their two divergence-free families have rough eigenvalue
`1-l(l+2)` at unit radius; adding Ric=2 gives the positive Hodge
eigenvalue `(l+1)^2` used above. Their complete vector expansion
includes the exact family and both coexact families. These are angular
unit-sphere spectra; they have not been substituted for the full
warped-cylinder connection Laplacian.

The warped background and radial weights commute with the round
isometries, Hodge decomposition, and the spectral projections. The
radial/exact and coexact adjoint tensor outputs are orthogonal after
angular integration. For example, their mixed one-form components
are exact versus coexact, and the angular tensor cross term reduces
by integration by parts to `<dE,alpha>=0`. Within a scalar block the
orthogonality of trace and trace-free Hessian gives the formulas below.
Thus Parseval sums their energies without omitted cross terms. An
invariant potential is exhausted by the monopole, all admitted
radial/exact blocks, and all admitted coexact blocks.

## 3. Exact all-mode adjoint formulas

For a scalar block set k=sqrt(lambda)/rho and
d_lambda=2(lambda-3)/3. Direct covariant differentiation gives

```
L^2:=<|P_m^*Z|^2>
 =4u_y^2+2(v_y-av+ku)^2
   +12(au-kv/3)^2+4d_lambda v^2/rho^2.                  (1)
```

This identity is valid for every lambda>=8, not just lambda=8,24.
It follows from `||Hess E||_L2^2=lambda(lambda-2)` and the orthogonal
trace-free Hessian square `lambda*d_lambda`. For the monopole,

```
L^2=4u_y^2+12a^2u^2.                                   (2)
```

For a coexact block the exact formula is

```
L^2=2(v_y-av)^2+2(mu-4)v^2/rho^2.                       (3)
```

To check (3), the two mixed tensor components contribute
`2(v_y-av)^2`. On the unit round S3, Weitzenbock gives
`||nabla alpha||_L2^2=(mu-2)||alpha||_L2^2`. Coclosedness and
integration by parts give
`integral nabla_i alpha_j nabla_j alpha_i=-2||alpha||_L2^2`.
Therefore `||nabla alpha+(nabla alpha)^T||_L2^2=2(mu-4)`.
This identity includes arbitrary multiplicity and both helicities.

For scalar blocks, the final square in (1) implies

```
k|v| <=sqrt[3lambda/(8(lambda-3))] L <=sqrt(3/5)L <(4/5)L,
|v| <=rho/(2sqrt(d_lambda)) L <36L.
```

The trace square then gives

```
|u| <=(1000/111)[1/sqrt(12)+sqrt(3/5)/3]L <5L,
|X|<37L.                                               (4)
```

For a purely rational replay it is enough to replace the first
inequality by `|u|<(1700/333)L<6L`; then
`36^2+6^2=1332<37^2` still proves the final bound. For the monopole,
`|u|<4L`. From (3), every remaining coexact block satisfies

```
|v| <=131/sqrt(10) L <42L.                              (5)
```

The delicate all-mode issue is the radial derivative. The old finite-
channel inequality `|X_y|<=L+(1/5)|X|` cannot be reused, because k
is unbounded. Instead keep the bilinear expression used by Hardy.
Writing j=v_y-av+ku, (1) gives

```
|u u_y+v v_y|
 <=[(1/2+sqrt(3/5))|u|+|v|/sqrt(2)]L+a v^2
 <(3/2)|X|L+a|X|^2.                                    (6)
```

Here the apparently dangerous product `k|u v|` is bounded using
`k|v|<=sqrt(3/5)L`, so no factor sqrt(lambda) remains. The elementary
inequality needed for the coefficient is
`(1/2+4/5)^2+1/2=219/100<9/4`. The monopole obeys the same bound.
For coexact blocks, (3) gives the stronger inequality

```
|v v_y| <=|v|L/sqrt(2)+a v^2
         <=(3/2)|v|L+a v^2.                            (7)
```

Equations (4)--(7) are independent of all angular eigenvalues and
provide the missing uniform tail input.

## 4. Endpoint-uniform Hardy and the all-mode L2 gap

Start with compactly supported smooth invariant potentials and finite
angular sums. On each endpoint tail set

```
t=1/r, sigma=sqrt(1-1/t), t>=2,
rho^3|dy|=h(t)dt, h=ell rho^4/(2t^2 sigma),
|q D_y|=[2kappa_w^2 sigma/(ell rho)]|D_t|.
```

The orientation of y differs between tails; all inequalities retain
absolute values. Since rho<=131 and sigma>=1/sqrt(2),

```
|D_t| <=10^-4 kappa_w^-2 |q D_y|.                       (8)
```

The actual geometric measure is essential. Its logarithmic derivative
is `h'/h=4a y_t-2/t-sigma'/sigma<0`, because
`|4a y_t|<0.000042/t^2<2/t`. Hence m=h exp(-2t) obeys
`m'<=-2m`. For a compact test F vanishing at t=2, integration by parts
gives the vector-valued Hardy inequality in its useful bilinear form

```
integral m |F|^2 <=integral m F dot F_t.                 (9)
```

No ordinary trace at infinity is assigned to the eventual completed
potential. Compact support removes that boundary term at this stage.

Choose a Lipschitz cutoff chi with chi=0 for t<=2, chi=1 for t>=3,
0<=chi<=1, and |chi_t|<=2. Smooth approximation preserves the
estimates. Sum (9) over both tails and all angular channels. Let

```
Y=||chi Z||_w,
I=||Z||_(w,t<=3), E=||q P_m^*Z||_w,
A_tail=10^-4 kappa_w^-2.
```

The cutoff derivative contributes at most `2Y I` by Cauchy--Schwarz:
it is the pairing of chi X with chi_t X, and the latter is supported
inside t<=3. Equations (6)--(8), Cauchy--Schwarz across all modes,
and q<=kappa_w^2/4 on the tail support yield

```
Y^2 <=2Y I+(3/2)A_tail E Y+delta Y^2,
delta=A_tail*(112/1000)*kappa_w^2/4=2.8*10^-6.
```

For Y>0 divide by Y; the case Y=0 is immediate. Since
`(1-delta)^-1<1.00001`,

```
Y<1.00001[2I+1.5*10^-4 kappa_w^-2 E].                  (10)
```

On the interior t<=3, q>=kappa_w^2/9. The common pointwise angular
estimate from (4)--(5) is |X|<42L, so

```
I<=378 kappa_w^-2 E.
```

The interior and the region where chi=1 cover the annulus. Consequently

```
||Z||_w <=I+Y
 <[378+1.00001(756+0.00015)] kappa_w^-2 E
 <=1200 kappa_w^-2 E=1.2*10^19 E.                      (11)
```

This includes arbitrarily large angular eigenvalues while t tends
to infinity. No relation between angular frequency and endpoint
distance was imposed. In particular it is not a fixed-interior
high-frequency estimate. The same constant passes to unrestricted
angular sums by Parseval and compact smooth approximation.

## 5. Full weighted H1 Korn estimate

The previous argument bounds the potential, not just its radial part.
A separate covariant integration identity also controls its full
weighted gradient. Put W=w^2q^2 and

```
J=||q nabla Z||_w, d_Z=||q div Z||_w,
delta_W=sup |q nabla log W|.
```

For compact one-forms the ordinary weighted Korn integration by parts
is

```
E^2=2J^2+2d_Z^2-2 integral W Ric_gamma0(Zsharp,Zsharp)
    +2 integral (nabla W dot Z) div Z
    -2 integral (nabla_i W) Z_j nabla_j Z_i.            (12)
```

This uses the symmetric gradient, not an injective map on arbitrary
unsymmetrized tensors. In particular it does not apply an invalid
pointwise lower bound from sym nabla Z to all of nabla Z.

The exact weight derivatives give

```
q nabla_y log W=2kappa_w^2(1+2r)r_y,
|r_y|<=2/(ell*129),
delta_W<=12kappa_w^2/(ell*129)<10^-11.
```

The retained background bounds also give `|a_y|=|a_z|/rho<10^-4`.
The radial/tangential Ricci eigenvalues of gamma0 are
`-3(a_y+a^2)` and `2/rho^2-a_y-3a^2`; hence their absolute values
are below 1. Denote Z0=||Z||_w. Cauchy--Schwarz in (12) gives

```
E^2 >=2J^2+2d_Z^2-2qmax^2 Z0^2
        -2delta_W Z0 d_Z-2delta_W Z0 J.
```

Minimizing the quadratic in d_Z and completing the square in J,

```
J^2 <=E^2+[(3/2)delta_W^2+2qmax^2]Z0^2,
J <=E+2*10^-11 Z0, qmax=kappa_w^2=10^-16.              (13)
```

Combine (11) and (13). With the derivative-square definition of H1,
or its larger literal-sum bound Z0+J, one obtains

```
||Z||_H1(q,w) <=2*10^19 ||q P_m^*Z||_w.               (14)
```

This is an all-mode Korn/coercivity bound with fixed evaluated
constant, uniform at both degenerating endpoints. All integrations
were for compact tests. Their inequalities extend by completion,
so no unproved endpoint boundary condition has entered the result.

## 6. Consequence for every source of the selected inverse

The actual selected potential for f in the invariant dotted H12
range belongs to the invariant dotted H14 domain by the retained
qualitative theorem. Approximate it in that norm by invariant compact
smooth potentials. Since E(Z)<=2||q nabla Z||_w, this approximation
also converges in the energy graph norm. The variational identity
therefore permits Z itself as a test and gives

```
E(Z)^2=<L_M Z,Z>_w=<f,Z>_w.
```

Using (11) and (14) proves, for every such source and the same
selected inverse,

```
E(Z) <=1.2*10^19 ||f||_w,
||Z||_w <=1.44*10^38 ||f||_w,
||Z||_H1(q,w) <=3*10^38 ||f||_w.                       (15)
```

Zero source gives zero potential. These bounds also identify the
energy-domain solution uniquely. They concern the actual inverse on
all invariant sources; no compact comparison tensor, altered right
inverse, or finite-support source assumption is used. Since
`||f||_w<=||f||_H12(q,w)`, (15) supplies a quantitative lower-order
part of the requested H14 estimate, but not its missing derivatives.

In fact the energy-domain weak solution exists for every invariant
weighted L2 source: complete invariant compact tests in E, use (11)
to bound the source functional, and apply the Hilbert-space Riesz
representation theorem to the E inner product. Equation (14) embeds
this completion continuously in the dotted H1 space. For sources in
the retained H12 range, uniqueness identifies this weak solution with
the selected H14 potential. This does not establish high-order
regularity for general L2 data.

## 7. Exact high-order obstruction to reusing the finite-channel proof

For completeness the all-mode transformed normal blocks can be written
without dropping their angular terms. Set

```
D=D_s=qD_y, U=w sqrt(qrho^3), V=U X,
b=D log U=kappa_w^2[(1+r)r_y+(3/2)ar^2],
c=q/rho, zeta=c sqrt(lambda).
```

For a radial/exact block the retained matrix formula remains exact at
every lambda>=8:

```
A_lambda=B D+C_lambda,
B=[[2,0],[0,sqrt(2)],[0,0],[0,0]],
C_lambda=[[-2b,0],
          [sqrt(2)zeta,-sqrt(2)(b+qa)],
          [sqrt(12)qa,-sqrt(12)zeta/3],
          [0,2c sqrt(2(lambda-3)/3)]],
U L_M U^-1=A_lambda^* A_lambda.
```

For a coexact block,

```
A_mu V=(sqrt(2)[D-b-qa]V, sqrt[2(mu-4)] c V),
U L_M U^-1=A_mu^* A_mu.
```

For the normalized source equation L_M Z=f, the transformed source
is Uf, not the transformed physical unnormalized canonical source
used in the old particular-source calculation. Unitarity is exact:
`||UZ||_L2(ds x S3)=||Z||_w`.

At any fixed interior point c>0, the matrix C_lambda contains entries
growing as sqrt(lambda), and its normal equation contains lambda c^2.
Thus the old estimate `||D^j K_i||<=10^-6` cannot hold uniformly in
angular eigenvalue. Its one-dimensional small-perturbation Fourier
argument is consequently unavailable for the all-mode operator. This
is a precise obstruction to that proof method, not a counterexample
to the retained weighted isomorphism.

Nor may angular high frequency simply be discarded in the endpoint
limit. Although c tends to zero, sequences with lambda tending to
infinity can keep zeta=c sqrt(lambda) fixed. In that joint limit the
principal radial/exact matrix is

```
L_pr=[[-4D^2+2zeta^2, 2zeta D],
      [-2zeta D-2zeta', -2D^2+4zeta^2]].                (16)
```

For frozen zeta its Fourier eigenvalues are exactly
`2(xi^2+zeta^2)` and `4(xi^2+zeta^2)`. Both transverse coexact
polarizations have principal eigenvalue `2(xi^2+zeta^2)`. Including
the nonzero exponential drift b at an endpoint replaces xi^2 by
xi^2+b^2 in the full frozen sym-gradient symbol; the corresponding
eigenvalues are again 2 and 4 times xi^2+b^2+zeta^2. This shows why
the full principal operator has to be retained in a uniform proof.

One useful exact positive-principal representation is

```
P=-D^2+zeta^2, d=(D,-zeta),
L_pr=2P I+2d^*d-2zeta' [[0,1],[1,0]].                  (17)
```

The last term has one fewer weighted derivative order because
`zeta'=(q_y-qa)zeta`. Quantifying this commutator structure, or a
full four-dimensional q-scale elliptic localization, is a possible
route to the missing estimate. It would still have to control all
mixed weighted derivatives through order fourteen and compare them
to the exact geometric dotted norms above. A frozen symbol alone
does not establish that regularity bound.

The remaining explicit target is a numerical C_reg with

```
||Z||_H14(q,w) <=C_reg[||L_M Z||_H12(q,w)+||Z||_w]
```

on the same whole invariant domain, with all endpoint, angular, and
norm-comparison constants charged. If supplied, (15) would give
`C_inv,M<=C_reg(1+1.44*10^38)`. Neither C_reg nor that full inverse
ceiling is supplied by the completed calculations in this note.

## 8. Scope

The completed result is a new uniform all-mode momentum gap and full
weighted H1 inverse estimate. It repairs the high-angular-frequency
failure in the previous tail argument, includes coexact sectors and
the correct reflection-kernel exclusion, and preserves the fixed
canonical density and dotted domain. It remains conditional on the
retained background and selected weighted existence theorem.

No evaluated H14 inverse constant, nonlinear initial ball, amplitude
radius, evolved tube, physical observable remainder, or publication
is claimed. All new bytes are private and the frozen checkpoints
remain unchanged.
