# Continuous full-angular coefficient and tensor-response certificate

All statements retain the frozen regular background, actual first master
and gauge, complete source identities, physical second initial prescription
and subsidiary premises. Parameter coefficients use the factorial convention
`g=g0+epsilon q+epsilon^2 k/2+...`, and similarly for the scalar.

## Angular and boundary spaces

The fixed coefficient ordering comprises K,L,N,J,w at degrees0,2,4;
v,u at degrees2,4; H at degrees2,4; and T+,T−. Every primitive operator
is block diagonal on the orthogonal groups of sizes5,8,8,1,1. The angular
Gram diagonal is the one authenticated from the full geometric operator.
Only L and v are Dirichlet. There are five Dirichlet and eighteen Robin
amplitudes; these count modal rows rather than spacetime dimensions.

The source and initial data for the registered state are invariant under
an orientation-reversing sphere isometry exchanging T+ with T−. Direct
all-point coefficient-map checks give equal TT sources, identical scalar
operators and zero TT initial data. Uniqueness of the inherited formal
problem gives T+=T−. Its 22-dimensional active representation has common
TT Gram weight1/3; the source is not divided by two. Both ambient slots
remain in all23 numerical and structural audits.

Writing `s=(-N+3J)/2`, `qtrace=N+J`, the Robin trace norm obeys
`N^2+3J^2=s^2+3qtrace^2/4`. The homogeneous physical Robin matrix is an
orthogonal direct sum of three identical3×3 blocks and nine zeros:

```
B3 = [[0, 2a, 2sigma1/3],
      [4a,-4a, 4sigma1/3],
      [0, p/2,-sigma2/2]],       M=-ell*rho*By.
```

Here `sigma1=-2p`, `sigma2=4phi0`. No inverse or coercivity of M is used.
The complete inhomogeneous shell source has an exact rank-eight packing;
that is a source representation, not a claim of an eight-dimensional
response. The inherited Codazzi/Bianchi subsidiary argument is required
in addition to algebraic angular closure for formal physical completion.

## Candidate and complete norm error

Use `ell=1e-6`, `x=t/ell`, `r=z/ell` and

```
D+={0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x},
X(u)=max(sup|u|H/3,sup|u_x|H,sup|u_r|H).
```

All fields/differences considered here vanish on `x+r<=1`. Hence
`sup|u|H<=3X(u)/4` on D+ and the shell-history trace is at most X(u)/4.
The inherited full-space contraction is `q=3409/8000000`.

Let U be the actual total smooth formal coefficient satisfying `U=b+TU`.
The exact rational registered displacement and a fixed continuous P1
velocity define B0 by the reflected d'Alembert formula. Reflection is
PN−PD. For feedback only, approximate the rational displacement and its
radial derivative independently by continuous P1 profiles. Their derivative
need not equal the derivative of the displacement P1 profile. With the
complete constant matrices, these give the ridge forcing and Robin lift B1.
The half-line source is reflected after its full coupled matrix action.
The resulting finite representation V=B0+B1 is C1 on D+.

Let deltaU,deltaUr,deltaV be actual-minus-rational initial H errors, and
epsU,epsUr,epsV the P1 errors. Set

```
efree=max((deltaU+(3/2)deltaV)/3,deltaUr+deltaV),
S=(3/2)*17387513156573+3851400=52162547172519/2,
e0=efree+epsV+S.
```

S is the entire retained bulk/shell forcing bound, charged once even when
its waveform is not inserted into B0. Nominal continuous H caps give

```
X0=max((nomU+(3/2)(nomV+epsV))/3,nomUr+nomV+epsV),
Ucap=(X0+e0)/(1-q).
```

For constant-matrix bounds nA,nB,nD,nM and actual-minus-constant bounds
dA,dB,dD,dM, put

```
df=(dA+dB+3dD/4)*X0+(nA+nB)*epsUr+nD*epsU,
dg=dM*X0/4+nM*epsU,
e1=(3/2)*df+dg.
```

Each matrix bound is the maximum of its five Gram-weighted block
Frobenius bounds. Off-block entries are checked exactly zero. The full
coefficient assembly retains angular Laplacian, connection and scalar
terms. In particular V,V',V'' have their registered normalization without
an extra half, `h=tanh(ell*x)/rho` is the time connection, and the boundary
coefficients are shell constants. All coefficient approximation errors are
charged, including the choices preserving the frozen mean matrices.

The exact identity and its measured bound are

```
U−V = T^2 U+(I+T)(b−B0)+(TB0−B1),
X(U−V) <= q^2 Ucap+(1+q)e0+e1 = E <1.1787044291e28.
```

This is a whole-domain estimate, not interpolation between sampled
space-time values. The six sample tables only illustrate the representation.
Ninety exact feedback entries agree with the unchanged five-mean candidate.
That global candidate is distinct from the sharper later mean-only
two-feedback calculation; its coarse error does not supersede those signs.

The full integral-equation residual `r=V−b−TV` satisfies

```
X(r) <= e0+(1+q)e1+q^2 X0 = R <1.1790566331e28,
X(U−V) <= R/(1−q) <1.1795592728e28.
```

The residual bound already includes the omitted iterates. No additional
tail is added to its resolvent bound. The two bounds obey the exact check
`R/(1−q)−E=2q e1/(1−q)`.

For the frozen homogeneous Robin matrix, the finite candidate obeys
`PN V_r−Mstar PN V=Mstar PN(uP1−u_rational−B1)` on the shell.
Its residual against the actual inhomogeneous law additionally includes
`(Mstar−M)PN V−G`. These coefficient, interpolation, source and iteration
defects are charged above. Its residual against the actual Robin law is
bounded in H by `(1+1/8000)E`, using the inherited `||M||<=1/2000` and
the error's shell zero-wedge estimate. This is an integral/C1 certificate,
not a numerical pointwise C2 Einstein residual.

## Angular derivative and pointwise scope

The retained selected-space evaluation factor is `sqrt(321/19)`.
It gives pointwise value error <3.6336399733e28 and scaled x/r first
derivative error <4.8448532977e28. These are frame-component bounds.

There is also intrinsic unit-sphere angular control. Saturating the
selected space by SO(4) gives the inherited287-dimensional section space W.
The equivariant derivative image Dkappa(W) has dimension at most287 and
constant evaluation-kernel trace under sphere rotations. Its evaluation
norm is at most sqrt287. The rough spectrum is at most24. Therefore

```
sup|Dkappa(U−V)|F <=(3/4)*sqrt(287*24)*E <(255/4)E.
```

This uses the intrinsic sphere connection and summed positive fibre norm.
A physical five-dimensional angular covariant derivative additionally
requires the warp/frame conversion. No second x/r derivative bound is
implied. The separate lemma provides the associated-bundle proof.

## Source-resolved total TT coefficient

The TT scalar equations have

```
A_T=-3ell*tanh(ell*x), B_T=-3ell*rho*a,
D_T=ell^2[-6/cosh(ell*x)^2
         +rho^2(-2a^2+2h^2−2bcurv−8dcurv+4V/3)], M_T=0.
```

They have zero initial data and no linear coupling to other amplitudes.
Take the exterior backward cone with apex (T,r0)=(1/4,5/4). It has
`1+x<=r<=3/2−x`, so it never reaches the shell or reflected region.
Consequently the free and shell terms are exactly zero there. The
quadratic source response is the total TT coefficient on this cone.

All78 terms of the complete actual geometric TT forcing are evaluated
with the authenticated actual first-field jet errors. Twenty clipped cells
at fixed N=4 cover the cone with exact total area1/16 and retarded factor1/2.
The Duhamel interval is widened by the entire scalar feedback tail:

```
qT=T*(aT+bT+T*dT/2)<1,
tail <=T^2*Fmax*qT/[2(1−qT)] <2.803e−8.
```

The resulting exact total coefficients have the outward enclosure

```
T+=T− in [−0.001311761536,−0.000001509027] <0.
```

The optional N=8 source evaluation was unnecessary. Independent geometry,
operator, source-API spot, interval and mutation reviews pass. The global
V has TT entries zero and cannot resolve this small value with its much
larger uniform error; the separate scalar calculation supplies the result.
This registered-gauge tensor coefficient is not a gravitational-wave
measurement, energy density or gauge-invariant tensor observable.

All numerical decisions use exact rationals and outward intervals. General
retarded integration, representation theory and norm arguments are standard
methods. The result is a model-specific verified construction; physical
novelty and observational relevance require separate work.
