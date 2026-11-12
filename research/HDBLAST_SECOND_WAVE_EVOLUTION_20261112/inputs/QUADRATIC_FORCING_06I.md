# Complete second-wave forcing and a finite H12 input estimate

6 November 2026. This proof is conditional on the retained smooth background, the exact
background-subtracted boundary chart, and the regularity explicitly stated
below. It does not bound an evolved H14 solution, a finite-amplitude radius,
or the physical observable remainder. No frozen predecessor is changed.

## 1. Exact second variation, before changing the equation sign

Use a fixed bundle connection D and its fixed product-base Hessian H. Write
the exact transformed equation as

```
q(x,V) : H V = c(x,V) + b(x,V) D V + d(x,V)(D V,D V).       (1)
```

Here d is symmetric in its combined derivative/fiber arguments. The
coefficients in (1) are the actual fixed-bundle coefficients after the
substitution described in section 2, not the coordinate coefficients with
partial derivatives relabeled as covariant derivatives. We suppress x in
state derivatives. The background is the zero section V0=0, identically
throughout the cylinder, and c(x,0)=0.

Put

```
u = partial_e V|0,  v = partial_e^2 V|0,
q0=q(x,0), A=c_V(x,0), B=b(x,0),
C2=c_VV(x,0), E2=b_V(x,0), D2=d(x,0), Q1=q_V(x,0),
Lplus w=q0:H w - B D w - A w.
```

The amplitude expansion uses V=e u+e^2 v/2+o(e^2), so

```
Lplus u=0,
Lplus v= R2[u],
R2[u]= C2(u,u)+2 E2[u]D u+2 D2(D u,D u)
                         -2 Q1[u]:H u.                   (2)
```

Every coefficient in (2) is fixed at V=0. In particular the last term is
present even though V0=0. The expansion of q(V):HV at second order is
q0:Hv+2q_V(0)[u]:Hu. Terms q_V(0)[v]:HV0 and
q_VV(0)[u,u]:HV0 vanish because HV0=0. Similarly d_V is absent because
DV0=0. This is exact algebra for a smooth parameter curve, not an argument
from a finite sample of perturbations.

The coefficient 2 in both the mixed b-variation term and the principal
variation term is mandatory. The 2 in the d-term follows because the
coefficient d itself was defined as one half of the first-jet Hessian.

## 2. Exact coordinate-to-bundle coefficients and connection signs

The PR133 coordinate equation is

```
q^{ij} partial_i partial_j V = c_c(V)+b_c(V)P+d_c(V)(P,P),
P_i=partial_i V.
```

Let D_i=partial_i+Omega_i on the fixed bundle and let gamma^k_ij be the
fixed product-base connection. Its Hessian is

```
H_ij V = partial_i partial_j V
       +Omega_i partial_j V+Omega_j partial_i V
       -gamma^k_ij partial_k V
       +(partial_i Omega_j+Omega_i Omega_j-gamma^k_ij Omega_k)V.
```

The PR133 correction formula is written in the old P=partial V variables.
To obtain (1), define Q_i=D_iV and substitute P_i=Q_i-Omega_iV in
the entire right-hand side, including the correction. The exact net identity
is

```
q:H V = f_c(V,Q-Omega V)
       +q^{ij}{Omega_i Q_j+Omega_j Q_i-gamma^k_ij Q_k
                       +(partial_i Omega_j-Omega_j Omega_i)V}. (3)
```

The two gamma Omega terms cancel. Keeping the old zero-order correction
after also replacing P by Q would double-count terms. Because q is symmetric,
the contracted Omega_j Omega_i can equivalently be written Omega_i Omega_j;
one should not commute the matrices before contraction.

With K(V)=(Omega_i V)_i, and

```
S(V)Q=q^{ij}(Omega_i Q_j+Omega_j Q_i-gamma^k_ij Q_k),
T(V)=q^{ij}(partial_i Omega_j-Omega_j Omega_i),
```

the exact coefficients in (1) are

```
c(V)=c_c(V)-b_c(V)K(V)+d_c(V)(K(V),K(V))+T(V)V,
b(V)Q=b_c(V)Q-2d_c(V)(K(V),Q)+S(V)Q,
d(V)=d_c(V).                                             (4)
```

All state derivatives of q,S,T in (4) are retained when forming A,B,C2,E2.
The full varying physical reference connection Gamma[g0] remains inside
the coordinate coefficients. The connection Omega here is a separate fixed
bundle connection; it does not replace that physical reference.

## 3. Ordered positive-energy convention

Let X_z=partial_z and X_A=sech(t)E_A in a round left-invariant orthonormal
frame with [E_A,E_B]=2 epsilon_ABC E_C. The fixed product base connection
has

```
nabla_t X_A=-tanh(t)X_A,  nabla_XA t=0,
nabla_XA X_B=sech(t) epsilon_ABC X_C.
```

Since D_t D_XA u=D_XA u_t-tanh(t)D_XA u, its time-angular Hessian is

```
H(t,X_A)u=D_t D_XA u-D_(nabla_t X_A)u=D_XA u_t,
H(X_A,t)u=D_XA u_t.
```

The moving-frame term and the product-base Hessian correction cancel.
The purely angular frame-Hessian correction also cancels after contraction
with symmetric q^{AB}, because epsilon_ABC is antisymmetric. Consequently,
in this exact convention,

```
q:H u = -Pord u,
Pord u=a u_tt-2 b_pr^i D_Xi u_t-c_pr^{ij}D_Xi D_Xj u,
a=-q^{tt}, b_pr^i=q^{ti}, c_pr^{ij}=q^{ij}.                 (5)
```

This ordering identity holds before specializing V to zero. It does not
discard the moving term in the energy derivative, where differentiating
the spatial first jet still gives -tanh(t)W_A.

At background the operator used by the positive-energy proof is

```
Lminus=Pord,0+B D+A=-Lplus.
```

Therefore its second-wave force is

```
F2[u]=-C2(u,u)-2E2[u]D u-2D2(D u,D u)+2Q1[u]:H u.          (6)
```

Equations (2) and (6) differ only by the declared equation-sign convention.
No normal derivative is given an invented homogeneous boundary condition.

## 4. A global tensor H3 to L-infinity constant

Let M=[0,Z]xS3, Z=6ell=6*10^-6, with fixed product spatial metric
q_sp=dz^2+kappa and full measure dz dmu_kappa. Norms are full covariant
derivative square sums. We prove for every tensor/direct-sum bundle section
T with at most 16 covariant angular-or-spatial slots that

```
||T||infty <= E3 ||T||H3(M),  E3=10^7.                    (7)
```

This estimate needs no boundary condition or extension across either end.
It will only be applied to field derivatives, never to high-rank coefficient
maps; their derivative bounds are supplied separately in section 6.

For a scalar, or a fixed-Hilbert-space-valued component array, the round
sphere has eigenvalues lambda_l=l(l+2), multiplicity (l+1)^2, and full
volume 2pi^2. The addition identity and Cauchy-Schwarz give

```
|f(omega)| <= [sum_l (l+1)^2/(2pi^2(1+lambda_l)^2)]^(1/2)
                   ||(1-Delta_comp)f||L2(S3)
             = (1/sqrt(12))||(1-Delta_comp)f||L2(S3).       (8)
```

Indeed 1+lambda_l=(l+1)^2 and sum_{n>=1}n^-2=pi^2/6.
The vector-valued estimate has the same constant, without a component-count
factor. Smooth sections suffice for the calculation, followed by density.

Use the global round left-invariant frame, with dt and dz added as parallel
slots. On a tensor with r<=16 slots the connection matrices Omega_A are
constant component matrices, and ||Omega_A||op<=r. Since
nabla_EA EA=0,

```
Delta_comp T=tr D_Omega^2 T-2 sum_A Omega_A D_A T
                                      +sum_A Omega_A^2 T.
```

It follows that

```
||(1-Delta_comp)T||2
 <= (1+3r^2)||T||2+2sqrt(3)r||D_Omega T||2
                                      +sqrt(3)||D_Omega^2 T||2
 <=1000 ||T||H2(S3).                                      (9)
```

For example the sum of the three displayed coefficients is below
769+64+2=835 when r<=16. Applying (8) only decreases the constant.
The Frobenius metric on symmetric tensors is used, so the induced
orthogonal representation has precisely these slotwise bounds.

For a Hilbert-valued h on [0,Z], averaging h and using the fundamental
theorem of calculus proves

```
sup_z ||h(z)|| <= Z^-1/2 ||h||L2_z+sqrt(Z)||h'||L2_z.
```

Apply it to h(z)=T(z) in the H2(S3) Hilbert space. The product connection
commutes with z differentiation, and the required derivative combinations
are included in the full H3(M) norm. The new radial derivative slot is
product-parallel and contributes no angular connection matrix. In particular
the angular argument is applied to T itself; the Hilbert-space differentiation
does not blindly reapply a rank-sixteen bound to a tensor with seventeen
rotating slots. Here
Z^-1/2+sqrt(Z)<1000: its square is 1/Z+2+Z<10^6.
Combining with (9) gives a constant less than 10^6, so the rounded
E3=10^7 in (7) is valid with ample reserve.

## 5. Finite H12 products without a W13,infinity hypothesis

For a fixed bilinear coefficient map K with
max_{a<=12}||D_sp^a K||infty,op<=C, the full covariant Leibniz rule gives

```
||K(U,W)||H12 <= 10^14 C ||U||H14 ||W||H12,                (10)
||K(G,H)||H12 <= 10^14 C ||G||H13 ||H||H13.                (11)
```

Here all input/output tensor slots are included in the norms, and the
operator norm of D_sp^a K controls its whole derivative-index tensor.
To prove (10), expand a derivative of order j<=12 with a+b+c=j assigned
to coefficient, U, and W. If b<=11, put D_sp^b U in L-infinity by (7),
since b+3<=14. If b=12, then a=c=0; instead put W in L-infinity using
its H3 norm, which is controlled by H12. In either case the other factor
is in L2. This is a finite regularity product estimate, not an assumption
that all twelve derivatives of U are pointwise bounded.

For (11), min(b,c)<=6. Put that derivative of its factor in L-infinity
using (7), and the other in L2. All required orders are at most 13.
The field-plus-derivative tensor rank is at most 16 in these applications;
the high-rank coefficients are never passed to (7).

For every order the multinomial sum is at most 3^12. Summing the 13
output derivative orders costs sqrt(13)<4. Thus the prefactor is at most
4*3^12*10^7<10^14, proving both displayed inequalities. No commutation of
normal derivatives with the boundary problem was used.

The companion linear multiplier estimate is

```
||K U||H12 <= 2^12 sqrt(13) C ||U||H12 <=10^6 C ||U||H12.   (12)
```

## 6. Actual coefficient envelope and derivative accounting

[ORDERED_BACKGROUND.md](ORDERED_BACKGROUND.md), section 4, proves on the actual
one-sided collar and |t|<=13ell/8, the common bound

```
C=10^152000
```

for every fixed-product spatial covariant derivative through order 12 of
the six actual background coefficient tensors

```
A=c_V(0), B=b(0), C2=c_VV(0), E2=b_V(0), D2=d(0), Q1=q_V(0).
```

These are Hilbert multilinear operator norms, with all base and bundle slots
charged. They are not derivative bounds for an evolved coefficient field.
The audit's finite accounting is as follows: PR133 supplies ordinary mixed
coefficient jets through total order 16 below 10^150000; the needed state
orders are at most two, leaving at least twelve base orders. The actual
q inverse jet is below 16!*10^17000<10^18000. The coordinate connection
jets through order 16 are below 10^60. Substitution (4) has at most two
Omega factors, at most 5^16 Leibniz terms and at most 16^2 state sums,
yielding coordinate coefficient partial-entry sums below 10^150200.
Expanding j<=12 covariant derivatives has branch factor below 10^4 and
at most j connection factors, each with derivatives below 10^60. Final
Frobenius/base-coordinate conversion costs less than 10^6. Thus

```
10^(150200+64*12+6)=10^150974 <10^151000.
```

The common 10^152000 also covers conversion to the moving ordered spatial
frame and its cosh/sech derivatives. The source pins and replay bind the complete ordered derivation; this
outline does not independently re-derive every coordinate coefficient entry. The bounds concern
the full maps A,B,C2,E2,D2,Q1, including connection and state-flow terms.

## 7. Eliminating the first-wave time acceleration at finite regularity

For any time where the first linear wave u solves Lplus u=0 and has

```
N14(t)=||u(t)||H14(M)+||u_t(t)||H13(M)<infinity,
```

its background principal coefficient is diagonal:
q0=rho^-2 diag(-1,1,sech(t)^2 kappa^-1) in the fixed product frame.
Consequently the exact acceleration equation is

```
u_tt=u_zz+sech(t)^2 Delta_kappa^D u
                         -rho^2(Au+B D u).               (13)
```

Delta_kappa^D is the bundle rough Laplacian, the trace of the fixed
angular Hessian. The trace is taken after differentiation; no reordering
of derivatives, false boundary row, or normal recovery lemma is inserted.
Since the fixed spatial principal tensor has Frobenius norm at most 2,
its H12 contribution is bounded by 2||u||H14.

The background source certificate gives
max_{j<=12}|partial_z^j rho^2|<10^210. All other spatial derivatives of
this scalar multiplier vanish. Hence the coefficient derivative envelope
for rho^2 A and rho^2 B is at most 2^12*10^210 C. Applying (12), and
using ||u||H12+||D u||H12<=2N14, yields

```
||u_tt||H12
 <=2N14+2*(2^24 sqrt(13))*10^210 C N14
 <=10^222 C N14.                                          (14)
```

For the arithmetic, 2^48*13<10^18 implies
2^24 sqrt(13)<10^9, and 2+2*10^219 C<10^222 C for C>=1.
This proves an H12 acceleration estimate from H14/H13 wave data. It does
not assert that a second time derivative belongs to H13.

Let J1=D u be the full spacetime first jet and J2=H u the full spacetime
Hessian, with positive fixed product base norms. Then

```
||J1||H13 <= N14,
||J2||H12 <= (2+10^222 C)N14.                             (15)
```

The 2 charges both mixed time-space Hessian entries and the purely spatial
Hessian; the time direction is product-parallel, so H_ti u=D_i u_t.

## 8. Complete forcing bound and the evaluated initial value

Apply (10)-(11) to every term of (6), retaining its displayed factors two.
Equations (14)-(15) give the following concrete conditional forcing map:

```
||F2[u](t)||H12
 <=10^14 C [5+2(2+10^222 C)] N14(t)^2
 <=10^240 C^2 N14(t)^2
 =10^304240 N14(t)^2.                                   (16)
```

The principal term Q1[u]:Hu, including u*u_tt and u*D_sp^2u, is present
in this estimate. The lower forcing terms do not replace it. All tensor
contractions are controlled by their full Hilbert operator norms, so there
is no hidden scalar-per-component estimate or omitted sixteen-field factor.

The preceding initial-transport capsule proves N14(0)<10^5200 in exactly
this fixed product norm. Therefore (13) evaluates the complete initial
time acceleration, and (16) yields

```
||F2[u](0)||H12 <10^(304240+10400)=10^314640 <10^315000.   (17)
```

In particular its L2 norm obeys the same ceiling. The margin also covers
two identical reflected copies; use each copy's positive inward coordinate.
The first initial packet vanishes in neighborhoods of both cylinder ends,
so its initial acceleration and the resulting forcing have that same
support property. No negative-physical-z background or cutoff derivatives
outside the identity collar are used.

## 9. Exact replay ledger and limitations

The exact replay checks the general second-amplitude polynomial
identity, all factors two, both equation-sign conventions, the substitution
P=Q-Omega V in the connection formula, and a noncommuting-matrix instance
that rejects a wrong Omega-order identity before symmetric contraction.
It also verifies the following exact rational/integer inequalities:

```
for r=0,...,16: (1+3*r*r)+4*r+2 < 1000,
1/Z+2+Z < 10^6,  Z=6/10^6,
4*3^12*10^7 < 10^14,
2^24*13 < 10^12,
2^48*13 < 10^18,
150200+64*12+6 < 151000 < 152000,
2+2*10^219 < 10^222,
10^14*(9+2*10^222) < 10^240,
304240+2*5200=314640 <315000.
```

The line bounding 2^24*13 supplies (2^12 sqrt(13))<10^6.
The acceleration and forcing inequalities with C>=1 follow monotonically
from the displayed integer checks after dividing by C and C^2 respectively.
The sphere spectral identity (8), tensor frame calculation (9), Leibniz
allocation argument, and source coefficient derivation remain mathematical
proofs; arithmetic replay is not a substitute for those arguments.

Result (17) is a numerical initial forcing bound. Result (16) is a valid
finite-Sobolev input map conditional on N14(t); it does not supply a
numerical evolved N14(t). Closing that wave norm still needs the retained
boundary-compatible high-order commutators and normal recovery. A force
bound in H12 does not by itself close the H14 energy hierarchy or produce
a nonzero nonlinear amplitude interval. The family patch/bootstrap and
proper-time, embedding, cut-area, denominator, and common-domain observable
estimates remain separate required steps.
