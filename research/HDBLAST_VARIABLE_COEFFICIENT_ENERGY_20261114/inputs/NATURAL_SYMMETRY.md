# Natural SO(4) symmetry, the actual first-wave quadrupole, and explicit angular constants

Conditional mathematical continuation for the retained smooth background and
unique background linear mixed-wall problem. SOURCE_PINS.json authenticates
the public inputs. The high-order theorem in HIGH_ORDER.md uses the generic
rotation identity below and does not require a Casimir-eight restriction.


## Result and source interface

The complete selected first wave, including the proper-y root, all sixteen fields and the nonzero temporal-radial gauge velocity, belongs to the **natural SO(4) Casimir-eight subspace**. The exact background linear evolution and both homogeneous boundaries preserve that subspace. This is a statement about natural tensor pullbacks, not scalar spherical harmonics of individual frame components.

For a section T in that subspace, covariant angular derivatives satisfy the particularly simple estimate

```
||D_Omega^j T||L2 <= 8^(j/2)||T||L2.
```

There is no tensor-rank factor in this inequality; the natural action on every derivative slot is essential. In particular,

```
||T||H14_Omega <= sqrt((8^15-1)/7)||T||L2 <= 3*10^6 ||T||L2,
||Delta_E T||L2 <= 8 sqrt(3)||T||L2 <= 14||T||L2.
```

On the fixed product collar, writing `||T||H_z^m²=sum_(k=0)^m ||partial_z^k T||L2²`,

```
||T||H_z^m <= ||T||H^m <= sqrt((9^(m+1)-1)/8)||T||H_z^m
                                      <=10^7||T||H_z^m, m<=14.
```

The last bound also applies to `partial_t T` and all existing mixed t,z derivatives of T. General homogeneous estimates are written with <= so zero fields are included. Strict coefficient comparisons above concern positive constants, not 0<0 assertions.

The directly inspected public sources are:

| Source | SHA256 | Role |
| --- | --- | --- |
| research/HDBLAST_INITIAL_TRANSPORT_20261106/INITIAL_TRANSPORT.md | 8a0cfa469037427f8636344c7635fcdd8249039564d16242d3f050ce912ed0bc | Complete first Cauchy values/velocities, fixed norms and exact initial chart differential |
| research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md | 1db76f9b9bd187cf264fcbc119107c76ceb66ff356feb7a979b7db83b01b4e5d | Literal frozen-shell auxiliary vector, constant shell normalizations and homogeneous boundary rows |
| research/HDBLAST_ORDERED_EVOLUTION_20261106/ORDERED_BACKGROUND.md | 0b85d3f5dcbdd2b04c8732687ce7691f5649e250fda8a27cec209fafd741dd62 | Exact product-bundle background operator, all lower maps and retained boundary conditions |

The pointwise Casimir/connection identity and resulting constants were also independently derived by a delegated internal reviewer. This is internal mathematical review, not external peer review or physical evidence.

## 1. Natural representation and normalization

Identify the unit round sphere with `{n in R4: |n|=1}`. For `1<=a<b<=4`, use exactly

```
K_ab = n_a partial_(n_b)-n_b partial_(n_a).
```

These are restrictions of unit skew generators in the ambient standard normalization; their scalar sum of squares is `sum K_ab²=Delta_kappa`. Equivalently the six ambient skew matrices are orthonormal for `-tr(AB)/2`.

The fixed bundle is

```
E=R_v + Sym²(R dt+T*S3)_H + R_w + (R dt+T*S3)_D.
```

For a sphere isometry R, extend it by the identity on t,z, and use its natural pullback on every tensor slot. The positive metric on E is the fixed direct-sum metric with full Frobenius norm on the symmetric block. It and the full sphere volume are preserved. The induced infinitesimal operators are

```
R_ab T = Lie_(K_ab) T,
C_E = -sum_(a<b) R_ab².
```

The symbols R_ab denote the natural Lie lift, including tensor-slot terms. They are not component differentiation by K_ab in a chosen moving frame. Time and radial covectors are invariant under the group action. Invariance of the positive fiber metric and sphere volume gives the exact skew-adjoint identity

```
<R_ab T,S>L2 = -<T,R_ab S>L2.
```

There is no angular boundary term. The same identity holds after integration over z because the Killing fields have no radial component.

For a covariant rank-r tensor,

```
R_ab T = D_(K_ab)T + rho_r(D K_ab)T,
```

where rho_r is the sum of the natural endomorphism action on its r slots. Parallel t,z slots contribute zero. The operator norm of `D K_ab` on the sphere is at most one, so `||rho_r(D K_ab)||<=r`. This conventional slot bound is useful, but the sharper orthogonal identity in section 4 removes it entirely from the Casimir-eight derivative estimate.

## 2. Exact commutation of the full operator and boundary domain

Every part of the retained nonlinear transformed system is natural under sphere isometries:

1. The actual background metric `dy²+rho(y)²[-dt²+cosh²(t)kappa]` and scalar `phi0(y)` are invariant. The Einstein–scalar equations are tensorial. The relative harmonic gauge uses the difference of the physical and background connections, a tensor. Retaining the full varying background connection is consistent with this invariance; removing selected connection terms is not justified by symmetry.
2. The raw state consists of a radial component square root, a tangential symmetric tensor, a scalar and a radial cross covector. The square root `sqrt(g_yy)` is a scalar under angular isometries. Constant rho_b normalizations and moving rho(z),cosh(t) frame factors are angularly invariant.
3. In the exact auxiliary vector, the tangential operations are `H^-1`, multiplication and `tr(H^-1 J4)` with `J4=-dt²+kappa`. The scalar tension functions and shell coefficients are invariant. Thus the vector is equivariant on all sixteen states. Uniqueness of its local autonomous ODE implies equivariance of both F_z and F_-z, not only their first differential.
4. The fixed product connection is preserved by sphere isometries. Consequently the exact coordinate-to-bundle substitution and the full transformed differential operator are equivariant. Differentiation at the identically zero section gives the background linear operator identity.

It follows that, for the complete operator including lower-order terms,

```
[L_bg,R_ab]=0,
[partial_t,R_ab]=[partial_z,R_ab]=0.
```

This does not assert that each coordinate coefficient matrix separately commutes with ordinary angular differentiation. Such component statements would omit frame transformations and connection terms.

The four-dimensional D block is the full radial cross covector `R dt+T*S3`; the twelve-dimensional N block contains the scalar root, the full symmetric tangential tensor and the scalar field. Their fixed projectors commute with the natural group action. Sphere rotations fix z=0 and z=Z and commute with D_z. Therefore they preserve

```
V_D=0, D_z V_N=0 at z=0;
V=0 at z=Z.
```

They preserve every homogeneous corner condition used by the retained linear problem as well. This argument applies to rotations and time derivatives; it does not impose homogeneous rows on normal derivatives of a solution.

## 3. The full first Cauchy jet and its transported Casimir

Write `Y=n1²-n2²`. It is the restriction of a homogeneous harmonic polynomial of ambient degree two, since the two ambient second derivatives cancel. The polar Laplacian identity in R4 gives

```
-Delta_kappa Y=2(2+2)Y=8Y,
-sum_(a<b) K_ab²Y=8Y.
```

The complete raw first data in the inspected initial-transport source are

```
g_tt,1=-2rho²uY, g_ti,1=0,
gamma1=uY(-4dy²+2gamma_angular), phi1=vY,
g_tt,t,1=-2rho²utY,
gamma_t,1=2utY(-2dy²+gamma_angular),
g_ty,t,1=-6rho²(u_y+4au)Y,
g_tA,t,1=0, phi_t,1=vtY.
```

Here u,v,ut,vt,rho,a are radial functions in this initial formula. Every nonzero term is Y times an invariant tensor: dt², dy², the angular metric, the symmetric dt–dy mixed tensor, or the scalar identity. Thus each is an equivariant linear image of the same scalar quadrupole module, and the entire raw value and velocity satisfy the natural tensor Casimir equation with eigenvalue eight. In particular the temporal-radial mixed velocity is nonzero and included; its invariant dt–dy tensor does not introduce a new angular representation. No angular mixed first velocity was silently dropped.

The positive raw-to-product fiber isometry I, the radial normalization map L(z), and the root differential R1=Dr(0) are natural maps. Their distinct time/angular coefficients remain invariant because angular rotations do not mix the temporal or radial lines with the sphere. The exact initial chart differential is

```
A(z)=DF_-z(0) R1 L(z).
```

Equivariance of F_-z implies that A intertwines the natural group actions, even when it mixes equivalent scalar/tensor copies. The transported initial first value is `V1(0)=A I(X1)`. The proved initial even-frame identity is `A_t(0)=0`, so the full first velocity is `V1,t(0)=A I(dot X1)`. Both therefore lie in the Casimir-eight subspace.

The retained unique linear solution operator commutes with sphere isometries because its full equation and boundary domain do. Equivalently, for a sufficiently regular solution u=V1, `(C_E-8)u` solves the same homogeneous linear mixed-wall problem with zero Cauchy data. Uniqueness gives

```
C_E u(t,z,.)=8u(t,z,.).
```

One may instead use the bounded compact-group isotypic projector onto the scalar quadrupole representation and uniqueness; this avoids imposing extra derivatives merely to commute a differential Casimir. At the retained H14/H13 regularity, the two-Casimir-derivative uniqueness argument also has ample low-energy regularity. All existing time and radial derivatives retain eigenvalue eight, since the group action is independent of t,z.

## 4. A sharp connection/Casimir identity including every derivative slot

Fix the north pole n=e4 and an orthonormal tangent basis e1,e2,e3. The six Killing fields satisfy

```
K_a4=-e_a,  D K_a4=0                 (a=1,2,3),
K_ab=0,    D_X K_ab=X_a e_b-X_b e_a (1<=a<b<=3).
```

The first three natural Lie lifts at this point are just the three covariant directional derivatives (with a sign), and the other three are pure algebraic slot rotations. Hence, exactly,

```
sum_(a<b)|R_ab T|²
  = |D_Omega T|² + sum_(1<=a<b<=3)|rho_r(A_ab)T|²
  >= |D_Omega T|².                                      (S1)
```

The sum over all six generators is invariant under SO(4): the generators transform by its orthogonal adjoint representation, and the tensor action preserves norms. Therefore (S1) holds at every sphere point with the corresponding isotropy endomorphisms. The orthogonality here removes potential cross terms exactly; it is stronger than a triangle inequality with rank-dependent errors.

Integrating and using skew-adjointness yields

```
||D_Omega T||² <= sum_(a<b)||R_ab T||² = <C_E T,T>.
```

Sphere isometries preserve the Levi-Civita/product connection, so covariant differentiation intertwines their **natural action on every old and new derivative slot**:

```
R_ab^(r+j)(D_Omega^jT)=D_Omega^j(R_ab^(r)T),
C_(r+j)(D_Omega^jT)=D_Omega^j(C_r T).                    (S2)
```

If C_r T=8T, the same is true of every derivative tensor under its own natural lifted Casimir. Iterating (S1) on those derivative tensors gives

```
||D_Omega^j T||L2² <= 8^j||T||L2².                      (S3)
```

There is no restriction on tensor rank in (S3). In particular it automatically includes the original rank-at-most-two field slots plus up to fourteen derivative slots. It neither treats these derivative slots as scalars nor relies on a componentwise harmonic band.

The calculation is valid for smooth sections and extends to the needed finite-Sobolev sections by density/closedness. If desired, group convolution smooths the angular variable while commuting with the Casimir, preserving the eigenvalue condition and the uniform estimates.

For completeness the Casimir is central in the enveloping algebra of the six generators in this normalization. Thus R_ab preserves its eigenvalue-eight subspace. For the sum over all ordered words of length k,

```
A_k(T)²=sum_(|I|=k)||R^I T||L2²=8^k||T||L2².            (S4)
```

This follows by applying the k=1 skew-adjoint identity to each word and induction; the generators need not commute with one another.

## 5. Explicit angular and full product norm constants

Summing (S3) through fourteen gives the exact squared angular constant

```
sum_(j=0)^14 8^j=(8^15-1)/7=5026338869833 <(3*10^6)².
```

At order fourteen alone, `8^(14/2)=2097152`. Taking the trace of the rank-two angular Hessian costs sqrt(3), hence

```
||Delta_E T||L2 <= sqrt(3)||D_Omega²T||L2
               <=8sqrt(3)||T||L2 <=14||T||L2,
```

since `3*8²=192<196=14²`. This includes tensor rough Laplacians; it does not identify them with a scalar component Laplacian or claim that `Delta_E T=-8T`.

For the full product derivative, angular and radial differentiation commute, including the fixed radial derivative slots. Decomposing the ordered derivative slots according to whether they are radial or angular gives exactly

```
||D_sp^s T||L2²
 =sum_(k=0)^s binom(s,k)||D_Omega^(s-k) partial_z^k T||L2²
 <=sum_(k=0)^s binom(s,k)8^(s-k)||partial_z^k T||L2².
```

All `partial_z^k T` remain in the same natural Casimir subspace. Summing s=0,...,m and bounding each coefficient by the total coefficient sum gives

```
||T||H^m² <= [sum_(s=0)^m 9^s] sum_(k=0)^m||partial_z^k T||L2²
           =[(9^(m+1)-1)/8]||T||H_z^m².
```

For m=14 the exact factor is `25736391511831<10^14`, proving the stated 10^7 norm conversion. The reverse inequality has constant one, since all-radial derivatives are among the full derivative components. These statements are unchanged for an identical reflected second copy, with its own positive inward radial coordinate.

## 6. Optional general rotation comparison

Even without restricting the angular spectrum, (S1) and the intertwining (S2) imply

```
||D_Omega^k T||L2 <= A_k(T).
```

Thus the conservative 40^k rotation-to-covariant reserve is conservative. For the opposite direction, (S1) and the at-most-three isotropy generators at a point give

```
sum_alpha ||R_alpha T||² <= ||D_Omega T||²+3r²||T||².
```

Induction using (S2), for original rank plus derivative count at most sixteen, yields

```
A_k(T) <= 28^k ||T||H_Omega^k,
```

because `1+3*16²=769<28²`. This also shows that the 60^k converse reserve is safe. These general comparisons are not needed for the sharper selected-first-wave Casimir route.

## 7. Consequence for the next high-order ledger and limits

For the selected first wave, a tangential energy sum may use only the fourteen time-commuted fields `partial_t^p u`, 0<=p<=13. The principal time-commutator angular term can use `||Delta_E partial_t^(p-k)u||2<=14||partial_t^(p-k)u||2`; it does not require an additional rotation word or an additional derivative of data. Radial normal recovery remains a separate finite triangular calculation. Once it controls the radial H14/H13 norms, section 5 supplies their conversion to full spatial N14 with constant 10^7. This observation simplifies the angular ledger but does not itself evaluate its remaining temporal/radial recurrence constants.

The existing general rotation-word route remains valid and can retain its conservative constants if already used. It does not need to be replaced to obtain a theorem.

The eigenvalue-eight restriction must not be transferred to the second wave or to a nonlinear family. Already the scalar product Y² has nonzero sphere average, so it cannot satisfy C(Y²)=8Y². Quadratic tensor products may populate additional natural representations; scalar component labels alone do not identify them. Likewise the simple Laplacian bound above is for the natural Casimir-eight field and is not an angular norm estimate for arbitrary coefficient maps.

This audit establishes symmetry and finite norm comparisons for the specified homogeneous background first wave. It supplies no finite nonlinear amplitude radius, actual-family metric-patch proof, evolved second-wave closure, proper-time observer/area remainder, external novelty claim, or evidence of a higher-dimensional cosmological origin.
