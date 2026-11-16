# Strong lower mixed limit of the actual frozen iteration

14 November 2026. Conditional continuation from the accepted conditional smooth Picard admission. This note proves the lower mixed limit once the evaluated lower-energy factorial Cauchy estimate is available. It does not assert continuous full X14, an X15 estimate, a parameter-Hessian bound, or an observable conclusion. Prior private roots and public capsules are read-only.

## 1. Inputs and scope

Use `Sigma=[0,6ell] x S3`, `ell=10^-6`, `0<=s<=T=13ell/8`, the retained positive sixteen-component product bundle

```
E=R_v + Sym^2(R dt+T*S3)_H + R_w + (R dt+T*S3)_D,
X_m(w;s)=sum_(j=0)^m ||partial_s^j w(s)||H^(m-j)(Sigma).
```

Spatial norms contain every ordered covariant derivative slot and full measure. `D` is the fixed, time-independent product connection. The discussion applies on one copy or the specified identical two-copy Hilbert union. All angular modes are retained.

Let `U^0` be the admitted compatible auxiliary seed and `U^(n+1)` the smooth frozen successors. The accepted admission provides

```
sup_s X14(U^n;s)<R/4,       R=10^-12,                    (1)
```

and common actual nonlinear initial jets `J0,...,J14`. Every successor satisfies the exact equation

```
U^(n+1)_ss = A(x,U^n)^i D_i U^(n+1)_s
           +B(x,U^n)^ij H_ij U^(n+1)
           +F(x,U^n,DU^n),                             (2)
F(x,v,Dv)=C(x,v)+E^mu(x,v)D_mu v
                         +G^mu nu(x,v)(D_mu v,D_nu v).
```

Here `D_mu` includes elapsed time, `H` is the full spatial product Hessian, and the normalized coefficient maps are the actual certified maps. In particular, their horizontal derivatives include the state-connection correction; their explicit base/time dependence and the physical reference-connection terms have not been discarded. The coefficient certificate is `Cstar=10^500000` through sixteen ordered horizontal/vertical derivatives; `C(x,0)=0`; tensor embedding has constant `E3=10^7` on all field arrays used here. Equation (2) is imposed only on successors; no equation is assigned to the auxiliary seed.

The lower-energy companion supplies a Cauchy estimate, uniform in time, for the actual pair:

```
sup_s [||U^n-U^m||H1+||U^n_s-U^m_s||L2] ->0
as n,m ->infinity.                                     (3)
```

Its evaluated factorial tail, not any constant below, is the quantitative convergence input. The qualitative Sobolev compactness constants used next enter neither the amplitude admission nor that tail.

Read-only analytic inputs are the accepted `ROTATION_MIXED_PRODUCTS.md`, `FROZEN_ENERGY_COMMUTATORS.md`, `CONDITIONAL_PICARD_ADMISSION.md`, and `NEXT_DIFFERENCE_AND_SMOOTHING_PLAN.md` in the preceding checkpoint's records, and the actual evolution source in the worktree at base `9b74d57b6caef415775abc069838a4bcea303831`. The new low-energy and top weak-continuity companions are joined separately.

## 2. Improve only the Cauchy pair by ordinary spatial compactness

On this compact smooth manifold with boundary, the ordinary bundle Sobolev embeddings

```
H14 -> H13 -> H1,       H13 -> H12 -> L2
```

have compact first arrows and continuous, injective second arrows. These are the usual Sobolev spaces of sections. No high mixed-boundary spectral domain or power of an auxiliary elliptic operator is identified with them.

For every `epsilon>0`, compactness gives finite constants `c_epsilon,d_epsilon` such that

```
||w||H13 <= epsilon ||w||H14+c_epsilon ||w||H1,
||z||H12 <= epsilon ||z||H13+d_epsilon ||z||L2.           (4)
```

For example, if the first inequality failed, normalizing a counterexample sequence in H13 would give bounded H14 norm, H1 norm tending to zero, and H13 norm one. Rellich compactness would produce an H13 limit of norm one whose H1 image is zero, a contradiction. The second inequality has the same proof. This argument is valid for arbitrary evolved Sobolev sections, irrespective of their higher normal traces. Standard local bundle trivializations establish the ordinary Rellich theorem; no numerical constant from that construction is imported into an evaluated estimate.

Apply (4) to differences of two iterates. Their H14 field and H13 velocity norms are each bounded by `R/2` by (1). Take the supremum in `s`, let `n,m` tend to infinity using (3), then let `epsilon` tend to zero. It follows that

```
U^n -> U_0       in C([0,T];H13),
U^n_s -> U_1    in C([0,T];H12).                        (5)
```

The limits are time-continuous because the smooth iterates are Cauchy in these Banach-valued continuous-function spaces. Passing the fundamental theorem of calculus in H12 gives `partial_s U_0=U_1`. This is the only interpolation/compactness improvement in the proof. Higher temporal jets will instead be recovered from (2).

In particular the state arguments converge uniformly. The limiting state remains strictly in the actual coefficient patch: the weaker bound (1) already gives

```
|U_0|star,infinity <= 2 E3 R/4=5*10^-6<10^-3.            (6)
```

## 3. Finite continuity lemma for differentiated rows

Write `U_l^n=partial_s^l U^n`. Fix `2<=j<=13`, set `r=j-2`, and aim at spatial regularity `m=13-j`. Assume inductively that

```
U_l^n -> U_l in C([0,T];H^(13-l)),   0<=l<=j-1.         (7)
```

Differentiate the actual equation (2) exactly `r` times in elapsed time and then up to `m` times covariantly in space, retaining the ordering. Every spatial derivative is a full array. The product connection is time-independent, so time commutes with the complete spatial Hessian; no two spatial derivatives need be commuted.

Replace the source's zero group by its exact identity

```
C(x,v)=Ctilde(x,v)[v],
Ctilde(x,v)=integral_0^1 C_U(x,theta v)dtheta.
```

As in the accepted coefficient certificate, radial state scaling commutes with horizontal transport; `Ctilde` inherits the needed coefficient derivatives. The five groups in (2) then start with at most two actual field factors. Each next ordered derivative has at most two coefficient-chain choices and one choice per existing field factor. Through `q` derivatives there are at most

```
5*(4*5*...*(q+3))=5*(q+3)!/3!
```

leaves, with repeated leaves retaining their multiplicity. For these rows, `q=r+k<=11`; hence there are at most `5*14!/6=72,648,576,000` leaves for any fixed ordered differentiated row. This finite count is used only to make the expansion exhaustive; no new numerical evolution constant is assigned.

Every leaf has the following properties:

* Total field derivative order is at most `r+k+2<=13`, counting elapsed and spatial derivatives equally.
* Its individual field factors have time indices at most `j-1`. Predecessor factors and successor factors keep their separate iterate labels.
* After a largest-order factor is selected for L2, every other field factor has total order at most six. Indeed two factors of order at least seven would already have total order fourteen.
* The L2 factor is controlled by one of the spaces in (7). Every other factor has three further spatial derivatives available in those same spaces and thus converges in L-infinity by the evaluated H3 embedding. Its tensor rank, including that reserve, is at most `2+6+3=11`.

These statements do not infer regularity of a high time jet from an unrelated low time jet. For a factor `D_sp^b U_l^n`, its L2 convergence follows from `l+b<=13` and (7) for that particular `l`. Its infinity convergence, when used, follows from `l+b+3<=9<=13` and that same jet's spatial Sobolev convergence.

The coefficient-array factors are bounded in operator-to-full-array norm by the accepted certificate. They converge uniformly after evaluation at `U^n`, because (5) gives uniform state convergence and one additional vertical derivative bounds the coefficient difference on the segment of states. With at most eleven row derivatives, the zero-group extraction and this difference derivative require at most thirteen derivatives of the original maps, within the reserve of sixteen. High-rank coefficient maps are not Sobolev-embedded.

For completeness, continuity of each leaf follows by a finite telescoping difference. When the changed factor is the designated L2 factor, use its strong L2 difference and bound all remaining factors in infinity norm. When a low factor changes, use its strong infinity difference, the bounded L2 norm of the designated factor, and the other bounded infinity norms. When the coefficient changes, use its uniform operator-array difference. Every convergence is uniform in time. These same estimates show time continuity of the limit product. Summing the finite leaves and the finitely many spatial derivative levels proves continuity of the entire row in `C([0,T];H^m)`.

## 4. Triangular recovery with the predecessor shift retained

Before taking the limit, the exact time-differentiated equation is

```
U_j^(n+1)
 = sum_(a=0)^r binom(r,a)
     [ partial_s^a A(x,U^n)^i D_i U_(r-a+1)^(n+1)
      +partial_s^a B(x,U^n)^ij H_ij U_(r-a)^(n+1) ]
   +partial_s^r F(x,U^n,DU^n),       r=j-2.             (8)
```

All coefficient time derivatives in (8) are derivatives of the composed coefficient along the predecessor, including explicit time dependence. The largest possible time indices on the right are:

| Term | Predecessor jet needed | Successor jet needed |
|---|---:|---:|
| Differentiated A coefficient | `r=j-2` | — |
| Mixed principal field | — | `r+1=j-1` |
| Differentiated B coefficient | `r=j-2` | — |
| Spatial-Hessian field | — | `r=j-2` |
| C, E, G coefficient derivatives | `r=j-2` | — |
| First-jet source factors | `r+1=j-1` | — |

Thus the row for `j` never calls on an unproved predecessor or successor jet of order `j`. In particular it is incorrect to suppress the predecessor index before the induction; (8) retains it.

Under (7), both `U_l^n` and `U_l^(n+1)` converge to the same `U_l` for every already established `l`. The finite continuity lemma therefore shows that the entire right side of (8) converges in `C H^(13-j)`. Define its limit as `U_j`. Since a one-step shift of a convergent tail has the same limit, this proves convergence of the entire sequence `U_j^n`, not only of a different subsequence for each row. Starting with (5), induction gives

```
U_j^n -> U_j in C([0,T];H^(13-j)),    j=0,...,13.        (9)
```

For `j<=12`, pass to the limit in

```
U_j^n(s)-U_j^n(0)=integral_0^s U_(j+1)^n(t)dt
```

in `H^(12-j)`. Consequently `U_j=partial_s^j U_0` as actual time derivatives in the indicated nested spaces and as spacetime distributions. Writing `U=U_0`, the conclusion is

```
sup_s X13(U^n-U;s) ->0,
partial_s^j U in C([0,T];H^(13-j)),   0<=j<=13.          (10)
```

There is no mixed-jet interpolation step in (8)--(10).

The same continuity argument is available for finite-data approximation: if a sequence of actual nonlinear solutions `V^nu` has the uniform top bound and a Cauchy field/velocity pair, replace both the predecessor and successor labels in (8) by `nu`. The induction is then unchanged. Common initial jets are not needed for convergence of the rows; their limiting values are obtained by evaluating the limiting triangular recursion at zero, starting from the limiting initial pair. This is a qualitative passage and supplies no lower-jet Lipschitz rate or smoothing error modulus.

## 5. Exact equation, initial values, boundary rows, and symmetry

The `j=2` row passes in `C H11`, so the limit satisfies the exact nonlinear equation

```
U_ss=A(x,U)^i D_i U_s+B(x,U)^ij H_ij U
      +C(x,U)+E^mu(x,U)D_mu U
                  +G^mu nu(x,U)(D_mu U,D_nu U)          (11)
```

in that topology. This includes the true quadratic source, the actual normalized principal coefficients, every direction-slot connection in `H`, and the full coefficient base dependence. Equation (11) follows from convergence of the individual displayed products; it is not a passage based only on weak convergence of the top derivatives. Since the normalized time coefficient remains positive by the admitted state margins, (11) is the certified normalized form of the original reduced equation.

Evaluation at zero in (9) gives

```
partial_s^j U(0)=J_j in H^(13-j),       0<=j<=13.        (12)
```

In particular the original initial Cauchy pair is unchanged. These are the actual nonlinear jets, because the limiting rows (8) are precisely their triangular defining recursion.

The field convergence `U^n->U` in `C H13` passes the four homogeneous Dirichlet shell rows, the twelve homogeneous Neumann shell rows, and the outer all-Dirichlet rows by the ordinary bounded trace maps. It also passes the corresponding rows of `U_s` using `C H12`, when those are desired. No normal derivative has been assigned a boundary condition in this proof. More generally, the lower topology alone permits differentiated Dirichlet traces through `j=12` and differentiated Neumann traces through `j=11`; no trace at a threshold where the ordinary trace map fails is asserted here. Any additional top-order weak trace passage belongs to the separate top-bound companion.

Each retained symmetry or identical-copy relation that holds for every iterate is a closed linear relation under these norms and hence survives (9). In particular, the prescribed Gamma invariance and the parallel elapsed-time/tensor reflection pass to the limit. This statement concerns the fixed symmetries of the actual iteration; it introduces neither a symmetry assumption on arbitrary rotation words nor a new physical matching assertion. All product arguments are intrinsic full-array arguments and remain valid on the full sphere and the retained two-copy union.

## 6. Negative-order handoff for the fourteenth jet

One more exact row is useful without claiming continuous full X14. Differentiate (2) twelve times in elapsed time and isolate only its two terms with no time derivative on the principal coefficient:

```
U_14^(n+1)
 = A(x,U^n)^i D_i U_13^(n+1)
  +B(x,U^n)^ij H_ij U_12^(n+1)+R_n.                    (13)
```

The remainder contains all `a>=1` principal coefficient terms with their binomial multiplicities and the complete `partial_s^12 F(x,U^n,DU^n)`. Its finite expansion has total field order at most fourteen and individual field order at most thirteen. Once a maximal factor is placed in L2, every other factor has order at most seven. Three spatial derivatives are available for those low factors in (9), so the same telescoping proof gives

```
R_n -> R_* in C([0,T];L2).                             (14)
```

There are at most `5*15!/6=1,089,728,640,000` leaves before grouping the two leading terms, and the original coefficient derivative requirement, including zero-state extraction and one difference derivative, is at most fourteen. Thus the accepted order-sixteen reserve suffices.

Use `H^-1=(H_0^1)^*` on the manifold with boundary. From (9), `U_13^(n+1)->U_13` in `C L2` and `U_12^(n+1)->U_12` in `C H1`. Their first derivative and Hessian, respectively, converge in `C H^-1`. From (5) and the coefficient chain rule, `A(x,U^n),B(x,U^n)` converge in `C W^(1,infinity)`. Such coefficients act continuously on H^-1 by duality with H_0^1; the relevant coefficient and coefficient-gradient terms are both bounded. Therefore (13) gives

```
U_14^n -> U_14 := A(x,U)^i D_i U_13
                     +B(x,U)^ij H_ij U_12+R_*
             in C([0,T];H^-1).                         (15)
```

Passing `partial_s U_13^n=U_14^n` distributionally identifies `U_14=partial_s^14 U`, and evaluation at zero gives `U_14(0)=J14` in H^-1. This is a genuine limit of the actual fourteenth jets, with the predecessor shift still present before passage. It does not use a fifteenth temporal row. Combining (15) with the admitted uniform L2 bound is a valid input for the companion's weak continuity and top-norm argument.

## 7. Delivered conclusion and open regularity interface

Given the evaluated low-energy Cauchy estimate and the admitted uniform X14 bound, the actual smooth iteration converges strongly in `C X13`; its fourteenth time jet also converges in `C H^-1`. The limit has the original nonlinear equation, Cauchy pair, lower initial jets, homogeneous mixed field rows, and retained symmetries. The proof uses ordinary spatial compactness for the pair, then the actual shifted equation for every higher time jet.

Weak retention of the top mixed bound, finite-data approximation, and uniqueness are joined in their companion. None of (4)--(15) establishes `C X14`, strong L2 continuity of the fourteenth jet, an X15 bound, or top-norm stability with respect to the data. Methods for those separate conclusions are not decided by this note. No qualitative compactness constant here changes the evaluated amplitude or the factorial tail.
