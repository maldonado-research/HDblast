# Selected 23-feature angular evaluation bound — unpublished

16 October 2026. This follow-on is outside the accepted frozen capsule. With
probability-normalized measure on the unit three-sphere, and the inherited
positive metric Frobenius norm, the explicit selected feature space satisfies

\[
 |U(n)|_F\leq \sqrt{321/19}\,\|U\|_{L^2(S^3)},\qquad
 \sqrt{321/19}=4.110320771\ldots .
\]

The constant is sharp for this 23-dimensional space. This is an angular
evaluation theorem, conditional on using precisely the selected space. It does
not assemble or certify the bulk, boundary, source, or initial projections.
The accepted full 287-dimensional certificate retains factor 17. A numerical
response certificate may use the smaller constant only after proving that its
state and the error being evaluated belong to this selected space. General
287-space errors still require their full-space bound.

## 1. Features, measure, and physical normalization

Write

\[
 S=n_1^2+n_2^2,\quad Y=n_1^2-n_2^2,\quad J=Y^2,
 \quad Z=S-\tfrac12,\quad F=J-\tfrac12S+\tfrac1{12}.
\]

The pointwise domain is exactly `0 <= S <= 1`, `0 <= J <= S^2`. The features
are the candidate's five scalar copies of `(1,Z,F)`, two mixed vector copies
of `(DZ,DF)`, and the angular trace-free tensor copy of
`(HZ,HF,T+,T-)`, where `Hf=(D^2 f)^TF`. Both TT directions are retained.

All scalar and tensor integrations use `integral(1)=1`. For a monomial,

\[
 \mathbb E[n_1^{2a_1}\cdots n_4^{2a_4}]
 =\frac{\prod_i(2a_i-1)!!}{\prod_{j=0}^{a_1+\cdots+a_4-1}(4+2j)},
\]

and a monomial odd in any coordinate integrates to zero. The exact norms are

| Copy | Unweighted squared norms |
| --- | --- |
| `1,Z,F` | `1, 1/12, 19/720` |
| `DZ,DF` | `2/3, 19/30` |
| `HZ,HF,T+,T-` | `20/9, 133/15, 1/6, 1/6` |

Each of these three lists is an orthogonal list in L2. In the scalar physical
slots `(k_yy,k_y0,k_00,angular trace,w)`, the Frobenius weights are
`(1,2,1,3,1)`. The vector slots `(k_ya,k_0a)` each have weight 2. The tensor
slot has weight 1. In particular the trace amplitude multiplies the unit sphere
metric, whose squared norm is 3; it does not multiply a unit tensor.

The checker constructs these physical scalar embeddings explicitly, including
both symmetric mixed entries and the tangent projector. It checks all their
weights and mutual orthogonality, the scalar/vector/tensor Gram matrices, and
the tangent, symmetric, trace-free identities of the tensor features. The
ambient representation uses the projector `P=I-n tensor n`; ambient component
contractions then equal intrinsic tangent Frobenius contractions.

For a field slot with constant physical weight `w`, the normalized feature is
`feature/sqrt(w*L2norm^2)`. Its pointwise squared physical norm is therefore
`|feature|^2/L2norm^2`. The weights cancel only after being included in both
the physical norm and the Gram matrix.

## 2. Orthogonal output blocks give the operator norm bound

The physical output fiber splits orthogonally into five scalar lines, two
mixed tangent-vector slots, and the angular trace-free symmetric tensor slot.
The corresponding coefficient spaces are orthogonal in L2. Thus the evaluation
operator is a direct sum of these blocks. For every point its squared operator
norm is bounded by

\[
 \max\{K_0(n),K_1(n),K_2(n)\},
\]

where each `K` is the trace of that block's positive evaluation kernel:

\[
\begin{aligned}
 K_0&=1+12Z^2+\tfrac{720}{19}F^2,\\
 K_1&=\tfrac32|DZ|^2+\tfrac{30}{19}|DF|^2,\\
 K_2&=\tfrac9{20}|HZ|^2+\tfrac{15}{133}|HF|^2
        +6|T_+|^2+6|T_-|^2.
\end{aligned}
\]

The scalar block maps to a single line, so its trace is its squared operator
norm. The vector and tensor blocks use the general fact that the largest
eigenvalue of a positive matrix is at most its trace. Repeated scalar/vector
slots do not add to the operator norm of a direct sum.

For comparison, the trace of the whole 23-feature evaluation kernel is
`5*K0+2*K1+K2`, with integral 23. It is nonconstant: its values at `e1` and
`e3` are respectively `59553/665` and `16353/665`. The homogeneous-kernel
argument that would identify this trace with the dimension does not apply.
The theorem here comes from explicit block estimates. The now-valid weaker
constant `sqrt(23)` also follows from the stronger theorem here; dimension
alone did not justify it.

## 3. Exact sphere polynomials

For a harmonic homogeneous ambient extension `p` of degree `l`, the checker
constructs

\[
 Dp=\nabla p-lpn,\qquad
 (D^2p)^{TF}=P(\nabla^2p)P+\tfrac{l(l-1)}3pP.
\]

Here `Z` extends as `n1^2+n2^2-|n|^2/2`, and `F` extends as
`Y^2-|n|^2*(n1^2+n2^2)/2+|n|^4/12`. Both are harmonic.
Direct contractions, identically on the sphere, give

\[
\begin{aligned}
 |DZ|^2&=4S-4S^2,\\
 |DF|^2&=-16J^2+24JS-8J-S^2+S,\\
 |HZ|^2&=\tfrac83(S^2-S+1),\\
 |HF|^2&=96J^2-240JS+128J+38S^2-6S+\tfrac23,\\
 |T_+|^2&=|T_-|^2=\tfrac16.
\end{aligned}
\]

The TT tensors use exactly the candidate's `C±=(C±*C)/2` and
`(T±)_ij=(C±)_ia jb n_a n_b`, including its curvature normalization.
The TT norms are constant individually, even though the full selected
evaluation kernel is not constant. These identities give

\[
\begin{aligned}
 K_0&=\tfrac3{19}(240J^2-240JS+40J+136S^2-96S+27),\\
 K_1&=\tfrac{48}{19}(-10J^2+15JS-5J-3S^2+3S),\\
 K_2&=\tfrac6{665}(1200J^2-3000JS+1600J+608S^2-208S+363).
\end{aligned}
\]

These are exact polynomial identities, not fits to angular samples. The checker
reduces the ambient norm polynomials by setting `n4^2=1-n1^2-n2^2-n3^2`,
then checks every `n3^2` coefficient cancels. Its certificate variable names are
`X=n1^2`, `Y=n2^2`, with `R=1-X-Y`; this certificate variable `Y` is different
from the quadrupole `Y=n1^2-n2^2` used above.

## 4. Global bounds and sharpness

The scalar estimate is elementary: `|Z|<=1/2` and
`-5/12<=F<=7/12`. The latter follows from `0<=J<=S^2`, since
`F>=1/12-S/2` and `F<=S^2-S/2+1/12`; the convex upper polynomial is largest
at an endpoint. Therefore

\[
 K_0\leq 1+3+\tfrac{720}{19}(\tfrac7{12})^2=\tfrac{321}{19}.
\]

For independent executable global bounds, the checker expands all three
polynomials in the triangle Bernstein basis

\[
 B^{d}_{ijk}(X,Y,R)=\frac{d!}{i!j!k!}X^iY^jR^k,
 \quad i+j+k=d,\quad X,Y,R\geq0,\quad X+Y+R=1.
\]

These nonnegative basis functions sum to 1. Thus a polynomial is bounded by
its largest Bernstein coefficient. Exact degree elevation produces:

| Block | Degree | Number of coefficients | Maximum coefficient |
| --- | ---: | ---: | ---: |
| scalar | 9 | 55 | `321/19` |
| vector | 6 | 28 | `320/19` |
| trace-free tensor | 5 | 21 | `7698/665` |

All are at most `321/19`. Negative Bernstein coefficients are harmless for
these upper bounds; no lower bound is inferred from them. The JSON report
contains every rational coefficient, and the checker independently expands the
Bernstein reconstruction back into ordinary monomials and checks equality.
No floating-point arithmetic, optimizer, angular mesh, or interval truncation
is used in this certificate.

At `n=e1`, `Z=1/2`, `F=7/12`, so `K0=321/19`. In any one of the scalar
slots the function

\[
 q(n)=1+6Z(n)+\tfrac{420}{19}F(n)
\]

satisfies `q(e1)=||q||_L2^2=321/19`; its evaluation ratio attains the claimed
constant. A physical scalar-slot weight multiplies its numerator and
denominator equally. Consequently the bound is sharp for the entire selected
23-feature space, without needing the exact vector or tensor block maxima.

## 5. Replay and wrong controls

Run the checker into a new output path:

```bash
python check_evaluation_bound.py --output NEW_NORMAL.json
python -O check_evaluation_bound.py --output NEW_OPTIMIZED.json
python replay_evaluation_bound.py --output NEW_REPLAY.json
```

The normal and optimized reports agree exactly. Nine targeted wrong controls
are rejected in each mode: wrong quartic L2 normalization, omitted mixed
Frobenius weight, omitted trace weight, ambient Hessian substituted for
intrinsic trace-free Hessian, one TT sector dropped, wrong sphere fourth
moment, dimension mistaken for a constant full kernel, unsafe constant 4,
and a corrupted Bernstein coefficient. The replay requires a nonzero exit,
the specific expected diagnostic, and no false certificate file.

This completes the evaluation estimate only. It changes no frozen files or
candidate/review files, and makes no numerical response claim.
