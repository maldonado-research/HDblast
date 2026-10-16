# Independent review of the selected-space evaluation constant

16 October 2026. The claimed sharp constant **sqrt(321/19)** is accepted for
the selected 23-dimensional feature space, with sphere probability measure
and the positive metric/scalar Frobenius norm. This review uses an independent
global argument, rather than relying on the producer's triangle Bernstein
certificate. The result applies to a response or approximation error only
after that object has been shown to lie in the selected space. It does not
replace the full 287-space constant 17 for arbitrary errors.

## Geometry, weights and blocks

The selected scalar templates are 1,Z,F, with
Z=S−1/2, F=J−S/2+1/12, S=n1²+n2² and J=(n1²−n2²)². Their squared norms are
1,1/12,19/720. Gradient norms are 2/3,19/30, trace-free Hessian norms are
20/9,133/15, and the two TT norms are 1/6 each with zero cross inner product.
These follow from exact sphere moments and the sphere commutation identities;
they agree with the independently reviewed quadratic-feature construction.

The physical output splits orthogonally into five scalar lines, two tangent
vector slots and one trace-free angular tensor slot. Scalar line weights are
1,2,1,3,1, because k_y0 occurs twice and the angular metric has squared norm
three; both mixed vector weights are two. For each block the same physical
weight multiplies both its pointwise norm and its L2 Gram matrix, and therefore
cancels in the normalized evaluation kernel. Repeated blocks contribute by
the maximum operator norm, not by adding their kernel traces. The block traces
are precisely

    K0 = 1+12Z²+(720/19)F²,
    K1 = (3/2)|DZ|²+(30/19)|DF|²,
    K2 = (9/20)|HZ|²+(15/133)|HF|²+6|T+|²+6|T−|².

The scalar block has rank-one output, so its trace is its squared evaluation
norm. For the other blocks the trace is a valid upper bound for the largest
kernel eigenvalue. The theorem does not assert that those trace maxima equal
the vector/tensor operator-norm maxima.

## Independent invariant contractions

Let Q=diag(1,−1,0,0), A=Q², P=I−n⊗n and v=PQn. Then Q³=Q,
tr Q=0, tr A=2, nᵀQn=Y and nᵀAn=S. In particular

    |DS|²=4S(1−S),  |DY|²=4(S−J),  DS·DY=4Y(1−S),
    DF=2Y DY−DS/2.

This gives the producer's gradient contractions without expanding coordinate
polynomials. For Hessians, HZ is the trace-free part of 2PAP and HF is the
trace-free part of

    B=8v⊗v+4Y PQP−PAP.

The omitted multiples of P have no trace-free part. The necessary contractions
are

    |v|²=S−J,                 tr(PQP)=−Y,       tr(PAP)=2−S,
    tr((PQP)²)=2−2S+J,       tr((PAP)²)=2−2S+S²,
    vᵀ(PQP)v=Y(1−2S+J),    vᵀ(PAP)v=S+J(S−2),
    tr(PQP PAP)=Y(S−2),      tr B=9S−12J−2.

Consequently |HF|²=tr(B²)−(tr B)²/3. The independent checker expands these
invariants and obtains, identically,

    |DZ|²=4S−4S²,
    |DF|²=−16J²+24JS−8J−S²+S,
    |HZ|²=(8/3)(S²−S+1),
    |HF|²=96J²−240JS+128J+38S²−6S+2/3.

For either chiral Weyl feature, |T±(n)|²=|C±|²/16=1/6 at every n.
One way to see the pointwise identity is to send n to a coordinate unit
vector by an orientation-preserving isometry. The electric restriction of a
self-dual or anti-self-dual Weyl tensor is a trace-free symmetric 3×3 matrix;
the chirality and curvature symmetries reconstruct its full tensor with
squared norm 16 times the electric squared norm. Both norms are invariant
under the isometry. With the candidate normalization |C±|²=8/3, the constant
is 1/6. The prior exact feature check also verifies this by polynomial
contraction. The scalar/tensor weight conventions are therefore consistent.

## A separate global proof

The physical parameter domain is 0≤S≤1, 0≤J≤S². These inequalities are
also sufficient, by choosing n1²=(S+sqrt(J))/2 and n2²=(S−sqrt(J))/2.

First |Z|≤1/2 and −5/12≤F≤7/12. Hence

    K0≤1+3+(720/19)(7/12)²=321/19.

Both maxima used in this bound occur simultaneously at n=e1.

For the vector block, factor 48/19 and write

    K1=(48/19)[−10J²+(15S−5)J+3S(1−S)].

When S≤1/3 the bracket decreases with J≥0 and is at most 2/3.
When S≥1/3 its vertex is J=(3S−1)/4, which lies in [0,S²]:
4S²−3S+1 is positive because its discriminant is negative. At the vertex
the bracket is (105S²−30S+25)/40, a convex quadratic in S. Its endpoint
maximum on [1/3,1] is 5/2. Thus the exact trace maximum is

    K1≤120/19 <321/19,

attained at S=1,J=1/2. This is substantially below the producer's valid
but more conservative triangle-Bernstein upper bound 320/19.

For the tensor block, factor 6/665. Its polynomial is convex in J, so
only J=0 and J=S² need consideration. At J=0 the remaining polynomial is
608S²−208S+363, whose maximum is 763 at S=1. At J=S² it is

    1200S⁴−3000S³+2208S²−208S+363.

Its ordinary degree-four Bernstein coefficients on [0,1] are exactly
(363,311,627,561,563). All are below 763. The checker reconstructs the
polynomial from these five coefficients, so this is an exact identity and
global bound. Therefore the exact tensor **trace** maximum is

    K2≤4578/665 <321/19,

attained at S=1,J=0. This again proves more than the producer's conservative
certificate needs, without computing the tensor kernel's largest eigenvalue.

The scalar block consequently controls the full evaluation map. The scalar
reproducing vector

    q=1+6Z+(420/19)F

satisfies q(e1)=||q||²=321/19. In any scalar physical slot its weight cancels
from the ratio. Thus the full selected-space constant sqrt(321/19) is sharp.

## Executable audit and limits

`check_independent_evaluation.py` uses only Fraction arithmetic. It derives
the contractions above, verifies the independent endpoint/vertex and ordinary
Bernstein argument, and compares the resulting invariant polynomials with
the authenticated producer report after substituting S=X+Y,J=(X−Y)².
It does not import the producer checker or use its global Bernstein bound.
The producer's triangle reconstruction was also read and found valid: its
monomial-to-Bernstein conversion is the falling-factorial identity, and its
back-expansion independently verifies that conversion.

Normal and optimized replay are identical. Five specific wrong-formula
controls reject missing Hessian trace subtraction, wrong quartic Gram entry,
one omitted TT contribution, the unsafe bound four, and a corrupted ordinary
Bernstein coefficient. Counts are software checks, not mathematical novelty.
Closure, actual source assembly and a validated numerical response are
separate obligations; this theorem alone proves none of them.
