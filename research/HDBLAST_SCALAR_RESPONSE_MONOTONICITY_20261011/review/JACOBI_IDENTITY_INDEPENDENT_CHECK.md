# Independent Jacobi and endpoint sign check

This note uses the supplied equations and is independent of numerical evidence. No accepted research files were edited.

Let k=ρ⁻², c=η(0), h=∂cη, q=∂cp, v=∂c logρ, and b=∂ca, with the radial coordinate held fixed. Primes on radial fields denote radial derivatives; primes on U and σ denote scalar derivatives.

## Exact variation identity

The linearized equations are

    h'=q,
    q'=U''h−4aq−4pb,
    v'=b,
    b'=2kv−(2/3)pq.

The linearized constraint is

    2ab=−2kv+pq/6−U'h/6.

Consequently j=b+ph/3 obeys exactly

    j'+4aj=−2kv,
    (ρ⁴j)'=−2ρ²v.

For a regular center with ρ'(0)=1,

    ρ=t−U(c)t³/36+O(t⁵),
    η=c+U'(c)t²/10+O(t⁴),
    p=U'(c)t/5+O(t³),
    h=1+U''(c)t²/10+O(t⁴),
    v=−U'(c)t²/36+O(t⁴),
    b=−U'(c)t/18+O(t³),
    j=U'(c)t/90+O(t³).

Thus the center contribution lim(t→0)ρ⁴j is zero, and

    j(T)=−2ρ(T)⁻⁴ ∫₀ᵀρ(t)²v(t)dt.

When β=η(T,c), h(T)≠0, and aβ denotes the fixed-T endpoint derivative, the normalized endpoint term is

    J:=aβ+p/3=j(T)/h(T).

The whole variation (including v) must be divided by h(T) if it is renormalized to h(T)=1. The unnormalized integral is divided by h(T), not merely its left side relabeled.

## Endpoint derivative along the scalar junction

Let F(T,β)=p(T,β)+σ'(β)/2. Along F=0, σ'=−2p. Let R=a−σ/6. Writing c=c(T),

    βT=p+h cT,
    RT=a'+b cT−σ'βT/6
      =−k+J(βT−p).

This verifies the proposed sign and all factors. Here RT is the total derivative along the junction, not the fixed-β derivative.

## Stronger bulk consequences of strictly convex U

For all positivity and monotonicity statements below, assume U''>0 throughout an interval containing the scalar range and assume ρ>0 on (0,T]. This is the strict convexity supplied by the primary theorem's contraction conditions. No sign assumption on a is required; this note makes no positivity claim for merely weak convexity.

For U'(c)≠0 put ε=sign U'(c). The scalar divergence equation (ρ⁴p)'=ρ⁴U'(η), together with strict convexity, gives εp>0 and εU'(η)>0 for t>0: η starts moving in direction ε, and a first failure contradicts the positive integral formula for εp. Put r=−εb. Differentiating a'=−a²−p²/4−U/6 gives b'=−2ab−pq/2−U'h/6, so the Jacobi equations take the damping form

    h'=q,
    (ρ⁴q)'=ρ⁴(U''h+4|p|r),
    (ρ²r)'=ρ²(|p|q/2+|U'|h/6).

The center expansions give q=U''(c)t/5+O(t³)>0 and r=|U'(c)|t/18+O(t³)>0 for all sufficiently small t>0, while h=1+U''(c)t²/10+O(t⁴)>1. Fix any sufficiently small t₀>0. On [t₀,T], the equations for (h,ρ⁴q,ρ²r) form an ordinary linear system with continuous coefficients and nonnegative off-diagonal entries; all three components have strictly positive initial values. The system therefore preserves positivity. Combining this with the strict center start gives h≥1, q>0, and r>0 on the full segment. No initial-value assertion is imposed at the singular endpoint. The weighted center terms also vanish: ρ⁴q=O(t⁵) and ρ²r=O(t³).

Since v(0)=0 and v'=b=−εr, one has εv<0. The exact integral identity for j then gives εj>0. If U'(c)>0, these signs yield b<0, v<0, j>0 and

    0<J<p/3,   aβ<0,   pβ=q/h>0.

If U'(c)<0, the same damping argument gives b>0, v>0, j<0 and

    p/3<J<0,   aβ>0,   pβ>0.

In either nonconstant case, pJ>0. The map c↦η(T,c) is locally and on each common connected existence interval globally injective because h>0. This does not prove global surjectivity or existence at arbitrary endpoints.

If U'(c)=0, uniqueness of the regular initial-value problem gives the constant-scalar background p=0; then b=v=j=J=0. The remaining system is (ρ⁴q)'=ρ⁴U''h and h'=q. Its strict center start and the same regular-interval argument give h≥1 and q>0 for t>0.

A further useful identity is obtained from d=hp'−pq:

    d'+4ad=4p(pj+kh).

Its regular-center boundary term ρ⁴d vanishes. Therefore d has the sign of p on a nonconstant strictly convex-U trajectory, and

    A:=p'−p pβ=d/h has the sign of p,
    p'/p>pβ>0.

In particular p' has the sign of p. This also follows from the first-crossing relation p''=p(U''+4k+4p²/3) at p'=0.

## Conditional branch monotonicity and the necessary limitation

Define the scalar-junction endpoint derivative

    D:=Fβ=pβ+σ''/2.

For a simple junction graph β(T),

    βT−p=−(p'+σ''p/2)/D=−p−A/D,
    RT=−k−pJ−JA/D.

Under the stated strict convexity, JA>0 on nonconstant backgrounds. Therefore D>0 is sufficient for

    RT<−k<0.

The stronger bound RT<−k−pJ also follows. The condition σ''≥0 implies D>0 because pβ=q/h>0 under the stated strict convexity. For a constant-scalar trajectory, RT=−k wherever a differentiable scalar-junction branch is defined, regardless of D.

For fixed β, the bulk identity separately yields

    (∂T R)β=−k−pJ<0

on nonconstant strictly convex-U trajectories, independent of σ''.

Even strict convexity of U alone does NOT determine RT along the scalar junction. At any nonconstant strictly convex-U endpoint choose σ'=−2p and choose σ''/2=−pβ−ε for ε>0. Then D=−ε and

    RT=−k−pJ+JA/ε>0

for sufficiently small ε. A smooth local quadratic σ realizes these prescribed derivatives; its value can independently be chosen as σ=6a to place this endpoint at R=0. This is a local counterexample to unconditional residual monotonicity or a uniqueness theorem based on bulk convexity alone.

When D>0 throughout a connected scalar-junction branch, strict residual monotonicity gives at most one R=0 point on that branch. Existence, connectedness of the whole admissible branch, and coverage of its domain are separate requirements for a global existence-and-uniqueness claim.
