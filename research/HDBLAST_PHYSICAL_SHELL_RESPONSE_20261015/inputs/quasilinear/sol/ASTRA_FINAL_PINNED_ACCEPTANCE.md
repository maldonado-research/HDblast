# Final internal acceptance of the Astra evolution proof

Reviewer: independent Sol derivation, 15 October 2026.

Reviewed artifact:

QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md in the sibling astra directory.

SHA-256:

    d6f80e6158226b253c56f2f60c8d9c6cfa2ee2940e79bde42dd8fe8e48a5d6c4

**Verdict: accepted for the stated prescribed-interval reduced mixed-wave
theorem at finite order \(m=14\), including strong top Sobolev continuity
and continuous initial-data dependence.** This is a completed internal
proof review of the pinned document, not acceptance of a future outline
or of claims outside its stated scope.

## Frozen linear construction

The fixed elliptic operator \(A_*\) has the asserted positive closed form
and compact resolvent. Its natural Neumann row agrees with \(u_{N,z}=0\)
because the principal coefficient is scalar, the connection preserves the
orthogonal splitting, and \(c^{zA}=0\) at the shell. There is no spatial
interface of boundary types. The outer all-Dirichlet boundary is harmless.

The full spatial eigenprojection is used solely as an \(L^2\) contraction.
No assertion that it commutes with a coefficient or is uniformly bounded
on arbitrary high mixed-boundary Sobolev spaces enters the proof. The
identity

\[
 A_*w_N=P_NA(t)w_N+P_N(A_*-A(t))w_N
\]

is valid because \(A_*\) preserves its eigenfunction subspace. The
principal coefficient difference multiplies \(D^2w_N\), and short-time
\(C^0\) smallness absorbs it into the elliptic estimate. Derivatives of the
coefficient difference multiply only \(Dw_N\); boundedness suffices.
This proves the uniform discrete \(H^2\) estimate.

The time-regularity argument is noncircular. To estimate the differentiated
energy at time level \(p\), its apparently unbounded \(D^2v^{(p-1)}\)
commutator is bounded by the discrete elliptic estimate in terms of
\(v^{(p+1)}\) in \(L^2\) and \(v^{(p)}\) in \(H^1\). Those are exactly the
two controlled parts of the same \(p\)-th energy. Lower \(H^2\) terms
descend in time level. Zero initial data and time-flat residual source make
all Galerkin initial energies zero.

After the limit, the elliptic regularity induction is valid: for each
spatial improvement, increase time level successively; derivatives of the
spatial operator multiply only already improved lower time levels.
All time levels are available. The argument supplies smooth terminal data.
Consequently the repeated short-interval construction is a legitimate
**smooth frozen** restart, with uniform short lengths on a compact smooth
coefficient slab. It does not require an unproved finite-Sobolev restart
theorem at an arbitrary nonzero shell corner.

## Nonlinear estimates and prescribed lifespan

The full variable energy flux and its vanishing under the four/twelve
split check. Tangential/time differentiation preserves those traces.
Normal derivatives are recovered from the equation and are never assigned
spurious homogeneous boundary data.

For \(m=14\), a two-factor principal commutator has total differential
order at most fifteen. At least one factor has order at most seven; its
\(L^\infty\) norm is controlled by \(X_{10}\), with three remaining spatial
derivatives in dimension four. In differentiated \(F(x,V,DV)\), two
factors of total order at least eight would consume at least fourteen
of the thirteen available differentiations. The stated composition
estimate therefore has at most one high factor. The derivative count is
valid for mixed time/spatial jets and does not invoke false fixed-time
interpolation between independent temporal jets.

The pointwise source cost in normal recovery is retained. The identity

\[
 X_{m-1}(w;t)\le X_{m-1}(w;0)+\int_0^tX_m(w;s)\,ds
\]

is the ordinary Banach-space fundamental theorem of calculus at each
displayed Sobolev order. Common exact initial nonlinear jets bound its
initial terms. It turns the pointwise lower terms into the displayed
Volterra bound without an unjustified small factor. The whole-interval
exponential majorant and the initial-size condition then give the claimed
fixed finite \(T\), including after the reference contact time when that
time lies inside the chosen slab.

The protected precontact seed supplies the correct finite jets and open
shell wedge. Its global \(X_m\) bound uses the protected Cauchy solution
and the fixed cutoff derivatives, rather than ratios of cutoff
derivatives or instantaneous spatial seed norms. The low-energy Picard
difference source is bounded in \(L^2\) because the common principal
coefficients depend on \(V\), without dependence on \(DV\), and the
uniform high bound controls \(D^2V\) in \(L^\infty\). Factorial convergence
on the full fixed interval and uniqueness follow.

## Strong top continuity

The Bona–Smith comparison has the correct asymmetric placement: the
coarser smooth solution carries the one additional derivative, while the
difference is controlled one derivative below the top. The derivative-loss
term is precisely a lower difference norm times the reference solution's
\(X_{m+1}\) norm. The lower difference estimate closes with the common
uniform \(X_m\) bound. Its initial nonlinear jet difference is controlled
by the initial data difference in the corresponding lower norm.

At \(m+1=15\), low index eleven controls the smaller factor of a
sixteen-order two-factor commutator. This index is bounded by the base
\(X_{14}\) estimate. Thus constants in the additional high estimate are
independent of the smoothing scale, even though
\(X_{15}(u^h)=O(h^{-1})\) need not be small. This is the explicit
interpretation of the document's low-index-eleven clause.

Compactly supported smoothing preserves the common open shell collar and
has top approximation error \(o(1)\) and lower error \(o(h)\). Their product
with \(O(h^{-1})\) tends to zero. This makes the smooth solutions Cauchy in
the strong top norm and proves the stated time continuity, rather than
inferring it from weak compactness. The same comparison gives continuous
data dependence; top-norm Lipschitz or smooth dependence is not asserted.

## Scope accepted and exclusions

The accepted hypotheses include a common fixed shell support gap and a
common protected precontact interval, a smooth background on the prescribed
finite slab, the scalar principal symbol, dependence of principal
coefficients on the state without its derivatives, and all stated positive
and boundary margins. The registered annular initial family and its
smoothing approximations have these common supports.

No blocking analytic gap was found in the pinned reduced-evolution proof.
The document expressly does not claim general finite-Sobolev arbitrary-
corner restart existence. Its a priori continuation bound and its direct
whole-interval small-data existence are distinct and correctly qualified.

Initial harmonic completion, subsidiary constraints and boundary identities,
both-copy matching, and the parameter-derivative budget are separate
interfaces. This acceptance does not silently certify them or identify a
second coefficient from a different initial correction. No external
novelty, spectral stability, Big Bang, or publication conclusion is covered.
