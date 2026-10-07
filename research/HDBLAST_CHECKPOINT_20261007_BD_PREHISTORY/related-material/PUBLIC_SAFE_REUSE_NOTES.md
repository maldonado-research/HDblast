# Transferable validation methods for HDBLAST

7 October 2026. This note audits selected public research methods. It does not
merge their physical models, establish their hypotheses, or enclose HDBLAST's
actual incoming state. No physical callback, saved quantum array or trajectory
was evaluated. The independent checks use exact rational algebra only.

## Preserve correlated contacts

The UTOE [sourced-EOM audit](https://github.com/maldonado-research/Unified-Theory-of-Everything/blob/c284afaa99c9a1e1fbde92a5547892f0ff9a2ce4/checkpoint_2026_10_01/SOURCED_STERILE_REDUCTION.md)
gives a useful bookkeeping control. For symmetric invertible M and a finite
shift A, set `Y'=Y+AᵀM` and `C'=C-YA-AᵀYᵀ-AᵀMA`. Then
`C'+Y'M⁻¹Y'ᵀ=C+YM⁻¹Yᵀ` by direct multiplication. The full action also carries
kinetic/derivative contacts. Five independent rational matrix examples passed;
15 omitted-contact controls failed. The transferable requirement is complete
source/contact accounting, not physical unification or amplitude cancellation.

## Separate consumption, mapping, moments and continuum coverage

The DMDE [source-operator gate](https://github.com/maldonado-research/dmde-research/blob/34c6347fab1018f9e4bd707c5270270bbabe5ba6/provider/docs/DMDE_v0920_SOURCE_OPERATOR_CLOSURE_GATE.md)
inspects the actual production source and its native-to-canonical map. For
`p=c q`, c>0, a density transforms as `S_q=c S_p(cq)`. Thirty-six independent
rational polynomial-moment cases passed; omitting c failed all corresponding
controls. These prove identities, not solver consumption.

Its [October 2 preflight](https://github.com/maldonado-research/dmde-research/blob/34c6347fab1018f9e4bd707c5270270bbabe5ba6/research/rounds/2026-10-02-source-adapter-preflight/README.md)
separately reports correct same-node consumption, failed moment/thermodynamic
gates, and nonmonotonic quadrature errors as a nonzero support endpoint crosses
grid nodes. A three-point apparent pass fails a denser phase scan. These are
the source's reported outcomes; its complete numerical pipeline was not rerun
for this note. For HDBLAST, input consumption, arithmetic, state normalization,
full contacts and quadrature must each retain their own evidence.

A small exact counterexample explains why finite moments cannot identify an
arbitrary shape. On [0,1], take `f_1=1` and `f_2=1+P/10000`, where

`P=1-42x+420x²-1680x³+3150x⁴-2772x⁵+924x⁶`.

The coefficient bound gives `f_2>=1011/10000>0`. Moments 0 through 5 agree
exactly; moment 6 differs by `1/120120000`. This is a classical
orthogonal-polynomial construction, not an external novelty claim. The analogous
need for a continuous spectral bound cannot be met by probe agreement alone.

## Keep signed retarded weights

AntiMatter's [linear-response derivation](https://github.com/maldonado-research/AntiMatter/blob/69c01f5e3fb6567933383cbaca8b1760371c34d8/research/rounds/2026-10-01_transport_gate/producer/DERIVATION.md)
keeps the retarded operator in the source integral. An exact scalar control
with `y_next=y/2+b/2`, zero initial state and biases +1 then -1 gives final
response -1/4 despite zero unweighted signed area. Reversing their order gives
+1/4. This supports signed pressure/source-work accounting, while supplying
no autonomous reservoir or heating mechanism.

## What would justify a spectral envelope

If a complex error e(k) has independently enclosed endpoint errors sigma_a,
sigma_b and a justified continuous bound `|e''|<=M` on [a,b], then

`|e(k)| <= ((b-k)sigma_a+(k-a)sigma_b)/(b-a) + M(k-a)(b-k)/2`.

The interpolation remainder is at most `M(b-a)²/8`. Alternatively, a justified
`|e'|<=L` gives `|e(k)|<=min_i(sigma_i+L|k-k_i|)`. These classical bounds are
conditional. Neither L, M nor endpoint sigma follows from finite-node agreement.
Applying them to HDBLAST's c(k), A(k) requires the unchanged target state and
independently verified derivative/endpoint premises.

For completeness, the nonnegative Dirichlet Green kernel is
`G(k,s)=(min(k,s)-a)(b-max(k,s))/(b-a)`. The difference from the endpoint
linear interpolant is `-integral G e'' ds`, while
`integral G ds=(k-a)(b-k)/2`; taking the complex norm proves the bound.

`check_reusable_algebra.py` reproduces the generic matrix, Jacobian,
moment-nonidentification and signed-response controls using Python's standard
library. Normal and optimized Python produced identical scientific receipts.
These checks are internal computational verification, not external peer review
or an enclosure of actual HDBLAST state, residual, contacts or continuum stress.
