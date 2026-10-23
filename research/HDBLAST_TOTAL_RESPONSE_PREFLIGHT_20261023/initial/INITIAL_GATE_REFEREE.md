# Independent review of the bounded initial/free-wave gate

23 October 2026. Internal mathematical/code review; this is not external peer
review or a replay of the inherited source/PDE campaigns.

**Decision:** the displayed construction supports the claimed coarse C1
free-wave error certificate and the failed provisional accuracy gate. It does
not certify a total-response approximation meeting the design budget and does
not prove that such accuracy is mathematically impossible.

I read `initial/certify_initial_representation.py`,
`initial/INITIAL_GATE_PROOF.md`, the registered `initial_frame_data` and
`pack_initial` implementations, the retained initial profile evaluator, and
the physical/harmonic initial proof notes. I independently derived the
free-wave norm estimate and the selected five-mean initial formulas.

The cutoff polynomial has value zero/one at its two ends and derivatives
one through four zero at both ends. Multiplying the rational packet by this
cutoff, then using its zero/uncut pieces, therefore gives one C4 function on
the full footprint. The registered scalar displacement uses at most three
packet derivatives and its scaled radial derivative uses at most four.
The initial velocity uses at most three. Smooth rational background factors
with nonvanishing denominators therefore yield a C1 initial displacement
and continuous initial velocity. This is sufficient for a C1 d'Alembert
wave; a C1 initial velocity is unnecessary. Exact zero near r=0 makes odd
reflection of L0 and even reflection of the other four slots compatible.

The 96 cells of width 1/128 cover [1,7/4], including the cutoff breakpoints
as partition endpoints. The actual packet on the initial tail where the
candidate is zero is enclosed by the inherited analytic tail bound, so it
is retained in the error. Both actual and candidate data are identically
zero on [0,1]. This covers every initial characteristic dependency in D+.
Cell interval subtraction is conservative despite overestimating error;
using a shared exact nominal rational function makes the derivative
interpretation legitimate.

The mean Gram weights (1,2,1,3,1) and universal ell conversion of ordinary
t velocities agree with the retained adapter. The lapse-square value,
mixed y0 velocity, second scalar value and second scalar velocity remain
present. Evaluating `Fmom=3K` and `Gmom-Fmom_y=p Pi` through adapter
arguments `Gmom=p Pi,Fmom_y=0` preserves the formula because this adapter
uses only that latter difference. It is an algebraically reduced evaluation
of the registered data; the individual physical Gmom profile has not been
redefined. The all-point physical formulas, rather than isolated checks,
remain the premise establishing those reductions.

For uniform H bounds delta0,deltar,deltav, d'Alembert's formula gives
value error at most delta0+x*deltav, with x<=3/2. Each first derivative
error is at most deltar+deltav. Reflection preserves H. Consequently

    e_free <= max((delta0+(3/2)deltav)/3, deltar+deltav).

This independently matches the producer and ledger. The exact ledger cap
is approximately 4.560670583840034e33; the safe outward display is
4.561e33. The gate 1e28 is not met by this upper bound. The fixed-point
oracle interval radii are evidence that the particular whole-enclosure
verification route needs sharper input/representation bounds; they are
not lower bounds on the actual error of the fixed solution. An interval
box may contain many points that no actual correlated background attains.

The producer's mathematical imports were moved after its complete
authentication call during review. Normal and optimized numerical reports
are recorded separately by the initial task. This review did not execute a
new full numerical initial campaign. No issue requires changing any frozen
physical input or completed auxiliary source integral.

A minor display caution was conveyed: rounded nearest `float()` numbers
are explanatory displays, not necessarily outward upper bounds or interval
endpoints. The accompanying manuscript uses safe outward caps. The exact
rational ledger controls every acceptance decision.
