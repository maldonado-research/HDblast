# Evaluated background first-wave H14 evolution

This continuation closes the numerical first-wave norm left open in PR #135.
On the retained background slab `0<=t<=13/8000000`, with the specified
homogeneous mixed shell and outer wall, the exact linear first wave obeys

```
N14(t)=||V1(t)||H14+||V1,t(t)||H13
     <10^4000000 exp(10^152017).
```

The norms use the fixed product connection, full sphere measure and the full
Frobenius metric of the sixteen-field bundle. All coefficient, multiplier,
time-commutator, initial-time-jet and radial-recovery constants are evaluated.
The proof assumes the retained actual background and compatible linear problem;
it does not estimate a finite-amplitude nonlinear family.

Substituting this bound into the complete previously proved quadratic-source
map now gives a numerical whole-slice estimate throughout that slab:

```
sup_t ||R2[V1](t)||H12 <10^8400000 exp(10^152018).
```

This includes the principal variation `V1 D²V1` and its factor two. It is a
spatial H12 estimate at each time; the full mixed time-source/corner hierarchy
and the evolved second-wave estimate remain separate requirements.

[HIGH_ORDER.md](HIGH_ORDER.md) gives the complete argument.
[INITIAL_NORMAL.md](INITIAL_NORMAL.md) supplies the literal initial and radial
recurrences. [SYMMETRY.md](SYMMETRY.md) proves exact natural SO(4) commutation,
all tensor-slot comparisons and preservation of the selected first wave's
Casimir-eight subspace. The high-order theorem uses the generic rotation
identity and does not require an angular spectral truncation. The special
first-wave subspace is not asserted for the second wave or a nonlinear family.

[REPRODUCE.md](REPRODUCE.md) describes 221 exact checks, identical normal and
optimized replays, and seven corruption controls. Eight original public
snapshots are authenticated in [SOURCE_PINS.json](SOURCE_PINS.json). The finite
algebra supports the written analytic proof; it is not independent formal
verification of the retained existence and uniqueness premises.

[LITERATURE.md](LITERATURE.md) compares four inspected primary methods papers
with the actual model. [CONTINUATION.md](CONTINUATION.md) records the next gaps.
These deliberately large constants use the fixed model normalization, have no
SI blast-energy calibration, and do not establish a cosmological origin,
external novelty, external peer review or a finite nonlinear observable result.
