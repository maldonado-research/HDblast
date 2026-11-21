# Occupation carried between births in collisionless canonical momentum

At the heavy witness, at least **0.00478% of total prescribed births** would require an already occupied state on both time grids. Resolving the angular sweep one birth at a time detects 0.002899%. The additional shortfall arises because emissions at two times can populate the same canonical momentum. The absolute fraction remains small. The calculation keeps the fluid history fixed and omits scattering, finite birth width, magnetic deflection and loss of particles.

| Daughter mass | Action | Native lower bound (% of total births) | Half-step lower bound (%) | Half-step bound using folded shells alone (%) |
|---|---|---:|---:|---:|
| 4 × 10¹⁰ | Published witness | 1.280617 | 1.281610 | 1.100083 |
| 4 × 10¹⁰ | Thermal drag floor | 0.0895158 | 0.0895310 | 0.0733998 |
| 4 × 10¹² | Five times thermal drag floor | 0.00478911 | 0.00478990 | 0.00403271 |
| 4 × 10¹² | Thermal drag floor | 0.000973495 | 0.000973730 | 0.000848520 |

Values in the table are rounded. Both grids use the central 99.998% of the referenced final comoving rest yield for each mass. Omitting the extreme tails preserves each lower bound on total collisionless excess. The underlying fluid applies a drag rate. The bounds test the limiting collisionless interpretation between births and cannot determine how many particles the fluid would actually lose after its production and transport laws change.

## Where two births use the same state

Reflecting the force axis for one charge sign makes the two signs equivalent. Measure comoving momenta in units of the daughter mass µ. Define the sphere center B=b/µ, shell radius K=k/µ, axial momentum z and squared transverse momentum y. Both z and the transverse momentum are in the same normalized units. Divide all counts by µ³. The density of states, including spins, charges and three equal species, is then 3 dz dy. Isotropy gives requested births P dt dz/(2K) wherever |z−B|≤K, with P equal to the total requested rate divided by µ³. Each newly born sphere is the curve

$$
y=g_z(B)=K(B)^2-(z-B)^2.
$$

On a monotone stretch in B, define J=P/(2K|dB/dt|) and Z(B)=B−K dK/dB. The derivative of the squared transverse momentum is 2[z−Z(B)]. At z=Z(B*) the sphere is tangent to the envelope of successive birth spheres. Births on both sides of B* can enter a single narrow interval in y.

## State-count inequality used by the code

A time block comprises consecutive accepted steps. The signed force d=dB/dt is linear in time on each step, the shell radius K is linear in time, and the rate P is constant. Write c=dK/dt. Where d has fixed nonzero sign, the ordinary derivative is

$$
Z'=1-c^2/d^2+Kc\dot d/d^3.
$$

Take a block in which these ordinary derivatives all have the same sign, with m≤|Z′|≤M, m>0. Let V be the sum of absolute jumps of Z at internal step boundaries and W the total travel in B. Let Jmin be a lower bound on J throughout the block. Select a half-width ℓ. For each B* lying at least ℓ inside the block, take z=Z(B*). For |B−B*|≤ℓ the graph g lies in a y interval of length at most Mℓ²+4Vℓ. To see this, integrate 2[Z(B*)−Z(B)] from B* to B. The smooth part is of fixed sign and has magnitude at most Mℓ². Jumps influence either edge of the interval by at most 2Vℓ. The code additionally requires ℓ≤min(K−|Kc/d|) using conservative bounds on each accepted step, ensuring |z−B|≤K throughout the selected neighborhood.

The set of z values selected this way has length at least [m(W−2ℓ)−V]₊. Along the interval trimmed by ℓ at either end, the continuous pieces move in the same direction. An opposite jump can repeat at most its magnitude in previously covered z. Equivalently, subtract forward gaps from the end-to-end change and bound the backward jumps by V. For each covered z one B* suffices, even if several exist. Counting only births in its neighborhood supplies at least 2ℓ Jmin per unit z. The available capacity of its single y interval is at most 3(Mℓ²+4Vℓ). Thus the additional occupation in the block is bounded below by

$$
[m(W-2\ell)-V]_+\,[2\ell(J_{\min}-6V)-3M\ell^2]_+.
$$

The code optimizes this expression over ℓ, including the shell accessibility restriction. Blocks that fail the sign or derivative tests contribute zero from this construction. The code also bounds the ordinary angular deficit on each accepted step using its maximum K, an exact logarithmic integral in the large-field region and lower rectangles in the remaining angular integral.

For nonnegative densities fᵢ, [Σᵢfᵢ−1]₊ ≥ Σᵢ[fᵢ−1]₊ pointwise. A verified count for a block may therefore be added to counts for disjoint time blocks even when their momentum regions overlap. Within one block the code takes the larger of the folded-shell count and the sum inherited from its two disjoint children. It never adds both counts for the same births. At a leaf the alternative is the angular deficit after integrating directions independently. The maximum dyadic block size is 8192 accepted steps.

The geometric inequality and the angular enclosure hold for the stated interpolant in exact arithmetic. The reported numbers are floating-point evaluations. Both half-step witness values agree to the displayed digits when the evaluation arrays use extended precision. The native and half-step comparisons test temporal resolution, without bounding error in the underlying continuum dynamics. The result gives a lower bound on the unblocked source's excess demand. A full occupation evolution could detect additional overlap across other times and momentum regions.

## Controls and reproduction

For a constant unit-radius sphere translating at unit speed over a distance W≥2, let a=P/12≤1/2. At |cos θ|=u the two branches occupy axial intervals of length W separated by 2u. Each contributes occupation a/u. Integrating excess over their overlapping and singly covered parts gives the exact fractional excess a−a²/W. For W=16 and a=0.2 the exact fraction is 0.1975, the angular deficit is 0.1 and the folded-shell certificate returns 0.1901313973. A separate scan over strip widths reproduces the optimized certificate. The checks also cover a field reversal and the effect of increasing the jump bound. The results file includes every dyadic refinement and separate angular and fold counts.

First run the directional replay with its `--slices` option. Set `capsule` to the directory containing upstream `kinetics/SWEEP.json` and `slices` to the replay output directory. From this directory run

```sh
python -B check.py
python -B reproduce.py --slices "$slices" --reference "$capsule/kinetics/SWEEP.json" --output /tmp/memory.json
```

[Machine-readable bounds](MEMORY.json) · [Canonical counting code](memory_bound.py) · [Analytic checks](check.py) · [Pinned upstream study](https://github.com/maldonado-research/HDblast/tree/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_COHERENT_DECAY_20261120)
