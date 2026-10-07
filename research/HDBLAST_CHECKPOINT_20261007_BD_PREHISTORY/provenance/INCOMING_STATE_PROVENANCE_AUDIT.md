# Actual incoming-state source audit

The upstream state prescription is explicit in the fixed first-order model.
The stored state at `a=-9/2` is an approximate response accumulated over
`[-5,-9/2]`, rather than an independently specified arbitrary state. This audit
used source/JSON text and streamed opaque hashes only: no physical source
callback, saved-array decode, trajectory, contact or momentum evaluation.

The actual error of the retained binary80 state remains **NOT_ENCLOSED**.
The full twelve-case pressure/contact certificate remains **UNRESOLVED**,
historical metric calibration **FAIL**, and a higher-dimensional Big Bang cause
**NOT_ESTABLISHED**.

## Exact prescription and normalization

Canonical producer:
`research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP/independent/forced_metric.py`.
It is byte-identical to the active-source package's text reference
`inputs/reference_code/forced_metric.py.txt`, SHA256
`4fa625f2b534d741e2a48d6f0bfe84a89ee17bf306805906f962a8552e005906`.
Line references below use that identical text reference.

* Lines 25–31 fix epsilon, initial/final times, observations, cutoffs, sources
  and coarse/fine grids. Lines 61–70 declare fixed mass squared 2, minimal
  coupling, H=1, incoming BD and no Wronskian projection.
* Lines 163–189 define `h=B` or `h=zB`, `z=t+4`, with
  `B=exp(1-1/(1-z^2))` inside `|z|<1` and zero outside, and
  `g=4L^2h-2Lh'-h''`, `L=-1/t`.
* Lines 253–260 define `v0=exp(-ikt)/sqrt(2k)`, `v0'=-ikv0`,
  `delta_v=v0*u`, `delta_v'=v0*(w-iku)` and the full metric operator direction.
  The frequency jets in lines 219–232 are subtraction inputs, **not** an
  adiabatic initial-state prescription.
* Line 368 sets `u=w=0` at `t=-6`. Lines 387–395 propagate the response.
  Lines 423–426 require exact zero response through `t=-5`.
* The update's real forcing contains epsilon (line 392). Thus stored `u,w`
  contain epsilon. Active-source `PROTOCOL.md` lines 94–95 define
  `U=u/epsilon`, `W=w/epsilon`.

For the exact prescribed linear target, therefore,

```
W(a,k) = - integral[-5,a] exp(2ik(a-s)) g(s) ds
U(a,k) = - integral[-5,a] Phi_k(a-s) g(s) ds
Phi_k(r) = (exp(2ikr)-1)/(2ik), Phi_0(r)=r
```

`U,W` are entire in `k` at zero. The singular homogeneous coordinate
`A=delta_W/(2ik)` must not be evaluated at zero. Real forcing and zero past
data give `Re U-Im W/(2k)=0` exactly for `k>0`; this is an analytic target
identity. It does not set the stored normalization defect to zero.

The represented epsilon is exactly
`3777893186295716171/37778931862957161709568`. The represented measure constant
Pi is `14488038916154245685/4611686018427387904`, with measure
`k^2 dk/(2 Pi^2)`. These pins appear in active-source `PROTOCOL.md` lines 48–53
and the input manifest. Substituting mathematical pi would change the inherited
numerical target.

For a stored node, the relevant differences are

```
c_stored = (Re u_1 - Im w_1/(2k))/epsilon
delta_W = w_1/epsilon - W_exact(a,k)
A_error = delta_W/(2ik)
```

`c_stored` is an exact rational function of the retained binary values once
decoded under a frozen reader; it needs no normalization projection. A bound
for `A_error` needs a rigorous enclosure of the incoming Duhamel moment.
Imaginary `delta_U` does not enter the inherited linear stress operators, but
must be retained when reporting mode error.

## Wronskian conventions

The producer's amplitude diagnostic (lines 397–400) is
`abs(2 Re u-Im w/k)/epsilon`; it is `2|c_stored|` in the exact algebra.
Its `1e-10` gate is a producer declaration, not a directed rounding enclosure.
The physical diagnostic (lines 305–308) uses

```
linear_W = delta_v*conj(v0') + v0*conj(delta_v')
           -delta_v'*conj(v0) -v0'*conj(delta_v)
```

and reports `max(abs(linear_W))/epsilon`. The canonical normalization used by
the later theorem is `i(conj(v)*v'-v*conj(v'))=1`; its first-order variation
is `2 epsilon c`. The producer's `linear_W` has the opposite conjugation
ordering and equals `2i epsilon c` in exact algebra. The absolute diagnostics
agree. Direct physical bilinears have additional native phase/amplitude
rounding, so those numerical maxima cannot certify an entire continuum state.

## All twelve cases and opaque identities

Active-source `PROTOCOL.md` lines 20–39 fixes two sources, two grids and
`K=64,128,256`, with anchors `u_1,w_1` at `-4.5`, midpoint `u_2,w_2` at `-4`
and endpoint `u_3,w_3` at `-3.5`. Producer lines 441–453 establish the archive
index-to-time mapping. The four capsules are the independent metric producer's
retained arrays; both diagnostic routes use these same inputs. They are not
four new physical states prepared at the diagnostic anchor.

| Source | Grid | K | Stored prefix nodes | Prehistory source cells |
|---|---|---:|---:|---:|
| positive_B | coarse | 64 | 2048 | 64 |
| positive_B | coarse | 128 | 4096 | 64 |
| positive_B | coarse | 256 | 8192 | 64 |
| positive_B | fine | 64 | 4096 | 128 |
| positive_B | fine | 128 | 8192 | 128 |
| positive_B | fine | 256 | 16384 | 128 |
| signed_uB | coarse | 64 | 2048 | 64 |
| signed_uB | coarse | 128 | 4096 | 64 |
| signed_uB | coarse | 256 | 8192 | 64 |
| signed_uB | fine | 64 | 4096 | 128 |
| signed_uB | fine | 128 | 8192 | 128 |
| signed_uB | fine | 256 | 16384 | 128 |

Coarse has time step `1/128`, momentum panels `1/2`; fine has `1/256` and
`1/4`. Momentum order is 16 and original source/time order 8. NumPy
`leggauss` computes their Gaussian constants in binary64 before promotion to
long double (producer lines 201–208 and 361–363).

`STATIC_PROVENANCE_AUDIT.json` records all twelve cases and full capsule,
original archive, selected `k`, weight, incoming `u_1,w_1`, and observation
member identities from the frozen manifest. Current whole-original and
capsule opaque SHA256 streams pass. Input manifest SHA256 is
`b70dea97bc7223a77f32c5e0964ee2a43812215d84bdac6393bedbb6601a291d`.
Historical public science lineage is
`71d00cc423e9049c8166b7ee7afbefbac28d8a18`, physical freeze
`19fde76912af6e2f30d6f55e066b27d88cc7e34b`, registration SHA256
`1c5bc9b21b34b1d036b7a7d12ccdb6766ac2ba7835ea04182dba868843866607`.

The `k.npy` byte hashes differ between sources within each grid. Binary80
padding could explain this, but no value equality is inferred. Preserve all
four capsules and verify any proposed equality through a registered exact
reader. No quantum payload values are supplied as text by the manifests;
individual retained input values cannot be recovered by metadata alone.

## Targets and a tractable next calculation

The existing anchored diagnostic treats retained binary values as exact and
holds them fixed. Its certified integration can therefore be conditional on
that exact discrete state; that does not assert that state matches the BD
prehistory target. A separate comparison can enclose the physical prescribed
linear target at each stored node. No continuum extension/interpolant for the
stored node values has been registered, so nodewise errors do not by themselves
yield a continuous stored-state envelope.

The analytic BD target has a continuous definition at every `k`. It is now
tractable to enclose that target before comparing any arrays. A narrow study
can retain the old nine exact momenta and prove its representation truncation
uniformly for `0<=k<=256`, while computing full arithmetic enclosures only at
those nine momenta. This computes the actual unique prescribed prehistory
state at nine probes, not the error of the original binary80 arrays and not
the full twelve-case certificate.

One independently proved flat-cap choice is `delta=1/128`, with
`d=2x-x^2`, `x=s+5`, `D=255/16384`, `lambda=128/639`. On the omitted left cap,
`|B'|<=2B/d^2`, `|B''|<=4B/d^4`. Monotonicity of
`d^-q exp(-1/d)` for `q=0,2,4`, together with `e>8/3`, gives a purely rational
majorant `B<=H=(3/8)^63`. Thus

```
G_B  = delta*H*(4 lambda^2 + 4 lambda/D^2 + 4/D^4)
G_zB = delta*H*(4 lambda^2 + 2 lambda + (4 lambda+4)/D^2 + 4/D^4)
```

bound the cap's `integral |g|`. The entire moment errors satisfy
`|Ucap|<=G/2`, `|Wcap|<=G`; consistently omitting that real source keeps
`c_cap=0`. Its `|A_cap|<=G/(2k)` is infrared integrable under the inherited
measure. On the entire subsequent interval, Pi>=3 and L<=2/7 give

```
Rcap <= G*(K^2/252+K/294)
Pcap <= G*(K^3/162+5K/1764)
|Icap| <= 2 Rcap
```

All six source/cutoff cap stress/work bounds are strictly below `1e-15` by
exact rational arithmetic in `EXACT_RATIONAL_FLAT_CAP_BOUND.json`. These are
analytic cap bounds, not evaluations or achieved total state accuracy. The
independent integrated-by-parts cap bound is tighter.

The complementary interior can use eleven exact geometric panels, degree112
Cauchy source jets, affine recentering onto subpanels of width at most `1/64`,
512-bit Arb coefficients/arithmetic and degree128 entire phase kernels. The
parent truncation remainder survives this subdivision. Fabricated whole-route
fixtures, controls, budgets, a complete registration, independent review and
public exact-byte readback must precede the 22 real source-jet callbacks.

The `2e-8` full-integral half-width remains an inherited proposed engineering
gate (next-study draft lines 23–25), not a current achieved result. Its full
acceptance still requires dense direct pressure, complete matched contacts and
source work, error propagation, exact original inputs if compared, momentum
and time integration, arithmetic and serialization. The narrow state study
must neither claim those terms vanish nor upgrade the metric failure.
