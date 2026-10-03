# Independent registered source-free ledger diagnosis

This route is staged before public freeze. It has not decoded a retained
physical array or evaluated a retained mode, profile or ledger. The permitted
input preflight reads only raw hashes, ZIP CRCs and NPY headers. The previous
metric calibration remains FAIL for every possible outcome of this diagnosis.

The one fixed interval is inclusive eta in[-2.5,-1.5], with129coarse and257fine
samples, both sources and K64/128/256. Each standalone route has one shared
900-second/262144-KiB envelope for all twelve cases and both80/100decimal-digit
passes. The independent route's envelope is separate from the primary route's;
no combined-runtime900-second claim is made.

## Independent mathematical route

Starting from the unprojected bare complex bilinears,

```
Rk=((2k²+3L²)Re(u)-k Im(w)-L Re(w))/(2k epsilon)
Pk=((2k²/3-L²)Re(u)-k Im(w)-L Re(w))/(2k epsilon).
```

With d=w_a/(i2k), c=u_a-d and E=exp(i2k(t-a)), the code retains Re(c).
The saved quadrature coefficient is A=weight*k/(4*pi_binary80²*epsilon_binary80).
It independently groups three Fourier moments:
Zr=sum Re(A*d*E), OZi=sum2k Im(A*d*E), K2Zr=sum k²Re(A*d*E),
with C0=sum A Re(c) and Ck2=sum A k²Re(c). Direct stresses and forcing are

```
R=2Ck2+3L²(C0+Zr)+L OZi
P=2Ck2/3-L²(C0+Zr)-4K2Zr/3+L OZi
F=6L³(C0+Zr)+4L K2Zr-2L²OZi.
```

The primitive is independently integrated from this F inventory as
G=3L²(C0+Zr)+L OZi. It does not use a measured final density or a
pressure defined through the Ward identity. The canonical I_ab sums its
per-node antiderivative using a fresh ordinary-mpmath exp(i2k(b-a)) from initial
retained modes. A second I_recurrence uses the independently accumulated
profile moments; |I_recurrence-I_ab|<=1e-12 is checked at EACH80/100precision.

The profile engine uses ordinary MPF nearest arithmetic, fixed64sample phase
refresh blocks and no integer quantization. Every precision pass starts from
the exact retained binary ratios. Only exact raw input arrays are reused;
no lower-precision coefficient, transcendental seed or recurrence is shared.
The refresh policy, interval, state and cutoffs never change after a result.

The mode endpoints are propagated with fresh direct MP phases. The actual
complex delta_u/delta_w, their signed density projection and triangle bound
are streamed per node. Their mathematical projection is a flow diagnostic,
not a bound on inherited initial-state or total physical error.

## Constants, inputs and Simpson arithmetic

The runtime must use Python3.12, NumPy2.2.6 and mpmath1.3.0 with native
little-endian16-byte x87 long double(nmant63,maxexp16384,minexp-16382) and
32-byte complex long double. Same-size quad-precision platforms are rejected.
Real/complex scalar guards reject binary64/Python floats BEFORE coercion.
Exact conversion uses a normalized integer mantissa and binary exponent;
roundtrips also preserve signed zero. Padding bytes remain checksum-protected
but are not arithmetic inputs.

```
epsilon=3777893186295716171/37778931862957161709568
pi_binary80=14488038916154245685/4611686018427387904.
```

These are the exact ratios of the producer's represented constants, not a
decimal epsilon substitution or exact mathematical pi. Analytic L=-1/eta uses
the exact represented dyadic eta at each MP precision. The original17members
of each capsule are read unchanged; their raw hashes, CRCs, dtype and shapes
are checked before payload decode. After GO, the saved source6jets,
forcing4jets and baseline contacts must be exactly zero across the selected
mask before the analytic free route is used.

S_ab reproduces the inherited long-double global-prefix Simpson function.
All three cutoff columns are summed together along axis0 in the original
arithmetic order. S_b-S_a rounds IN LONG DOUBLE before exact-ratio MP
conversion. S_direct_reset is a separate interval-slice diagnostic. Fine
double-step Simpson uses the whole global F[::2],2dt and corresponding global
endpoint prefixes, again subtracting in long double. DeltaR instead subtracts
separately exact-ratio converted stored endpoints in MP, as registered.

## Output, guards and interpretation

The CLI is independent/diagnostic_independent.py with --checkpoint-root,
--registration-sha256,--freeze-commit and --output-dir. A fresh output directory
outside the immutable checkpoint is required. The source and configuration
must be protected by independent/MANIFEST.json and FULL_REGISTRATION.json.
The explicit registration SHA, own/imported source hashes and shared integrity
helper SHA are checked before helper import; shared verify_frozen then checks
the freeze receipt, all registered files and safe paths before any np.load.
The same frozen guard is repeated after the calculation.

diagnostic.json contains four ordered records(source,setting), each with
levels80/100 and three canonical K rows. All scientific scalars/profiles are
decimal strings. It retains every sampled analytic/stored R/P/F and signed
profile differences; progress.jsonl and streamed per-node flow files are
separate. Every non-K row scalar/profile and every non-node_index streamed
decimal is in the <=1e-12 precision-gap universe. Serialization compares
exact Fraction(Decimal(string)) with the exact MPF dyadic, also at1e-12.
No scientific gate casts a value through binary64.

The signed decomposition D_S=D_cont+E_Q must close within1e-12 at EACHprecision.
R/P profiles, D_cont and signed E_flow use the original2e-7 consistency scale.
These scientific negative outcomes return exit0 and classification
CONSISTENCY_FAILURE. Precision/source/schema/resource/closure failures return
nonzero. All consistency checks must pass before gate-scale attribution can
be LEDGER_ERROR_DEMONSTRATED; otherwise the classification is
NO_GATE_SCALE_ATTRIBUTION or CONSISTENCY_FAILURE. old_metric_status is always
FAIL. The diagnostic is an empirical arithmetic and saved-flow check; it is
not an interval certificate or a new physical theory verification.

## Synthetic preparation evidence and resource plan

Ten synthetic math/conversion/Simpson/serialization checks pass normally and
with optimization, including nonzero Re(c), native LD/CD1+2^-60 witnesses,
extreme/subnormal values, direct raw bilinear profiles and independent MP
quadrature. Eleven additional guards check input schema, source/contact
nonzero mutations, exact-decimal parsing, precision gaps and fabricated
capsule/registration hashes. Eight separately authored exact symbolic
identities and three rejected mutations pass normally and with optimization.

SYNTHETIC_PIPELINE_TIMING_FINAL.json uses1024fabricated nodes and257samples,
including exact conversion, direct phases, streamed flow I/O, serialization and
80/100comparison. It projects273.5seconds for both precision passes over the
full10,534,912momentum/time pairs per pass and peaks37492KiB in the synthetic
process. This is feasibility evidence, not a retained-data resource proof.
Only one capsule is resident at a time; Fourier accumulators are3*257MPF
values plus small scalar sums, per-node flow is synchronous streamed I/O,
and only scalar/profile strings remain between capsules. Expected retained
working memory is below80MiB. The external replay driver records authoritative
whole-process wall/RSS/exit; the route's own receipt is supplemental.
