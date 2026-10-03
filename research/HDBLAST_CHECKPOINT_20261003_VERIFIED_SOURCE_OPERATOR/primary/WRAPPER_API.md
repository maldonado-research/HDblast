# Reviewed primary producer interface

`primary_route.py` is intended to be copied to `primary/route.py` with sibling
`primary/verified_moments.py` in the future checkpoint. Import it as a package
module or include its parent directory in sys.path for explicit path loading.

Production interface: `run(auth) -> (payload, method_evidence)`.
The root driver must finish authentication of the exact reviewed public freeze,
all sources/schema/registrations and runtime before this call. The auth callable
receives only the128 `real_source_taylor` events, one per source/panel, before
each actual source evaluation. There is no additional producer event. Root owns
the capability's complete-state check and source/center membership checks.

Fabricated interface: `run_fabricated() -> (payload, method_evidence)`.
This passes an explicit manufactured coefficient provider and its own dummy
authorization callable. That dummy rejects any real-source callback event.
It does not call root's real authorization object or create real-source events.
Root still authenticates the complete fixture/sourcefile contract before import.

Payload shape exactly matches protocol schema1 with four arrays:
1152 panel moment rows,18 whole moment rows,128 panel Lg work rows and2 whole
Lg work rows. All rational endpoint/interval/momentum/radius strings are
canonical reduced n/d, including integers and zero as denominator1.
Complex moment radius is the exact sum of real and imaginary outward endpoint
half-widths. The real moments M0 and Lg always have exact imaginary zero.
Mexp/Mu at k=0 are projected to exact imaginary zero using the registered
real-source target identity. This does not project generic incoming complex
mode responses from the reusable engine.

Method evidence is a separate artifact. It records the exact source and kernel
analytic models, source/phase positive disk contributions, coefficient/finite
polynomial arithmetic baseline enclosure widths, and every complete output
radius. The baseline with source/local model tails suppressed is not an
additive partition of complete interval width or a bound on every full-run
rounding operation after adding those tails. Complete radii determine the
registered gate. No signed residual subtraction isolates a favorable error.

Final wrapper SHA256:
334a1cdee4992493231c4646ea2dcf242196ca67542560054a3ddd45424b456f.
Reusable engine SHA256:
de995d2edeca3944400fd90b892521f05c941b6d8a1f75c2a820740b82d995c0.

Final full-entry manufactured fixtures pass normal and optimized Python with
identical numerical payload and evidence bytes. Payload SHA256:
186055875b3504d8513217fdb6662d68d2c4e03cbb87f0ed2b810bd02f37b22c.
Evidence SHA256:
266dcd893c93af0122afc997a3d13a92b2e4cdad832ee92a166f63e64ac9335f.
The standalone protocol validator passes. `run(None)` rejects before source
evaluation and the fabricated root alias exactly reproduces both fixtures.

No physical source callback or checkpoint array decoder has run in preparing
this producer. The complete bootstrap, output envelope, external resource
receipt, independent route, prospective public freeze/readback and scientific
publication remain root's responsibility.
