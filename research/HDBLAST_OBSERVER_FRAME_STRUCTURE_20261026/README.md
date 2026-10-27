# Same-event observer-frame structure of the first Weyl perturbation

This research note extends the restricted scalar observer-tilt calculation of
the [screen magnetic curvature checkpoint](../HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/).
For the exact registered first metric, the three-screen magnetic Weyl
tensor's second parameter derivative has the same local components
in two consistently oriented frames that agree at the same background frame
and event, for **arbitrary first-order relative frame variations**. The result
is conditional on the same metric family and exact correlated first-sector
representation in [PROOF.md](PROOF.md); no field equation is used in this
geometrical lemma. It does not construct a second-order solution.

The proof converts the master first metric to a lapse-only representative by
an infinitesimal conformal rescaling. Constant-y slices have constant intrinsic
curvature and are umbilic. Gauss–Codazzi and five-dimensional Weyl trace
subtraction then yield a vanishing one-normal block and a tangential
Kulkarni–Nomizu form. A pure first-order Lie-gauge block does not contribute
because the background Weyl tensor vanishes. These structures imply cancellation
for all ten infinitesimal Lorentz generators; H1=0 alone would not suffice.

The earlier bounded coefficient remains conditional on its inherited inputs
and hypotheses. This note makes no new numerical evaluation of those inputs.
It does not compare different background observers, identify events across gauge
changes, construct a physical observer or clock, establish general second-order
gauge invariance, cover the entire five-dimensional Weyl tensor, or prove a
finite-amplitude spacetime. It gives no evidence of a higher-dimensional origin
of the Big Bang. A coupled bulk and shell boundary solution remains open.

## Reproduction

Python 3 and its standard library are sufficient. From this directory run:

```sh
python checks/check_gauss_codazzi.py
python checks/check_frame.py
python checks/frame_curvature_control.py
```

Each may also be run with `python -O`; the explicit checks remain active. The
first uses a complete 15-case algebraic Gauss–Codazzi basis, verifies 960 mixed
and 3,840 tangential Weyl components, and checks the normal-electric trace.
The second checks nine tracefree tangential tensors against ten Lorentz
generators, including a generic Weyl counterexample with H1=0 and nonzero
frame variation. The third is a separate implementation, retaining a nonzero
second-order Weyl tensor and using exact polynomial interpolation to check the
factor and full frame variation; it has a separate generic counterexample.
These tests check finite tensor algebra and signs; the geometric proof supplies
the all-point model-to-structure argument. They are internal AI-assisted
checks, not a proof assistant or external specialist review.

## Related literature and next step

[Mukohyama](https://arxiv.org/abs/hep-th/0202100), section 2.2, distinguishes
bulk and brane gauge transformations at linear order. It motivates careful
identification of the shell, observer and event; its setting does not prove the
present second-order statement. [Maartens](https://arxiv.org/abs/gr-qc/0304089)
explains why brane perturbations generally need bulk information. His
brane-projected electric Weyl tensor is not the present screen projection.
These are focused comparisons with established primary literature, not an
exhaustive novelty review.

A useful next step is to specify marked shell/clock data and event transport,
then work out their combined transformation with the curvature tensor while
solving the compatible second-order bulk and shell equations.
