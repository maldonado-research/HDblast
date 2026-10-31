# A conditional curvature projection for proper-time geodesic observers

For the previously selected exact Einstein–scalar family, a new calculation
encloses the **complete second parameter derivative** of a specified magnetic
screen-curvature projection:

```
h_geo in [-0.404312161091, -0.014602053466].
```

The interval is strictly negative under the retained analytic and actual-field
premises. It includes the boosted time/radial tensor derivative, the full
scalar quadratic curvature, and the geodesic observer and screen tilts. This
is a new calculation for the defined proper-time observation. The older
coordinate-observer sign is not assumed to hold for this observer.

Read [the assembled theorem and proof](DERIVATION.md), the
[background and clock construction](background/DERIVATION.md), and the
[complete numerical enclosure](numeric/NORMAL.json).

## Observation and scope

Launch future unit-normal geodesics from the initial physical hypersurface,
retaining its radial and angular labels. Elapsed time along each trajectory
is proper time. Within each proper-time slice, choose the outward normal to
the initial-radial-label level sets. Normalize the angular connecting vectors
using their positive inverse Gram square root. These are Gaussian-label
screen axes; they are not assumed to be parallel transported.

The exact launch radius and clock are defined implicitly so the background
observation sphere is `(x,r)=(5/4,1/8)`. Their certified intervals enclose
those exact quantities; a rounded clock value is not substituted into the
definition. All numbers use the model normalization, with no SI calibration.

The scalar functional is

```
J_geo(epsilon)=3 mean_S3 <Hhat_geo(epsilon),Tplus-Tminus>_kappa,
J_geo(0)=J_geo'(0)=0,   J_geo''(0)=h_geo.
```

Its templates, orientation and round probability measure belong to the
initial label data. Transporting the metric, hypersurface and anchor data
together transports the geodesics, events and frames by ODE uniqueness,
giving a covariant anchored functional. Changing the observer, labels,
screen convention, orientation or nonlinear completion can change the
quantity. The result is not observer independent.

There exists an uncomputed `delta>0` such that

```
0<|epsilon|<delta  implies
J_geo(epsilon)<=-0.0036505133665 epsilon^2<0.
```

This follows from C2 regularity and the negative second coefficient. It gives
no usable numerical amplitude window or sign at a specified nonzero
amplitude. Negativity of the projection does not assert pointwise negativity
at every angle, a physical force, detector strain or radiation energy.

The observation is at positive inward radius, on the interior side of the
shell. It evaluates no distributional shell curvature and does not describe
an observer confined to the shell.

The analytic localized-constraint, mixed-boundary evolution, parameter-
regularity and actual-first-field premises remain explicit. The exact family
is the specific weighted-range construction identified in the
[preceding checkpoint](../HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/).
The result is not asserted for every family with its first perturbation.

The auxiliary background certificate covers proper time up to `1/5000`.
The actual family is used only on the smaller observer tube
`tau<=1/5500`, whose coordinate time is strictly below `3ell/2`.
The [domain join](background/FAMILY_DOMAIN_JOIN.md) establishes this
restriction and the exact zero initial observer variations. The larger
background overcover is not an extension of the exact-family theorem.

## Reproduce

From a complete repository checkout, with Python 3.12.14:

```sh
python -m pip install -r research/HDBLAST_PROPER_CLOCK_SCREEN_CURVATURE_20261031/requirements.txt
python -I -S -B research/HDBLAST_PROPER_CLOCK_SCREEN_CURVATURE_20261031/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-proper-clock-fresh
```

The output directory must not exist and must be outside the source checkout.
The runner authenticates the closed capsule and inherited inputs before
scientific imports and again afterward. It reproduces the new canonical
reports normally and under optimized Python, and checks the stated
corrupted-execution controls. Authentication of retained upstream results
is distinguished from replaying their proofs or numerical campaigns.

The full coefficient uses 64 frozen whole time cells for actual first-field
sources and one apex two-jet call. The separate TT calculation reconstructs
the reflected endpoint integrals, inhomogeneous shell history and complete
feedback bound from authenticated retained source ranges. Its source-oracle
campaign is not repeated. The independent numerical audit reconstructs
aggregation from the recorded actual-jet boxes, sharing the inherited oracle
premises rather than independently re-solving that first PDE.

The coefficient enclosure and its arithmetic checks use exact rational
outward intervals. All-point symbolic identities and a separate rational
curvature fixture support the geometric formula. An independent coordinate
fixture uses floating-point tolerance checks as a diagnostic; it supplies
neither the outward enclosure nor an analytic proof. Its dependencies and
tolerances are explicit. Internal AI-assisted reviews and reproducible tests
do not constitute external specialist or proof-assistant verification.

The Big Bang mechanism and cosmological interpretation remain unestablished.
No Zenodo version or journal submission accompanies this checkpoint. Earlier
capsules, manuscripts, candidates and websites are preserved.
