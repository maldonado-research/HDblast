# Smooth candidate evaluated at an active shell point

This checkpoint encloses all 23 components of the tighter smooth candidate
`V=B0+B1`, including reflected bulk first feedback and its Robin lift, at
`(x,r)=(9/8,0)`. This is an active shell point: its characteristic cone meets
the nonzero initial data beyond the zero collar. The previous complete sample
was before shell contact.

The inherited conditional C1 error estimate gives a positive registered
second-order scalar mean and scaled-time derivative here. Using the exact
actual Robin law also establishes a positive inward radial derivative, whose
sign the general C1 derivative-error interval alone leaves undecided. The
candidate itself has a provably nonzero boundary defect. It therefore still
requires a boundary correction.

## One conditional formal shell sample

For the scalar-mean coefficient `w:0`, the detailed rational enclosures imply
the following deliberately outward-rounded strict bounds:

| Quantity | Lower bound | Upper bound |
| --- | ---: | ---: |
| `U_w:0` | `8.287e30` | `8.294e30` |
| `partial_x U_w:0` | `1.9744e32` | `1.9749e32` |
| `partial_r U_w:0`, using actual Robin law | `2.153e27` | `2.155e27` |
| Candidate actual-law scalar-mean boundary defect | `5.84e22` | `6.67e22` |

Here `x=t/ell`, `ell=1e-6`, and `r` is the registered inward scaled radial
coordinate. These are parameter-derivative coefficients in a fixed formal
decomposition and gauge. Their signs do not give the sign of the full scalar
field at every angle, a physical force, an invariant observable or evidence
for a higher-dimensional cause of the Big Bang. Large coefficients are not
physical signal or energy scales. The broad all-mode table does not replace
earlier sharper tensor-mode results.

The [detailed method](theory/NUMERIC_NOTE.md) uses whole-panel interval ranges
for two displacement integrals, exact smoothed P1 velocity primitives and
proved convolution remainders. The [independent coarse method](theory/COARSE_NOTE.md)
uses local cone norms and wider Lipschitz velocity-convolution bounds; it
imports neither the new feedback evaluator nor its output and confirms all
three primary signs. Complete exact endpoints are in
[reference/SHELL.json](reference/SHELL.json) and
[reference/COARSE.json](reference/COARSE.json).

## Boundary defect and coefficient uncertainty

The [reflection derivation](theory/CONTACT_FORMULAS.md) couples the physical
source first, then reflects its output. No commutation of the bulk matrices
with the parity projector is assumed. The candidate uses frozen rational
Robin coefficients `M0`; the actual formal boundary operator uses the
interval-enclosed `M`. With inward derivative convention, its full defect is

```
R = PN V_r - M PN V - G
  = -M0 PN B1 - (M-M0) PN V - G.
```

Both coefficient uncertainty and the full affine shell-source enclosure are
included. The [source check](theory/SOURCE_NOTE.md) derives uniform magnitude
bounds covering the sample from authenticated earlier reports; it does not
newly resolve a point value or sign of `G`. The positive scalar-mean defect
diagnoses the specified approximation's truncation. It does not imply an
inconsistency in the exact boundary equation or contradict the inherited error
budget. All five homogeneous Dirichlet candidate traces vanish exactly;
candidate TT zeros do not assert zero actual TT response.

The conditional formal actual radial interval follows from
`PN U_r = M PN U + G`, intersected with the transferred C1 derivative interval.
This use of a known boundary equation sharpens a first derivative estimate;
it supplies neither a second derivative error nor a new existence theorem.

## Reproduction

Python 3 standard library is sufficient. From the repository root:

```sh
C=research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030
python "$C/checks/check_contact_identities.py"
python "$C/checks/check_affine_shell_source.py" --repo-root . > /tmp/shell-source.json
python "$C/checks/check_shell_contact.py" --repo-root . > /tmp/shell-detailed.json
python "$C/checks/check_coarse_shell.py" --repo-root . > /tmp/shell-coarse.json
python "$C/checks/compare_shell_enclosures.py" --detailed /tmp/shell-detailed.json --coarse /tmp/shell-coarse.json --source /tmp/shell-source.json
```

Repeat under `-O`; outputs must agree byte for byte and numerical reports must
match the three reference files. The workflow runs both modes. The independent
polynomial audit has 34 named exact groups, including noncommuting coupling,
physical reflected integration, Robin-lift signs and actual-versus-frozen
coefficient corrections. The comparison checks all 23 rows across both
methods. Hashes authenticate inputs, while the underlying field certificates
and analytic transfer proof remain inherited conditional premises. The old
nominal producer authenticates its upstream roots before scientific imports
and reauthenticates them after evaluation. These internal AI-assisted checks
are not external specialist review or proof-assistant verification.

[Focused primary literature](RELATED.md) explains the classical delayed Robin
convolution and why recent scalar reconstruction results do not supply the
missing coupled-system theorem. This checkpoint makes no originality claim.

The next reconstruction can solve the frozen Robin trace Volterra equation
without truncating its boundary feedback, then include the affine shell source
and reassess the bulk and subsidiary residuals. That requires validated trace
integration and a fresh transfer ledger. The complete second-order boundary
solution, numerical C2 error, pointwise Einstein/constraint residual and exact
nonlinear family remain unfinished.
