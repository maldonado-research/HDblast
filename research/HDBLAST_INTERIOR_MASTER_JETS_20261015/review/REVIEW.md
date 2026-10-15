# Independent internal mathematical review

**Accepted within the stated conditional scope.** The new potential, first-master
C3, and complete longitudinal bulk-source envelope calculations are consistent.
This is an independent internal AI-assisted review, not external peer review or
a formal verification of the inherited analytic PDE hypotheses.

## Mathematical assessment

The actual-background identity is anchored to the retained branch, polynomial
potential, shell box, and wider first-exit rectangle. I checked the signs and
product rules in `hzz`, `Qzz`, and
`Wzz=(6 rho^4 a^2-2 rho^3 a_z)Q-4 rho^3 a Q_z+rho^2 Qzz`.
The independent checker authenticates six retained potential inputs using fixed
SHA256 pins. Its separate interval Taylor propagation of the background ODE,
without importing a scientific program or using the handwritten Wzz expression,
proves `|Wzz|<3e10` on all sixteen closed Cartesian subrectangles covering the
entire outer background box. An unsplit automatic-differentiation box is too
wide for that bound; subdivision closes it. These are interval enclosures,
not sampled trajectories. Twenty-seven exact point-jet witnesses supplement
that continuous calculation.

The reflected C2/C3 identities retain the correct shell-ray terms: `-B,+B`
at order two and `+B1,-B1` at order three. The retained flat initial collar
gives `f(0,0)=f_x(0,0)=0` and the zero initial boundary jets needed for the
separate retarded operator and boundary lift to join through C3. This is a
sufficient condition for this state; separate flatness is not necessary for
every possible compatible nonflat problem. The incoming identity for Hprime
retains its moving endpoint `-f(x,0)`. Its use with first interior jets closes
the boundary third derivative before the interior third-jet step, so there is
no circular dependence or differentiation of independent H interval choices.

The bump derivatives and envelopes are correct, including the smooth zero
extension: `Q<=1`, `|Q'|<=16`, `|Q''|<=256`, and
`|Q'''|<=46656/7`. The positive-coefficient inequalities after `s=u+4` and
`exp(2)>7` establish these bounds over the whole bump interval.

I independently reconstructed all 4,353 H cells by propagating boxes for the
first and second integrated primitives. This uses a different primitive
representation from the scientific checker. Seven boundary/interior scalar
fractions and all four C0-C3 order norms agree exactly with the final report.
The powers of `ell^-j` converting scaled to physical derivatives are retained.

The full longitudinal Einstein formulas, scalar Hessian, and on-shell scalar
cancellation are consistent with the factorial second-variation convention.
The reconstruction coefficients satisfy the correct identities for
`h=U'/p`, including `k=rho(5a/2-h)`, and `partial_y=-rho^-1 partial_z`.
I checked the sphere gradient/Hessian bounds and the slice connection terms.
Independent rational reconstruction reproduces all ten field-jet envelopes
and all four full-source envelopes exactly. The coarse continuous bounds are
`|S^E_AB|<1.2e21` for every individual component in the adapted background
orthonormal frame and `|S^phi|<1e13` on the full triangle times S3. The frame
qualification matters; arbitrary Lorentz boosts are not included.

The bounds concern the actual **first longitudinal representative**, matching
the retained coordinate shell-source representative. They do not supply a
first relative-harmonic gauge vector or the transformed harmonic source. The
full scalar source is a coefficient times Y squared on the actual first
solution. Generic jet checks alone would not establish that actual-solution
statement; the analytical first-field equations and authenticated norm
certificates are used separately.

## Replayed checks

`fresh-replay/REPLAY_ACCEPTANCE.json` records the final portable replay.
Normal and optimized scientific results agree, allowing only the declared
optimization field. Five jet/bulk/control receipts match their saved goldens
byte for byte. The potential controls also require complete report-byte
identity in the two modes.

- Potential: two positive runs, seven formal Laurent identities per run,
  and eight intended-failure subprocesses.
- Jets: two positive runs; the checker reports nine named groups per run,
  252 reflected derivative witnesses, 168 characteristic traces, 73 compatible
  joins, 48 shell traces, twelve incoming-integral witnesses and six lift-sign
  checks. Twenty-two mathematical-error subprocesses fail as intended.
- Bulk: two positive runs with 51 generic exact component/cancellation checks
  each, plus eight intended-failure subprocesses. Rejections are checked for
  their specific mathematical reason, including the positive-mode physical
  scale identity. An earlier unconditional scale-control rejection was fixed
  before this acceptance.
- Independent referee: two positive runs, each with the sixteen-subrectangle
  potential proof, 27 point-jet witnesses, 192 reflected C2/C3 witnesses,
  twelve C3 characteristic joins, two incompatible-corner counterexamples,
  150 manufactured H/Hprime identities, and the independent actual-bound
  reconstructions described above.

These counts describe different activities. They are not summed into a count
of mathematical theorems or an exhaustive proof over all smooth forcings.
The jet error controls include deliberately incorrect local mathematical
calculations; they are not all edits to the production source program.

## Limits and reproduction

The reviewed claims remain conditional on the retained actual regular branch
and smooth compatible first solution. This work supplies coefficient and
forcing envelopes. It does not provide pointwise interior tables, physical
source signs, a localized second initial constraint correction, a solved first
gauge bridge, compatible transformed boundary data for a common quantitative
second evolution, a computed U2, or a nonlinear remainder bound. The existing
conditional finite-angular construction remains a specification for that next
work. No cosmological observation or finite-amplitude nonlinear family follows.

Run the complete review from the capsule root, with a fresh external output
directory:

```sh
python -I -S -B review/run_independent_acceptance.py --root . --output-dir /tmp/hdblast-independent-review
```

The independent verifier and replay runner use only Python's standard library,
read scientific inputs relative to the capsule, and write only to the selected
output directory. `ACCEPTANCE.json` records the exact reviewed source hashes and
the corresponding final replay artifacts. Subsequent source changes require
review against those pins.
