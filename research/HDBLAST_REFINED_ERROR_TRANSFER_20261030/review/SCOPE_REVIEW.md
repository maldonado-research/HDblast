# Independent publication scope review

30 October 2026. This review checks the refined error-transfer checkpoint for
unchanged `W1_quadratic`, its mathematical input roles, displayed numerical
claims and limitations. It is an internal mathematical and implementation
review, not external specialist endorsement, a novelty determination or a
formal proof-assistant result.

## Assessment

The staged derivations and authenticated mathematical reports pass this
review. The new interval subdivisions tighten upper certificates for the
same functions; they do not change the candidate, target, original sources,
initial data or true error. The final README wording review also passes; its reviewed bytes are recorded
below.

The following advertised strict bounds are supported by the exact reports:

| Quantity | Strict bound |
| --- | ---: |
| Uniform initial displacement-curvature error H cap | <6.495e26 |
| Each evolved scaled coordinate second-error H cap | <6.643e26 |
| Global initial rr H-sup ratio | <2e-7 |
| Global evolved rr H-sup ratio | <2.045e-7 |
| Global maximum-second-derivative H-sup seminorm ratio | <2.045e-7 |
| Initial local rr H-norm relative error at (0,65/64) | >0.999 |
| Conditional trace-reversed Einstein second-parameter residual | <1.476e35 |
| Conditional scalar second-parameter residual | <2.118e30 |

The final two rows require additional field identities specified below.
They are large absolute comparison-norm bounds and establish no physical
smallness. They have a different normalization from the coordinate ratios.

## Initial enclosure and error transfer

The explicit nominal third derivative is enclosed using ordinary scaled
`6 ell^3 c[3]` jets of the registered nominal expressions. This does not
assume an actual sixth derivative or differentiate an unknown actual Taylor
remainder. Complete nominal coverage, forced branches and authenticated
parent minima preserve the same nominal derivative. The paired calculation
keeps all three telescoping actual-minus-nominal channels; child minima use
only fresh child caps and authenticated old seed caps. There is no hidden
unverified internal-parent bound.

The unchanged isolated curvature-jump allowance and convolution first
absolute moment give component caps

```
s_i=(35 delta/128) M3_i+j_i/2.
```

All 23 positive Gram weights are retained. Independently recombining the
paired channel suprema and the packed smoothing vector reproduces the
combined initial error cap exactly. The complete original velocity,
absolute candidate, source, C1 and reduced-defect ledgers remain fixed.
The smaller initial curvature-error cap is not substituted for an absolute
candidate curvature bound in the differentiated forcing calculation.

The independent arithmetic review recomputes the new differentiated initial
majorant `J_new` from both branches of its maximum and verifies all three
new C2 caps by the exact identity

```
new_cap_ab=PR103_cap_ab+(J_new-J_old)/(1-q),
q=3409/8000000, ab in {xx,xr,rr}.
```

This avoids subtraction of rounded displayed bounds. The separate derivative
review recomputes the 319 new derivative leaves and full covers; this scope
review does not claim a further independent replay of those expensive leaf
calculations.

## Global normalization and its local limitation

The denominator is supplied by actual target initial curvature, not by the
old affine candidate's curvature error or a nominal value. Both authenticated
actual endpoint methods give `u0_rr,w0(7/4)<-3.249e33`, and `Gram(w:0)=1`.
The retained classical target initial trace and physical-side C2 regularity
imply `U_rr(0,7/4)=u0_rr(7/4)`. That corner belongs to the declared closed
`D+`; continuity supplies the same supremum lower bound if the interior is
used instead. The relevant global target scales are therefore finite and
strictly positive.

The rr H-sup quotient and the specifically defined maximum of the three
second-derivative H suprema are valid global coordinate ratios. The latter
is a seminorm, not a full C2 norm. A lower scale for rr does not separately
normalize xx or xr. These bounds do not give error divided by the target at
every point, coordinate-independent accuracy or a physical reference scale.

The local calibration makes this limitation quantitative for this very
candidate. At `(0,65/64)`, all 23 authenticated actual component intervals
provide an upper squared target H norm `A2`, and the existing strict local
H-error lower bound `Lerr` satisfies the exact rational inequality

```
Lerr^2>(999/1000)^2 A2.
```

The positive scalar-mean actual interval makes the local denominator
nonzero. Thus the same local H-norm relative rr error exceeds 99.9%, even
while its ratio of global suprema is small. This is a clipped-tail numerical
approximation mismatch in a tiny local target relative to the outer-endpoint
scale. It is not a physical model obstruction or an observed signal.

## Conditional natural-equation bridge

The natural equations are trace-reversed Einstein–scalar and scalar
`E=Ric-dphi*dphi-(2/3)Vg`, `S=Box phi-Vprime`. They are not the Einstein
tensor G itself. With the registered `lambda^2/2` perturbation convention,
the reported quantities are second parameter derivatives; the coefficient
of `lambda^2` is one half of each quantity.

Under the actual first harmonic identity on a neighbourhood, the retained
20-row source identity is `QE-QEred=sym(nabla QC)`, where sym includes 1/2.
The additional formal-target completion must hold as the C1 covector
identity `DC[U]+QC=0` over the physical domain. Pointwise interval values or
C2 closeness alone do not establish it. With `e=W1q-U`, exact shared-source
cancellation gives

```
E2[W1q]=R_metric[W1q]/(2ell^2rho^2)+sym(nabla DC[e]),
S2[W1q]=-R_scalar[W1q]/(ell^2rho^2).
```

The metric/scalar normalization factors, one-half symmetrization and error
orientation are consistent. The primitive scalar assembly retains its
natural harmonic coupling; it is not added a second time. The retained
reduced-defect cap has the opposite defect sign, which does not affect its
norm. No filtered-candidate jet, source averaging discrepancy or `2F_H`
charge enters the new route.

The harmonic derivative estimate is applied to error jets, with the safe
value bound `3E1` from the X norm. It requires no stronger zero-wedge value
estimate. The positive covector output Gram carries no mixed-tensor factor
two. The physical conversion retains the radial rho product derivative,
angular derivatives and connection charges. Independently rebuilding those
charges and expanding the full numerator reproduces both natural residual
caps exactly.

These remain conditional bounds in a fixed positive background-orthonormal
frame angular L2 comparison norm. The explicit CLI assumption flags record
the theorem's antecedents; they do not prove complete physical initial,
shell or subsidiary constraints. No Lorentz-invariant scalar magnitude,
small physical residual, finite-amplitude nonlinear spacetime or cosmological
prediction follows.

## Source authentication and review evidence

The six principal report SHA256 values checked in this review are:

```
M3: b1cbb3d50b49da04313feba77fab2eec1e6e4f8b7778d36eec89a9ce0ea0223d
PAIRED: 458f1fcddd06d14ac02bb69e609e37d431fc0e2a3c2232179fb89a9c51b3ad1c
COMBINED: 3e4fb0174df7e88c566346a1b68aedbe4e731fbd062de11eafa28a2cae0a4989
C2: 85bf4133d931d93401b339dc508553ae488bf91f43f75bd7e736432cc856a923
NORMALIZATION: e63711829cad2e3dea6ac7e133f754d0e25d04fc40c460c681c7ed1ef3a67a5a
NATURAL: e66489946ad67db59ecb8eab635a49a4de199e1b1a49de89e1079429075eed21
```

The old actual endpoint, independent endpoint, join, C2 and background-time
reports were authenticated independently. The private compact arithmetic
checker imports no scientific producer code. Its 32 explicit report/arithmetic
checks pass in isolated normal and optimized Python with byte-identical
receipt SHA256 `ba1245b10349126346c550284750cf024b99564219588af0a79ab0249ee00767`.
These 32 checks are not presented as32 new mathematical theorems. The source
closure, derivative-leaf replay and broader negative controls are documented
by the dedicated producer and audit checkers.

No new model assumptions beyond the stated conditional antecedents are
silently promoted to proved facts. Frozen manuscripts and prior candidates
remain separate. No external journal decision, originality, physical
accuracy, nonlinear realization or Big Bang-origin claim is inferred.

## Final public wording and packaging review

The final capsule README and the new top-level README paragraph agree with
the authenticated reports and the restricted mathematical conclusions. The
README explicitly distinguishes a global coordinate seminorm ratio from the
same-H local error, states the natural transfer's two field-identity premises
before its numerical result, and preserves the parameter convention, positive
frame comparison norm and lack of physical smallness. Radial coordinate
amplitude derivatives are described as second derivatives, avoiding an
unintended claim about physical Riemann curvature. The initial candidate
trace is written as `u_W(r)=W(0,r)` to avoid confusion with the distinct W0
finite iterate.

All relative Markdown links in the staged capsule resolve. All 13 reference
reports have successful statuses. The reproduction description accurately
counts 13 normal/optimized reference replays, including the retained
independent endpoint report, and 43 corrupted executions in each mode:
ten refinement injections plus eight C2, seven normalization and eighteen
natural-transfer scope/premise cases. The runtime check-entry count from the
leaf audit is not presented as that many distinct mathematical results.

The only tracked pre-existing file modified in the reviewed worktree is the
top-level README; the workflow and research capsule are additions. Earlier
manuscripts, candidates and source checkpoints remain separate.

Final reviewed README SHA256 values:

```
Capsule README: e610cf4a23ea5a3300c08fb61bb34935a9798e42c8a6c602fbbbcc5b665f24db
Top-level README: d593dc477b4539442af473b0bdfab6f3c01951b4f4d467b6f2a7a70929416535
```

**Final assessment: PASS, with no outstanding scientific wording or scope
blocker.** Publication CI and remote-file readback are separate operational
checks performed by the parent workflow.
