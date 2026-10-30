# One-cell refinement of the unchanged full23 initial-data certificate

This calculation tightens the continuous actual-minus-nominal initial-data
upper bounds for the **same** frozen candidate used by the preceding Robin
reconstruction. It does not change its initial displacement, P1 velocity,
smoothing widths, frozen matrices, shell law or source. A smaller upper
certificate does not measure a change in the actual data or actual error.

The inherited initial report partitions `[1,7/4]` into 192 intervals of width
`1/256`. Its channel named `rational_packet_cutoff_and_E_approximation`
dominates the old radial-derivative bound in cell7,
`[263/256,33/32]`, at the end of the cutoff transition. The other191 cell
certificates are retained unchanged. Both actual and nominal initial data are
zero on `[0,1]`, under the inherited zero-collar premise.

## Enclosing exactly the same difference on smaller intervals

For each of three difference channels `c`, each of three kinds
`k=(U,partial_r U,partial_x U)`, and each of 23 registered modes `i`, let
`D_c,k,i(r)` denote the difference channel represented by the old paired
algebra. The fixed candidate and this difference decomposition are not
modified. The old cell certificate establishes

```
|D_c,k,i(r)| <= B_old,c,k,i       on the complete parent interval.
```

On a dyadic subinterval `I`, the same authenticated interval evaluator gives a
new cap `B_eval,I,c,k,i`. Both apply to the same mathematical function on `I`,
so the implemented intersection is the valid cap

```
B_I,c,k,i = min(B_old,c,k,i, B_eval,I,c,k,i).
```

This intersection makes the new cell certificates componentwise
nonincreasing. It is valid even if a differently ordered outward evaluation
has a last-bit excess over its parent bound. The unclipped evaluator result
and effective cap are both saved for every final leaf.

The final leaves are an exact dyadic cover of the parent cell. The binary
path records the subinterval: path `p` of length `d`, interpreted as integer
`j`, corresponds to

```
[left_old+(right_old-left_old)*j/2^d,
 left_old+(right_old-left_old)*(j+1)/2^d].
```

The checker sorts the leaf intervals and verifies continuity, nonempty
interiors, the exact parent endpoints and the path formula. All leaves lie
inside the inherited analytic cutoff branch; none crosses a branch join.
Together with the191 untouched intervals they cover the complete inherited
nonzero footprint. Neither endpoints nor an interior gap are discarded.

## Component maxima precede Gram norms

For every `c,k,i`, the global component cap is the maximum of that component
over every final leaf and every inherited untouched cell. With the frozen
positive Gram weights `H_i`, the reconstructed channel norm and total
difference caps are

```
C_c,k,i = max_all_cells B_cell,c,k,i,
N_c,k = sqrt_upper(sum_i H_i C_c,k,i^2),
E_k = sum_c N_c,k.
```

The maxima occur **before** the Gram norm. Taking a maximum of independently
computed cell norms instead would be a different bound and is not used here.
The sum over channels is the triangle inequality for the inherited
telescoping paired difference identity. The checker first reconstructs the
old component maxima, channel norms and total caps exactly from the pinned
old rows, then applies the same aggregation to the refined cover.

Every square-root upper bound uses exact rational arithmetic and a112-bit
outward grid. No floating-point value determines selection, stopping,
coverage or a published inequality.

## Deterministic stopping and review record

The algorithm starts with the inherited parent certificate and bisects the
leaf with the largest channel1 radial Gram-square, breaking ties by binary
path order. It stops when the global total radial-derivative H bound is
strictly below `2e24`. It has hard budgets of128 leaves and depth8 and fails
explicitly if either is exhausted. It also checks a strict100-fold
improvement over the old radial upper certificate and that every global
component cap is nonincreasing.

The machine report records the complete192-row cover description, every
evaluated final leaf, raw and intersected3-by-3-by-23 caps, split history,
final component maxima and norm arithmetic. The191 unchanged rows point to
their exact indices and intervals in the authenticated inherited report.

The default run finishes with10 leaves of maximum depth5 in the refined
parent, giving201 continuous intervals in total. Its exact fractions imply
these outward-rounded strict bounds for the continuous actual-minus-nominal
initial profiles:

| Kind | Refined H upper bound |
| --- | ---: |
| Displacement `U` | `1.121e21` |
| Radial derivative `partial_r U` | `1.890e24` |
| Scaled-time derivative `partial_x U` | `2.683e13` |

The old radial upper certificate is about `1.5334379593e27`; the new exact
bound is about `1.8898899191e24`, over800 times smaller. These are upper
certificates for unchanged mathematical functions. Smoothing and P1 velocity
interpolation errors are separate and remain charged when these refined
profile bounds are propagated to the Robin reconstruction.

## Authentication and limitations

The checker pins `full23_initial.py` and `INITIAL_REPORT.json` by SHA-256
before importing the producer. Its `Context` authenticates all seven complete
upstream scientific roots before scientific imports and reauthenticates them
after calculation; the producer and report are pinned again at completion.
It does not change or regenerate any frozen file. Normal and optimized Python
must produce identical JSON.

These checks authenticate and refine inherited conditional mathematical
certificates. They do not independently establish the underlying first-field
or background equations, nonlinear existence, a numerical C2 error, a physical
observable or a higher-dimensional origin of the Big Bang. The smaller caps
can be substituted in the direct C1 error ledger for the unchanged Robin
reconstruction; its remaining bulk forcing is still frozen.
