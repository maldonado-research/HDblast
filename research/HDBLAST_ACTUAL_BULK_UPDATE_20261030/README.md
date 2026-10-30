# Boundary-preserving coupled bulk update and local enclosures

The [previous reconstruction](../HDBLAST_ROBIN_RECONSTRUCTION_20261030/)
satisfied its affine shell boundary law while retaining a frozen bulk source.
This checkpoint constructs a first coupled bulk update `W1`, using the
retained actual variable coefficients and affine bulk source acting on `W0`.
The shell law remains exact for the integral-defined update. Its bulk forcing
is evaluated at `W0`, so it is still a finite iteration rather than a
self-consistent solution at `W1`.

Three conditional results are established:

* The old reconstruction has a nonzero reduced forcing defect at the selected
  shell point. The new integral-defined candidate has a uniform defect upper
  bound smaller than the old point's proved lower bound by a factor over 2.5.
* The whole-cone error upper certificate for the new candidate is
  `X(U-W1)<3.541e25`, over 5.7 times smaller than PR96's certificate for `W0`.
  This compares upper bounds for different candidates; it does not measure a
  reduction in their unknown true errors.
* Local causal estimates narrow 53 of the 69 best target point boxes, with all
  69 contained in the preceding boxes. The scalar-mean widths shrink by about
  1.81 for value/radial derivative and 2.68 for time derivative.

The actual bulk convolution is **analytically enclosed**, not numerically
integrated or assigned a sign. Here “actual” refers to the fixed retained
coefficient/source functions under the stated model assumptions, not to
observationally established quantities. No full coupled Einstein solution,
nonlinear family, physical observable or Big Bang origin is established.

## A concrete forcing defect

With smoothed free wave `Bs`, frozen bulk operator `L0` and actual reduced
operator `La Y=A(x,r)Y_x+B(x,r)Y_r+D(x,r)Y`, define

```
f0=L0 Bs,
r0=La W0+F-f0.
```

The sign convention is actual right-hand side minus principal forcing.
At `(x,r)=(9/8,0)`, the [point calculation](theory/POINT_RESIDUAL_NOTE.md)
encloses all 23 components and proves 17 are nonzero. For the scalar mean,

```
-9.417e22 < (r0)_w:0 < -9.416e22,
||r0(9/8,0)||_H > 9e22.
```

The actual coefficient matrices are assembled at the point. The affine bulk
source is enclosed using genuine uniform component bounds covering it; its
exact point value and sign are not newly evaluated. Ninety-three source Gram
entries, field multiplicities, the packing map and wave-forcing normalization
are verified. The response bound `X(KF+LG)` is never substituted as a pointwise
source bound. An [independent arithmetic audit](theory/POINT_AUDIT_NOTE.md)
uses a third residual decomposition and reproduces the same nonzero rows.
Both computations share the inherited coefficient/source premises.

This is a reduced-evolution forcing defect computed from first jets and the
known principal forcing of `W0`. It is not a full pointwise Einstein or
constraint residual, a physical force, or the sign of the causal correction.

## One coupled update with the boundary law preserved

Let `K` be the mixed reflected wave propagator, `L` the causal boundary lift,
`PN` the Robin projector and `PD` the Dirichlet projector. Put

```
Bdry V=L(M PN V),       P=(I-Bdry)^(-1),
Psi(Y)=P[Bs+K(La Y+F)+LG],
W1=Psi(W0),            C=W1-W0=P K r0.
```

The [derivation](theory/ACTUAL_BULK_UPDATE.md) keeps the same fixed actual
`M,F,G` throughout. It forms the physical coupled source before reflecting
its output and assumes no matrix/parity commutation. The boundary inverse is
a Volterra/Neumann construction, requiring no inverse or positive spectrum
of the matrix `M`. The correction obeys

```
PN C_r=M PN C,       PD C=0.
```

With `k` the bulk operator cap and `b` the boundary cap, the update has
contraction majorant `beta=k/(1-b)`. The retained actual-initial-data target
`U` differs from the smoothed-data iterates. Its initial mismatch remains in
the [global estimate](theory/BULK_BOUND_NOTE.md):

```
E(n+1) <= e_initial/(1-b)+beta E(n),
X(U-W1) <= 3.540335165343e25 < 3.541e25.
```

The norm is `X(e)=max(sup||e||_H/3,sup||e_x||_H,sup||e_r||_H)` on
`0<=x<=3/2, r>=max(0,x-5/4), x+r<=7/4`. Coefficient uncertainty, initial
enclosures, smoothing, velocity interpolation and affine-source response
units are retained. The auxiliary smoothed-data fixed point is distinct from
the target's physical initial prescription; smoothing is not assumed to
preserve its physical constraints. The limiting initial-data majorant is an
upper estimate, not a positive lower bound on the unknown true error.

The remaining forcing defect, in the same sign convention, is

```
r1=La W1+F-(La W0+F)=La C,
sup||r1||_H <= (nA+nB+(3/4)nD) X(C) < 3.574e22.
```

The [cross-report arithmetic](reference/SUMMARY.json) verifies that this new
uniform upper bound is below the old sample's lower bound divided by 2.5.
It is a reduction of a registered reduced forcing defect, rather than a
measured field-error reduction or a numerical C2/Einstein residual certificate.

## Additional local information narrows the target boxes

The global norm-only [baseline](theory/COARSE_NOTE.md), with consistently
evaluated exact majorants, obeys `D0+E1=E0`: adding its correction and error
radii leaves the corresponding total norm radius unchanged. The
[local calculation](theory/LOCAL_NOTE.md) supplies genuinely additional
information on the causal triangle feeding the selected point,

```
Omega={(s,y):s>=0,y>=0,s+y<=9/8}.
```

Inherited exact zero-wedge support leaves active physical area `17/128`.
Each physical source point has at most two reflected copies and the wave
prefactor is one half. This bounds the integral continuously without moving
coefficients through parity. Local data, actual coefficient and genuine
affine-source bounds give

```
sup_Omega ||r0||_H <= 1.905589378958e25,
sup_Omega ||C||_H  <= 2.530947563205e24.
```

Combined with the new global error and the shared actual Robin law, these
give strict target scalar-mean bounds at the shell sample:

| Registered formal target | Lower | Upper |
| --- | ---: | ---: |
| `U_w:0` | `8.29036e30` | `8.29042e30` |
| `partial_x U_w:0` | `1.974621e32` | `1.974623e32` |
| `partial_r U_w:0` | `2.154067e27` | `2.154080e27` |

The sample is in scaled time `x=t/ell`, `ell=1e-6`, and inward scaled radial
coordinate `r`. These parameter-derivative modal coefficients are gauge- and
normalization-dependent; they do not give the full angular scalar sign,
proper-time observable, physical force or finite-amplitude signal. Broad
all-mode bounds do not replace earlier specialized tensor certificates.
Zero-M rows retain their source-only radial enclosures.

## Reproduction and remaining work

Python 3 standard library is sufficient. From the repository root:

```sh
C=research/HDBLAST_ACTUAL_BULK_UPDATE_20261030
python "$C/checks/check_bulk_update_identities.py"
python "$C/checks/check_bulk_update_bound.py" --repo-root . > /tmp/bulk-bound.json
python "$C/checks/check_bulk_coarse.py" --repo-root .
python "$C/checks/check_bulk_point_residual.py" --repo-root . > /tmp/bulk-point.json
python "$C/checks/check_point_residual_audit.py" --repo-root . --residual-report /tmp/bulk-point.json
python "$C/checks/check_local_bulk_enclosure.py" --repo-root . --bound-report /tmp/bulk-bound.json > /tmp/bulk-local.json
python "$C/checks/check_result_consistency.py" --repo-root . --bound-report /tmp/bulk-bound.json --point-report /tmp/bulk-point.json --local-report /tmp/bulk-local.json
```

Normal and optimized output must match byte for byte and reproduce all six
[reference reports](reference/). The workflow checks both modes. The exact
polynomial/geometry audit has 39 named groups, including 90 active-source
polygon/edge fixtures. Input hashes authenticate the recorded dependencies;
the coefficient driver's scientific roots are authenticated before imports
and reauthenticated after evaluation. Both point-coefficient programs share
that authenticated assembler and scientific-root checks. The remaining norm
and consistency programs use authenticated JSON without scientific imports.
Internal AI-assisted reviews and finite controls are
not external specialist review, proof-assistant verification, or a novelty
claim for the iteration method.

The next numerical task is to resolve the causal bulk convolution and its
derivatives, followed by quantitative second-derivative and subsidiary
constraint checks. Source support, source projection, actual background,
first-field and formal-response premises remain inherited. The finite update
does not complete the coupled second-order solution, prove a nonlinear
spacetime family, or establish a higher-dimensional origin of the Big Bang.
