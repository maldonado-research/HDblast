# Independent full23 reflected-kernel review

26 October 2026. This checks the new `compute/angular_ridge.py` with source
SHA256 `d6e935bd958512cb032e1ef3d9503d6fa159418bcf9b28051de6d2f43a43c02f`.
It supplies exact bounded kernel/representation controls. It does not by
itself evaluate the actual HDBLAST source or certify a continuous response
approximation error.

## Authentication and independence

The reviewer authenticates the complete seven upstream closures and the
latest shell precision closure before importing any frozen scientific code.
It separately checks the latest authentication helper pin, reauthenticates
every closure at exit, and checks the new kernel source pin before and after
execution. It reconstructs the new source byte for byte from the frozen
five-mean source using exactly the ten recorded replacements in
`compute/KERNEL_TRANSFORMATION.json`.

The expected integral values are computed by a separate standard-library
rational implementation. For the physical retarded diamond, coordinates
are `(s,v)` and its vertices are `(0,r-x),(x,r),(0,r+x)`. The independent
oracle clips first at `v=0` and then at every relevant piecewise-polynomial
ridge boundary. It integrates monomials by Green's formula

`integral s^i*v^j ds dv = contour s^(i+1)*v^j dv/(i+1)`.

First x and r derivatives are independent direct integrals along the upper
and lower characteristic edges, with the lower edge split where it crosses
the shell. Neither oracle uses a producer kernel or its retarded-integral
formula to obtain an expected value. Normalized profile primitives are
integrated segment by segment and joined by their constants of integration.

## Coverage and exact findings

The complete 23 labels, positive angular Gram weights and five reflection
signs agree with the authenticated full angular operator. Every slot is
active in the degree-four polynomial tests. Separate dense, noncommuting
23-by-23 A, B and D fixtures couple every input/output amplitude. The
piecewise-linear displacement and its independent piecewise-linear radial
derivative approximant are kept distinct, just as in the producer.

The checks cover:

- full23 degree-four forcing integrals on exterior and shell-crossing diamonds;
- both first derivatives by direct characteristic-edge integration;
- all normalized piecewise-linear profile primitives of orders 0 through 5;
- all coupled F/G primitives and the independently expanded physical source
  `A*free_x+B*free_r+D*free_value`;
- five Dirichlet and eighteen Neumann/Robin-lift modal rows;
- zero initial feedback and the exact initial zero wedge;
- nontrivial C1 joins across `r=x`, `r+x=1`, `r+x=5/4`, `r-x=1`,
  and `x-r=1`;
- all six authenticated old mean feedback samples, embedded into the 23 slots,
  with every inactive slot zero;
- all five individually malformed Dirichlet output rows, each tested with
  both a Robin input column and a Dirichlet input column.

The C1 interface controls use exact one-sided polynomial continuation. On
each sampled side no other profile knot is crossed. Seven rational sample
values determine the constant term of a polynomial of degree at most six;
the full value and both derivative limits coincide. The generic degree-four
source has integral degree at most six. The coupled piecewise-linear source
has degree at most two in each cell and integral degree at most four.

The physical source must be formed on `v>=0`, including every A/B/D coupling,
and only then reflected by the output projector `J=PN-PD`. Blindly evaluating
the original ridge decomposition at negative v gives a different integral.
The review retains four exact nonzero rational witnesses. No commutation of
J through A, B or D is assumed.

Fresh normal and optimized subprocesses also test four actual in-memory
producer mutations: all-even reflection, omission of the `v4` Dirichlet row
guard, blind extension of the coupled source, and omission of the PN
projection before M. Each must fail at its specified reason and write no
report. These mutate disposable in-memory source text; the root kernel and
all frozen files remain unchanged.

## Finite-candidate Robin law

The producer's first feedback satisfies at the shell

`PN (B1)_r = M PN (uP1 + primitive(vP1))`.

It is an affine lift of the free trace. It does not impose exact Robin data
on the completed finite candidate `V=Free(u0,vP1)+B1`. The exact defect is

`PN V_r-M PN V = M PN(uP1-u0-B1)`.

Thus even if `uP1=u0`, the remaining defect is `-M PN B1`. The review tests
both this coherent-input identity and a changed exact displacement with the
same P1 feedback, producing the additional `M PN(uP1-u0)` term. The root's
actual free displacement is the registered rational u0, while uP1 and duP1
occur only inside the forcing approximant. The review does not identify duP1
with the derivative of uP1 or assert that a piecewise-linear free displacement
is a globally C1 candidate.

## Complete two-term error identity

Let `U=b+T U` be the exact full23 integral equation. For any reconstructed
`B0,B1`, define `E0=b-B0`, `E1=T B0-B1`. Direct substitution gives the exact
identity

`U-B0-B1=T^2 U+(I+T)E0+E1`.

Consequently, in one common full23 X norm,

`||U-B0-B1|| <= q^2*Ucap+(1+q)*e0+e1`,

provided `||T||<=q`, `||U||<=Ucap`, `||E0||<=e0` and `||E1||<=e1`.
The reviewer checks the vector identity with a dense rational 23-by-23 T;
the inequality follows from the triangle and operator norm inequalities.

Here b is the complete affine term: free initial evolution plus the actual
quadratic bulk and shell lifts. Therefore e0 must include any omitted or
approximated affine source contribution as well as the initial approximation.
Likewise e1 must include coefficient errors, independent uP1/duP1 errors and
the shell displacement mismatch against the actual B0. Mean-only bounds or
an unweighted component norm cannot be inserted without a proved conversion
to this same full23 norm.

An alternative complete residual certificate is
`||U-V|| <= ||b+T V-V||/(1-q)`. For `V=B0+B1`, its residual is
`E0+E1+T B1`, giving the bound
`(e0+e1+q*||B1||)/(1-q)` when those quantities are available. Neither
identity requires M to be invertible or positive. The inherited physical
subsidiary theorem identifies the exact coefficient U; these kernel controls
and a separately complete numerical ledger then bound the reconstruction's
distance to U.

## Artifacts and replay

`check_angular_ridge.py` contains the independent integrator and bounded
fixtures. `verify_kernel_review.py` runs normal/optimized positive checks and
the eight specified rejection subprocesses. The canonical report and receipt
are under `verified-final/`.

Only new bounded checks were executed. No old numerical campaign, mesh or
source certificate was rerun, and no root kernel change was made.
