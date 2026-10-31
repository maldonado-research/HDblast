# Exact packet Cauchy data and conditional response bounds

31 October 2026. This checkpoint separates two results: tighter certificates
for the unchanged `W1_exact_packet`, and a distinct approximation
`W1_exact_Cauchy_packet` with a new initial velocity. The actual formal
second-response target U, original F/F_x/G/M, operators, physical-side
domain and positive all23 H norm are unchanged. Earlier candidates and
frozen manuscripts are preserved.

## A genuine local velocity improvement

At the initial point (x,r)=(0,65/64), V=U_x(0,r), and the actual target
velocity derivative has positive H norm near 2.2502771243e-12. The new
candidate satisfies

```
||V_r - Vnew_r||_H < 6.556e-25,
||V_r - Vnew_r||_H / ||V_r||_H < 2.914e-13.
```

A fresh all23 calculation proves the old smoothed P1 velocity's error at
that same point exceeds 6.194e-6. The new error upper is below half that
old error lower. This proves a true local improvement in the stated
coordinate norm, not just comparison of two upper allowances. Since the
finite response shares its first Cauchy jet, this is also its initial mixed
second-derivative error. It is not a physical or observational relative error.

The old derivative is evaluated exactly: the convolution window isolates
one P1 node, so an even normalized kernel gives half the adjacent secant
slopes. Actual target and new velocity are evaluated independently with
exact packet jets, positive Gram weights and a separately proved local
denominator. A global normalization is not substituted for that denominator.

## Complete conditional evolved certificates

The original source and shell law are retained. Complete absolute candidate
jets and paired target-error jets have different roles and are recomputed
where the candidate changes.

| Certificate | Published `W1_exact_packet` | Same candidate, refined bounds | New `W1_exact_Cauchy_packet` |
| --- | ---: | ---: | ---: |
| Evolved xx/xr errors, each H supremum | <2.135e25 | <1.330e25 | <3.623e24 |
| Evolved rr error, H supremum | <2.144e25 | <1.338e25 | <3.700e24 |
| Complete C1 error, X norm | <1.603e23 | <1.365e23 | <1.304e23 |

The new candidate's second-derivative upper certificates are more than
5.7 times tighter than the published exact-packet counterparts. Absolute
global bounds remain large. Global factors compare upper certificates and
do not prove corresponding reductions of unknown global true errors.

For the unchanged candidate, three bisections refine a complete 49-leaf
cover to 52 leaves. Its absolute initial displacement rr cap falls from
3.177e35 to 9.399e34, over 3.38 times tighter. Parent/child intersections
of component magnitude allowances are valid by subset inclusion, and all
23 Gram contributions are reaggregated. The function and its true error
are unchanged. Separate arithmetic and a direct all23 expression table
also establish complete independent bounds.

## New Cauchy candidate definition

The displacement remains the completed PR109 nominal displacement, using
q=exp(4-1/((r-1)(2-r))) on 1<r<=7/4, q=0 for r<=1, and E=exp(-ell*r)
with ell=10^-6. All five mean overrides and all23 slots are retained.

The new velocity is the completed registered expression with the same
exact q/E and nominal background polynomial. It uses neither the old P1
interpolant nor its convolution. At r=7/4 each completed velocity component
has an affine exterior with its own exact value and first radial derivative.
The displacement's completed-component quadratic exterior is unchanged.
The data are C2/C1 as required for the physical-side C2 response argument;
global C3 displacement or C2 velocity is not claimed. Extending only a raw
packet is not this construction. Exact flatness preserves the support
collar, parity, compatible corner and retained causal wedge.

Shared q and E differences vanish by analytic identity. Background
coefficient differences and separately validated ODE derivative remainder
channels remain. Enclosure midpoints never define the functions.

The finite causal construction is unchanged in form:

```
W0=P[Bs+K(L0 Bs)+LG],
W1=P[Bs+K(La W0+F)+LG].
```

Fresh absolute v/v_r allowances can exceed their old P1 counterparts and
are used in the correction ledger. Fresh velocity-error allowances replace
the old charges only for the new candidate. The old candidate retains its
complete published initial-error ledger. Each response proof keeps the
W0-Bs cancellation, nonzero correction initial accelerations, original
unfiltered F_x, full commutator and Robin/G_x terms. The exact r1(0)=0
follows from the common first Cauchy jet; it does not set the accelerations
or the initial actual source to zero.

## Verification, assumptions and limits

Use a complete checkout and Python 3.12.14. From the repository root, choose
a fresh output directory outside the checkout:

```sh
python -I -S -B research/HDBLAST_EXACT_CAUCHY_PACKET_20261031/checks/run_checkpoint.py \
  --repo-root . --output-dir /tmp/hdblast-exact-cauchy-fresh
```

The closed source catalog and replay plan specify 17 complete canonical
reports, each checked in normal and optimized execution, and 340 deliberate
corruption subprocesses. Counts describe executed tests, not independent
theorems. Runtime wrapper copies receive documented interpreter-isolation
flags; mathematical source files and canonical reports remain unchanged.

The package supplies frozen definitions, complete interval covers, exact
rational reports, independent internal evaluations, effective deliberate
corruption tests and portable normal/optimized replay. Input authentication
does not claim to rerun every inherited analytic proof. Internal AI-assisted
reviews do not replace external specialist review or establish originality.

Native checks replay the new subdivision nodes. Separate arithmetic
evaluators establish their own complete enclosures; the wider independent
velocity boxes are not numerical reproductions of every tighter native box.
Report audits check coverage and all23 aggregation, while independent exact
response checks rebuild the ledger from the authenticated input bounds.

The inherited domain and norm are

```
D+={0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x},
X(V)=max(sup ||V||_H/3,sup ||V_x||_H,sup ||V_r||_H).
```

The actual target/reconstruction, source regularity, physical shell,
corner and subsidiary assumptions remain explicit analytic premises.
These nominal auxiliary data are not assumed to satisfy exact physical
constraints. No old point or quadrature table is inherited after changing
velocity, and no new natural-equation residual or invariant observable is
computed. Earlier conditional exact initial-data/precontact families retain
their own hypotheses. Full finite-amplitude nonlinear evolution through
shell contact and a higher-dimensional Big Bang origin remain unestablished.
No Zenodo version or journal submission accompanies this checkpoint.
