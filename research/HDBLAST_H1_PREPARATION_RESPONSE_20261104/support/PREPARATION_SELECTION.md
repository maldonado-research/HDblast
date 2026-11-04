# Preparation selection: first-order energy and initial cut charges

The main [H1 response estimate](../DERIVATION.md) is a conditional bound on
the response to a specified **second preparation**. Its mathematical energy
does not supply a physical rule choosing that preparation. This note explains
why two different candidate selection rules—an energy depending only on the
shared first tangent, or charges depending only on shared initial boundary
jets—cannot supply the missing information. It does not rule out other
energy/state principles or weaker-norm response estimates.

Retain the background, physical initial section, scalar-label calibration,
parameter normalization and exact-family premises of the
[preparation-dependence theorem](../../HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/DERIVATION.md).
In a common initial-data chart,

```
x_lambda(e)=x0+e v+(e^2/2) w_lambda+o(e^2),
w_lambda-w_base=lambda(H_h,0,0,0).
```

Here `H_h` is the admitted compact metric-mean increment; the other canonical
second sectors are fixed. Its support avoids the initial shell collar. Every
`lambda` is fixed and finite, with its own possibly different sufficiently
small parameter interval. The comparison fixes the supplied physical initial
section; it neither selects that section nor excludes a separate rule comparing
different physical sections.

Write `D` for the difference of actual second areas of the continued scalar
levels, divided by background area. Its radial derivative uses the fixed
background scalar-label parametrization, as in the main estimate.

## A quadratic of the fixed first tangent

For a specified quadratic functional `Q(v)=B(v,v)`, all preparations have the
same value because their complete first tangent is the same. The integration
domain, background structures, generator, corner convention and normalization
must also be shared. This includes an expression `Omega_0(v,A v)` when the
symplectic form and linear operator are actually defined. A threshold or
extremization using only `Q(v)` cannot distinguish the second preparations.
First-tangent dependence of canonical energy is already explicit in
Hollands–Wald; it is not a new general principle.[^HW]

This statement must not be extended to every nonlinear energy. For a C2
functional of the full data,

```
(F o x_lambda)''(0)=DF_0[w_lambda]+D^2F_0[v,v],
Delta F''(0)=DF_0[lambda(H_h,0,0,0)].
```

The last term can be nonzero. Stationarity on the allowed increment must be
proved before removing it. Field-dependent generators or reference choices
also require their variation terms.

## What the retained relative Noether construction defines

The [shell Noether calculation](../../HDBLAST_SHELL_NOETHER_20261013/derivation/FULL_NONLINEAR_SHELL_NOETHER_IDENTITY.md),
especially equations (21)–(29), includes two compact interiors, their matched
internal shell, both GHY terms, one tension and the specified relative corner.
For fixed common shell-tangent generators it defines

```
H_xi=sum_i integral_Ci C_xi,i + integral_S B_xi,
B_xi=-i_(I dot xi) epsilon_Sigma.
```

`I` is the Israel residual. This relative constraint generator vanishes on
exact matched solutions; its zero is not a positive physical energy. Under
that construction's **full classical matched-domain hypotheses**, a reducible
background Lorentz generator with fixed compensator gives the necessary test

```
Omega_rel,0(v,Kcal v)=D^2H_K|0[v,v]=0, with DH_K|0=0.
```

Those full-domain hypotheses are not established by the local five-mean
response problem. The first-tangent non-selection observation above does not
require this test. No asymptotic charge boundary, positive conserved energy
or globally timelike Killing generator is assumed. In particular the
background angular metric contains `cosh^2(t)`, so `partial_t` is not a
stationary Killing field merely because it labels the shell clock.

## Initial finite-boundary-jet lemma

Suppose a specified initial charge factors as `J(e)=G(b(e))`, with `G` and
the boundary-input map C2, where `b` comprises finite boundary jets of the
complete initial data and all auxiliary inputs. The latter include the cut
and embedding, lapse/shift and needed generator jets, reference subtraction,
corner convention and normalization. Require these auxiliary inputs to agree
through second order too. Fix a trace order for which the retained regularity
justifies this map; lower C4 control does not license arbitrary higher jets.

The common first tangent and compact interior second difference give equal
`b0,b1,b2` on the shared initial collar. Consequently

```
J2=DG_b0[b2]+D^2G_b0[b1,b1]
```

and the zeroth, first and second charge coefficients agree. The common-first
quadratic terms cancel in the difference; they were not discarded. This is
the standard chain rule and locality, not a new charge identity. A bare
Noether cut integral is not automatically an integrable Hamiltonian: boundary
conditions, symplectic potential and shell/corner completion remain essential.

The lemma does not cover a cut through the packet, a future cut after
propagation, or an interior-dependent choice of generator, reference or
embedding whose boundary jets change. Second-order agreement alone does not
prove equality of full nonlinear charges. With only C2 dependence, no cubic
or quartic coefficient, or quartic sufficiency, follows.

## Consequence and compatibility with the H1 estimate

Any rule using only these shared inputs, and admitting the base, also admits
the previously certified `lambda=-10^21` preparation with opposite proper
volume and expansion coefficients. The
[initial area-profile certificate](../../HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_20261104/DERIVATION.md)
also gives `|D_lambda(T)|>3.7621e33`. More generally, on `I=[1,5/4]`,

```
D_lambda=lambda D_h,
R_A(lambda h)=|lambda| R_A(h),
||D_lambda||_H1=|lambda| ||D_h||_H1,
||D||_H1^2=integral_I (D^2+D_r^2) dr.
```

The existing nonzero profile has positive finite H1 norm. Thus both norms are
unbounded along the admitted fixed-finite-`lambda` line while `Q(v)` and the
specified initial charge coefficients stay fixed. No finite bound on either
norm can depend only on those shared inputs over this whole line. This uses
no `lambda(e)` limit or uniform nonlinear amplitude interval.

The main tilted wave energy instead depends on the homogeneous **second-field
difference** and its data. It can therefore vary with the preparation and
support an H1-to-observable estimate under its stated operator and regularity
hypotheses. It is not identified with a positive conserved physical Noether
energy. Failure of a first-tangent selection rule, or of conversion from H1
control to a derivative supremum, does not prohibit that direct response
estimate. Selecting a physical H1 class and radius remains additional work;
a uniform nonlinear parameter window requires separate remainder control.

[^HW]: S. Hollands and R. M. Wald, *Stability of Black Holes and Black Branes*,
    [arXiv:1201.0463v4](https://arxiv.org/abs/1201.0463v4), equation (46) and the
    paragraph after (47); Commun. Math. Phys. **321** (2013) 629–680. Their
    stationary black-hole and boundary hypotheses are not imported here.
    See the [primary-source comparison](../REFERENCES.md). The contribution
    of this diagnostic is its application to the retained preparation witness.
