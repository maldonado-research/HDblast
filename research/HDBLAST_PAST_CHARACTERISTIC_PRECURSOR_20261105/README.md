# A local characteristic precursor for the registered shell packet

5 November 2026. Conditional first-order Einstein–scalar research checkpoint.

This continuation constructs a short-past linear precursor for the previously
registered compact quadrupole packet. It specifies two compatible characteristic
traces by a convergent integral construction and bounds their values. It does
not specify a remote blast source or a nonlinear cosmological preparation.

## What is established

The actual normalized master satisfies

```
w_tt - w_zz + mathcalW(t,z) w = 0,
w(0,z) = q_ell(z),  w_t(0,z) = exp(-z) q_ell(z),  ell = 10^-6.
```

The shell is at `z=0`; depth increases away from it. With `x=t/ell`,
`r=z/ell`, and `U(x,r)=w(ell*x,ell*r)`, the supplied data have support
`1<=r<=2`. A fresh rational first-exit calculation checks the static collar
through `z=8ell`, conditional on the retained actual branch and shell box.
The exact potential formula then gives, directly on
`-ell<=t<=0`, `0<=z<=8ell`,

```
5500 < mathcalW(t,z) < 103000,
|ell^2 mathcalW| < 103/10^9.
```

This uses the explicit static background and time dependence. It does not
transfer future five-mean operator or evolved-field estimates to the past.

**Both principal wave directions persist before possible past shell contact.**
For every closed interval `-a<=x<=0` with `a<1`, the backward solution has
support in `1-|x|<=r<=2+|x|` and an open gap from the shell. Its characteristic
profiles `J_plus=U_x+U_r` and `J_minus=U_x-U_r` have L2 norms within a relative
`10^-6` of `||Q_r||_L2`. Each normalized squared principal
weight therefore has share greater than `0.499999` on `-1<x<=0`.
These are diagnostics in the specified master normalization; they are not
physical energy or radiation flux.

**Two nontrivial past characteristic traces reproduce the current packet.**
In null coordinates `u=x+r`, `v=x-r`, take

```
D = {u>=1, v>=-2, u+v<=0}.
```

Its past corner is the interior point `(x,r)=(-1/2,3/2)`; its future base is
the entire current packet interval `[1,2]`. The restrictions to its two past
faces share the same actual corner value. The convergent Cauchy construction
and its uniform remainder give:

| Master value | Strict bound |
|---|---:|
| Common interior past corner | less than `-4*10^-8` |
| Midpoint of either past face | greater than `0.499999` |

Consequently each face trace has an interior zero and is nonconstant. Their
two-face Goursat integral operator has norm at most `103/(16*10^9)<1`, so
these compatible traces determine the local master solution and reproduce
both supplied Cauchy data on the base. Setting either whole-face trace to
the common corner constant cannot reproduce this packet on this diamond.

The sign changes and zeros survive a positive smooth conjugation, including
the factor relating this normalized master to the longitudinal metric master.
The numerical principal-weight bounds do not remain unchanged under arbitrary
conjugations or gauge choices.

## Interpretation and limits

The retained reconstruction identities give the coupled first linearized
Einstein–scalar fields, rather than a separate test scalar. Their first shell
junction rows are trivial on this past interval because of the open support
gap. Identification with a specified harmonic representative retains the full
initial gauge completion and linear uniqueness premises.

The backward master is constructed analytically by a convergent Volterra
series. No new numerical PDE evolution, source-oracle sample or background
trajectory is used. Finite rational checks verify numerical inequalities and
input identities; they do not replace the continuous proofs or the retained
branch and reconstruction hypotheses.

An earlier incoming history may still supply these two traces after scattering
or boundary interaction. This result neither selects that history nor rules
out all incoming histories. It does not complete second-order constraints,
control a nonlinear family on the past domain, specify a radiation-producing
matter sector, give a common finite perturbation interval, or establish a
higher-dimensional origin of the Big Bang.

The standard Volterra, energy and characteristic methods are established
mathematics. The contribution here is a quantitative construction and data
constraint for this registered packet; broader originality requires specialist
assessment. Earlier papers and capsules remain unchanged. No new Zenodo edition
or journal submission accompanies this checkpoint.

## Proof and reproduction

Read the [primary proof](primary/PROOF.md), the
[independent characteristic construction](characteristic_review/PROOF.md),
the [independent coefficient and remainder derivation](independent/PROOF.md),
the [primary-source comparison](REFERENCES.md), and the [next tasks](NEXT_TASKS.md).
The working paper is available as [PDF](manuscript/PAST_CHARACTERISTIC_PRECURSOR.pdf)
and [editable TeX](manuscript/PAST_CHARACTERISTIC_PRECURSOR.tex) for author and
specialist review.

The compact capsule includes all selected runtime input documents. It does not
need the earlier programs or complete earlier campaigns. With Python 3.12.14,
from this capsule directory and using a fresh external output directory:

```sh
python -I -S -B checks/run_checkpoint.py --output-dir /tmp/hdblast-past-precursor-fresh
```

The driver authenticates the complete inventory before running new arithmetic
and afterward, reproduces canonical reports with normal and optimized Python,
and requires identical report bytes. Source authentication is also available
with `--authenticate-only`. These finite checks verify inputs and quantitative
guards; the continuous proof and retained analytic premises remain explicit.
