# Independent review of the natural-equation residual identity

30 October 2026. This review concerns the unchanged original-source
`W1_quadratic` candidate. It checks the analytic identity and the premises
needed to use its C1/C2 error bounds. It does not independently certify the
additional physical constraint identity assumed below.

## Assessment and necessary premises

The proposed identity has the correct signs and factors. Let `U` denote the
same formal target and let `e=W1_quadratic-U`. Use the background covariant
derivative throughout, and define the scaled reduced defect by

```
R[V]=V_xx-V_rr-La V-F.
```

The natural equations are the trace-reversed Einstein–scalar equation
`E=Ric-dphi tensor dphi-(2/3)V(phi)g` and the scalar equation
`S=Box_g phi-Vprime(phi)`. Expansions use the factorial convention
`g=g0+lambda q+(lambda^2/2)W`, and similarly for the scalar. Thus `E2` and
`S2` are second parameter derivatives; their `lambda^2` coefficients are
one half of the displayed derivatives.

The retained complete original-source identities give, for any suitable
second field `V`,

```
E2[V]=R_metric[V]/(2 ell^2 rho^2)
       + nabla_(a DC[V]_b) + (QE-QEred),
S2[V]=-R_scalar[V]/(ell^2 rho^2).
```

Parentheses denote symmetrization with factor `1/2`. These signs follow
from the registered metric forcing `-2ell^2 rho^2 QEred` and scalar forcing
`+ell^2 rho^2 Qphi`. The scalar row already retains its linear `DC[V]`
coupling; adding a second scalar constraint-gradient charge would double
count that coupling.

The complete source-map identity

```
QE-QEred=nabla_(a QC_b)
```

is the retained actual harmonic specialization, requiring the actual first
harmonic constraint to vanish on a neighbourhood. Its arbitrary off-shell
version has further first-constraint terms and cannot be substituted here.
The old exact 20-row source audit applies because the first field and
original sources have not changed.

An additional second-target premise is essential:

```
DC[U]+QC=0
```

as a `C1` covector identity on the physical domain, with the appropriate
one-sided traces. Merely bounding the reduced target, its C2 regularity or
its distance to a numerical candidate does not establish this identity.
Pointwise values without derivative regularity would also be insufficient
to differentiate it. Under this explicit premise, linearity of `DC` gives

```
E2[W1q]=R_metric[W1q]/(2 ell^2 rho^2)
        + nabla_(a DC[e]_b),
S2[W1q]=-R_scalar[W1q]/(ell^2 rho^2).
```

The error convention matters: if error is instead written `U-W1q`, the
symmetrized-gradient term has the opposite sign. Norm bounds are unchanged.
The older reduced forcing-defect convention is opposite to `R` in some
reports, likewise leaving its absolute norm cap unchanged.

## Applying the retained harmonic-gradient map

The PR101 map and connection estimate apply to any sufficiently regular
full23 field, and therefore to `e`. With

```
P=ell rho DC[e]=Rx e_x+Rr e_r+Zmat e,
```

the retained component and coefficient operator bounds imply

```
P0 <= (15/8) max(Vx,Vr)+z0 V0,
Px <= (15/8) max(Vxx,Vxr)+z0 Vx+zx V0,
Pr <= (15/8) max(Vxr,Vrr)+z0 Vr+zr V0.
```

`V0=3E1` and `Vx=Vr=E1` follow directly from the published X norm. This
safe value bound requires no extra zero-wedge improvement. The new C2
error estimates supply the three second-jet caps. The retained scaled
radial background identity `rho_r=-alpha rho` gives

```
e_y DC[e]=-(P_r+alpha P)/(ell^2 rho^2),
e_0 DC[e]=P_x/(ell^2 rho^2).
```

The `+alpha P` term in the radial numerator and the physical-frame
connection terms remain necessary. The complete PR101 estimate is

```
G=(Px+Pr+[2 alpha_abs+5ell+(7/4)(alpha_abs+eta_abs)]P0)
    /(ell^2 rho_lower^2).
```

It bounds the positive angular L2 norm of the full background covariant
gradient. Symmetrization is an orthogonal contraction in the positive
tensor Frobenius norm, so its norm is no larger. The metric input norm
counts mixed tensor components twice; the covector output norm does not.
Those two different Gram normalizations must retain their authenticated
map conventions.

If `F0` bounds the full reduced defect `R[W1q]`, positivity of its
metric/scalar direct-sum norm bounds either block by `F0`. Consequently

```
||E2[W1q]|| <= F0/(2 ell^2 rho_lower^2)+G,
||S2[W1q]|| <= F0/(ell^2 rho_lower^2).
```

No causal-filter discrepancy or `2F_H` charge belongs in this unchanged
original-source candidate. The prior filtered candidate has a different
defect ledger and must remain separate.

## Limits of the result

This is a conditional natural-equation coefficient-residual bound with an
explicit second-target harmonic/subsidiary premise. It is not an independent
verification of all constraints. Numerical error bounds cannot remove that
premise by themselves. A finite positive background-frame/angular L2 norm
bound does not establish physical smallness, a Lorentz-invariant scalar
residual, an observational prediction, a finite-amplitude nonlinear remainder
or a higher-dimensional Big Bang origin. The note checks the analytic
reduction; the independently reviewed numerical checker must authenticate
all retained map, background, source and new-error inputs separately.
