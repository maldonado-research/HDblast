# Independent review of the differentiated fourth-gauge-jet symbol

This review concerns the highest-jet part of the complete reduced Einstein
source. It does not provide a bound for the actual differentiated source,
prove that the actual source derivative is nonzero, or prove that all possible
data-specific identities require an independently certified C4 norm.

## Independent indexed reduction

The complete source convention is the one in
`research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/geometry/GEOMETRIC_SOURCE_DERIVATION.md`:
`QEred=QE-nabla_(a QC_b)` on the actual harmonic first field, with full
metric-lowered `QC` and one-half symmetrization. The affine second variation
is used, not an epsilon-squared Taylor coefficient.

At a background-normal event, write

```
P_abc=(d_b q_ca+d_c q_ba-d_a q_bc)/2,
J^c_ab=-2 q^(cd) P_dab,
QC_b=-2 q^(cd) P_bcd.
```

Keep only the terms that change when a third metric derivative changes while
all metric jets through order two stay fixed. Then

```
delta d_t QRic_ab
 = -2 q^(cd) [delta d_t d_c P_dab-delta d_t d_b P_dac],
delta d_t nabla_(a QC_b)
 = -q^(cd) [delta d_t d_a P_bcd+delta d_t d_b P_acd].
```

Expanding P and canceling the paired mixed derivatives gives

```
delta d_t QEred_ab = q^(cd) delta d_t d_c d_d q_ab.       (1)
```

The other Einstein/scalar matter terms have fewer derivatives of the metric
and do not change under this particular highest-gauge-jet variation. The
general off-harmonic correction `2 P^d_ab C1_d` also has only metric first
jets and therefore cannot change this principal symbol. Background curvature,
connections and coefficient derivatives multiply lower gauge jets, which
are held fixed. Thus this is the highest-jet part of the **complete** source,
not an isolated Ricci or potential summand.

`check_independent_principal_symbol.py` computes the left-hand side from the
indexed P/J/QC expressions for generic rational symmetric q values in all
five frame slots. It imports neither the geometric coefficient map nor the
parent source-symbol checker. Eight generic signature/sign fixtures agree
with (1) exactly. Three controls reject Euclidean time-index raising, omitted
QC gradient and incorrect symmetrization.

## Radial gauge fourth jet

The registered first-gauge transformation has

```
A=-4u+2f_y,  B=f_t/rho-rho b_y,
C=f/(rho c)+rho c d_y,  D=-2(u+a f+b_t),
E=c d_t-b/c,  L=2(u+a f+theta b),  M=2d.
```

Vary only the fourth derivatives of f. Put
`zeta_y=1`, `zeta_t=sigma*rho`, `sigma=+1 or -1`, and

```
delta f_(y^i t^j) = zeta_y^i zeta_t^j tau,   i+j=4.
```

All gauge jets of total order at most three are unchanged. The physical
frame covector is `(1,sigma,0,0,0)`, which is null in signature
`(+,-,+,+,+)`. Consequently every pair trace of the quartic jet is zero.
This variation survives the **highest-order** differentiated wave/gauge
trace conditions. It must not be called a separately constructed global
solution of the full curved registered gauge problem.

The only third ordinary metric derivatives affected are

```
delta (A_yy)_t = 2 sigma rho tau,
delta (A_yt)_t = 2 rho^2 tau,
delta (A_tt)_t = 2 sigma rho^3 tau,
delta (B_yy)_t = rho tau,
delta (B_yt)_t = sigma rho^2 tau,
delta (B_tt)_t = rho^3 tau.
```

The lower metric null contraction is

```
q^(cd) zeta_c zeta_d=(A+D-2 sigma B)Y.
```

Since the metric forcing is `Fg=-2 ell^2 rho^2 QEred` and
`d_x=ell d_t`, (1) gives

```
delta (Fg_yy)_x=-4 sigma ell^3 rho^3
                (A+D-2 sigma B) Y^2 tau,
delta (Fg_y0)_x=-2 ell^3 rho^3
                (A+D-2 sigma B) Y^2 tau.                 (2)
```

All other metric sectors have zero symbol for this particular radial-gauge
variation. In the registered scalar basis,
`Y^2=1/6+(1/2)Z+F`. Thus at the purely algebraic specialization
`A=1`, `B=D=0`, `ell=rho=c=1`, the yy coefficient triple is
`sigma*(-2/3,-2,-4)` and the y0 triple is `(-1/3,-1,-2)`, agreeing with the
complete-map extraction. This specialization is a coefficient witness; it
is not claimed to be the actual HDBLAST background or first solution.

## What follows, and what does not

Equation (2) is not identically zero. Merely replacing wave-trace combinations
of fourth gauge derivatives by lower derivatives does not universally remove
the differentiated source's highest gauge dependence. The actual lower-field
combination `A+D-2 sigma B` and the actual fourth gauge jet may have special
relations, but those would require separate evidence. A formal nonzero
coefficient does not prove a nonzero actual source derivative.

The calculation sharpens the remaining obligation: bound the relevant fourth
gauge combinations, or prove an equivalent cancellation/estimate for the fixed
actual data and complete source. It is invalid to differentiate a C0 magnitude
cap, treat a C3 jet box as a C4 box, or infer an actual derivative value from
this principal-symbol witness. No quantitative C2 error, full Einstein
residual, nonlinear family or physical Big Bang mechanism is established.
