# Quantitative persistence of the selected formal scalar mean

This is a consequence of the retained conditional full causal construction and
its accepted refined norm bound. It gives signed enclosures for an evolved
formal coefficient without evaluating its waveform. The coefficient remains
in the registered gauge and second initial prescription.

Write `ell=10^-6`, `x=t/ell`, `r=z/ell`, and let `U=(k,w)` be the actual
integral-defined formal second parameter coefficient on
`D={0<=x<=3, r>=0, x+r<=8}`. The retained solution belongs to `C1(D;H)` and
satisfies

    ||U||_X=max(sup_D |U|_H/3, sup_D |U_x|_H, sup_D |U_r|_H)<1.03e33.

Here H integrates the positive metric Frobenius norm plus the scalar square
against sphere probability measure. Define `m(x,r)=integral_S3 w(x,r,omega)dmu`.
Because `W0=1` and `integral dmu=1`, Cauchy–Schwarz gives

    |m_x| <= ||w_x||_L2(dmu) <= |U_x|_H < 1.03e33.

The mean functional has norm exactly one: the constant scalar attains it.
The selected angular point-evaluation constant `sqrt(321/19)` is irrelevant
to this mean. The derivative has unit weight in X; the factor three applies
only to values. Since U is C1 and each vertical segment lies in D,

    |m(x,r)-m(0,r)| <= integral_0^x |m_x(s,r)|ds < x*1.03e33.

The accepted initial certificate bounds the actual initial mean over two
complete closed radius cells:

| Sign | Closed scaled radius cell | Initial strict magnitude lower bound |
| --- | --- | ---: |
| positive | [73/64,585/512] | 7e30 |
| negative | [951/512,119/64] | 7e30 |

For the closed positive time interval `1/1000<=x<=1/200`, the drift is strictly
less than `5.15e30`. Consequently, on the full time/radius rectangles,

    m(x,r)>1.85e30 on [1/1000,1/200] x [73/64,585/512],
    m(x,r)<-1.85e30 on [1/1000,1/200] x [951/512,119/64].

All endpoints are included, and every rectangle lies within D. The physical
time interval is exactly `10^-9<=t<=5*10^-9`. The checker also retains the
unrounded rational enclosures, which are slightly stronger than these rounded
strict bounds.

For every fixed time in this interval, radial continuity gives at least one
zero of m in the open radius bracket `(585/512,951/512)`. No uniqueness,
continuous zero curve, shell zero, or pointwise angular sign follows. The two
signed means do imply that the scalar field takes the corresponding sign at
some sphere point in each cell, but the certificate's stated object is the
sphere mean.

The premises retain the complete bulk source, dynamic boundary rows, Robin
lift and all metric/scalar couplings. This proof uses the resulting global
derivative bound; it makes no local free-wave replacement. In particular it
does not discard a source or boundary term simply because no waveform was
computed.

`check_evolved_scalar_signs.py` authenticates frozen local copies with
`INPUT_PINS.json`, checks that the witness and refined bound use the same actual
initial report, reconstructs the cheap exact norm join, and checks the new
projection, time conversion, domain and sign inequalities with rational
arithmetic. It does not replay the expensive parent profile or all-point
operator proofs; those remain explicit accepted conditional premises. The
original locations in the pin file are provenance metadata, not runtime
dependencies. The scripts require only the Python standard library.

`run_sign_controls.py` checks normal/optimized byte equality and seven
reason-specific invalid variants: an underestimated mean functional norm,
an incorrect derivative weight, an omitted drift integral, wrong physical
time scaling, an invented initial margin, expanded unauthenticated radius
coverage, and a time interval with insufficient sign room. No Python assert
statement controls acceptance.

This is quantitative persistence for the selected formal coefficient. It is
not an independent waveform solve, a complete numerical residual certificate,
a gauge-invariant observable or force, a finite-amplitude nonlinear family or
remainder estimate, or a claim of external novelty.

For comparison, the full fixed point also gives
`|mean_w(U-b)|<=3q||U||_X<6.51681e30`, where
`b=U_free+KF+LG`, `q=2109/1000000`. The three comes from the X value weight,
and the mean still has norm one. To transfer initial signs by that weaker
route one must separately bound the initial-velocity and bulk Duhamel terms,
keep both characteristic feet inside a witness cell, and use `LG=0` only in
the strict region `r>x`. None of those additional local claims is needed for
the derivative proof above.
