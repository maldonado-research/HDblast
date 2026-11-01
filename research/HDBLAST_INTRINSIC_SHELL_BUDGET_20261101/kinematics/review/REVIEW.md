# Independent review of the frozen actual shell trace aggregation

The independent exact Fraction reconstruction passes. It confirms strict
positivity of the intrinsic curvature-defect second derivative, expansion
variance second derivative, and Buchert coefficient on the complete frozen
scaled-time interval `[31/25,63/50]` and at the original point `x=5/4`.
There is no strictly positive certified lower bound for squared shear.

This review treats the frozen actual scalar jets, background boxes, and
inherited continuous source-error boxes as a **shared scientific premise**.
It does not rerun the oracle or independently prove its PDE/error theorem.
It imports no producer, oracle or third-party module. The normalized input
projection pins all raw exact boxes and inherited input hashes/byte counts;
private historical origin strings and execution metadata are excluded.
This permits the same review of a sanitized portable report containing
identical raw scientific inputs. The original full report hash is recorded
for provenance in `SOURCE_PINS.json`.

The reviewer uses its own unrounded rational interval arithmetic. It
reconstructs the trace coefficients, all 64 forcing cells, global feedback
bounds and 1,728 exact kernel weights and terms. Its binomial polynomial
integration is algebraically independent of the producer antiderivative.
Every reconstructed target field and continuous hull lies inside the
producer's outward dyadic enclosure. The smaller independent upper shear
bound reflects avoiding dyadic rounding near zero; it supplies no strict
lower sign. In particular this check verifies enclosures, not merely
agreement of sign flags or overlap of reported intervals.

For the actual shell equation

    u_xx+a(x)u_x+b(x)u=f(x), u(0)=u_x(0)=0,
    a=3ell*tanh(t), b=ell²[8/cosh²(t)−4],
    f=ell²rho_b²F'(phi_b)*psi/6,
    F'=intersection(2p_b B_b, 2V'_b−8a_b p_b+sigma''_b p_b),

put T=63/50, amax=sup|a|, bmax=sup|b| and F1=sum(cellwidth*sup|f|).
Integrating the equation once, using |u(x)|<=x sup|u_x|, gives

    Vcap=F1/[1−amax*T−bmax*T²/2], Ucap=T*Vcap,
    G=amax*Vcap+bmax*Ucap.

The denominator is strictly positive. Whole-history integration from zero
retains precontact scalar errors. The remainders after integrating the
forcing alone are bounded by xG, x²G/2 and x³G/6 for u_x, u and integral u.
Positive polynomial kernel weights are integrated exactly over every full
or truncated source cell; target-cell endpoint weights bound every
intermediate target time. No point quadrature or precontact zero shortcut
is used.

The geometric transfer uses proper time tau=rho_b*t, P=rho_b*ell*integral u dx,
A=rho_b*cosh(t), H=tanh(t)/rho_b, H'=1/A², u'=u_x/(rho_b*ell),
E=u−HP, E'=u'−Hu−H'P, L=P/A² and theta1=3E'+8L. The verified formulas are

    C''=(10/3)E²/A²,
    <|shear|²>''=(80/9)L²,
    Var(theta)''=theta1²/3,
    Q_B''=(2/9)[theta1²−40L²].

All suffixes `second` denote actual second parameter derivatives, with
coefficient epsilon²/2. Conditional exact-family existence, actual-first
identification and proper-clock history remain premises. The coefficient
signs do not give a numerical finite-amplitude radius, a shell energy,
hot matter, reheating, or cosmological-origin conclusion.

Run from this directory with an authenticated report containing the frozen
raw inputs:

```sh
python -I -B review_trace.py --raw-report /path/to/trace-report.json --output /tmp/review-normal.json
python -I -O -B review_trace.py --raw-report /path/to/trace-report.json --output /tmp/review-optimized.json
python -I -B run_controls.py --raw-report /path/to/trace-report.json --output /tmp/review-controls.json
```

Normal and optimized science reports are byte-identical. Eleven deliberate
mutations reject wrong kernel powers, point rectangle weights, a missing
integration half, improper time factors, an epsilon²/second-derivative
factor, missing feedback, a history gap, dropped actual-error metadata,
precontact zero shortcuts and a shrunk raw source box. Five independent
polynomial fixtures cover full, truncated and inactive cells, and the
interval square guard checks a box crossing zero.
