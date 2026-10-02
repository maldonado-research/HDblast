# Post hoc interpretation of the Simpson ledger's frequency response

This explanatory note was prepared after the public freeze, while the
registered source-free diagnosis was running and after a provisional primary
result had been reported. It adds no criterion, case, tolerance or numerical
algorithm. No retained physical array was read or rerun. Its numerical checks
use fabricated exponentials; its conclusions about those exponentials do not
constitute a quantitative proof of the actual ledger error. These are standard
quadrature and sampling identities, with no novelty claim.

Let f(t)=A exp(i Omega(t-t_a)), sampled at t_j=t_a+j Delta, with N=2M>0
subintervals, theta=Omega Delta and q=exp(i theta). In exact arithmetic,
composite Simpson integration and the exact integral are

```
Q = (A Delta/3)*(1+4q+q^2)*sum(q^(2m),m=0,...,M-1)
I = A Delta*(q^N-1)/(i theta),                 theta != 0.
```

If q^2!=1, the geometric sum gives

```
Q = (A Delta/3)*(1+4q+q^2)*(q^N-1)/(q^2-1).
```

Where I!=0 and sin(theta)!=0, cancellation yields the phase-independent
normalized response

```
Q/I = T(theta) = theta*(2+cos(theta))/(3*sin(theta)).
```

The response is even in theta and independent of the number of complete
two-cell blocks. The same multiplier applies to a constant-amplitude sine,
cosine or derivative harmonic by linearity. A zero real integral can still
make its literal real-valued ratio undefined.

The exceptions must be distinguished:

- At theta=0, f is constant and Q=I=A N Delta. The removable response is T(0)=1.
- At a nonzero global cancellation frequency q^N=1 with sin(theta)!=0,
  both Q and I vanish. Literal Q/I is0/0; T gives its finite analytic extension.
- At theta=m pi!=0, N even makes I=0. For even m, q=1 and
  Q=A N Delta; for odd m, q=-1 and Q=-A N Delta/3. These are genuine alias
  failures of the sampled integral. T has a pole, so a normalized ratio is
  inappropriate. No claim is made that an actual retained node lies at a pole.

The small-frequency expansion is

```
T(theta)=1+theta^4/180+theta^6/1512+theta^8/14400
           +17 theta^10/2395008+O(theta^12).
```

Thus the leading relative error of this particular pure harmonic is fourth
order. Halving Delta reduces that leading term by16 only in the small-theta
regime at fixed Omega. The expansion does not justify a factor15 Richardson
certificate for the actual ledger. Its nearest nonremovable poles are at
theta=+-pi, so the Taylor expansion about zero is unsuitable for theta=4.

For the canonical mode frequency Omega=2k, the fixed cutoff bounds give the
following illustrations. Actual quadrature nodes satisfy k<K, so the listed
theta is a nominal upper bound, not a saved node; T is evaluated on the
abstract pure exponential at that bound.

| Setting | Delta | K | Nominal theta=2K Delta | Pure-exponential T at the bound |
| --- | --- | ---: | ---: | ---: |
| Coarse | 1/128 | 64 | 1 | 1.006294276 |
| Coarse | 1/128 | 128 | 2 | 1.161228524 |
| Coarse | 1/128 | 256 | 4 | -2.372008351 |
| Fine | 1/256 | 64 | 0.5 | 1.000357835 |
| Fine | 1/256 | 128 | 1 | 1.006294276 |
| Fine | 1/256 | 256 | 2 | 1.161228524 |

The coarse K=256 band crosses the canonical Nyquist point theta=pi,
k=pi/(2 Delta)=64 pi. The fine band is below that point but extends to
theta=2, which is not uniformly a small-theta regime. Being below Nyquist
therefore does not establish accurate Simpson integration. Applying double-step
Simpson to the same fine trajectory doubles each theta; it reaches the same
nominal1,2,4 bounds without changing that trajectory.

The source-free stress ledger is not a single constant-amplitude exponential.
Its oscillatory terms have rational L(t)=-1/t envelopes and its unrestricted
Re(c) term also varies with L. A frozen-envelope frequency response is therefore
an interpretive heuristic. Since L'=L^2, the fractional envelope rate for
L^n is nL: at high k this can be slow compared with2k and motivate a local
frozen-envelope picture. That scale separation need not cover small k and
does not turn the global transfer into an exact formula. Envelope derivatives, endpoint phases, momentum
superposition and signed cancellation determine the actual absolute error.
They prevent assigning the values in the table as relative errors or bounds
for a stored case. A small exact integral can also make a relative error
misleading; the registered signed absolute decomposition remains the test.

In exact arithmetic, subtracting the two even global Simpson prefixes equals
Simpson integration on the selected interval. The registered S_ab instead
subtracts long-double prefixes in long double, which can introduce additional
rounding effects. Its separately retained direct-reset diagnostic covers that
distinction. The pure-exponential transfer does not quantify those effects.

The quantitative evidence, when the registered checks finish, must come from
the unchanged signed D_S=D_cont+E_Q decomposition, all sampled R/P consistency
checks, endpoint flow projections and fixed80/100 arithmetic checks. This note
does not establish mode accuracy inherited before the anchor, the active-source
ledger, an old calibration PASS or a new physical result.

The accompanying NORMAL.json and OPTIMIZED.json each verify31 exact SymPy
identities/series statements and10 fabricated exponential cases at100decimal
digits. The proof script's SHA256 is recorded in both receipts. The .py script
belongs outside the frozen checkpoint; only this report and its JSON evidence
are eligible for the declared post-freeze report directories.
